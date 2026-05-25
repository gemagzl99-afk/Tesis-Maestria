
# Guía para identificar TSS, región promotora y motivos palindrómicos

## Objetivo

Identificar el posible inicio de transcripción (TSS), la región promotora y posibles motivos palindrómicos en la región upstream del gen `NGR_RS27720`.

---

## 1. Ubicarse en la carpeta de trabajo

```bash
cd /LUSTRE/bioinformatica_data/virolab/gemag/HISAT2_BAMs

En esta carpeta se encuentran los archivos de alineamiento en formato .bam.
```

## 2. Identificar las coordenadas del gen

Buscar el gen de interés en el archivo .gff:

```bash
grep "NGR_RS27720" /LUSTRE/bioinformatica_data/virolab/gemag/GenomaFredii/GCF_000018545.1_ASM1854v1_genomic.gff
```

Gen localizado en la cadena negativa

## 3. Localizar el programa samtools

```bash
find /home/gema.garcia/miniconda3/envs /home/gema.garcia/.conda/envs -type f -name "samtools" 2>/dev/null
```
Ejemplo

```bash
/home/gema.garcia/miniconda3/envs/bio-utils/bin/samtools
```

## 4. Revisar la cobertura cerca del extremo 5’ del gen

Se uso la muestra AWT (WT) como referencia:

```bash
/home/gema.garcia/miniconda3/envs/bio-utils/bin/samtools depth -r "NC_012587.1:3499870-3500100" AWT.dedup.bam | less
```
La cobertura seguía siendo alta, por lo que todavía no se había localizado claramente el inicio de transcripción. 

## 5. Ampliar la búsqueda de cobertura

```bash
/home/gema.garcia/miniconda3/envs/bio-utils/bin/samtools depth -r "NC_012587.1:3500100-3500500" AWT.dedup.bam | less
```
En esta región la cobertura volvió a aumentar, lo que sugirió que se estaba entrando en la transcripción de un gen vecino.

## 6. Revisar los genes vecinos

```bash
grep -C 15 "NGR_RS27720" /LUSTRE/bioinformatica_data/virolab/gemag/GenomaFredii/GCF_000018545.1_ASM1854v1_genomic.gff | grep -w "gene"
```

Se observó que el gen vecino NGR_RS27725 se encuentra en la cadena positiva +

```bash
NC_012587.1 RefSeq gene 3500380 3501414 . + .
```

Por lo tanto, NGR_RS27720 y NGR_RS27725 son genes divergentes.


## 7. Analizar la región intergénica

```bash
/home/gema.garcia/miniconda3/envs/bio-utils/bin/samtools depth -r "NC_012587.1:3499870-3500380" AWT.dedup.bam | less
```

Esta región corresponde al espacio donde probablemente se encuentran los promotores de ambos genes.

## 8. Encontrar el punto con menor cobertura

Para evitar revisar manualmente todas las líneas, se ordena la cobertura de menor a mayor

```bash
/home/gema.garcia/miniconda3/envs/bio-utils/bin/samtools depth -r "NC_012587.1:3499870-3500380" AWT.dedup.bam | sort -k3,3n | head -n 10
```

Resultado
```bash
NC_012587.1 3500094 36
NC_012587.1 3500095 36
NC_012587.1 3500096 36
NC_012587.1 3500097 36
NC_012587.1 3500098 36
NC_012587.1 3500099 36
NC_012587.1 3500100 36
```

El punto de menor cobertura se localizó aproximadamente entre:

```bash
3500094-3500100
```





# Búsqueda de genes específicos en análisis transcriptómicos (RStudio)

## Objetivo

Buscar genes específicos dentro de los resultados de expresión diferencial obtenidos en RNA-Seq.

---

## Código

```bash
busqueda_gen <- results_df %>%
  filter(gene_label == "NGR_RS19460") %>%
  select(gene_label, baseMean, log2FoldChange, padj, is_significant)

print(busqueda_gen)
```




# Guía para buscar motivos palindrómicos en RStudio

## Objetivo

Buscar motivos palindrómicos en una secuencia promotora con una estructura tipo:

```text
7 pb - 2 pb - 7 pb
```


## 1. Cargar la librería
library(Biostrings)


## 2. Cargar la secuencia

```bash
mi_secuencia <- DNAString("GTTCTAAGCTCTCTTTCTGTCTCTCCTCGTCGGGTGGAGGCGGAGCCGCGCGGACCGTCCCGCGAAGATCCTGTCCCCGGCCGTACCGGGGAGAGAGCGGAGGGCCAGAGACGTCAGACGGTTTTCGTCGTCGTTTTGGCCGTTGTTTTGCTATTGAGCGAGAAAGAGCGTGCCTGCCCGGCGTGGGCCGGCAAGTTCAGCGCAGCGATGGACGAATTGTTCGGTCGCATCTCGTTTCCTCGTTCGCCGCTGCTCTTACCGGCTTTTCCACAAGAGTCAAGCCGC")
```

## 3. Buscar palíndromos tipo 7-2-7

```bash
palindromos_encontrados <- findPalindromes(mi_secuencia, 
                                           min.armlength = 7, 
                                           max.looplength = 2, 
                                           max.mismatch = 2)
```


## 4. Convertir el resultado a tabla

```bash
tabla_resultado <- as.data.frame(palindromos_encontrados)
```

## 5. Mostrar resultados

```bash
cat("\n--- Palíndromos Encontrados ---\n")
print(tabla_resultado)
```









