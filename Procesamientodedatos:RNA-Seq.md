Existen dos modalidades para acceder al servidor de Ómica del CICESE: de manera local o remota. Para la conexión en la institución, basta con utilizar una conexión vía Ethernet e iniciar sesión en la terminal ingresando el nombre de usuario asignado seguido de la contraseña correspondiente. 

```
gema.garcia - -> Usuario proporcionado
```

Por otro lado, para el acceso remoto es necesario instalar y configurar el software de VPN institucional siguiendo los manuales del CICESE. Una vez establecida la conexión, se debe acceder a la terminal y ejecutar el siguiente comando:

```
ssh  gema.garcia@omica
```
Cabe señalar que, por motivos de seguridad en sistemas Unix/Linux, los caracteres de la contraseña no se visualizan en pantalla durante su ingreso. 


Con el fin de mantener la organización institucional del CICESE, el trabajo se realizó dentro del directorio de almacenamiento correspondiente al grupo de investigación asignado. Para acceder directamente a este directorio de trabajo, se ejecutó el siguiente comando:

```
cd /LUSTRE/bioinformatica_data/virolab/
```

Pasar de mi PC los archivos de los datos de Novogene o bien los archivos con terminación fq.gz a la terminal. Para ello, se debe abrir una pestaña nueva en la terminal local, donde se visualice el nombre de usuario de la computadora. Después, ejecutar el comando de transferencia hacia la ruta correspondiente en el servidor Ómica.
Ejemplo

```
scp -r "G:\Proyecto de Maestria\Novogene\Datos Novogene" gema.garcia@omica:/LUSTRE/bioinformatica_data/virolab/gemag
scp [ruta_del_archivo_en_pc] [usuario]@[servidor]:[ruta_destino_en_servidor]
```

Posteriormente, se verificó la integridad de la transferencia mediante el archivo MD5.txt, el cual contiene códigos hash MD5 utilizados para confirmar que los archivos no sufrieron daños ni modificaciones durante la transferencia. Cada archivo original posee una "huella digital" única generada mediante el algoritmo MD5; por lo tanto, cualquier cambio en el archivo, incluso de solo un byte, produce una modificación en dicho código. Para realizar esta verificación se ejecutó el siguiente comando.

```
md5sum -c MD5.txt
```
La confirmación “Ok” en los archivos indica que estos se encuentran en perfecto estado. 
Es importante mencionar que los datos obtenidos de secuenciación constan de dos archivos por muestra; las lecturas Forward (1) y Reverse (2). Esta configuración es característica de Paired-End, la cual mejora la resolución del análisis transcriptómico al proporcionar información posicional de ambos extremos de la biblioteca de insertos, resultando en un alineamiento más confiable contra el genoma de referencia


## Limpieza de lecturas

Para la etapa de limpieza de lecturas, se utilizó fastp, con el objetivo de eliminar lecturas duplicadas (--dedup), identificar y eliminar adaptadores presentes en los extremos de lecturas paired-end (--detect_adapter_for_pe), generar archivos de lecturas limpias y mejorar la calidad general de las secuencias mediante el filtrado y recorte de bases de baja calidad. 

Se procedió a la creación de un archivo mediante Vim que sirve para crear y editar scripts. Para iniciar la edición del script, se ejecutó el siguiente comando:

```
vim
```
Comandos de función para Vim.
```
:w    Guardar el archivo
:q     Sale del editor 
:wq   Guardar y salir del vim
Los scripts se guardan con terminación .slrm
```


## Fastp

Script usado en Vim.

```
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=fastp
#SBATCH -o fastp.%N.%j.out
#SBATCH -e fastp.%N.%j.err
#SBATCH -N 1
#SBACTH -n 24
#SBATCH --mem=100GB

#ruta a fastp

FASTP="/LUSTRE/bioinformatica_data/virolab/gemag/fastp"

#directorio con archivos fq.gz

DATA="/LUSTRE/bioinformatica_data/virolab/gemag/DatosNovogene/01.RawData"
OUTDIR="/LUSTRE/bioinformatica_data/virolab/gemag/FastClean"

#ejecutar fastp

$FASTP/./fastp \
 -i $DATA/AWT/AWT_1.fq.gz -I $DATA/AWT/AWT_2.fq.gz \
 -o $OUTDIR/AWT_1_clean.fq.gz -O $OUTDIR/AWT_2_clean.fq.gz \
 --dedup --detect_adapter_for_pe \
 --thread 16

$FASTP/./fastp \
 -i $DATA/BWT3/BWT3_1.fq.gz -I $DATA/BWT3/BWT3_2.fq.gz \
 -o $OUTDIR/BWT3_1_clean.fq.gz -O $OUTDIR/BWT3_2_clean.fq.gz \
 --dedup --detect_adapter_for_pe \
 --thread 16

$FASTP/./fastp \
 -i $DATA/CMUT/CMUT_1.fq.gz -I $DATA/BWT3/CMUT_2.fq.gz \
 -o $OUTDIR/CMUT_1_clean.fq.gz -O $OUTDIR/CMUT_2_clean.fq.gz \
 --dedup --detect_adapter_for_pe \
 --thread 16


$FASTP/./fastp \
 -i $DATA/DMUT2/DMUT2_1.fq.gz -I $DATA/DMUT/DMUT2_2.fq.gz \
 -o $OUTDIR/DMUT2_1_clean.fq.gz -O $OUTDIR/DMUT2_2_clean.fq.gz \
 --dedup --detect_adapter_for_pe \
 --thread 16
```
Revisar si se está ejecutando el script utilizando el siguiente comando:
```
squeue
```
Con el objetivo de optimizar la inspección de los parámetros de calidad y obtener una visión global del set de datos, se emplea MultiQC. Esta herramienta integra los reportes individuales en un único informe interactivo en formato HTML, facilitando la identificación de tendencias o anomalías en el conjunto de las muestras. Script usado en vim.

```
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=multiqc
#SBATCH --output=mul-%N.%j.log
#SBATCH --error=mul-%N.%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

export PATH=/LUSTRE/apps/Anaconda/2023/miniconda3/bin:$PATH
source activate qiime2-2023.2

multiqc /LUSTRE/bioinformatica_data/virolab/gemag/FastQC_Reports/*.html -o /LUSTRE/bioinformatica_data/virolab/gemag/FastQC_Reports --data-format json --export
```

Para comparar la calidad entre los datos crudos y las lecturas después de la limpieza, se realizó un análisis con FastQC de los datos originales y, posteriormente, se integraron los reportes mediante MultiQC para obtener una evaluación global y facilitar la comparación entre ambos conjuntos de datos. 

## FastQC

Script usado en Vim. 
```
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc-%j.log
#SBATCH --error=fastqc-%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 6-00:00:00

export PATH=$PATH:/LUSTRE/apps/bioinformatica/FastQC_v0.11.7/:$PATH
module load gcc-7.2.0

fastqc /LUSTRE/bioinformatica_data/virolab/gemagi/DatosNovogene/01.RawData/*.fq.gz -t 8 -o /LUSTRE/bioinformatica_data/virolab/gemag/DatosNovogene/01.RawData/fastqc_results
```

## MultiQC de datos originales

Script usado en Vim.

```
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=multiqc
#SBATCH --output=mul-%j.log
#SBATCH --error=mul-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

export PATH=/LUSTRE/apps/Anaconda/2023/miniconda3/bin:$PATH
source activate qiime2-2023.2

multiqc /LUSTRE/bioinformatica_data/virolab/gemag/DatosNovogene/01.RawData/fastqc_results/*.zip -o /LUSTRE/bioinformatica_data/virolab/gemag/DatosNovogene/01.RawData/fastqc_results --data-format json --export
```


## Descarga y transferencia del genoma de referencia

Se descargó el genoma completo de Sinorhizobium fredii NGR234 desde la base de datos de NCBI en formato Download Package. Posteriormente, se creó la carpeta GenomaFredii para almacenar los archivos del genoma de referencia y, finalmente, estos archivos fueron transferidos desde la computadora local hacia la terminal de trabajo en Linux.


## Indexación y alineamiento

Script usado en Vim.

```
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=HISAT2_Full_Pipeline
#SBATCH --output=hisat2_pipeline-%j.log
#SBATCH --error=hisat2_pipeline-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB           
#SBATCH --ntasks-per-node=24  
#SBATCH -t 06-00:00:00

module load conda-2023

conda activate hisat2_env

module load trinityrnaseq-v2.15.1

# Directorio donde está el archivo .fna del genoma
GENOME_DIR="/LUSTRE/bioinformatica_data/virolab/gemag/GenomaFredii"
GENOME_FILE="GCF_000018545.1_ASM1854v1_genomic.fna"
INDEX_PREFIX="S_fredii_index"

# Directorio donde están tus archivos FASTQ limpios (salidas de Fastp)

# Directorio de salida para los archivos BAM
OUTPUT_DIR="/LUSTRE/bioinformatica_data/virolab/gemag/HISAT2_BAMs"
mkdir -p "$OUTPUT_DIR"

# ---  Indexación del Genoma ---
echo "--- Iniciando Indexación del Genoma ---"
# Se ejecuta la indexación en el directorio del genoma.
cd $GENOME_DIR
hisat2-build $GENOME_FILE $INDEX_PREFIX
echo "--- Indexación Completa ---"
cd - # Volver al directorio de trabajo anterior


# --- Lista de Muestras ---
SAMPLES="AWT BWT3 CMUT DMUT2"

# Bucle para Mapeo, Alineamiento y Ordenamiento
for SAMPLE in $SAMPLES
do
    echo "--- Procesando muestra: $SAMPLE ---"

 # Construcción de las rutas completas

 R1_READS="${INPUT_DIR}/${SAMPLE}_1_clean.fq.gz" # Resultado: AWT_1_clean.fq.gz (CORRECTO)
    R2_READS="${INPUT_DIR}/${SAMPLE}_2_clean.fq.gz" # Resultado: AWT_2_clean.fq.gz (CORRECTO)

    # --- A. Alineamiento con HISAT2 ---
    hisat2 -x $GENOME_DIR/$INDEX_PREFIX \
    -1 $R1_READS \
    -2 $R2_READS \
    -S $OUTPUT_DIR/$SAMPLE.sam \
    --threads 12

    # --- B. Conversión SAM a BAM y Ordenamiento (SAMtools) ---
    samtools sort -o $OUTPUT_DIR/$SAMPLE.sorted.bam $OUTPUT_DIR/$SAMPLE.sam -@ 12

    # --- C. Limpieza ---
    rm $OUTPUT_DIR/$SAMPLE.sam

    echo "Finalizado: $SAMPLE.sorted.bam creado en $OUTPUT_DIR"
done

echo "--- Pipeline de Alineamiento Completo ---"

```

## Eliminación de duplicados

Los datos obtenidos tras el alineamiento se exportaron en formato BAM (Binary Alignment Map) y posteriormente se procesaron para eliminar duplicados de PCR generados durante la construcción de las bibliotecas, utilizando el comando markdup. Script usado en Vim. 


```
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=Eliminacion_duplicado
#SBATCH --output=eliminacion_duplicado-%j.log
#SBATCH --error=eliminacion_duplicado-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=12
#SBATCH -t 06-00:00:00


SAMTOOLS="/home/gema.garcia/miniconda3/envs/bioinfo/bin/samtools"


# Lista de archivos BAM ordenados (ARRAY)
BAM_FILES=("AWT.sorted.bam" "BWT3.sorted.bam" "CMUT.sorted.bam" "DMUT2.sorted.bam")

echo "Iniciando pre-procesamiento y eliminación de duplicados"

    echo "--- Procesando $BAM ---"

    # Reemplaza .sorted.bam por el nuevo sufijo
    NAME_SORTED="${BAM/sorted.bam/name_sorted.bam}"
    FIXMATE_BAM="${BAM/sorted.bam/fixmate.bam}"
    DEDUP_BAM="${BAM/sorted.bam/dedup.bam}"

    #  Ordenar por nombre
    echo "1. Ordenando por nombre: $NAME_SORTED"
    $SAMTOOLS sort -n -@ 12 -o "$NAME_SORTED" "$BAM"

    #  Arreglar mates y añadir etiquetas MC/MD
    echo "2. Arreglando mates (fixmate): $FIXMATE_BAM"
    $SAMTOOLS fixmate -m -@ 12 "$NAME_SORTED" "$FIXMATE_BAM"

    # Reordenar por posición 
    echo "3. Reordenando por posición"
    $SAMTOOLS sort -@ 12 -o "temp.$BAM" "$FIXMATE_BAM"

    # Eliminar duplicados
    echo "4. Eliminando duplicados: $DEDUP_BAM"
    $SAMTOOLS markdup -r -@ 12 "temp.$BAM" "$DEDUP_BAM"

    #  Indexación del nuevo BAM deduplicado
    echo "5. Indexando $DEDUP_BAM..."
    $SAMTOOLS index "$DEDUP_BAM"
#  Limpieza de archivos para ahorrar espacio
    echo "6. Limpiando archivos intermedios..."
    rm "$NAME_SORTED" "$FIXMATE_BAM" "temp.$BAM"

    echo "Proceso completado para $BAM."
done
echo "Todos los archivos BAM han sido procesados y deduplicados."

```


Posteriormente, se verificaron los archivos de salida generados tras el procesamiento de las lecturas, con el fin de confirmar que la ejecución se realizó correctamente y que los archivos fueron generados sin errores. Además, se comparó el tamaño y la cantidad de lecturas entre los datos originales y los datos procesados, permitiendo evaluar la limpieza.


```
ls -l *.dedup.bam
ls -lh AWT.sorted.bam AWT.dedup.bam
```


## Descarga y preparación del archivo de anotación (GFF)

Con el objetivo de identificar la localización exacta de los genes dentro del genoma de Sinorhizobium fredii NGR234, se descargó el archivo de anotación genómica en formato GFF (General Feature Format). Este archivo contiene información sobre la posición y organización de los elementos genéticos presentes en el genoma, incluyendo genes codificantes, regiones de ARN y otras características funcionales.


```
wget -P /LUSTRE/bioinformatica_data/virolab/gemag/GenomaFredii https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/545/GCF_000018545.1_ASM1854v1/GCF_000018545.1_ASM1854v1_genomic.gff.gz
gunzip /LUSTRE/bioinformatica_data/virolab/gemag/GenomaFredii/GCF_000018545.1_ASM1854v1_genomic.gff.gz

```


## Cuantificación de la expresión con FeatureCounts

Finalmente, se realizó la cuantificación de la expresión génica mediante la herramienta FeatureCounts, la cual permite asignar y contabilizar las lecturas alineadas a cada gen anotado en el genoma de referencia. Script usado en Vim. 

```
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=feature_counts
#SBATCH --output=$OUTPUT_DIR/featurecounts-%j.log
#SBATCH --error=$OUTPUT_DIR/featurecounts-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=12
#SBATCH -t 06-00:00:00

export PATH="$CONDA_ENV_BASE/bin:$PATH"

ANNOTATION_FILE="/LUSTRE/bioinformatica_data/virolab/gemag/GenomaFredii/GCF_000018545.1_ASM1854v1_genomic.gff"

# Lista de archivos BAM deduplicados para conteo
DEDUP_BAM_FILES=(
    "AWT.dedup.bam"
    "BWT3.dedup.bam"
    "CMUT.dedup.bam"
    "DMUT2.dedup.bam"
)
# Ruta completa para el archivo de salida de la matriz de conteo
OUTPUT_FILE="$OUTPUT_DIR/counts_matrix.txt"

# EJECUCIÓN DE FEATURECOUNTS
echo "Iniciando conteo de reads con featureCounts..."
echo "Los resultados se guardarán en la carpeta: $OUTPUT_DIR"
# featureCounts, con la corrección -g Parent
featureCounts \
-p \
-t exon \
-g Parent \
-a "$ANNOTATION_FILE" \
-o "$OUTPUT_FILE" \
-T 8 \
"${DEDUP_BAM_FILES[@]}"

echo "Conteo completado. Matriz de conteo guardada en $OUTPUT_FILE."
```


