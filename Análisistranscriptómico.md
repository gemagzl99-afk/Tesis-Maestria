MUTANTES VS WT

Cargar librerías requeridas

```
library(dplyr)
library(tidyr)
install.packages("ggplot2")
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(pheatmap)
```

Establecer el directorio de trabajo 

```
setwd("C:/Users/gemag/OneDrive/Documentos/Proyecto de Maestria/Novogene/FeatureCounts")
```

Importar matriz de conteos
```
counts_data <- read.table(
  "counts_matrix.txt",
  header = TRUE,
  sep = "\t",
  skip = 1
)
```

Convertir Geneid a rownames y eliminar metadatos

```
count_matrix <- counts_data %>%
  column_to_rownames(var = "Geneid") %>%
  select(-(1:5))

head(count_matrix)
```

Metadatos de muestras

```
sample_names <- colnames(count_matrix)

conditions <- factor(c("WildType", "WildType", "Mutant", "Mutant"))

col_data <- data.frame(
  sample = sample_names,
  condition = conditions
)

rownames(col_data) <- sample_names

col_data

# Verificar orden
all(rownames(col_data) == colnames(count_matrix))
```
DESEQ2

```
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition
)

dds <- DESeq(dds)

# Comparación Mutant vs WildType
results_raw <- results(
  dds,
  contrast = c("condition", "Mutant", "WildType")
)
```

Procesar resultados

```
results_df <- as.data.frame(results_raw) %>%
  na.omit() %>%
  arrange(padj)

head(results_df, 10)

# Filtrar genes con expresión diferencial significativa

significant_genes <- results_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1.5)

print(paste("Número de genes DEGs significativos:", nrow(significant_genes)))

write.csv(
  significant_genes,
  "significant_genes_list.csv",
  row.names = TRUE
)
```

Preparar datos para Volcano Plot

```
# Convertir nombres de fila a columna
results_df <- results_df %>%
  rownames_to_column(var = "gene")

# Crear etiquetas sin "gene-"
results_df <- results_df %>%
  mutate(
    gene_label = gsub("^gene-", "", gene)
  )

# Clasificar genes según su significancia
results_df <- results_df %>%
  mutate(
    is_significant = case_when(
      padj < 0.05 & log2FoldChange >= 1.5 ~ "Up-regulated",
      padj < 0.05 & log2FoldChange <= -1.5 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Seleccionar los 10 genes más significativos
top_genes <- results_df %>%
  arrange(padj) %>%
  slice_head(n = 10)
```

VOLCANO PLOT

```
volcano_plot <- ggplot(
  results_df,
  aes(x = log2FoldChange, y = -log10(padj))
) +
  geom_point(
    aes(color = is_significant),
    size = 2.2,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = gene_label),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = 0.6,
    point.padding = 0.5,
    segment.color = "black",
    max.overlaps = Inf
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.25))   
  ) +
  scale_color_manual(
   
    values = c(
      "Down-regulated" = "#6184d8",    
      "Up-regulated" = "#e8384f",      
      "Not significant" = "#4d4d4d"    
    ),

    labels = c(
      "Down-regulated" = "Down",
      "Up-regulated" = "Up",
      "Not significant" = "Not significant"
    )
  ) +
  labs(
    x = expression(Log[2]~"(fold change)"),
    y = expression(-Log[10]~italic(P)),
    color = NULL
  ) +
  
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 22,
      color = "#0B4C8A"
    ),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.line = element_line(color = "black", linewidth = 1.2),
    axis.ticks = element_line(color = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    
    legend.position = "top",                 
    legend.direction = "horizontal",         
    legend.background = element_blank(),     
    legend.key = element_blank(),            
    legend.text = element_text(size = 16)
  )

print(volcano_plot)
```

Limpiar nombres de muestras en dds

```
colnames(dds) <- gsub("\\.dedup\\.bam$", "", colnames(dds))

# Transformación rlog

# Usa el objeto dds 
rld <- rlog(dds, blind = FALSE)

# Resultados de DESeq2
res <- results(dds)

# Eliminar genes con padj = NA
res <- res[!is.na(res$padj), ]

# Ordenar por significancia (padj)
res <- res[order(res$padj), ]

# Seleccionar TOP 50 Genes

top_genes_names <- rownames(res)[1:50]
```

Extraer datos para HEATMAP

```
heatmap_data <- assay(rld)[top_genes_names, ]

#  Escalar (Z-score por gen)
scaled_heatmap_data <- t(scale(t(heatmap_data)))

# Anotación de muestras

annotation_col <- data.frame(
  Condition = colData(dds)$condition
)

rownames(annotation_col) <- colnames(scaled_heatmap_data)

ann_colors <- list(
  Condition = c(
    WildType = "#1f78b4",
    Mutant   = "#e31a1c"
  )
)
```
HEATMAP

```
# Quitar "gene-"

tiff(
  filename = "heatmap_S_fredii.tiff", 
  width = 10,           
  height = 12,         
  units = "in", 
  res = 300,           
  compression = "lzw"
)


rownames(scaled_heatmap_data) <- gsub("gene-", "", rownames(scaled_heatmap_data))


ann_colors <- list(
  Condition = c(WildType = "#74c476", Mutant = "#d95f02") 
)


pheatmap(
  scaled_heatmap_data,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 14,       
  fontsize_row = 10,    
  fontsize_col = 14,   
  cellheight = 10,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
```

PCA

```
pca_data <- plotPCA(
  rld,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(100 * attr(pca_data, "percentVar"))

#  PCA con etiquetas de las muestras
ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  
  geom_point(size = 4) +
  
  geom_text_repel(
    aes(label = name),
    size = 5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  
  
  scale_color_manual(
    values = c(
      WildType = "#1f78b4",
      Mutant   = "#e31a1c"
    )
  ) +
  
  theme_minimal(base_size = 10) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 16),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.2
    ),
    
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

```



MUTANTES VS AWT



Cargar librerías requeridas

```
library(dplyr)
library(tidyr)
install.packages("ggplot2")
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(pheatmap)
```

Establecer el directorio de trabajo 

```
setwd("C:/Users/gemag/OneDrive/Documentos/Proyecto de Maestria/Novogene/FeatureCounts")
```

Importar matriz de conteos
```
counts_data <- read.table(
  "counts_matrix.txt",
  header = TRUE,
  sep = "\t",
  skip = 1
)

count_matrix <- counts_data %>%
  column_to_rownames(var = "Geneid") %>%
  select(-(1:5))
```

Excluir muestra BWT3

```
count_matrix <- count_matrix %>%
  select(-contains("BWT3"))

colnames(count_matrix)
```

Metadatos de muestras

```
sample_names <- colnames(count_matrix)

conditions <- factor(
  c("WildType","Mutant","Mutant"),
  levels = c("WildType","Mutant")
)

col_data <- data.frame(
  sample = sample_names,
  condition = conditions
)

rownames(col_data) <- sample_names

# verificar
all(rownames(col_data) == colnames(count_matrix))
```
DESEQ2

```
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition
)

dds <- DESeq(dds)
```

Procesar resultados

```
res <- results(dds, contrast = c("condition","Mutant","WildType"))

results_df <- as.data.frame(res) %>%
  na.omit() %>%
  arrange(padj)

head(results_df)

significant_genes <- results_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1.5)

print(paste("Número de DEGs:", nrow(significant_genes)))

write.csv(
  significant_genes,
  "significant_genes_MUT_vs_WT_sin_BWT3.csv",
  row.names = TRUE
)

top50_genes <- results_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 50)

print(top50_genes)

write.csv(
  top50_genes,
  "Top50_Genes_Significativos_Sin_BWT3.csv",
  row.names = TRUE
)
```

Preparar datos para volcano plot

```
results_df <- results_df %>%
  rownames_to_column(var = "gene")

# quitar prefijo gene-
results_df <- results_df %>%
  mutate(gene_label = gsub("^gene-","",gene))

# clasificación de genes
results_df <- results_df %>%
  mutate(
    is_significant = case_when(
      padj < 0.05 & log2FoldChange >= 1.5 ~ "Up-regulated",
      padj < 0.05 & log2FoldChange <= -1.5 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Selección de los 10 genes más significativos

top10_genes <- results_df %>%
  arrange(padj) %>%
  slice_head(n = 10)

print(top10_genes$gene_label)
```

Volcano Plot

```
volcano_plot <- ggplot(results_df,
                       aes(x = log2FoldChange,
                           y = -log10(padj))) +
  
  geom_point(
    aes(color = is_significant),
    size = 2.2,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    data = top10_genes,
    aes(label = gene_label),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = 0.6,
    point.padding = 0.5,
    segment.color = "black",
    max.overlaps = Inf
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  
  geom_vline(
    xintercept = c(-1.0, 1.0),
    linetype = "dashed",
    color = "black"
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.25))   
  ) +
  
  scale_color_manual(
    values = c(
      "Down-regulated" = "#6184d8", 
      "Up-regulated" = "#e8384f",   
      "Not significant" = "#4d4d4d" 
    )
  ) +
  
  labs(
    x = expression(Log[2]~"(fold change)"),
    y = expression(-Log[10]~italic(P)),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 22, color = "#0B4C8A"
    ),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", linewidth = 1.2),
    axis.ticks = element_line(color = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    

    legend.position = "top",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 16)
  )

print(volcano_plot)

```

Limpiar nombres de muestras en dds

```
# Limpiar nombres de muestras en dds
colnames(dds) <- gsub("\\.dedup\\.bam$", "", colnames(dds))

# Transformation rlog
rld <- rlog(dds, blind = FALSE)

# Resultados de DESEQ2
res <- results(dds)

# Eliminar genes con padj = NA
res <- res[!is.na(res$padj), ]

# Ordenar por significancia (padj)
res <- res[order(res$padj), ]

# Selección top 50 genes
top_genes_names <- rownames(res)[1:50]
```

Extraer datos para HEATMAP

```
heatmap_data <- assay(rld)[top_genes_names, ]

scaled_heatmap_data <- t(scale(t(heatmap_data)))

annotation_col <- data.frame(
  Condition = colData(dds)$condition
)

rownames(annotation_col) <- colnames(scaled_heatmap_data)

ann_colors <- list(
  Condition = c(
    WildType = "#1f78b4",
    Mutant   = "#e31a1c"
  )
)
```

HEATMAP

```
# LIMPIAR NOMBRES DE GENES (Quitar "gene-")

tiff(
  filename = "heatmap_S_fredii.tiff", 
  width = 10,
  height = 12,          
  units = "in", 
  res = 300,            
  compression = "lzw"
)

rownames(scaled_heatmap_data) <- gsub("gene-", "", rownames(scaled_heatmap_data))

ann_colors <- list(
  Condition = c(WildType = "#74c476", Mutant = "#d95f02") 
)

pheatmap(
  scaled_heatmap_data,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 14,        
  fontsize_row = 10,   
  fontsize_col = 14,   
  cellheight = 10,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
```

PCA

```
pca_data <- plotPCA(
  rld,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(100 * attr(pca_data, "percentVar"))

# PCA con etiquetas de muestra

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  
  geom_point(size = 4) +
  
  geom_text_repel(
    aes(label = name),
    size = 5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  
  
  scale_color_manual(
    values = c(
      WildType = "#1f78b4",
      Mutant   = "#e31a1c"
    )
  ) +
  
  theme_minimal(base_size = 10) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 16),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.2
    ),
    
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )
```














