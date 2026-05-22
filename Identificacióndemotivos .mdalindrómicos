
Para encontrar el inicio de la transcripción de un gen RStudio

mapa_genoma <- import("C:/Users/gemag/OneDrive/Documentos/Proyecto de Maestria/Genoma Fredii/ncbi_dataset/ncbi_dataset/data/GCF_000018545.1/genomic.gtf") %>%
  as.data.frame()

# 3. Buscar las coordenadas de un gen específico (ej. NGR_RS19460)
inicio_transcripcion <- mapa_genoma %>%
  # Filtramos para que solo busque la línea principal del gen
  filter(type == "gene") %>% 
  # Buscamos tu gen por su código
  filter(grepl("NGR_RS19460", gene_id)) %>% 
  # Seleccionamos solo las columnas que nos importan
  select(seqnames, start, end, strand, gene_id)

# 4. Ver el resultado
print(inicio_transcripcion)
