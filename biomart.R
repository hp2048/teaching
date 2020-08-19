library(biomaRt)
library(karyoploteR)
library(GenomicRanges)
library(pheatmap)
library(tidyverse)

##list datasets availables for use from ensembl biomart
emart <- useMart("ensembl")
listDatasets(emart)

##define the required connections to the human gene related information database
hsmart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
##list what information is available for human genes
hsattr <- listAttributes(hsmart) 
##check data
head(hsattr)
##check what homology information is available for human genes with the chicken genome
hsattr %>% filter(page=="homologs" & grepl("gallus", name))

##get chicken information column names
chicken_columns <- hsattr %>% filter(page=="homologs" & grepl("gallus", name)) %>% pull(name)

chickenhomologs <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", chicken_columns), mart = hsmart)
colnames(chickenhomologs)

##individual steps
chickenhomologs <- filter(chickenhomologs, ggallus_homolog_orthology_type == "ortholog_one2one")
hs2gg <- table(chickenhomologs$chromosome_name, chickenhomologs$ggallus_homolog_chromosome)
pheatmap(hs2gg, cluster_rows = F, cluster_cols = F)

##all in one
chickenhomologs %>% 
filter(ggallus_homolog_orthology_type=="ortholog_one2one") %>% 
group_by(chromosome_name, ggallus_homolog_chromosome) %>% 
dplyr::summarise(genecounts = n()) %>% 
ggplot(aes(x=chromosome_name, y = ggallus_homolog_chromosome, fill = genecounts)) + 
geom_tile(stat = "identity") + 
scale_fill_viridis_b() + 
theme_bw()

##ggmart <- useMart("ensembl", dataset="ggallus_gene_ensembl")

hsproteingenes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position"),
      filters = "transcript_biotype",
      values = c("protein_coding"),
      mart = hsmart)
hsproteingenes$chromosome_name <- paste("chr",hsproteingenes$chromosome_name,sep="")

hsgenesranges <- makeGRangesFromDataFrame(hsproteingenes, start.field = "start_position", end.field = "start_position", seqnames.field = "chromosome_name")

kp <- plotKaryotype(genome = "hg38")
kpPlotDensity(kp, data=hsgenesranges, col = "brown", window.size = 1000000)
kp <- plotKaryotype(genome = "hg38",chromosomes = c("chr1"))
kpPlotDensity(kp, data=hsgenesranges, col = "brown", window.size = 1000000)


trackNames(ucscTableQuery(session))

