library(biomaRt)
library(karyoploteR)
library(GenomicRanges)


hsmart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

chickenhomologs <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "ggallus_homolog_ensembl_gene", "ggallus_homolog_chromosome", "ggallus_homolog_orthology_type"), mart = hsmart)
chickenhomologs <- chickenhomologs[chickenhomologs$ggallus_homolog_orthology_type=="ortholog_one2one",]

hs2gg <- table(chickenhomologs$chromosome_name, chickenhomologs$ggallus_homolog_chromosome)
pheatmap(hs2gg, cluster_rows = F, cluster_cols = F)

ggmart <- useMart("ensembl", dataset="ggallus_gene_ensembl")

hsproteingenes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position"),
      filters = "transcript_biotype",
      values = c("protein_coding"),
      mart = hsmart)
hsproteingenes$chromosome_name <- paste("chr",hsproteingenes$chromosome_name,sep="")

hsgenesranges <- makeGRangesFromDataFrame(hsproteingenes, start.field = "start_position", end.field = "start_position", seqnames.field = "chromosome_name")

kp <- plotKaryotype(genome = "hg38")
kpPlotDensity(kp, data=hsgenesranges, col = "brown", window.size = 200000, r0=0.2, r1=1)

