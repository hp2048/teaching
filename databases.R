.libPaths(c(.libPaths(), "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))

library(biomaRt)
library(karyoploteR)
library(GenomicRanges)
library(pheatmap)
library(tidyverse)

##ensembl genome browser
##Study the evolution of chromosomes in the human genome compared to the chicken genome
##With X chromosome as an example
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

##with fractions
filter(chickenhomologs, ggallus_homolog_orthology_type == "ortholog_one2one") %>% 
group_by(chromosome_name, ggallus_homolog_chromosome) %>% 
summarise(homologs = n()) %>% 
group_by(chromosome_name) %>% 
summarise(chromosome_name, ggallus_homolog_chromosome, homologs, homologsF = homologs/sum(homologs)) %>% 
ggplot(aes(x=ggallus_homolog_chromosome, y = chromosome_name)) + 
geom_tile(aes(fill = homologsF)) + 
geom_text(aes(label=round(homologsF, 1))) + 
scale_fill_viridis_b() + 
theme_bw()


##ggmart <- useMart("ensembl", dataset="ggallus_gene_ensembl")
##learn about gene density in different regions of the human genome
##this shows us that genome is not random and it has a specific structure
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

###UCSC genome browser
###study gene expression data from human samples
###which genes are generally highly expressed
###which genes are low expressed
library(rtracklayer)
library(GenomicRanges)
library(gprofiler2)
session <- browserSession("UCSC")
genome(session) <- "hg38"
##learn about tracks available for use
trackNames(ucscTableQuery(session))
##learn about tables available from that track
tableNames(ucscTableQuery(session, track="gtexGeneV8"))
##get expression values for all genes from chr22
geneexp <- getTable(ucscTableQuery(session, track="gtexGeneV8", range=seqinfo(session), table="gtexGeneV8"))
##inspect returned information
colnames(geneexp)
head(geneexp)
##check what types of genes are expressed and their expression score distribution
geneexp %>% filter(score>0) %>% ggplot(aes(x=score,fill=geneType)) + geom_density(alpha=0.5) + facet_wrap(.~geneType) + theme(legend.position='none')
##subset the data for protein coding genes only
pcoding <- geneexp %>% filter(geneType == "protein_coding")
##see how gene expression scores are distributed
hist(pcoding$score, breaks = 100)
##identify four quantiles
quantile(pcoding$score)
##get gene ontology enrichment for low expressed genes
lowexp <- gost(filter(pcoding, score<327) %>% pull(name), organism = "hsapiens", custom_bg = pcoding$name)
gostplot(lowexp)
lowexp <- lowexp$result[,c(1:12)]
filter(lowexp, source == "GO:BP")
##get gene ontology enrichment for low expressed genes
highexp <- gost(filter(pcoding, score>492) %>% pull(name), organism = "hsapiens", custom_bg = pcoding$name)
gostplot(highexp)
highexp <- highexp$result[,c(1:12)]
filter(highexp, source == "GO:BP")

###NCBI data
##will study divergence between different corona virus groups
library(Biostrings)
library(philentropy)
library(xml2)
library(XML)
library(tidyverse)

fastafiles <- list.files("~/ncbi-genomes-2022-07-24/", full.names = T)
fastafiles <- fastafiles[grep("_genomic.fna.gz", fastafiles)]
fastafiles <- tibble(fnames = fastafiles, asmid = str_match(fastafiles, pattern = "\\S+\\/(\\S+_\\d+\\.\\d)_")[,2])

metadata <- read_xml("~/ncbi-genomes-2022-07-24/assembly_result.xml")
metadata <- xmlParse(metadata)
metadata <- bind_cols(xmlToDataFrame(nodes = getNodeSet(metadata, "//AssemblyAccession")), xmlToDataFrame(nodes = getNodeSet(metadata, "//SpeciesName")))
colnames(metadata) <- c("asmid", "organism")
metadata <- left_join(metadata, fastafiles)

gnames <- list.files("~/ncbi-genomes-2022-07-24/")
gnames <- gnames[grep("_genomic.fna.gz", gnames)]
gnames <- gsub("_genomic.fna.gz", "", gnames)
gnames <- gsub("GCF_009858895.2_ASM985889v3", "SARS-CoV-2", gnames)

getkmerfreq <- function(x) {
      genomeseq <- readDNAStringSet(x)
      kfreq <- oligonucleotideFrequency(genomeseq,9) ##9 mers obtained
      return(kfreq)
}
metadata %>% filter(!is.na(fnames))

kmerfreq <- sapply(pull(metadata, fnames), getkmerfreq)
kmat <- t(as.matrix(kmerfreq))
covjsd <- JSD(kmat, test.na = TRUE, unit = "log2", est.prob = "empirical")
row.names(covjsd) <- pull(metadata, organism)
colnames(covjsd) <- pull(metadata, organism)
pheatmap(covjsd)



