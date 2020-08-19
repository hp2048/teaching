##include central libraries in your paths
.libPaths(c(.libPaths(), "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library")) 

##load genome sequence library
library(biovizBase)
library(BSgeome)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(GenomicRanges)

hchrinfo <- data.frame(chr = Hsapiens@seqinfo@seqnames, start = 1, end = Hsapiens@seqinfo@seqlengths)
myslidingwindows <- GenomicRanges::slidingWindows(GRanges(hchrinfo), width = 250000, step = 250000)
myslidingwindows



library(Biostrings)
library(ggplot2)

getGCfraction <- function(x){
	genome <- readDNAStringSet(x)
	af <- alphabetFrequency(genome)
	bc <- colSums(af)
	return(sum(bc[c("C","G")])/sum(genome@ranges@width))
}

fastafiles <- list.files("/mnt/data/", pattern = ".fa.gz", full.names=T)
gc <- lapply(fastafiles, getGCfraction)
gc <- data.frame(species=sapply(strsplit(fastafiles, "/"), "[", 5), gc)
ggplot(gc, aes(x=species, y=gc)) + geom_bar(stat="identity")
