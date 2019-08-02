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
