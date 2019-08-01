library(rentrez)
library(XML)
library(xml2)
library(RCurl)
library(ggplot2)

unescape_html <- function(str){
  xml2::xml_text(xml2::read_html(paste0("<x>", str, "</x>")))
}

getAssemblyInfo <- function(x) {
	return(entrez_search("assembly", paste(x, "[Organism]", sep = "") , retmax=2000))
}

getSummaryInfo <- function(x, group){
	suminfo <- entrez_summary(db="assembly", id=x)
	a <- xmlToList(xmlParse(paste("<doc>",unescape_html(suminfo$meta),"</doc>"), asText = T))
	b <- a[["Stats"]]
	df <- as.data.frame(t(matrix(unlist(b), nrow = 3)))
	colnames(df) <- c("value", "attribute", "type")
	df$value <- as.numeric(as.character(df$value))
	df <- cbind(df, ncbitaxid = suminfo$taxid)
	df <- cbind(df, organismname = suminfo$organism)
	df <- cbind(df, group = group)
	return(df)
}

getGroupInfo <- function(group, maxvalue){
	searchresults <- getAssemblyInfo(group)
	maxvalue <- ifelse(length(searchresults$ids) > maxvalue, maxvalue, length(searchresults$ids))
	return(as.data.frame(do.call(rbind, lapply(sample(searchresults$ids, maxvalue), getSummaryInfo, group))))
}

listofgroups <- as.list(c("mammalia", "aves", "amphibia", "teleostei"))

assemblySummary <- do.call(rbind,lapply(listofgroups, getGroupInfo, 10))

ggplot(assemblySummary[assemblySummary$attribute=="ungapped_length",], aes(x = ncbitaxid, y=value, fill=group)) + geom_bar(position = "dodge", stat = "identity")
