#' Convert a vector of raw LGC data to FAPS readable format.
vectorise_locus <- function(locus, sep=":"){
  loclist <- strsplit(as.character(locus), split = sep)
  #alleleA <- sort(names(table(unlist(loclist))))[1] # assign minor allele
  #alleleB <- 
  alleles <- sort(names(table(unlist(loclist)))) # assign major allele
  alleles <- alleles[alleles %in% c("A", "C", "G", "T")]
  
  loclist <- sapply(loclist, function(x) paste(x, collapse = ""))
  loclist[loclist == paste(c(alleles[1], alleles[1]), collapse="")] <- 0
  loclist[loclist == paste(c(alleles[1], alleles[2]), collapse="")] <- 1
  loclist[loclist == paste(c(alleles[2], alleles[1]), collapse="")] <- 1
  loclist[loclist == paste(c(alleles[2], alleles[2]), collapse="")] <- 2
  
  as.integer(loclist)
}

#' Convert a dataframe from LGC to FAPS format.
#' 
#' Convert a dataframe from LGC to FAPS format.
#' 
#' Data should have a column for genotype IDs, followed by columns of
#' string data. These data should be "G:G", "C:A" etc
lgc2faps <- function(lgc_raw, sep=":"){
  faps <- apply(lgc_raw, 2, vectorise_locus, sep)
  faps <- data.frame(faps)
  faps[,1] <- lgc_raw[,1]
  faps
}

