# USAGE:
# Place structure file and this script in a directory, and set the R
# current working directory to be the same. Then type the following 
# commands, where "file" is the name of your structure file:
# source("convert_new_hybrids.r")
# convert_snapp("file") 


convert_nh <- function(file) {

	x <- read.table(file, sep = "\t")
	# remove the population identification column or not
	if (is.na(x[1, 2]) == TRUE) {
		x <- cbind(x[,1], x[,3:ncol(x)])
	}
	
	x <- x[2:nrow(x),2:ncol(x)]
	
	n.inds <- nrow(x) / 2
	loci <- c()
	for(c in 1:ncol(x)) {
		x.rep <- x[,c]
		x.rep <- x.rep[x.rep>0]
		x.rep <- sort(as.vector(table(x.rep)))
		x.min <- sum(x.rep) * 0.3
		if(x.rep[1]>x.min) {loci <- c(loci, c)}
	}
	
	if(length(loci)>500) {
		loci <- loci[1:500]
	}
	n.loci <- length(loci)
	
	write(paste("NumIndivs",n.inds),file="new_hybrids_input.txt",ncolumns=2)
	write(paste("NumLoci",n.loci),file="new_hybrids_input.txt",ncolumns=2,append=T)
	write(paste("Digits 1"),file="new_hybrids_input.txt",ncolumns=1,append=T)
	write(paste("Format Lumped"),file="new_hybrids_input.txt",ncolumns=2,append=T)
	write(paste(""),file="new_hybrids_input.txt",ncolumns=1,append=T)
	
	for(a in 1:n.inds) {
		rep2 <- a * 2
		rep1 <- rep2 - 1
		
		for(b in loci) {
			snp.rep1 <- x[rep1,b]
			snp.rep2 <- x[rep2,b]
			snp.rep <- paste(snp.rep1,snp.rep2,sep="")
			if(snp.rep == "00") {
				snp.rep <- 0
			}
			if(b == loci[1]) {
				snps <- paste(a, snp.rep)
			} else {
				snps <- paste(snps, snp.rep)
			}
		}
		write(snps,file="new_hybrids_input.txt",ncolumns=1,append=T)
	}
	
}
