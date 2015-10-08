convert_lfmm <- function(file) {
	file <- file
	
	x <- read.table(file, sep = "\t")
	# remove the population identification column or not
	if (is.na(x[1, 2]) == TRUE) {
		x <- cbind(x[,1], x[,3:ncol(x)])
	}
	
	
	snp.names <- as.character(x[1,2:ncol(x)])	
	x <- x[2:nrow(x),]
		
	bi_allelic <- c()
	for(a in 2:ncol(x)) {			#determine which SNPs are bi-allelic
		
		bi <- x[,a]
		bi <- unique(bi)
		bi <- bi[ bi != 0 ]
		bi <- length(bi)
		bi_allelic <- c(bi_allelic, bi)
		
	}
	
	bi_allelic <- bi_allelic == 2
	
	x2 <- x[,bi_allelic]
	snp.names <- snp.names[bi_allelic]
	snp.names <- as.numeric(snp.names)
	snp.names <- cbind(snp.names,seq(1:length(snp.names)))
	
	write.table(snp.names,file="xxx.lfmm.snp.names.txt",quote=FALSE,col.names=FALSE,		row.names=FALSE, sep="\t")
	###determine possible alleles of each SNP
	allele.a <- c()
	allele.b <- c()
	for(b in 1:ncol(x2)) {
		
		alleles <- unique(x2[, b])
		alleles <- alleles[ alleles != 0 ]
		allele.a <- c(allele.a, alleles[1])
		allele.b <- c(allele.b, alleles[2])
		
	}
	c <- c()
	d <- c()
	n.inds <- nrow(x2) / 2
	
	file <- c()
	individual <- c()
	for(c in 1:n.inds) {
		
		indiv.a <- c*2 - 1
		indiv.b <- c*2
		individual <- c(individual, as.character(x2[indiv.a, 1]))
		indiv2 <- c()
		
			if(x2[indiv.a,2] == x2[indiv.b, 2] & x2[indiv.a,2] != 0 & x2[indiv.a,2] == allele.a[2]){indivi <- 0}
			if(x2[indiv.a,2] == x2[indiv.b, 2] & x2[indiv.a,2] != 0 & x2[indiv.a,2] == allele.b[2]){indivi <- 2}
			if(x2[indiv.a,2] != x2[indiv.b, 2]){indivi <- 1}
			if(x2[indiv.a,2] == 0){indivi <- 9}
		
		for(d in 3:ncol(x2)){
			
			if(x2[indiv.a,d] == x2[indiv.b, d] & x2[indiv.a,d] != 0 & x2[indiv.a,d] == allele.a[d]){indiv <- 0}
			if(x2[indiv.a,d] == x2[indiv.b, d] & x2[indiv.a,d] != 0 & x2[indiv.a,d] == allele.b[d]){indiv <- 2}
			if(x2[indiv.a,d] != x2[indiv.b, d]){indiv <- 1}
			if(x2[indiv.a,d] == 0){indiv <- 9}
			indivi <- paste(indivi, indiv, sep=" ")
			
		}
		
		file <- c(file, indivi)
		
	}
	
	for(a in 1:length(file)) {
		if(a==1) {
			write(file[a], file="xxx.lfmm",ncolumns=1)
		} else {
			write(file[a], file="xxx.lfmm",ncolumns=1, append=TRUE)
		}
	}
			
}




