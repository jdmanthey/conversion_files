# add structure file and popmap

convert_treemix <- function(file, popmap) {
	
	x <- read.table(file, sep = "\t")
	# remove the population identification column or not
	if (is.na(x[1, 2]) == TRUE) {
		x <- cbind(x[,1], x[,3:ncol(x)])
	}
	
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
	
	###determine possible alleles of each SNP
	allele.a <- c()
	allele.b <- c()
	for(b in 2:ncol(x2)) {
		
		alleles <- unique(x2[, b])
		alleles <- alleles[ alleles != 0 ]
		allele.a <- c(allele.a, alleles[1])
		allele.b <- c(allele.b, alleles[2])
		
	}
	
	# match individuals to populations
	pops <- read.table(popmap,sep="\t")
	inds <- as.character(pops[,1])
	species <- as.character(pops[,2])
	pops <- cbind(inds, species)
	unique.species <- unique(species)
	
	
	# population allele frequencies loop
	allele.freqs <- list()
	for(a in 1:length(unique.species)) {
		z.a <- c()
		z.b <- c()
		inds.rep <- as.vector(pops[pops[,2]==unique.species[a],1])
		inds.rep <- match(inds.rep,x2[,1])
		inds.rep <- c(inds.rep, inds.rep+1)
		inds.rep <- x2[inds.rep,2:ncol(x2)]
		for(b in 1:ncol(inds.rep)) {
			snp.rep <- inds.rep[,b]
			snp.rep <- snp.rep[snp.rep != 0]
			z.a <- c(z.a, length(snp.rep[snp.rep == allele.a[b]]))
			z.b <- c(z.b, length(snp.rep[snp.rep == allele.b[b]]))
		}
		allele.freqs[[a]] <- paste(z.a, ",", z.b, sep="")
	}
	
	# put all lists together and add column names
	for(a in 1:length(allele.freqs)) {
		if(a == 1) {
			final <- allele.freqs[[a]]
		} else {
			final <- cbind(final, allele.freqs[[a]])
		}
	}
	colnames(final) <- unique.species
					
	write.table(final, file = "treemix", sep = " ", row.names = FALSE, quote=FALSE)
}
