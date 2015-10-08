convert_bayenv <- function(structure, popmap) {
	
	x <- read.table(structure, sep = "\t")
	pops <- read.table(popmap, sep = "\t")
	inds <- as.character(x[2:nrow(x), 1])
	n.pops <- length(unique(pops[,2]))
	
	# remove the population identification column or not
	if (is.na(x[1, 2]) == TRUE) {
		x <- cbind(x[,1], x[,3:ncol(x)])
	}
	
	snp.names <- as.vector(as.character(x[1,2:ncol(x)]))
	
	
	#######################################################################
	#
	# determine that a SNP is in each population, remove SNPs that are not
	for (a in 1:ncol(x)) {
		if(a == 1) {
			x2 <- as.character(x[,1])
		} else {
			x.test <- x[2:nrow(x),a]
			skip <- 0
			for(b in 1:n.pops) {
				test <- match(inds,as.character(pops[pops[,2]==b,1]))
				test[is.na(test)] <- 0
				test <- x.test[test > 0]
				test <- test[test > 0]
				if(length(test) == 0) {
					skip <- skip + 1
				}
			}
			if(skip == 0) {
				x.test <- c(snp.names[a-1], x.test)
				x2 <- cbind(x2,x.test)
			}
		}
	}
	# remove items 
	x <- x2
	rm(a,b,skip,test,x.test,x2,popmap,structure)
	snp.names <- x[1,2:ncol(x)]
	
	
	######################################################################
	#
	# determine if locus is biallelic, remove if not
	for(a in 1:ncol(x)) {
		if(a == 1) {
			x2 <- x[,1]
		} else {
			x.test <- x[2:nrow(x),a]
			x.test <- x.test[x.test > 0]
			if(length(unique(x.test))==2) {
				x.test <- c(snp.names[a-1], x[2:nrow(x),a])
				x2 <- cbind(x2, x.test)
			}
		}
	}
	# remove items
	x <- x2
	rm(x2, a, x.test)
	snp.names <- x[1,2:ncol(x)]
	
	
	######################################################################
	#
	# create allele frequency and sample size matrices
	for(a in 1:ncol(x)) {
		if(a == 1) {
		} else {
			x.rep <- as.vector(x[2:nrow(x),a])
			allele.a <- x.rep[x.rep > 0]
			allele.a <- allele.a[1]
			alleles <- c()
			samples <- c()
			for(b in 1:n.pops) {
				test <- match(inds,as.character(pops[pops[,2]==b,1]))
				test[is.na(test)] <- 0
				test <- x.rep[test > 0]
				alleles <- c(alleles, length(na.omit(match(test, allele.a))))
				samples <- c(samples, length(test[test>0]))
			}
			if(a == 2) {
				allele.matrix <- alleles
				sample.size.matrix <- samples
			} else {
				allele.matrix <- cbind(allele.matrix, alleles)
				sample.size.matrix <- cbind(sample.size.matrix, samples)
			}
		}
	}
	# remove items and give row/column names
	colnames(allele.matrix) <- snp.names
	rownames(allele.matrix) <- 1:n.pops
	colnames(sample.size.matrix) <- snp.names
	rownames(sample.size.matrix) <- 1:n.pops
	allele.b.matrix <- sample.size.matrix - allele.matrix
	rm(sampl.size.matrix)
	
	
	######################################################################
	#
	# Create a combined file that has allele 1 on line 1, and allele 2 on 
	# line 2. Also write files to a directory of each locus
	
	dir.create("SNP_files")
	for (a in 1:ncol(allele.matrix)) {
		temp.a <- allele.matrix[,a]
		temp.b <- allele.b.matrix[,a]
		if (a == 1) {
			output.matrix <- rbind(temp.a, temp.b)
		} else {
			output.matrix <- rbind(output.matrix, temp.a, temp.b)
		}
		snp.num <- paste("SNP_files/", a, sep="")
		snp.output <- write.table(rbind(temp.a, temp.b), file = snp.num, 
			quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
	}	
	
	output.matrix <- cbind(output.matrix, rep("", nrow(output.matrix)))
	write.table(output.matrix, file = "output.matrix", quote=FALSE,
		col.names=FALSE, row.names=FALSE, sep="\t")
	snp.names <- cbind(snp.names, 1:length(snp.names))
	write.table(snp.names, file="SNP.identifiers", quote=FALSE, 
		col.names=FALSE, row.names=FALSE, sep="\t")
	
}

