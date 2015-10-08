convert_bayescan <- function(structure, popmap) {
	# Read input files
	file.struc <- read.table(structure, sep="\t")
	file.popma <- read.table(popmap, sep="\t")
	
	# Remove populations columns from structure file if it exists
	if (is.na(file.struc[1, 2]) == TRUE) {
		file.struc <- cbind(file.struc[,1], file.struc[,3:ncol(file.struc)])
	}
	
	# Create list of loci from structure file
	loci <- as.list(as.numeric(file.struc[1,2:ncol(file.struc)]))

	# Loop for each population to check for missing loci (due to method of 
	# Stacks output)
	pops <- sort(unique(file.popma[,2]))
	test1 <- c()
	for(a in 1:length(pops)) {
		# Extract sample names from popmap file
		names.inds <- file.popma[file.popma[,2]==pops[a],1]
		
		# Extract structure data for population
		pop.struc <- sort(c(match(names.inds, file.struc[,1]),
			match(names.inds, file.struc[,1]) + 1))
		pop.struc <- file.struc[pop.struc,2:ncol(file.struc)]
		
		# Check for missing data for a population (Stacks misrepresentation)
		for(c in 1:ncol(pop.struc)) {
			test <- pop.struc[,c]
			test[test==0] <- NA
			if(length(na.omit(test)) == 0) { test1 <- c(test1, c)}
		}
	}
	
	# check if biallelic, add bad loci to test1 for blacklist
	file.struc2 <- file.struc[2:nrow(file.struc), 2:ncol(file.struc)]
	for(e in 1:ncol(file.struc2)) {
		if(!is.na(match(0,file.struc2[,e]))) {
			num.alleles <- length(as.vector(table(file.struc2[,e])))
			if(num.alleles != 3) {
				test1 <- c(test1,e)
			}
		} else {
			num.alleles <- length(as.vector(table(file.struc2[,e])))
			if(num.alleles != 2) {
				test1 <- c(test1,e)
			}
		}
	}
	
	# test1 is now a blacklist
	# create a file of all loci (test2) to remove the blacklist
	test2 <- c()
	for(a in 1:ncol(pop.struc)) {
		test2 <- c(test2, a)
	}
	if(length(test1)!=0) {
		test1 <- sort(unique(test1))
		test3 <- match(test2, test1)
		test <- test2[is.na(test3)]
	} else {
		test <- test2
	}
	
	# Write names of snps to file
	snp.names.file <- as.vector(loci[test])
	length.loci <- seq(1:length(loci))
	loci.numbers <- length.loci[test]
	write(c("locus", "snp_column_from_structure"),
			file="bayescan.SNP.names.txt", ncolumns=2)
	for(a in 1:length(snp.names.file)) {
		write(c(snp.names.file[[a]], loci.numbers[a]),
			file="bayescan.SNP.names.txt", append=TRUE)
	}
	
	# Write header of bayesan file
	pops <- sort(unique(file.popma[,2]))
	header <- list()
	header[[1]] <- paste("[loci]=", length(test), sep="")
	header[[2]] <- ""
	header[[3]] <- paste("[populations]=", length(pops), sep="")
	header[[4]] <- ""
	write(header[[1]], file="bayescan.infile.txt", sep=" ", ncolumns=10)
	write(header[[2]], file="bayescan.infile.txt", sep=" ", ncolumns=10,
		append=TRUE)
	write(header[[3]], file="bayescan.infile.txt", sep=" ", ncolumns=10,
		append=TRUE)
	write(header[[4]], file="bayescan.infile.txt", sep=" ", ncolumns=10,
		append=TRUE)
	
	# remove blacklist loci from file.struc2
	file.struc2 <- file.struc2[,test]
	
	# obtain the two SNPs for each locus
	allele.a <- c()
	allele.b <- c()
	for(f in 1:ncol(file.struc2)) {
		if(!is.na(match(0,file.struc2[,f]))) {
			allele <- sort(unique(file.struc2[,f]))
			allele.b <- c(allele.b, allele[3])
			allele.a <- c(allele.a, allele[2])
		} else {
			allele <- sort(unique(file.struc2[,f]))
			allele.b <- c(allele.b, allele[2])
			allele.a <- c(allele.a, allele[1])
		}
	}
	
	# Loop for each population
	a <- 1
	for (a in 1:length(pops)) {
		header[[1]] <- paste("[pop]=", a, sep="")
		write(header[[1]], file="bayescan.infile.txt", sep=" ", ncolumns=10,
			append=TRUE)
		
		# Extract structure data for population
		names.inds <- file.popma[file.popma[,2]==pops[a],1]
		pop.struc <- sort(c(match(names.inds, file.struc[,1]),
			match(names.inds, file.struc[,1]) - 1))
		pop.struc <- file.struc2[pop.struc,1:ncol(file.struc2)]
		pop.struc[pop.struc==0] <- NA
		# loop for each snp
		b <- 1
		for(b in 1:ncol(file.struc2)) {
			count.rep <- length(na.omit(pop.struc[,b]))
			unique.rep <- unique(na.omit(pop.struc[,b]))
			if(length(unique.rep)==2) {
				allele.a.rep <- as.vector(table(na.omit(pop.struc[,b]))[1])
				allele.b.rep <- as.vector(table(na.omit(pop.struc[,b]))[2])
			} else {
				homo.rep <- unique(na.omit(pop.struc[,b]))
				length.homo.rep <- length(na.omit(pop.struc[,b]))
				if(is.na(match(homo.rep, allele.a[b]))) {
					allele.a.rep <- 0
					allele.b.rep <- length.homo.rep
				} else {
					allele.a.rep <- length.homo.rep
					allele.b.rep <- 0
				}
			}
			line.rep <- c(b, count.rep, 2, allele.a.rep, allele.b.rep)
			write(line.rep, file="bayescan.infile.txt", sep=" ", 
				ncolumns=10, append=TRUE)
		}
		if(a != length(pops)) {
			write(header[[2]], file="bayescan.infile.txt", sep=" ", ncolumns=10,
				append=TRUE)	
		}
	}
}
