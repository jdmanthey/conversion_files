# Converts a structure file to SNAPP format using the first SNP from 
# each locus. Output is named SNAPP.nex and the mutation rates (U and
# V) are printed.

# USAGE:
# Place structure file and this script in a directory, and set the R
# current working directory to be the same. Then type the following 
# commands, where "file" is the name of your structure file:
# source("convert_snapp.r")
# convert_snapp("file") 


convert_snapp <- function(file) {
	file <- file
	
	x <- read.table(file, sep = "\t")
	# remove the population identification column or not
	if (is.na(x[1, 2]) == TRUE) {
		x <- cbind(x[,1], x[,3:ncol(x)])
	}
	
	# take first SNP from each locus
	list.loci <- unique(as.numeric(x[1,2:ncol(x)]))
	matches.loci <- c(1, match(list.loci, x[1,]))
	x <- x[,matches.loci]
		
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
			if(x2[indiv.a,2] == 0){indivi <- "-"}
		
		for(d in 3:ncol(x2)){
			
			if(x2[indiv.a,d] == x2[indiv.b, d] & x2[indiv.a,d] != 0 & x2[indiv.a,d] == allele.a[d]){indiv <- 0}
			if(x2[indiv.a,d] == x2[indiv.b, d] & x2[indiv.a,d] != 0 & x2[indiv.a,d] == allele.b[d]){indiv <- 2}
			if(x2[indiv.a,d] != x2[indiv.b, d]){indiv <- 1}
			if(x2[indiv.a,d] == 0){indiv <- "-"}
			indivi <- paste(indivi, indiv, sep="")
			
		}
		
		file <- c(file, indivi)
		
	}
	output <- cbind(individual, file)
	header <- list()
	header[[1]] <- "#NEXUS"
	header[[2]] <- ""
	header[[3]] <- "Begin data;"
	header[[4]] <- paste("	Dimensions ntax=", length(individual), " nchar=", nchar(file[1]), ";", sep="")
	header[[5]] <- paste("	Format datatype=standard symbols=", '"', "012", '"',
		" missing=-;", sep="")
	header[[6]] <- "	Matrix"
	header[[7]] <- "	;"
	header[[8]] <- "End;"
	write(header[[1]], file="SNAPP.nex",ncolumns=1)
	write(header[[2]], file="SNAPP.nex",ncolumns=1, append=TRUE)
	write(header[[3]], file="SNAPP.nex",ncolumns=1, append=TRUE)
	write(header[[4]], file="SNAPP.nex",ncolumns=1, append=TRUE)
	write(header[[5]], file="SNAPP.nex",ncolumns=1, append=TRUE)
	write(header[[6]], file="SNAPP.nex",ncolumns=1, append=TRUE)
	for(a in 1:nrow(output)) {
		write(output[a,], file="SNAPP.nex", ncolumns=2, sep="\t",append=TRUE)
	}
	write(header[[7]], file="SNAPP.nex",ncolumns=1, append=TRUE)
	write(header[[8]], file="SNAPP.nex",append=TRUE)
	
	mutate("SNAPP.nex")	
		
		
		
		
}








# input is snapp file
mutate <- function(x) {
	#obtain number of rows of data
	y <- read.table(x, skip=3, nrows=1)
	y <- as.character(y[1,2])
	y <- strsplit(y, "=")
	y <- y[[1]]
	y <- y[2]
	y <- strsplit(y, ";")
	y <- y[[1]]
	y <- y[1]
	y <- as.numeric(y)
	
	x <- read.table(x, skip=6, nrows=y)
	x <- as.matrix(x[,2])
	x <- as.vector(x[,1])
	for(a in 1:length(x)) {
		if(a == 1) { y <- x[a] }
		if(a != 1) { y <- paste(y, x[a], sep="")}	
		}
	y <- strsplit(y, "")
	y <- y[[1]]
	
	x1 <- as.vector(table(y)[2])
	x2 <- as.vector(table(y)[3])
	x3 <- as.vector(table(y)[4])
	u <- x1 + 0.5 * x2
	v <- x3 + 0.5 * x2
	w <- u + v
	u <- u / w
	v <- v / w
	u2 <- 1 / (2 * u)
	v2 <- 1 / (2 * v)
	print(paste("U = ", u2, sep=""))
	print(paste("V = ", v2, sep=""))
}
