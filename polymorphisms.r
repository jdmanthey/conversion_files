# Function to identify private, shared, and fixed polymorphisms
# for all populations defined in a popmap and snp data in a structure
# formatted file.
# The structure file is a tab delimited matrix with all SNP data in columns
# and two rows for each diploid individual, and including a column
# for both the individual names and their a priori locality (not used in function).
# Missing data in the structure file should be 0 or a negative integer. Each base
# pair should be a unique positive integer.
# The popmap is a tab delimited file with the first column
# indicating the sample names, and the second column a unique 
# number for each population.
# 
# With this script in current working directory: source("polymorphisms.r")
polymorphisms <- function(file, popmap) {
	x <- read.table(file,sep="\t")
	inds <- read.table(popmap, sep="\t")
	pops <- unique(as.numeric(inds[,2]))
	
	# loop to calculate for each population
	for(a in pops) {
		
		# define population 1 and all others
		inds1 <- as.character(inds[inds[,2]==pops[a],1])
		inds2 <- as.character(inds[inds[,2]!=pops[a],1])
		
		# extract all the lines from the structure file
		# pertaining to each group
		names <- as.character(x[2:nrow(x),1])
		inds1 <- match(inds1,names)
		inds1 <- sort(c(inds1, (inds1 + 1)))
		inds2 <- match(inds2,names)
		inds2 <- sort(c(inds2, (inds2 + 1)))
		x.new <- x[2:nrow(x), 3:ncol(x)]
		inds1 <- x.new[inds1,]
		inds2 <- x.new[inds2,]
		
		# reset all counters
		shared <- 0
		fixed <- 0
		private <- 0
		mono <- 0
		nbi <- 0
		
		# loop for each snp to identify if fixed, shared, etc.
		for(b in 1:ncol(inds1)) {
			i1.rep <- inds1[,b]
			i2.rep <- inds2[,b]
			i1.rep <- i1.rep[i1.rep>0]
			i2.rep <- i2.rep[i2.rep>0]
			if(length(unique(c(i1.rep,i2.rep)))==2) {
				if(length(unique(i1.rep))==1 & length(unique(i2.rep))==1) {
					fixed <- fixed + 1
				}
				if(length(unique(i1.rep))>1 & length(unique(i2.rep))==1) {
					private <- private + 1
				}
				if(length(unique(i1.rep))>1 & length(unique(i2.rep))>1) {
					shared <- shared + 1
				}
			}
			if(length(unique(c(i1.rep,i2.rep)))==1) {
				mono <- mono + 1
			}
			if(length(unique(c(i1.rep,i2.rep)))>2) {
				nbi <- nbi + 1
			}
		}
		
		# print output to screen for each population
		writeLines(" ")
		writeLines(paste("Population", pops[a],sep=" "))
		writeLines(paste("Fixed Differences:", fixed))
		writeLines(paste("Shared Polymorphisms:", shared))
		writeLines(paste("Private Polymorphisms:", private))
		writeLines(paste("Monomorphic:", mono))
		writeLines(paste("More than two alleles:", nbi))
		
	}
}
