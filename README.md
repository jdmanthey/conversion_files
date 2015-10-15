# conversion_files

Included are many R scripts for conversions between file types for different genetic and ecological analyses. 

.

Generally, these files can be sourced into the current iteration of R followed by their usage. 

Conversion files for use with a STRUCTURE file and/or population map (for use with SNP data):

BayEnv2: convert_bayenv.r

BayeScan: convert_bayescan.r

BEDASSLE: convert_bedassle.r

Latent Factor Mixed Models: convert_lfmm.r

NewHybrids: convert_new_hybrids.r

SNAPP: convert_snapp.r

TreeMix: convert_treemix.r


Here, the STRUCTURE file is assumed to have diploid organisms, two rows for each SNP, and have a header column containing the names of each locus. The population map has the name of each individual and their associated population/locality (numeric) in a two column tab-delimited file.



Script to calculate shared, fixed, and private polymorphisms (annotated and explained within script): polymorphisms.r
