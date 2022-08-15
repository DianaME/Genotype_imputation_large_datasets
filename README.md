# Genotype_imputation_large_datasets
This imputation was performed based on conditional probabilities of marker genotypes, given the estimated recoombination rates with the nearest nonmissing flanking markers. 
The information for the imputation was obtained from the next two papers
Wu, R., C. Ma, and G. Casella. 2007. Statistical genetics of quantitative traits: Linkage, maps and QTL. Springer, New York.
Jacobson A, Lian L, Zhong S, Bernardo R. Marker Imputation Before Genomewide Selection in Biparental Maize Populations. Plant Genome. 2015 Jul;8(2):eplantgenome2014.10.0078. doi: 10.3835/plantgenome2014.10.0078. PMID: 33228308.

In this repository there are three files a bash script for submision of the file to the server; an R script file to perform the genotype imputation by families, which allows to use 40 node at the time and speed up the computation time (server_imputation), and ohter R script (pre_files.R) to estimate the recombination rates among markers which is needed for the imputation
