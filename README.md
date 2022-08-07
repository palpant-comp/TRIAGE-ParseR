# pca_gmm
TRIAGE gene clustering by Principal component analysis - Gaussian Mixture model (PCA-GMM)

> python3 pca_gmm_new.py -i YOUR_INPUT_TEXT_FILE -o OUTPUT_DIRECTORY 

options

-i, input file

-r, pre-calculated H3K27me3 principal components (Default = ./data/pca_x)

-p, Number of PCs to use (Default = 67)

-g, Number of top genes to use (Default = 100)

-o, Output directory (Default = ./results)

-e, Whether to perform GO enrichment analysis (1: Yes, default, 0: No)

-a, input type (option: table or list, default = list)

For example of the input type, see included example text files. 
