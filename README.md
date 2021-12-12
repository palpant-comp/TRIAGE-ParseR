# pca_gmm
TRIAGE gene clustering using Principal component analysis - Gaussian Mixture model (PCA-GMM)

> python3 pca_gmm.py -i YOUR_INPUT_DS_TABLE 

options
-i, input discordance table
-r, pre-calculated H3K27me3 principal components (Default = ./data/pca_x)
-p, Number of PCs to use (Default = 67)
-g, Number of top genes to use (Default = 100)
-o, Output directory (Default = ./results)
-e, Whether to perform GO enrichment analysis (1: Yes, default, 0: No)
