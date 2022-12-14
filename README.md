# pca_gmm
TRIAGE gene clustering by Principal component analysis - Gaussian Mixture model (PCA-GMM)

> python3 TRIAGE_parser.py -i YOUR_INPUT_TEXT_FILE -o OUTPUT_DIRECTORY 

options

-i, input file

-r, pre-calculated H3K27me3 principal components (Default = ./data/pca_x)

-p, Number of PCs to use (Default = 10)

-g, Number of top genes to use (Default = 100, applicable only if the input is a table)

-o, Output directory (Default = ./results)

-e, Whether to perform GO enrichment analysis (1: Yes, default, 0: No)

-a, input type (option: table or list, default = list)

-v, verbose level (option: 1 or 0, default = 1)

-w, max. number of clusters allowed (default=10)

-j, gene order, gene sort direction by value (option: ascending or descending, default = descending)

For example of the input type, see included example text files. 

If you have any questions, please contact Woo Jun (Chris) Shim, w.shim@uq.edu.au 

For the citation, please check our paper in BioRxiv https://doi.org/10.1101/2022.10.12.512003 
