# Single Cell Deconvolution and Survival Analysis
Glioblastoma multiforme (GBM) is an aggressive brain cancer with a median survival time of 15-16 months.  Traditional treatment for the average 14k cases of GBM per year includes surgery, radiation, and chemotherapy focused on extending survival time.  Survival time may be associated with increased immune cell presence in these solid tumors (requires citation).  Deconvolution of bulk RNA data from GBM tumors by immune single cell RNA profiles to determine immune cell content may provide a glimpse into estimating patient survival.

## related links:
- https://github.com/xzhoulab/DECComparison/tree/master/data_download
- https://github.com/bnwolford/BDSI

# Background:

- Glioblastoma multiforme (GBM) is the most aggressive type of brain cancer with a median survival time of 15-16 months and an average of 14k cases in the U.S. each year. 
- Traditional treatments (surgery, radiation, and chemotherapy) and on-going research concentrate on understanding and estimating patient survival. 
- Immune cell presence in tumors has been found to be associated with patient survival, which can be estimated with the application of RNA deconvolution methods.  


# RNA sequencing methods:
## Single Cell RNA Sequencing: 
- Provides the gene expression profile of individual cells
- Reveals expression variability of different cell types and subpopulations
## Bulk RNA Sequencing:
- Provides the average gene expression level of all cells in a tissue 
- Masks the cell type heterogeneity in the tissue
- Relatively cheap 

# General outline for deconvolution:
## Input:
- Gene signatures (single cell): Expression level of m genes in k cells
- Mixture (bulk):  Expression level of m genes in n tissue samples
## Output:
Cell type composition in the n tissue samples

# References(MLA)
Aran, Dvir, Zicheng Hu, and Atul J. Butte. "xCell: digitally portraying the tissue cellular heterogeneity landscape." Genome biology 18.1 (2017): 220.

Azizi, Elham, et al. "Single-cell map of diverse immune phenotypes in the breast tumor microenvironment." Cell 174.5 (2018): 1293-1308.

Frishberg, Amit, et al. "Cell composition analysis of bulk genomics using single-cell data." Nature methods 16.4 (2019): 327.

Gong, Ting, and Joseph D. Szustakowski. "DeconRNASeq: a statistical framework for deconvolution of heterogeneous tissue samples based on mRNA-Seq data." Bioinformatics29.8 (2013): 1083-1085.

Li, Bo, et al. "Comprehensive analyses of tumor immunity: implications for cancer immunotherapy." Genome biology 17.1 (2016): 174.

Newman, Aaron M., et al. "Robust enumeration of cell subsets from tissue expression profiles." Nature methods 12.5 (2015): 453.

Racle, Julien, et al. "Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data." Elife 6 (2017): e26476.

Wang, Xuran, et al. "Bulk tissue cell type deconvolution with multi-subject single-cell expression reference." Nature communications 10.1 (2019): 380.
