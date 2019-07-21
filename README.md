# Single-Cell-Deconvolution-and-Survival-Analysis
Glioblastoma multiforme (GBM) is an aggressive brain cancer with a median survival time of 15-16 months.  Traditional treatment for the average 14k cases of GBM per year includes surgery, radiation, and chemotherapy focused on extending survival time.  Survival time may be associated with increased immune cell presence in these solid tumors (requires citation).  Deconvolution of bulk RNA data from GBM tumors by immune single cell RNA profiles to determine immune cell content may provide a glimpse into estimating patient survival.

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
# Output:
Cell type composition in the n tissue samples
