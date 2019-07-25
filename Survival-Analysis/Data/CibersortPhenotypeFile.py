import pandas as pd
import csv
import math

#####################################################################
# Making of cibersort phenotype class matrix
#####################################################################
'''
Read cluster name, genes, and single cell survival data
'''
clusterNamePath = '/home/stephen/Desktop/clusterName.xlsx' # 
clusterNameDF = pd.read_excel(clusterNamePath, header=None)

clusterGenePath = '/home/stephen/Desktop/clusterGene.xlsx' # 
clusterGeneDF = pd.read_excel(clusterGenePath, header=None)

sc_counts_path = '/home/stephen/survival_sc_counts.txt'# 
sc_counts_matrix = pd.read_csv(sc_counts_path, sep=" ")

'''
Rename cell types
'''
cellTypes = clusterNameDF[1].tolist() # list of cell types

cellNewTypesDict = {
    'B:':'B',
    'MACROPHAGE:':'MACROPHAGE',
    'NEUTROPHIL:':'NEUTROPHIL',
    'MONOCYTE:precursor':'MONOCYTE',
    'MONOCYTE:':'MONOCYTE',

    'MAST:':'MAST',
    
    'NK:CD56-16+3-':'NK',
    'NK:CD56+16+3-':'NK',
    'NK:CD56+16+3+NKT':'NK',
    'NKT':'NK',

    'T:Reg':'T:Reg',
    
    'T:CD8+EM':'T:CD8',
    'T:CD8+CM':'T:CD8',
    'T:CD8+NAIVE':'T:CD8',
    
    'T:CD4+EM':'T:CD4',    
    'T:CD4+CM':'T:CD4',
    'T:CD4+NAIVE':'T:CD4',

    'mDC:':'DENDRITIC',
    'pDC:':'DENDRITIC'
    }

renamedCellTypes = [] # ordered list of celltypes
for cellType in cellTypes: 
    try:
        renamedCellTypes.append(cellNewTypesDict[cellType])
    except:
        renamedCellTypes.append(cellType)

# filter and create a cluster name dictionary:
#   key: cluster
#   value: cell type
cluster = 1
clusterNameDict = {}
for cellName in renamedCellTypes:
    if type(cellName) == str: # removes the nan values
        clusterNameDict.update({cluster:cellName})
    else:
        pass
    cluster += 1


#####################################################################
'''
Output: geneID_cellType_Dict
    A dictionary of key: geneID
                    value: Cell type
'''
geneNames = clusterGeneDF[0].tolist()
geneCluster = clusterGeneDF[1].tolist()

geneID_cellType_Dict = {}

position = 0
for geneID in geneNames:
    try:
        geneID_cellType_Dict.update({geneID:clusterNameDict[(geneCluster[position])]})
    except:
        pass
    position += 1

    
#####################################################################
# Making of max genes matrix
#####################################################################
'''
immunceCellGenes and nonImmuneCellGenes in a list
'''
immuneCellGenes = list(geneID_cellType_Dict.keys())
nonImmuneCellGenes = [geneName for geneName in geneNames if geneName not in immuneCellGenes]
nonImmuneCellGenes_iterative = nonImmuneCellGenes[:]

'''
Remove non-immune cell types 
'''
for gene in nonImmuneCellGenes_iterative:
    if gene not in sc_counts_matrix.index.values:
        nonImmuneCellGenes.remove(gene)
        
sc_counts_matrix_filtered1 = sc_counts_matrix.drop(nonImmuneCellGenes)
rowNames = sc_counts_matrix_filtered1.index.values.tolist()

dropping = []
for rowName in rowNames:
    if rowName in nonImmuneCellGenes:
        dropping.append(rowName)
    elif rowName not in geneID_cellType_Dict.keys():
        dropping.append(rowName)

sc_counts_matrix_filtered2 = sc_counts_matrix_filtered1.drop(dropping)


'''
changing row names:
    gapminder.rename(index={0:'zero',1:'one'}, inplace=True)
    print(gapminder.head(4))
'''
sc_counts_matrix_filtered = sc_counts_matrix_filtered2.rename(index={rowName:geneID_cellType_Dict[rowName] for rowName in rowNames if rowName in geneID_cellType_Dict.keys()})


'''
Preparing data for phenotype matrix
    make matrix with row and celltype items:
        [1, 'MONOCYTE']
'''
max_cell_type_list = []

columnNumber = 1
for column in sc_counts_matrix_filtered.columns.values.tolist():
    if sc_counts_matrix[column].max() == 0:
        max_cell_type_list.append([columnNumber, 0]) # problem when all values in column is 0
    else:
        max_cell_type_list.append([columnNumber, sc_counts_matrix_filtered[column].idxmax()]) # problem when all values in column is 0
    columnNumber += 1


'''
Generating and saving phenotype matrix where rows are samples and values are as follows:
    0: ignore comparison
    1: membership in class
    2: class that sample will be compared against
'''
cellTypeRowOrder = [
    'B', 
    'DENDRITIC', 
    'MACROPHAGE', 
    'MONOCYTE', 
    'NEUTROPHIL',
    'MAST', 
    'T:CD4', 
    'T:CD8', 
    'T:Reg',
    'NK'
]


totalColumnsValues = []

for column_cellType in max_cell_type_list:
    columnsValues = [2] * len(cellTypeRowOrder) # makes column filled with 2's
    try:
        columnsValues[cellTypeRowOrder.index(column_cellType[1])] = 1 # switch 2 to 1 where is correct cellType
    except:
        pass
    totalColumnsValues.append([column_cellType[0], columnsValues])
   

phenotype_classes = pd.DataFrame({
    'B':[], 
    'DENDRITIC':[], 
    'MACROPHAGE':[], 
    'MONOCYTE':[], 
    'NEUTROPHIL':[],
    'MAST':[], 
    'T:CD4':[], 
    'T:CD8':[], 
    'T:Reg':[],
    'NK':[]
}).transpose()

for column in totalColumnsValues:
    phenotype_classes[column[0]] = column[1]

    
# saving GBM phenotype class file
phenotype_classes.to_csv('GBM_phenotype_class.txt', header=False, index=True, sep='\t', mode='a')
