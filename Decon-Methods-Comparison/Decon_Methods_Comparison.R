library(tidyverse)
library(ggplot2)

############################################################
##Comparing Deconvolution Methods
############################################################

####Data the deconvolution methods use###

##single cell RNAseq data is obtained from brain tissue.
## Single cell data set from Gene Expression Omnibus (GEO)
## The paper describing the single cell RNAseq data set:
## Darmanis, S. et al. A survey of human brain transcriptome diversity at the single cell level.
## Proc. Natl. Acad. Sci. U. S. A. 112, 7285–7290 (2015)

## Bulk RNAseq data is obtained from the GTEx study.
## A description of the GTEx study is available in the following papers and website:
## The Genotype-Tissue Expression (GTEx) project, Nature Genetics, 2013
## The Genotype-Tissue Expression (GTEx) pilot analysis: Multitissue gene regulation in humans, Science, 2015
## https://www.gtexportal.org/home/

#########################CIBERSORT####################
#load Data
Cibersort_output = read.table("/home/juejuew/project/CIBERSORT/Cibersort_output.txt", header = T,fill=T)
truth = read.table("/home/juejuew/project/Music/truth.txt")

#select cell type columnes 
Cibersort = Cibersort_output[,2:10]
names(truth)
names(Cibersort)

#Prepare data frame
new_truth = truth %>% mutate(neurons = sum(glutamatergic_neurons_from_the_PFC, 
                                           granule_neurons_from_the_hip_dentate_gyrus_region,
                                           pyramidal_neurons_from_th_hip_CA_region)) %>%
  dplyr::select(astrocytes,endothelial,OPC,oligodendrocytes,neurons) %>% mutate(row_num = row_number()) %>%
  gather(cell_type, proportion, -row_num)
new_Cibersort= Cibersort %>% mutate(row_num = row_number()) %>% 
  select(row_num,  astrocytes:OPC) %>% gather(cell_type, proportion, -row_num) %>% 
  filter(cell_type %in% c("astrocytes","endothelial","OPC","oligodendrocytes", "neurons")) #%>% 

#Draw graph
new_Cibersort  %>% ggplot + 
  geom_boxplot(aes(x = cell_type, y =proportion)) + 
  geom_boxplot(data = new_truth,aes(x = cell_type, y = proportion), color = "red") + 
  labs(title = "CIBERSORT") +
  scale_y_continuous(breaks = seq(0, 1, by=0.25), limits = c(0,1))
new_Cibersort  %>% ggplot + 
  geom_boxplot(aes(x = cell_type, y =proportion)) +
  geom_jitter(aes(x = cell_type, y =proportion,color = cell_type), size = 0.5, alpha = 0.5) +
  geom_boxplot(data = new_truth,aes(x = cell_type, y = proportion), color = "red") + 
  labs(title = "CIBERSORT") +
  scale_y_continuous(breaks = seq(0, 1, by=0.25), limits = c(0,1)) +
  scale_color_manual(values=c('#7fc97f','#beaed4','#fb9a99','#1f78b4','#fdc086'))

#########################Epic####################
#prepare data for epic
new_truth_epic = new_truth %>% filter(cell_type == "endothelial") %>% mutate(cell_type = "Endothelial")
epic_estimate = read.table("/home/juejuew/project/epic/cellFractions.txt", header = T)
new_epic= epic_estimate %>% mutate(row_num = row_number()) %>% 
  select(row_num,  Bcells:otherCells) %>% gather(cell_type, proportion, -row_num)

#Draw graph
new_epic %>% ggplot(aes(x = cell_type, y =proportion)) + 
  geom_boxplot() + coord_flip() + 
  geom_boxplot(data = new_truth1,color = "red") + 
  labs(title = "EPIC")
new_epic  %>% ggplot(aes(x = cell_type, y =proportion)) + 
  geom_violin() +
  geom_jitter(size = 0.001) + coord_flip() + 
  geom_boxplot(data = new_truth1, color = "red") + labs(title = "EPIC")

#########################MuSIC####################
#prepare data for MuSIC
Music_estimate = read.table("/home/juejuew/project/Music/estimate.txt", header = T)
new_Music_estimate= Music_estimate %>% mutate(row_num = row_number()) %>% 
  select(row_num,  astrocytes:OPC) %>% gather(cell_type, proportion, -row_num) %>% 
  filter(cell_type %in% c("astrocytes","endothelial","OPC","oligodendrocytes", "neurons")) 

#Draw graph
new_Music_estimate  %>% ggplot(aes(x = cell_type, y =proportion)) + 
  geom_boxplot() + 
  geom_boxplot(data = new_truth,color = "red") + labs(title = "MuSIC")
new_Music_estimate  %>% ggplot + 
  geom_boxplot(aes(x = cell_type, y =proportion)) +
  geom_jitter(aes(x = cell_type, y =proportion, color = cell_type), size = 0.5, alpha = 0.5)+
  geom_boxplot(data = new_truth,aes(x = cell_type, y = proportion), color = "red") + 
  labs(title = "Music") +
  scale_y_continuous(breaks = seq(0, 1, by=0.25), limits = c(0,1)) +
  scale_color_manual(values=c('#7fc97f','#beaed4','#fb9a99','#1f78b4','#fdc086'))

#####################DeconRNASeq########################
#Prepare data for DeconRNASeq
DeconRNASeq_estimate = read.table("/home/juejuew/project/DeconRNASeq/DeconRNASeq_estimate.txt", header = T)
new_DeconRNASeq_estimate=DeconRNASeq_estimate %>% mutate(row_num = row_number()) %>% 
  dplyr::select(row_num,  astrocytes:OPC) %>% gather(cell_type, proportion, -row_num) %>% 
  filter(cell_type %in% c("astrocytes","endothelial","OPC","oligodendrocytes", "neurons"))

#Draw graph
new_DeconRNASeq_estimate  %>% ggplot + 
  geom_boxplot(aes(x = cell_type, y =proportion)) + 
  geom_boxplot(data = new_truth,aes(x = cell_type, y = proportion), color = "red") + 
  labs(title = "DeconRNASeq") +
  scale_y_continuous(breaks = seq(0, 1, by=0.25), limits = c(0,1))
new_DeconRNASeq_estimate  %>% ggplot + 
  geom_boxplot(aes(x = cell_type, y =proportion)) +
  geom_jitter(aes(x = cell_type, y =proportion, color = cell_type), size = 0.05,alpha = 0.5)+
  geom_boxplot(data = new_truth,aes(x = cell_type, y = proportion), color = "red") + 
  labs(title = "DeconRNASeq")+
  scale_color_manual(values=c('#7fc97f','#beaed4','#fb9a99','#1f78b4','#fdc086')) +
  scale_y_continuous(breaks = seq(0, 1, by=0.25), limits = c(0,1))


###################Calculate RMSE for Comparison####################
#1)CIBERSORT
#install.packages("mefa")
library(mefa)
new_truth_RMSE = truth %>% mutate(neurons = sum(glutamatergic_neurons_from_the_PFC, 
                                                granule_neurons_from_the_hip_dentate_gyrus_region,
                                                pyramidal_neurons_from_th_hip_CA_region)) %>%
  dplyr::select(astrocytes,endothelial,OPC,oligodendrocytes,neurons)
tool_truth = rep(new_truth,nrow(new_Cibersort)) 
new_Cibersort_RMSE= Cibersort %>% as.matrix() %>% .[,names(new_truth_RMSE)]

#install.packages("hydroGOF")
library(hydroGOF)
CIBERSORT_test_RMSE = rmse(new_Cibersort_RMSE, tool_truth)
sum_CIBERSORT_test_RMSE = sum(CIBERSORT_test_RMSE )

#2)Music
new_Music_estimate= Music_estimate %>% 
  dplyr::select(-label) %>% 
  as.matrix() %>% 
  .[,names(new_truth)] 
Music_test_RMSE = rmse(new_Music_estimate, tool_truth)
sum_Music_test_RMSE = sum(Music_test_RMSE)

#3)Decon
new_DeconRNASeq_estimate=DeconRNASeq_estimate %>% 
  as.matrix() %>% 
  .[,names(new_truth)]
Decon_test_RMSE = rmse(new_DeconRNASeq_estimate, tool_truth)
sum_Decon_test_RMSE = sum(Decon_test_RMSE)

#Draw graph to compare three models
RMSE_three_model = tibble(
  Music = Music_test_RMSE,
  CIBERSORT = CIBERSORT_test_RMSE,
  Decon = Decon_test_RMSE,
)
cell_type = c(names(new_truth_RMSE), "sum")
RMSE_sum = c(sum_Music_test_RMSE, sum_CIBERSORT_test_RMSE, sum_Decon_test_RMSE)
RMSE_three_model = rbind(RMSE_three_model, RMSE_sum)
RMSE_three_model = RMSE_three_model %>% cbind( cell_type,.) %>%
  gather(Method, RMSE, -cell_type)

#Draw graph
RMSE_three_model %>% data.frame(.) %>%  
  ggplot(aes(x = fct_reorder(Method, RMSE), y = RMSE, color = cell_type, group = cell_type)) +
  geom_point(size = 2) +
  labs(title = "Compare Three Deconvolution Method based on RMSE") +
  scale_color_manual(values=c('#7fc97f','#beaed4','#fb9a99','#1f78b4','#fdc086', '#e41a1c'))