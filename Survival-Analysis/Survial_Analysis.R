library(tidyverse)
library(stringr)
library(broom)
require(survival)
setwd("/home/juejuew/project")
load("tcga.gbm.RData")
############################################################
#################Prepare data for analysis##################
############################################################
#load survival data
surv_dat <- download_tcga[[3]]
#read in cell proportion data from the output of Timer deconvolution result
cell_composition <- read.table("/home/yujiew/Project/timer_composition.txt",header=T)
cell_composition <- data.frame(cell_composition)

##Select the samples who contain >=4 immune cell types(totally 6 immune cell types) from cell_composition
cell_composition2 <- cell_composition[apply(cell_composition[,-1] == 0, 1, sum) <= 2, ] %>%
  mutate(sampleID = as.character(sampleID))
##Select the samples who contain >=4 immune cell types(totally 6 immune cell types) from survival data
r <- anti_join(surv_dat %>% mutate(barcode = str_replace_all(barcode, "-", ".")),
               cell_composition2, by=c("barcode" = "sampleID"))
throw_names <- r$barcode
surv_dat_new <- surv_dat %>% mutate(barcode = str_replace_all(barcode, "-", ".")) %>% 
  filter(!(barcode %in% throw_names))

###filter out samples missing information
#idx without time information
surv_dat <- as.data.frame(surv_dat_new)
idx_time <- which(is.na(as.numeric(surv_dat$days_to_last_follow_up)) & is.na(as.numeric(surv_dat$days_to_death)))
#idx without survival information
idx_sur <- which(is.na(surv_dat$vital_status))
# idx with wrong survival information
time <- ifelse(surv_dat$vital_status %in% c("Alive", "alive"), surv_dat$days_to_last_follow_up, surv_dat$days_to_death)
censor <- ifelse(surv_dat$vital_status %in% c("Alive", "alive"), 0, 1)
idx_err <- which(is.na(time) | is.na(censor))
# idx without follwoing up information
idx_f <- which(time==0)
idx <- unique(c(idx_time, idx_sur, idx_err, idx_f))
s <- Surv(time, censor)[-idx]

cell_composition <- cell_composition2[-idx,]

#final cell_composition after filtering out invalid samples
cell_composition[1:5,1:5]
dim(cell_composition)#126 rows(samples), 7 columns(6 cell types)

#########################################################
################optimize number of clusters##############
#########################################################
###Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
data <- cell_composition[,-1]
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=10000,iter.max = 100000)$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

###Log rank test for finding the optimize number of clusters
vec <- vector("numeric")
i = 2
while(i <=15) {
  test_group <- kmeans(cell_composition[,-1], i, nstart=10000, iter.max=10000000)$cluster
  chisq_val <- survdiff(s ~ test_group)$chisq
  vec[i - 1] = logrank <- 1 - pchisq(chisq_val, i-1) 
  i = i+1
}
length(vec)
vec
compare = tibble(
  number_of_clusters = 2:15,
  logrank = vec
)

compare %>% ggplot +
  geom_line(aes(x = number_of_clusters, y = logrank)) +
  scale_x_continuous(limits = c(2,15), breaks = c(2:15)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
  geom_hline(yintercept=vec[which.min(vec)], linetype="dashed", color = "blue") 

#optimization, choosing the cluster number that has the least logrank
num_clust <-  which.min(vec) + 1

#####################################################
###################k-means clustering################
#####################################################

#subset the data into 5 clusters
group <- kmeans(cell_composition[,-1], num_clust, nstart=10000, iter.max=10000000)$cluster

#filter out samples missing information from survival data
set.seed(123)
test_surv = data.frame(surv_dat)
test_surv[1:5,1:5]
test_s <- test_surv[-idx,]

#final survival data after filtering out invalid samples
test_s[1:5,1:5]
dim(test_s)#126 rows(samples), 115 columns(features)

#Add cluster numbers to survival data
group = data.frame(group) %>% mutate(r = row_number())
test_s = test_s %>% mutate(r = row_number()) %>% left_join(group, by = "r") %>%
  select(sample, patient, group, everything())
test_s[1:5,1:5]
test_s %>% count(group)
#Create new columns for survival time and censor status
test_s = test_s %>% 
  mutate(time = ifelse(vital_status %in% c("Alive","alive"), 
                       days_to_last_follow_up, days_to_death)) %>%
  mutate(censor = ifelse(vital_status %in% c("Alive", "alive"), 0, 1)) %>% 
  select(sample, patient, group, time, censor,everything())
test_s <- test_s[,-117]

#Combined survival data(used for survial analysis, K-M curve)
test_s[1:5,1:5]
dim(test_s)

################draw graph for k-means clustering#####################
#install.packages("Rtsne")
#devtools::install_github("jkrijthe/Rtsne")

###Basic graph
library(Rtsne)
class(cell_composition)
cell_composition[1:5,1:5]
timer_output_unique <- unique(cell_composition)
set.seed(50)
tsne_out <- Rtsne(as.matrix(timer_output_unique[,-1]))
d_tsne_1 = as.data.frame(tsne_out$Y)
ggplot(d_tsne_1, aes(x=V1, y=V2)) +  
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")

d_tsne_1_original = d_tsne_1

###Colorful t-SNE graph by groups
fit_cluster_kmeans=kmeans(scale(d_tsne_1), num_clust)  
library(factoextra)
fviz_cluster(fit_cluster_kmeans, cell_composition[,-1])
group = fit_cluster_kmeans$cluster
cell_composition_with_clusters = cell_composition %>% mutate(cluster = group) %>% count(cluster) %>% print()
sum(cell_composition_with_clusters$n)
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

## setting num_clust clusters as output
d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=num_clust))  
#Plotting the cluster models onto t-SNE output
#Now time to plot the result of each cluster model, based on the t-SNE map.

plot_cluster=function(data, var_cluster, palette)  
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=4) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") #+ 
  #scale_colour_brewer(palette = palette) 
}


plot_k=plot_cluster(d_tsne_1_original, "cl_kmeans", "Accent")  
plot_h=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")

## and finally: putting the plots side by side with gridExtra lib...
library(gridExtra)  
grid.arrange(plot_k, plot_h,  ncol=2)

##############################################################
###################PCA-plot to locate outliers################
##############################################################

###combine survival data with cell type composition from Timer result
copy <- test_s
paste <- cell_composition[,-1]
names(paste)
bind <- cbind(copy,paste)
dim(bind)
bind[1:5,1:5]
bind[,3] = as.factor(bind[,3])

###PCA
comp.pca <- prcomp(bind[,c(117:122,4,5)],scale=TRUE)
comp.pca

#install.packages("ggfortify")
library(ggfortify)
pca.plot <- autoplot(comp.pca, data = bind,colour="group")
pca.plot

###identify outliers
#install.packages("factoextra")
library(factoextra)
fviz_pca_var(comp.pca)


fviz_pca_var(comp.pca, col.var="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.5) + theme_minimal()

fviz_pca_ind(comp.pca, col.ind="cos2") +
  scale_color_gradient2(low="cyan", mid="blue", 
                        high="dark blue", midpoint=0.50) + theme_minimal()

p2 <- fviz_pca_ind(comp.pca, label="none", habillage=bind$group, 
                   addEllipses=TRUE, ellipse.level=0.95) %>% print


fviz_pca_ind(comp.pca) +
  labs(title ="PCA", x = "PC1", y = "PC2")

##############################################################
###################Cox Regression Model#######################
##############################################################
library("survival")
library("survminer")

###Analyze various factors seperately
#group
res.cox <- coxph(Surv(time, censor) ~ factor(group), data = bind)
summary(res.cox)

#immune cell types
covariates <- c("B_cell", "T_cell.CD4", "T_cell.CD8", "Neutrophil","Macrophage","DC")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, censor)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = bind)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res

############################################################################
###################Cox beta comparison for prediction#######################
############################################################################

###Compare beta among all 6 immune cells(seperately)
covariates <- c("B_cell", "T_cell.CD4", "T_cell.CD8", "Neutrophil","Macrophage","DC")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, censor)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = bind)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         y = x %>% tidy()
                         estimate=y$estimate
                         conf.low=y$conf.low
                         conf.high=y$conf.high
                         res<-c(estimate, conf.low, conf.high)
                         names(res)<-c("estimate", "low.conf", 
                                       "high.conf")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res = as.data.frame(univ_results)
res = t(res)
new_res = data.frame(res) %>% mutate(cell_types = rownames(res)) %>% select(cell_types, everything())
new_res

#Draw graph
ggplot(new_res, aes(x = cell_types, y = estimate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = high.conf, ymin = low.conf)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  coord_flip() + labs(y = "estimate") +
  labs(y = "beta_estimate", title = "Seperately")

###T_cell.CD8 + Macrophage, beta interval
res.cox2 <- coxph(Surv(time, censor) ~ T_cell.CD8 + Macrophage, data =  bind)
summary(res.cox2)
conf = res.cox2 %>% tidy() %>%
  mutate(
    estimate=estimate,
    conf.low=conf.low,
    conf.high=conf.high
  ) %>%
  select(term, estimate, starts_with("conf"))

#Draw graph
require(ggplot2)
ggplot(conf, aes(x = term, y = estimate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = conf.high, ymin = conf.low)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  coord_flip() + labs(title = "jointly T_cell.CD8 + Macrophage")