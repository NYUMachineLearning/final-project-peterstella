library(caret)
library(ROCR)
library(pROC)
library(MASS)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggfortify)
library(glmnet)
library(e1071)
library(mlbench)
library(tidyverse)
#library(ISLR)
library(tree)
library(randomForest)
library(gbm)
library(ranger)
library(Hmisc)
library(FSelector)


#Load data and modify column names
load('/gpfs/scratch/pas376/HTA20_RMA.RData')
ano=read.csv("anoSC1_v11_nokey.csv",stringsAsFactors = FALSE)
data_df <- as.data.frame(eset_HTA20)
annotations <- ano #%>% column_to_rownames(., var= "SampleID")
data_transposed <- as.data.frame(t(data_df))
data_transposed <- rownames_to_column(data_transposed)
data_transposed <- dplyr::rename(data_transposed, SampleID = rowname)
data_annotated <- inner_join(annotations,data_transposed)
data_2rownames <- data_annotated %>% column_to_rownames(.,var = "SampleID")
data_train <- filter(data_annotated, Train == 1)  
data_test <- filter(data_annotated, Train == 0)
data_train_rownames <- filter(data_2rownames, Train == 1)  
data_test_rownames <- filter(data_2rownames, Train == 0)  

# select only data from first batch
#data_train1 <- data_train[1:335,]
data_train1 <- data_train
data_train1 <- data_train1 %>% dplyr::select(-"Batch",-"Set",-"Train", -"Platform")
data_train1 <- data_train1 %>% mutate(term = ifelse(GA<37,0,1))
data_train1 <- data_train1 %>% dplyr::select("SampleID", "GA", "term", everything())
data_train_term <- data_train1 %>% dplyr::select(-"GA")
data_train_GA <- data_train1 %>% dplyr::select(-"term")

#data_train_rownames1 <- data_train_rownames[1:335,]
data_train_rownames1 <- data_train_rownames
data_train_rownames1 <- data_train_rownames1 %>% dplyr::select(-"Batch",-"Set",-"Train", -"Platform")
data_train_rownames1 <- data_train_rownames1 %>% mutate(term = ifelse(GA<37,0,1))
data_train_rownames1 <- data_train_rownames1 %>% dplyr::select( "GA", "term", everything())
data_train_rownames_term <- data_train_rownames1 %>% dplyr::select(-"GA")
data_train_rownames_GA <- data_train_rownames1 %>% dplyr::select(-"term")
data_test_rownames1 <- data_test_rownames %>% dplyr::select(-"Batch",-"Set",-"Train", -"Platform")
data_test_rownames1 <- data_test_rownames1 %>% mutate(term = ifelse(GA<37,0,1))
data_test_rownames1 <- data_test_rownames1 %>% dplyr::select( "GA", "term", everything())
data_test_rownames_term <- data_test_rownames1 %>% dplyr::select(-"GA")
data_test_rownames_GA <- data_test_rownames1 %>% dplyr::select(-"term")

data_train_rownames_GA_label <- data_train_rownames_GA %>% select("GA")
data_train_rownames_term_label <-data_train_rownames_term %>% select("term")
data_train_rownames_GA_no_label <- data_train_rownames_GA %>% select(-"GA")
data_train_rownames_term_no_label <- data_train_rownames_term %>% select(-"term")
GA_no_lab_mat <- data.matrix(data_train_rownames_GA_no_label, rownames.force = TRUE)
GA_lab_mat <- data.matrix(data_train_rownames_GA_label, rownames.force = TRUE)
GA_no_lab_mat2 <- GA_no_lab_mat[,1:32800]

data_test_rownames_GA_label <- data_test_rownames_GA %>% select("GA")
data_test_rownames_term_label <-data_test_rownames_term %>% select("term")
data_test_rownames_GA_no_label <- data_test_rownames_GA %>% select(-"GA")
data_test_rownames_term_no_label <- data_test_rownames_term %>% select(-"term")
GA_test_no_lab_mat <- data.matrix(data_test_rownames_GA_no_label, rownames.force = TRUE)
GA_test_lab_mat <- data.matrix(data_test_rownames_GA_label, rownames.force = TRUE)

#convert numeric scoring to factor for term
data_train_term <- data_train_term %>% mutate(term = ifelse(term==1, "yes", "no"))
data_train_term$term <- as.factor(data_train_term$term)
data_train_rownames_term <- data_train_rownames_term %>% mutate(term = ifelse(term==1, "yes", "no"))
data_train_rownames_term$term <- as.factor(data_train_rownames_term$term)


#PCA
#data_pca_GA <- data_train1 %>% dplyr::select(-"SampleID", -"GA", -"term")
#pca_GA <- prcomp(data_pca_GA,center = TRUE, scale. = TRUE)

#FSelector

lct <- linear.correlation(GA~., data_train_rownames_GA)
chi <-chi.squared(term~., data_train_rownames_term)


#Build datasets from feature selection 

lct <- rownames_to_column(lct)
lct <- rename(lct, gene = rowname)
lctest_ordered <- lct %>% arrange(desc(attr_importance))
lct_top100 <- lctest_ordered[1:100,]
lct_top200 <- lctest_ordered[1:200,]
lct_top500 <- lctest_ordered[1:500,]
lct_top1000 <- lctest_ordered[1:1000,]

#ordered lists of genes, subsets of database by tops 

lct_100 <- lct_top100$gene
lct_200 <- lct_top200$gene
lct_500 <- lct_top500$gene
lct_1000 <- lct_top1000$gene

data_train_GA_lc_top100 <- subset(data_train_GA, select= lct_100)
data_train_GA_lc_top200 <- subset(data_train_GA, select= lct_200)
data_train_GA_lc_top500 <- subset(data_train_GA, select= lct_500)
data_train_GA_lc_top1000 <- subset(data_train_GA, select= lct_1000)

chi <- rownames_to_column(chi)
chi <- rename(chi, gene = rowname)
chitest_ordered <- chi %>% arrange(desc(attr_importance))
chi_top100 <- chitest_ordered[1:100,]
chi_top200 <- chitest_ordered[1:200,]
chi_top500 <- chitest_ordered[1:500,]
chi_top1000 <- chitest_ordered[1:1000,]

chi_100 <- chi_top100$gene
chi_200 <- chi_top200$gene
chi_500 <- chi_top500$gene
chi_1000 <- chi_top1000$gene

data_train_term_chi_top100 <- subset(data_train_term, select= chi_100)
data_train_term_chi_top200 <- subset(data_train_term, select= chi_200)
data_train_term_chi_top500 <- subset(data_train_term, select= chi_500)
data_train_term_chi_top1000 <- subset(data_train_term, select= chi_1000)


save.image(file="univar2.RData")
saveRDS(lct, file = "lct2.rds")
saveRDS(chi, file = "chi2.rds")
saveRDS(pca_GA, file = "pca_GA2.rds")
save(lctest,chitest, rf_san_GA, rf_san_term, ranger_san_GA, ranger_san_term, file = "sancheck2.RData")

