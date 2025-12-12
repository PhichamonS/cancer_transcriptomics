#-------------------------------------------------------------------------------
# R Script: Comprehensive Omics Data Analysis for Kidney Cancer (TCGA)
#-------------------------------------------------------------------------------

rm(list = ls())

# 1. Setup and Library Loading

# Utility
library(tidyverse)
library(ggplot2)
library(cowplot)
library(here)

# Bioconductor/Omics Packages
library(BiocManager)
library(biomaRt)
library(limma)
library(edgeR)
library(DESeq2)
library(tidyverse)
library(tximport)
library(rhdf5)

# Survival and ROC Analysis
library(survival)
library(survminer)
library(pROC)

# Classification
library(randomForest) 
library(caret)

# Note:
# biomaRt: Accessing Ensembl annotation
# limma:  Linear Models for Microarray and RNA-Seq Data
# edgeR:  differential analysis of sequence read count data
# DESeq2: Analyzing RNA-seq data
# survual: Survival Analysis

rm(list = ls())

# Define working directory
data_path <- here("data")
setwd(data_path)

# 2. Data Import and Quality Control -------------------------

## 2.1. Import Read Count Data
readCount <- read.table("TCGA_KIDNEY_Gene_ReadCounts.txt", header = TRUE, row.names = 1)
colnames(readCount) <- gsub('\\.', '-', colnames(readCount), perl = TRUE)
readCount <- readCount[, sort(colnames(readCount))]

## 2.2. Import Kallisto-derived Transcript Abundance

# 2.2.1. Retrieve Gene Annotation Map (t2g)
ensembl_mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get transcript id and gene name
t2g <- getBM(
  attributes = c("ensembl_transcript_id_version", "external_gene_name"), 
  mart = ensembl_mart
)

# Use "tximport" to import transcript abundance and convert transcript id to gene name
kallisto_count <- tximport("./output/abundance.h5", type = "kallisto",tx2gene = t2g)
kallisto_df    <- as.data.frame(kallisto_count$counts[-c(1),])
colnames(kallisto_df) <- "kallisto"

# Subset a sample for comparison
sample_count <- readCount %>% dplyr::select("TCGA-KN-8435-11A")
colnames(sample_count) = "sample"

# Merge outer join
combined = merge(sample_count, kallisto_df, by=0, all = TRUE)

ggplot(combined,aes(x=kallisto,y = sample))+
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color='red', linewidth=2)+
  scale_x_log10() + scale_y_log10() + 
  theme_minimal(base_size = 14) +
  labs(x = "Kallisto Counts (Log Scale)", 
       y = "Sample Counts (Log Scale)", 
       title = "Comparison of Kallisto vs. Sample Counts")

# Commend on how it follow one-end or two-end , 
# the count will be double of the given table if it is one-end.
# There is more error at the lower counts.

# 3. Pre-processing, Filtering, and Normalization -----------

## 3.1. Filter Low-Count Genes (QC/Pre-processing)
source("plot_templates.R", local = knitr::knit_global())
readCount <- readCount[!rowSums(is.na(readCount)) == ncol(readCount),]
readCount <- readCount[, colSums(is.na(readCount)) < ncol(readCount)]

# Calculate total read count per gene and log-transform
readCountSum <- as.data.frame(log10(rowSums(readCount, na.rm = TRUE)))
colnames(readCountSum) <- 'sumCount'

# Define a count threshold based on the distribution
count_threshold <- 4.5

ggplot(readCountSum, aes(x=sumCount)) + geom_density()+
  labs(title='Distribution of Gene Counts',
       x='Log10(Sum of Gene Counts)',y='Density')+
  geom_vline(xintercept = count_threshold, linetype='dashed')+ 
  theme_bw()

# Filter out small count
filtcounts <- subset(readCount, readCountSum$sumCount > count_threshold)
nrow(filtcounts)

readCountSum_new <- as.data.frame(log10(rowSums(filtcounts, na.rm = TRUE)))
colnames(readCountSum_new) <- 'sumCount'

ggplot(readCountSum_new, aes(x=sumCount)) + geom_density()+
  labs(title='Distribution of Gene Counts (Before Filtering)',
       x='Log10(Sum of Gene Counts)',y='Density')+
  theme_bw()

## 3.2. Normalization: TMM and Voom
dge <- DGEList(counts = filtcounts)
dge <- calcNormFactors(dge, method= "TMM")

# Voom transformation (Variance modeling at the observational level)
voom_results <- voom(dge, plot = T)
lognormcounts <- as.data.frame(voom_results$E)
normcounts <- as.data.frame(2^lognormcounts)

dge_2 <- DGEList(readCount)
dge_2 <- calcNormFactors(dge_2, method= "TMM")
voom_results_2 <- voom(dge_2, plot=T)

lognormcounts_2 <- as.data.frame(voom_results_2$E)
normcounts_2 <- as.data.frame(2^lognormcounts_2)

# 4. Dimensionality Reduction (PCA) and Clinical Annotation ---------------------

## 4.1. Import Clinical Data
clinicalAnnotation = read.table('TCGA_KIDNEY_ClinicalAnnotation.txt',sep = '\t',
                                header = 1)

## 4.2. Principal Component Analysis (PCA)
PCA <- prcomp(t(lognormcounts), scale=FALSE, center=TRUE)

# Extract PC scores
PCACounts <- PCA$x
PCACounts <- as.data.frame(PCACounts)[,1:50]

# Create new column with patient ID
PCACounts$patientID <- row.names(PCACounts)
PCACounts$patientID <- substr(PCACounts$patientID,1,nchar(PCACounts$patientID)-4)
PCACounts$sample <- row.names(PCACounts)

# Merge PCA count with clinical annotation
pcadata <- merge(PCACounts, clinicalAnnotation, by="patientID", all.x=TRUE)
pcadata$tissuetype <- str_sub(pcadata$sample, start=-3)

## 4.3. PCA Visualization

# 4.3.1. PCA by Cancer Type
ggplot(pcadata, aes(x=PC1, y=PC2, color = Type))+geom_point()+
  theme_bw()+
  labs(title= "PCA of Gene Expression Data by Cancer Types", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)")

# 4.3.2. PCA by Tissue Type (Tumor vs. Normal)
ggplot(pcadata, aes(x=PC1, y=PC2, color = tissuetype))+geom_point()+
  theme_bw()+
  labs(title= "PCA of Gene Expression Data by Tissue Types", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)")


# 5. Survival Analysis -------------------------------------

## 5.1. Prepare Survival Data
data_surv <- clinicalAnnotation %>% 
  dplyr::select(days_to_death.demographic,vital_status.demographic, 
          days_to_last_follow_up.diagnoses,laterality,Type) %>% 
  mutate(event = ifelse(vital_status.demographic=='Alive',0,1),
         time = ifelse(event == 0, days_to_last_follow_up.diagnoses, days_to_death.demographic)) %>% 
  dplyr::select(time, event, laterality, Type)

## 5.2. Perform and Plot Kaplan-Meier Curves
# Survival comparison by Cancer Type pairs
surv1 =  survfit(Surv(time, event) ~ Type, data = data_surv %>% filter(Type %in% c("KICH","KIRC")))
ggsurvplot(surv1, data = data_surv, pval = TRUE)

surv2 =  survfit(Surv(time, event) ~ Type, data = data_surv %>% filter(Type %in% c("KICH","KIRP")))
ggsurvplot(surv2, data = data_surv, pval = TRUE)

surv3 =  survfit(Surv(time, event) ~ Type, data = data_surv %>% filter(Type %in% c("KIRP","KIRC")))
ggsurvplot(surv3, data = data_surv, pval = TRUE)

# Survival comparison by Laterality
surv4 =  survfit(Surv(time, event) ~ laterality, data = data_surv %>% filter(laterality %in% c("Left","Right")))
ggsurvplot(surv4, data = data_surv, pval = TRUE)

# 6. Gene Expression vs. Clinical Labs (ROC/Association) ------------------

## 6.1. Prepare Gene and Lab Data
normcounts_t = as.data.frame(t(normcounts))
normcounts_t$patientID <- row.names(normcounts_t)
normcounts_t$patientID <- substr(normcounts_t$patientID,1,nchar(normcounts_t$patientID)-4)
normcounts_t$sample <- row.names(normcounts_t)

# Merge gene counts with clinical/lab data
genedata = merge(normcounts_t, clinicalAnnotation, by="patientID", all.x=TRUE) %>%
  dplyr::select(VHL, HIF1A, VEGFA,
         hemoglobin_result, platelet_qualitative_result,
         serum_calcium_result, white_cell_count_result)


# Define genes and lab tests of interest
predictor_genes <- c("VHL", "HIF1A", "VEGFA")
lab_tests <- c("hemoglobin_result", "platelet_qualitative_result",
               "serum_calcium_result", "white_cell_count_result")

## 6.2. Boxplot Visualization of Associations
gene_lab_long = genedata %>% 
  pivot_longer(all_of(lab_tests), names_to = "test_name", values_to = "test_result") %>% 
  pivot_longer(all_of(predictor_genes), names_to = "gene_name", values_to = "gene_counts") %>% 
  filter(test_result %in% c("Low","Elevated")) %>%
  mutate(test_result = factor(test_result, levels = c("Low", "Elevated")))

ggplot(gene_lab_long,
       aes(x= gene_counts, y = test_result))+
  geom_boxplot()+ geom_point()+
  facet_wrap(gene_name~test_name, scales = "free")+
  theme_bw()

## 6.3. ROC Analysis for Predictive Power

# divide dataset into training and test set
set.seed(1)
sample <- sample(c(TRUE, FALSE), nrow(genedata), replace=TRUE, prob=c(0.8,0.2))
train <- genedata[sample, ]
test <- genedata[!sample, ]

train = train %>% mutate(hemoglobin_result = ifelse(hemoglobin_result == "Low",0,1),
                         platelet_qualitative_result = ifelse(platelet_qualitative_result == "Low",0,1),
                         serum_calcium_result = ifelse(serum_calcium_result == "Low",0,1),
                         white_cell_count_result= ifelse(white_cell_count_result == "Low",0,1))

test = test %>% mutate(hemoglobin_result = ifelse(hemoglobin_result == "Low",0,1),
                        platelet_qualitative_result = ifelse(platelet_qualitative_result == "Low",0,1),
                        serum_calcium_result = ifelse(serum_calcium_result == "Low",0,1),
                        white_cell_count_result= ifelse(white_cell_count_result == "Low",0,1))



#>>> hemoglobin

#fit logistic regression model to training set
model_1 <- glm(hemoglobin_result~VHL, family="binomial", data=train)
predicted_1 <- predict(model_1, test, type="response")
rocobj_1 <- roc(test$hemoglobin_result, predicted_1)

model_2 <- glm(hemoglobin_result~HIF1A, family="binomial", data=train)
predicted_2 <- predict(model_2, test, type="response")
rocobj_2 <- roc(test$hemoglobin_result, predicted_2)

model_3 <- glm(hemoglobin_result~VEGFA, family="binomial", data=train)
predicted_3 <- predict(model_3, test, type="response")
rocobj_3 <- roc(test$hemoglobin_result, predicted_3)

auc_1 <- auc(rocobj_1)
auc_2 <- auc(rocobj_2)
auc_3 <- auc(rocobj_3)

ggroc(list(VHL = rocobj_1, HIF1A = rocobj_2, VEGFA = rocobj_3)) +
  ggtitle("Hemoglobin Test: ROC Curves") +
  theme_minimal() +
  labs(color = "Predictors") +
  scale_color_manual(
    values = c("VHL" = "blue", "HIF1A" = "green", "VEGFA" = "red"),
    labels = c(
      paste0("VHL (AUC = ", round(auc_1, 2), ")"),
      paste0("HIF1A (AUC = ", round(auc_2, 2), ")"),
      paste0("VEGFA (AUC = ", round(auc_3, 2), ")")
    )
  )

#>>> platelet_qualitative_result
model_1 <- glm(platelet_qualitative_result~VHL, family="binomial", data=train)
predicted_1 <- predict(model_1, test, type="response")
rocobj_1 <- roc(test$platelet_qualitative_result, predicted_1)

model_2 <- glm(platelet_qualitative_result~HIF1A, family="binomial", data=train)
predicted_2 <- predict(model_2, test, type="response")
rocobj_2 <- roc(test$platelet_qualitative_result, predicted_2)

model_3 <- glm(platelet_qualitative_result~VEGFA, family="binomial", data=train)
predicted_3 <- predict(model_3, test, type="response")
rocobj_3 <- roc(test$platelet_qualitative_result, predicted_3)

auc_1 <- auc(rocobj_1)
auc_2 <- auc(rocobj_2)
auc_3 <- auc(rocobj_3)

ggroc(list(VHL = rocobj_1, HIF1A = rocobj_2, VEGFA = rocobj_3)) +
  ggtitle("Platelet Qualitative Test: ROC Curves") +
  theme_minimal() +
  labs(color = "Predictors") +
  scale_color_manual(
    values = c("VHL" = "blue", "HIF1A" = "green", "VEGFA" = "red"),
    labels = c(
      paste0("VHL (AUC = ", round(auc_1, 2), ")"),
      paste0("HIF1A (AUC = ", round(auc_2, 2), ")"),
      paste0("VEGFA (AUC = ", round(auc_3, 2), ")")
    )
  )


#>>> serum_calcium

#fit logistic regression model to training set
model_1 <- glm(serum_calcium_result~VHL, family="binomial", data=train)
predicted_1 <- predict(model_1, test, type="response")
rocobj_1 <- roc(test$serum_calcium_result, predicted_1)

model_2 <- glm(serum_calcium_result~HIF1A, family="binomial", data=train)
predicted_2 <- predict(model_2, test, type="response")
rocobj_2 <- roc(test$serum_calcium_result, predicted_2)

model_3 <- glm(serum_calcium_result~VEGFA, family="binomial", data=train)
predicted_3 <- predict(model_3, test, type="response")
rocobj_3 <- roc(test$serum_calcium_result, predicted_3)

auc_1 <- auc(rocobj_1)
auc_2 <- auc(rocobj_2)
auc_3 <- auc(rocobj_3)

ggroc(list(VHL = rocobj_1, HIF1A = rocobj_2, VEGFA = rocobj_3)) +
  ggtitle("Serum Calcium Test: ROC Curves") +
  theme_minimal() +
  labs(color = "Predictors") +
  scale_color_manual(
    values = c("VHL" = "blue", "HIF1A" = "green", "VEGFA" = "red"),
    labels = c(
      paste0("VHL (AUC = ", round(auc_1, 2), ")"),
      paste0("HIF1A (AUC = ", round(auc_2, 2), ")"),
      paste0("VEGFA (AUC = ", round(auc_3, 2), ")")
    )
  )

#>>> white_cell_count
# Fit logistic regression models for each predictor
model_1 <- glm(white_cell_count_result ~ VHL, family = "binomial", data = train)
predicted_1 <- predict(model_1, test, type = "response")
rocobj_1 <- roc(test$white_cell_count_result, predicted_1)
auc_1 <- auc(rocobj_1)

model_2 <- glm(white_cell_count_result ~ HIF1A, family = "binomial", data = train)
predicted_2 <- predict(model_2, test, type = "response")
rocobj_2 <- roc(test$white_cell_count_result, predicted_2)
auc_2 <- auc(rocobj_2)

model_3 <- glm(white_cell_count_result ~ VEGFA, family = "binomial", data = train)
predicted_3 <- predict(model_3, test, type = "response")
rocobj_3 <- roc(test$white_cell_count_result, predicted_3)
auc_3 <- auc(rocobj_3)

# Create a ggplot of ROC curves with AUC values in the legend
ggroc(list(VHL = rocobj_1, HIF1A = rocobj_2, VEGFA = rocobj_3)) +
  ggtitle("White Cell Count Test: ROC Curves") +
  theme_minimal() +
  labs(color = "Predictors") +
  scale_color_manual(
    values = c("VHL" = "blue", "HIF1A" = "green", "VEGFA" = "red"),
    labels = c(
      paste0("VHL (AUC = ", round(auc_1, 2), ")"),
      paste0("HIF1A (AUC = ", round(auc_2, 2), ")"),
      paste0("VEGFA (AUC = ", round(auc_3, 2), ")")
    )
  )

# 7. Cancer Subtype Classification (Random Forest) -----------------

## 7.1. Prepare Classification Data
# Filter clinical data to include only relevant columns (cancer type: KIRC, KIRP, KICH)
cancer_samples = merge(normcounts_t %>% dplyr::select(-sample), 
                       clinicalAnnotation %>% dplyr::select(patientID, Type), 
                       by="patientID", all.x=TRUE) 

colnames(cancer_samples) <- gsub('-','_',colnames(cancer_samples))

# Partition the dataset into training and test sets
set.seed(123)  # Ensure reproducibility
cancer_samples <-  cancer_samples[complete.cases(cancer_samples), ]
cancer_samples$Type <- as.factor(cancer_samples$Type)


# 7.1.1. Train/Test Split
trainIndex <- createDataPartition(cancer_samples$Type, p = .7, list = FALSE, times = 1)
train_data <- cancer_samples[ trainIndex,]
test_data <- cancer_samples[-trainIndex,]

## 7.2. Feature Selection: Top Variance Genes 
gene_variance <- apply(train_data[, -ncol(train_data)], 2, var)
gene_variance_norm <-gene_variance/sum(gene_variance, na.rm=T)

selected_genes <- data.frame(
  genes = names(sort(gene_variance_norm, decreasing = TRUE)[1:30]),
  variance = sort(gene_variance_norm, decreasing = TRUE)[1:30]
)

# Plot genes variances
ggplot(selected_genes, aes(x= reorder(genes, -variance),y=variance))+geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = "Genes", y = "Variance (normalized)", title = "Top 30 Genes by Variance")

reshape2::melt(cancer_samples %>% dplyr::select(selected_genes$genes,"Type") , id = "Type") %>%
  ggplot(aes(x=variable,y=value, fill=Type))+
  geom_boxplot()

# 7.2.1. Standardization
processCenter <- preProcess(train_data[,selected_genes$genes], method = c("center"))
train_data <- predict(processCenter,train_data)
test_data <- predict(processCenter,test_data)

## 7.3. Model Training (Random Forest)
# Model training with Random Forest on selected genes
rf_model <- randomForest(Type ~ ., data = train_data[, c(selected_genes$genes, "Type")], 
                         importance = TRUE, ntree = 100)

## 7.4. Performance Evaluation and Interpretation
# Predictions on the test set
predictions <- predict(rf_model, test_data[, selected_genes$genes])
conf_matrix <- confusionMatrix(predictions, test_data$Type)

cat("\n--- Random Forest Classification Performance ---\n")
print(conf_matrix)

# Feature Importance
gene_importance <- as.data.frame(importance(rf_model, type = 1)) %>% # Type 1: Mean Decrease Accuracy
  rownames_to_column("Gene") %>%
  rename(MeanDecreaseAccuracy = 'MeanDecreaseAccuracy') %>%
  arrange(desc(MeanDecreaseAccuracy))

cat("\n--- Top 10 Most Important Genes in RF Model ---\n")
print(head(gene_importance, 10))

plot_importance <- ggplot(head(gene_importance, 10), aes(x = reorder(Gene, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "Random Forest: Top 10 Gene Importance",
       x = "Gene", y = "Mean Decrease Accuracy")
print(plot_importance)

# Plot the Confusion Matrix
conf_df <- as.data.frame(conf_matrix$table)

# Plot the confusion matrix using ggplot2
ggplot(conf_df, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Freq), color = "white") +
  geom_text(aes(label = Freq), vjust = 1) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  ggtitle("Confusion Matrix") +
  theme_minimal()

## 7.5. Cross-reference with Biological Signature
# Cross-reference genes with VEGF pro-angiogenic gene signature
vegf_genes_up <- readLines("VEGF_A_UP.V1_UP.v2024.1.Hs.txt")
vegf_genes_dn <- readLines("VEGF_A_UP.V1_DN.v2024.1.Hs.txt")

# Genes in the signature and importance in the model
signature_genes <- selected_genes$genes[selected_genes$genes %in% c(vegf_genes_up)]
cat("Number of pro-angiogenic genes in the selected signature: ", length(signature_genes), "\n")
cat("Pro-angiogenic genes in the classifier's signature:\n", signature_genes)


signature_genes <- selected_genes$genes[selected_genes$genes %in% c(vegf_genes_dn)]
cat("Number of pro-angiogenic genes in the selected signature: ", length(signature_genes), "\n")
cat("Pro-angiogenic genes in the classifier's signature:\n", signature_genes)


# 8. Prediction of Unclassified Samples ---------------------------------------

## 8.1. Prepare Unclassified Samples Data
unclassified_data = merge(normcounts_t %>% dplyr::select(-sample), 
                          clinicalAnnotation %>% dplyr::select(patientID, Type), 
                          by="patientID", 
                          all.x=TRUE) 

colnames(unclassified_data) <- gsub('-','_',colnames(cancer_samples))
unclassified_data <- unclassified_data[is.na(unclassified_data$Type),]

# 8.1.1. Standardize Unclassified Data
# Apply the *same* scaling (center) as derived from the training data
unclassified_data[,selected_genes$genes]=predict(processCenter,unclassified_data[,selected_genes$genes])

## 8.2. Predict Cancer Subtype
pred_NA <- predict(rf_model, unclassified_data[, selected_genes$genes])
unclassified_data$Type <- paste0(pred_NA,"_pred")

# Create a data frame for predicted samples
pcadata_new <-  pcadata
pcadata_new$Type[is.na(pcadata_new$Type)] <- unclassified_data$Type[match(pcadata_new$patientID[is.na(pcadata_new$Type)], unclassified_data$patientID)]

# Plot predicted label in the same plot with reference
custom_colors <- c("KICH" = "red", "KIRC" = "blue", "KIRP" = "green",
                   "KICH_pred" = "darkred", "KIRC_pred" = "darkblue", "KIRP_pred" = "darkgreen")
ggplot() +
  geom_point(data = pcadata, aes(x = PC1, y = PC2, color = Type), alpha=0.3) +
  geom_point(
    data = pcadata_new %>% filter(Type %in% c("KICH_pred", "KIRC_pred", "KIRP_pred")),
    aes(x = PC1, y = PC2, color = Type),
    alpha = 0.8, shape = 17, size = 2 ) +
  theme_bw() +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "PCA of Gene Expression Data by Cancer Types",
    x = "Principal Component 1 (PC1)",
    y = "Principal Component 2 (PC2)"
  )
