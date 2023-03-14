setwd("/Users/kojiro/desktop/qbio490/qbio_490_KainingFeng/Midterm_project_Feng")
dir.create("outputs")
setwd("outputs")
clinical <- read.csv("/Users/kojiro/Desktop/qbio490/qbio_490_KainingFeng/analysis_data/brca_clinical_data.csv")
#Install and load all libraries
if(!require(maftools)){
  install.packages("maftools")
}
library(maftools)
if(!require(TCGAbiolinks)){
  install.packages("TCGAbiolinks")
}
library(TCGAbiolinks)
if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)
if(!require(survival)){
  install.packages("survival")
}
library(survival)
if(!require(survminer)){
  install.packages("survminer")
}
library(survminer)
if (!require("SummarizedExperiment", quietly = TRUE))
  BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
# query and prepare the clinical data from GDCdata
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")
clinic <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
#query and prepare the mutation data from GDC data frame
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(maf_query)

maf <- GDCprepare(maf_query) # as long as it runs, ignore any errors

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)
GDCdownload(rna_query)
#calculate the mdium age among the patients and determine the "old" and "young" patients
median_age <- median(as.numeric(maf_object@clinical.data$age_at_initial_pathologic_diagnosis))
print(median_age)

maf_object@clinical.data$age_category <- ifelse(maf_object@clinical.data$age_at_initial_pathologic_diagnosis > median_age,'Old', 'Young')
#create young mask and old mask
young_mask <-ifelse(maf_object@clinical.data$age_category == "Young",T,F)
old_mask <- ifelse(maf_object@clinical.data$age_category == "Old",T,F)

#apply the mask to maf_object
young_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[young_mask]
young_maf <- subsetMaf(maf = maf_object,
                       tsb = young_patient_barcodes)
old_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[old_mask]
old_maf <- subsetMaf(maf = maf_object,
                     tsb = old_patient_barcodes)
#find the top10 mutation genes in the old and young groups using Oncoplots
oncoplot(maf=old_maf,top=10)
oncoplot(maf=young_maf,top=10)
coOncoplot(m1 = old_maf, 
           m2 = young_maf,
           genes = c("PIK3CA", "TP53","GATA3","TTN","CDH1","MUC16", "KMT2C", "MAP3K1", 
                     "FLG","RYR2", "PTEN"),
           m1Name = "Old Patients", 
           m2Name = "Young Patients", 
           borderCol = NA)
maf_object@clinical.data$Overall_Survival_Status <-
  ifelse(maf_object@clinical.data$vital_status == "Alive",
         maf_object@clinical.data$Overall_Survival_Status <- TRUE,
         maf_object@clinical.data$Overall_Survival_Status <- FALSE)
# draw the KM plot to compare the survival rates of mutant and wildtype gene
par(mar=c(1,1,1,1))
mafSurvival(maf = maf_object,
            genes = "GATA3",
            time = "days_to_last_followup", 
            Status = "Overall_Survival_Status", 
            isTCGA = TRUE)



# query and prepare rna sequencing data frames
rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
rna_se <- GDCprepare(rna_query,directory = "/Users/kojiro/GDCdata")
#separate the coldata from rna_se
rna_clinical <- rna_se@colData
rna_clinical <- as.data.frame(rna_clinical)
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments", F,T)
rna_clinical <- rna_clinical[,treatments_mask]
#separate row ranges from the rna_se
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
#separate the unstrand DNA from the rna_se and change its column names to barcodes and rownames to gene_id. 
rna_counts <- rna_se@assays@data$unstranded
rna_counts <- as.data.frame(rna_counts)
row.names(rna_clinical) <- rna_clinical$barcode
rna_clinical$age_category <- ifelse(rna_clinical$age_at_index > 58, "Old", "Young")
row.names(rna_genes) <- rna_genes$gene_id
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id
unique(rna_clinical$definition)
sample_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal",F,T)
#apply the solid tissue sample mask to the rna_clinical and rna_counts
rna_clinical <- rna_clinical[sample_mask,]
rna_counts <- rna_counts[,sample_mask]

#make a GATA3 mask and apply it to the rna_counts
GATA3_mask <- ifelse(rna_genes$gene_name == "GATA3",T,F)
GATA3_ensembl <- rna_genes$gene_id[GATA3_mask]
GATA3_counts <- rna_counts[GATA3_ensembl,]
#make a old and young mask to the clinical data
new_old_mask <- ifelse(rna_clinical$age_category == "Old",T,F)
new_young_mask <- ifelse(rna_clinical$age_category == "Young",T,F)
old_ensembl <- rownames(rna_clinical[new_old_mask,])
young_ensembl <- rownames(rna_clinical[new_young_mask,])
#combine the age category mask and GATA3 mask in the rna_counts data frame
#to compare the GATA3 gene expression in old and young groups
GATA3_old_mask <- ifelse(colnames(GATA3_counts) %in% old_ensembl,T,F)
GATA3_old_counts <- GATA3_counts[,GATA3_old_mask]
GATA3_young_mask <- ifelse(colnames(GATA3_counts) %in% young_ensembl,T,F)
GATA3_young_counts <- GATA3_counts[,GATA3_young_mask]
GATA3_old_counts <- as.integer(GATA3_old_counts)
GATA3_young_counts <- as.integer(GATA3_young_counts)
summary(GATA3_old_counts)
summary(GATA3_young_counts)
age_counts <- data.frame(GATA3_old_counts[1:605],GATA3_young_counts) 
# make boxplots to compare the GATA3 gene expression in old and young groups
par(mar=c(1,1,1,1))
boxplot(age_counts,
        xlab = "Age_categories",
        ylab = "Counts")


BiocManager::install("DESeq2",force = TRUE)
library(DESeq2)
#renew the rna_counts, rna_clinic, and rna_gene
rna_counts <- rna_se@assays@data$unstranded[,!is.na(rna_se@colData$age_at_index)]
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
row.names(rna_genes) <- rna_genes$gene_id
rna_clinical <- rna_se@colData[!is.na(rna_se@colData$age_at_index),]
rna_clinical <- as.data.frame(rna_clinical)
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments", F,T)
rna_clinical <- rna_clinical[,treatments_mask]
# use rowSums() to create a list with the total number of counts of each gene
row_sums <- rowSums(rna_counts)

# create a boolean mask where genes with < 10 total counts are FALSE, and genes with >= 10 total counts are TRUE
low_counts_mask <- ifelse(row_sums<10,F,T)

# rewrite the rna_counts df, subsetting for only genes with >= 10 total counts
rna_counts <- rna_counts[low_counts_mask,]

#update rna_genes with the low_counts_mas
rna_genes <- rna_genes[low_counts_mask,]
#set the covariables and main varible to factors
rna_clinical$gender <- factor(rna_clinical$gender, levels = c("male","female"))
rna_clinical$age_category <- ifelse(rna_clinical$age_at_index > 58, "Old", "Young")
rna_clinical$age_category <- factor(rna_clinical$age_category, levels = c("Young","Old"))
rna_clinical$ajcc_pathologic_stage <- factor(rna_clinical$ajcc_pathologic_stage, levels = c("Stage IA", "Stage IB","Stage I","Stage IIA","Stage IIB","Stage IIIA", "Stage IIIB","Stage IV", "Stage X"))
#eliminate any NA values in rna_counts and rna_clinical
na_mask <- ifelse(is.na(rna_clinical$ajcc_pathologic_stage),F,T)
rna_clinical <- rna_clinical[na_mask,] 
rna_counts <- rna_counts[,na_mask]
rownames(rna_clinical) <- colnames(rna_counts)
#apply RNA sequencing with gender and ajcc pathologic stage as covariables and age categories as main variable
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~ gender + ajcc_pathologic_stage + age_category)
dds_obj <- DESeq(dds)
resultsNames(dds_obj) 
rna_clinical$age_category
# get the young vs. old comparison and set the young group as control 
results <- results(dds_obj, format = "DataFrame", contrast = c("age_category", "Old", "Young")) 
results <- data.frame(rna_genes[(rownames(results) %in% rna_genes$gene_id),"gene_name"], 
                      rownames(results), results$log2FoldChange, results$pvalue,results$padj,-log10(results$padj))
colnames(results) <- c("gene_name", "gene_id", "log2foldchange", "pvalue","padj", "-log10(padj)")
sig_results <- results[ifelse(results$padj <0.05, T,F),]
log2foldchange_order = order(results$log2foldchange) # order by column "y"
log2foldchange_order
# this rewrites the df based on the sorted rows
results = results[log2foldchange_order, ] 
# order the genes based on the log2 fold change to find the most up/down regulated gene
mask_log2foldchange <- ifelse(results$log2foldchange >1, T,F)
up_results <- results[mask_log2foldchange,]
log2foldchange_order_larger1 <- order(up_results$log2foldchange) 
up_reg_results <- up_results[nrow(up_results)-log2foldchange_order_larger1 +1, ]
up_reg_results

mask_log2foldchange <- ifelse(results$log2foldchange < -1, T,F)
down_results <- results[mask_log2foldchange,]
log2foldchange_order_smaller <- order(down_results$log2foldchange) 
down_reg_results <- down_results[log2foldchange_order_smaller, ]
down_reg_results
# modify the rownames of results
rownames(up_reg_results) <- up_results[rownames(up_reg_results),"gene_id"]
rownames(down_reg_results) <- down_reg_results[rownames(down_reg_results),"gene_id"]
row.names(results) = results[rownames(results),"gene_id"]
#make a volcano plot
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(results,
                lab = results$gene_name,
                pCutoff = 0.05,
                x = 'log2foldchange',
                y = 'padj',  
                pointSize = 1.0,
                labSize = 3.0,
                col=c('black', 'green', 'blue', 'red3'),
                colAlpha = 1)





