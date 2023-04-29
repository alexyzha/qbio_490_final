dir.create("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/final_project_group1/outputs")
setwd("/Users/lilly/Desktop/QBIO490/qbio_490_lilly/final_project_group1/outputs")


if(!require("BiocManager"))  
  install.packages("BiocManager")
if(!require("TCGAbiolinks"))  
  BiocManager::install("TCGAbiolinks")
if(!require("dplyr"))  
  install.packages("dplyr")
if(!require("maftools"))  
  install.packages("maftools")
if(!require("SummarizedExperiment"))  
  install.packages("SummarizedExperiment")
if(!require("ggplot2"))  
  install.packages("ggplot2")
if(!require("DESeq2"))  
  install.packages("DESeq2")
if(!require("tidyr"))  
  install.packages("tidyr")
install.packages(c("data.table", "Biobase", "BiocGenerics", "GenomicRanges"))
if(!require("EnhancedVolcano"))
  install.packages("EnhancedVolcano")
install.packages("ggrepel")
BiocManager::install("apeglm")

library(dplyr)
library(ggplot2)
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(DESeq2)
library(tidyr)
library(data.table)
library(Biobase)
library(BiocGenerics)
library(GenomicRanges)
library(EnhancedVolcano)
library(ggrepel)
library(apeglm)

# load in clinical data from TCGA 
clin_query <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Clinical",
                       file.type = "xml")

GDCdownload(clin_query)

clinic <- GDCprepare_clinic(query = clin_query,
                            clinical.info = "patient")

# filter out the blank values in the stage_event_pathologic_stage column
clinic_mask <- ifelse(clinic$stage_event_pathologic_stage == "", FALSE, TRUE)
clinic <- clinic[clinic_mask, ]

# Generalize each substage (A/B/C) into a single numerical stage
clinic$stage_event_pathologic_stage <- as.character(clinic$stage_event_pathologic_stage)

StageI_mask <- ifelse(clinic$stage_event_pathologic_stage == "Stage IA" | clinic$stage_event_pathologic_stage == "Stage IB", TRUE, FALSE)
clinic$stage_event_pathologic_stage[StageI_mask] <- "Stage I"

StageII_mask <- ifelse(clinic$stage_event_pathologic_stage == "Stage IIA" | clinic$stage_event_pathologic_stage == "Stage IIB", TRUE, FALSE)
clinic$stage_event_pathologic_stage[StageII_mask] <- "Stage II"

StageIII_mask <- ifelse(clinic$stage_event_pathologic_stage == "Stage IIIA" | clinic$stage_event_pathologic_stage == "Stage IIIB", TRUE, FALSE)
clinic$stage_event_pathologic_stage[StageIII_mask] <- "Stage III"

clinic$stage_event_pathologic_stage <- as.factor(clinic$stage_event_pathologic_stage)


# load in the mutation data from TCGA (need to include a Tumor_Sample_Barcode column in the clinical dataframe)
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_query <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(maf_query)

maf <- GDCprepare(maf_query) 

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinic,
                       isTCGA = TRUE)
mutations <- maf_object@data

# oncoplot of the top 10 most mutated genes
oncoplot(maf = maf_object,
         top = 10,
         borderCol = NA)
jpeg("/top10_oncoplot.jpg")
dev.off()


# co-oncoplot comparing the mutation rates of the top 10 most mutated genes between Stage I and Stage II patients 
S1_mask <- ifelse(maf_object@clinical.data$stage_event_pathologic_stage == "Stage I", T, F)

S1_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[S1_mask]

S1_maf <- subsetMaf(maf = maf_object,
                    tsb = S1_barcodes)

S2_mask <- ifelse(maf_object@clinical.data$stage_event_pathologic_stage == "Stage II", T, F)

S2_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[S2_mask]

S2_maf <- subsetMaf(maf = maf_object,
                    tsb = S2_barcodes)

coOncoplot(m1 = S1_maf, 
           m2 = S2_maf, 
           m1Name = "Stage I Patients", 
           m2Name = "Stage II Patients", 
           genes = c("TP53", "TTN", "MUC16", "CSMD3", "RYR2", "LRP1B", "ZFHX4", "USH2A", "KRAS", "XIRP2"),
           borderCol = NA)

jpeg("/StageI_II_cooncoplot.jpg")
dev.off()

# co-oncoplot comparing the mutation rates of the top 10 most mutated genes between Stage II and Stage III patients
S3_mask <- ifelse(maf_object@clinical.data$stage_event_pathologic_stage == "Stage III", T, F)

S3_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[S3_mask]

S3_maf <- subsetMaf(maf = maf_object,
                    tsb = S3_barcodes)

coOncoplot(m1 = S2_maf, 
           m2 = S3_maf, 
           m1Name = "Stage II Patients", 
           m2Name = "Stage III Patients", 
           genes = c("TP53", "TTN", "MUC16", "CSMD3", "RYR2", "LRP1B", "ZFHX4", "USH2A", "KRAS", "XIRP2"),
           borderCol = NA)

jpeg("/StageII_III_cooncoplot.jpg")
dev.off()

# co-oncoplot comparing the mutation rates of the top 10 most mutated genes between Stage III and Stage IV patients
S4_mask <- ifelse(maf_object@clinical.data$stage_event_pathologic_stage == "Stage IV", T, F)

S4_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[S4_mask]

S4_maf <- subsetMaf(maf = maf_object,
                    tsb = S4_barcodes)

coOncoplot(m1 = S3_maf, 
           m2 = S4_maf, 
           m1Name = "Stage III Patients", 
           m2Name = "Stage IV Patients", 
           genes = c("TP53", "TTN", "MUC16", "CSMD3", "RYR2", "LRP1B", "ZFHX4", "USH2A", "KRAS", "XIRP2"),
           borderCol = NA)

jpeg("/StageIII_IV_cooncoplot.jpg")
dev.off()

# load in transcriptomic data from TCGA
rna_query <- GDCquery(project ="TCGA-LUAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

GDCdownload(rna_query)

rna_se <- GDCprepare(rna_query)

rna_clinical <- rna_se@colData[!is.na(rna_se@colData$age_at_index), ]
rna_clinical <- as.data.frame(rna_clinical)
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
rna_counts <- rna_se@assays@data$unstranded[, !is.na(rna_se@colData$age_at_index)]
rna_counts <- as.data.frame(rna_counts)

# filter out columns with weird data in rna_clinical (in lecture slides)
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments", F, T) 
rna_clinical <- rna_clinical[, treatments_mask]
primary_site_mask <- ifelse(colnames(rna_clinical) == "primary_site", F, T)
rna_clinical <- rna_clinical[, primary_site_mask]
disease_type_mask <- ifelse(colnames(rna_clinical) == "disease_type", F, T)
rna_clinical <- rna_clinical[, disease_type_mask]

# label the rows and columns with more descriptive names 
rownames(rna_clinical) <- rna_clinical$barcode
rownames(rna_genes) <- rna_genes$gene_id
rownames(rna_counts) <- rownames(rna_genes)
colnames(rna_counts) <- rownames(rna_clinical)

# filter out the na values in the stage column, the control patients, and the very low rna counts 
na_mask <- ifelse(is.na(rna_clinical$ajcc_pathologic_stage), F, T)
rna_clinical <- rna_clinical[na_mask, ]
rna_counts <- rna_counts[ ,na_mask] 

normal_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", F, T)
rna_clinical <- rna_clinical[normal_mask, ]
rna_counts <- rna_counts[ ,normal_mask]

row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums >= 10, T, F)
rna_counts <- rna_counts[low_counts_mask, ]
rna_genes <- rna_genes[low_counts_mask, ]

# generalize each substage (A/B/C) into a single numerical stage 
rna_clinical$ajcc_pathologic_stage <- as.character(rna_clinical$ajcc_pathologic_stage)

StageI_mask <- ifelse(rna_clinical$ajcc_pathologic_stage == "Stage IA" | rna_clinical$ajcc_pathologic_stage == "Stage IB", TRUE, FALSE)
rna_clinical$ajcc_pathologic_stage[StageI_mask] <- "Stage I"

StageII_mask <- ifelse(rna_clinical$ajcc_pathologic_stage == "Stage IIA" | rna_clinical$ajcc_pathologic_stage == "Stage IIB", TRUE, FALSE)
rna_clinical$ajcc_pathologic_stage[StageII_mask] <- "Stage II"

StageIII_mask <- ifelse(rna_clinical$ajcc_pathologic_stage == "Stage IIIA" | rna_clinical$ajcc_pathologic_stage == "Stage IIIB" | rna_clinical$ajcc_pathologic_stage == "Stage IIIC", TRUE, FALSE)
rna_clinical$ajcc_pathologic_stage[StageIII_mask] <- "Stage III"

# filter out the Stage X patients
StageX_mask <- ifelse(rna_clinical$ajcc_pathologic_stage == "Stage X", FALSE, TRUE)
rna_clinical <- rna_clinical[StageX_mask, ]
rna_counts <- rna_counts[ ,StageX_mask] 

# convert the variables of interest to factors and turn them into categorical variables 
rna_clinical$ajcc_pathologic_stage <- as.factor(rna_clinical$ajcc_pathologic_stage)
rna_clinical$ajcc_pathologic_stage <- factor(rna_clinical$ajcc_pathologic_stage)

rna_clinical$age_category <- ifelse(rna_clinical$age_at_index <= 58, "young", "old")
rna_clinical$age_category <- factor(rna_clinical$age_category)
rna_clinical$gender <- factor(rna_clinical$gender)
sum(is.na(rna_clinical$age_category)) # check for NA values
sum(is.na(rna_clinical$gender))

# run DESeq to analyze RNA expression 
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~age_category + gender + ajcc_pathologic_stage)

dds_obj <- DESeq(dds)
resultsNames(dds_obj)

# differential RNA expression between Stage I and Stage II
resLFC1_2 <- lfcShrink(dds_obj, coef="ajcc_pathologic_stage_Stage.II_vs_Stage.I", type="normal") 
gene_name_mask <- ifelse(rna_genes$gene_id %in% resLFC1_2@rownames, T, F)
resLFC1_2 <- data.frame(rna_genes$gene_name[gene_name_mask], resLFC1_2@rownames, resLFC1_2$log2FoldChange, resLFC1_2$pvalue, resLFC1_2$padj, -log10(resLFC1_2$padj))
colnames(resLFC1_2) <- c("gene_name", "gene_id", "log2FoldChange", "pvalue", "padj", "-log10(padj)")

# volcano plot showing the differential RNA expression between Stage I and Stage II
EnhancedVolcano(resLFC1_2,
                lab = resLFC1_2$gene_name,
                x = "log2FoldChange",
                y = "pvalue",
                title = "Stage II vs. Stage I",
                legendLabels = c("Not sig.", "Log2FC", "pvalue", "pvalue & Log2FC"),
                legendPosition = "right")
jpeg("/StageI_II_volcano_plot.jpg")
dev.off()

# differential RNA expression between Stage I and Stage III + volcano plot
resLFC1_3 <- lfcShrink(dds_obj, coef="ajcc_pathologic_stage_Stage.III_vs_Stage.I", type="normal")
gene_name_mask <- ifelse(rna_genes$gene_id %in% resLFC1_3@rownames, T, F)
resLFC1_3 <- data.frame(rna_genes$gene_name[gene_name_mask], resLFC1_3@rownames, resLFC1_3$log2FoldChange, resLFC1_3$pvalue, resLFC1_3$padj, -log10(resLFC1_3$padj))
colnames(resLFC1_3) <- c("gene_name", "gene_id", "log2FoldChange", "pvalue", "padj", "-log10(padj)")
EnhancedVolcano(resLFC1_3,
                lab = resLFC1_3$gene_name,
                x = "log2FoldChange",
                y = "pvalue",
                title = "Stage III vs. Stage I",
                legendLabels = c("Not sig.", "Log2FC", "pvalue", "pvalue & Log2FC"),
                legendPosition = "right")
jpeg("/StageI_III_volcano_plot.jpg")
dev.off()

# differential RNA expression between Stage I and Stage IV + volcano plot
resLFC1_4 <- lfcShrink(dds_obj, coef="ajcc_pathologic_stage_Stage.IV_vs_Stage.I", type="normal")
gene_name_mask <- ifelse(rna_genes$gene_id %in% resLFC1_4@rownames, T, F)
resLFC1_4 <- data.frame(rna_genes$gene_name[gene_name_mask], resLFC1_4@rownames, resLFC1_4$log2FoldChange, resLFC1_4$pvalue, resLFC1_4$padj, -log10(resLFC1_4$padj))
colnames(resLFC1_4) <- c("gene_name", "gene_id", "log2FoldChange", "pvalue", "padj", "-log10(padj)")
EnhancedVolcano(resLFC1_4,
                lab = resLFC1_4$gene_name,
                x = "log2FoldChange",
                y = "pvalue",
                title = "Stage IV vs. Stage I",
                legendLabels = c("Not sig.", "Log2FC", "pvalue", "pvalue & Log2FC"),
                legendPosition = "right")
jpeg("/StageI_IV_volcano_plot.jpg")
dev.off()

# filter out only the gene expressions with a log2foldchange value of >= 2 or <= -2 for each results dataframe
log_mask <- ifelse(resLFC1_2$log2FoldChange >= 2 | resLFC1_2$log2FoldChange <= -2, TRUE, FALSE)
resLFC1_2_sig <- resLFC1_2[log_mask, ]
log_mask <- ifelse(resLFC1_3$log2FoldChange >= 2 | resLFC1_3$log2FoldChange <= -2, TRUE, FALSE)
resLFC1_3_sig <- resLFC1_3[log_mask, ]
log_mask <- ifelse(resLFC1_4$log2FoldChange >= 2 | resLFC1_4$log2FoldChange <= -2, TRUE, FALSE)
resLFC1_4_sig <- resLFC1_4[log_mask, ]

# combine all the significant RNA into one dataframe 
gene_name <- c(resLFC1_2_sig$gene_name, resLFC1_3_sig$gene_name, resLFC1_4_sig$gene_name)
gene_id <- c(resLFC1_2_sig$gene_id, resLFC1_3_sig$gene_id, resLFC1_4_sig$gene_id)
log2FoldChange <- c(resLFC1_2_sig$log2FoldChange, resLFC1_3_sig$log2FoldChange, resLFC1_4_sig$log2FoldChange)
sig_RNA <- data.frame(gene_name, gene_id, log2FoldChange)

# export the dataframes to be used in the final machine learning dataframe
write.csv(sig_RNA, "/sig_RNA.csv")
write.csv(rna_clinical, "/rna_clinical.csv")
write.csv(rna_genes, "/rna_genes.csv")
write.csv(rna_counts, "/rna_counts.csv")
write.csv(mutations, "/mutations.csv")