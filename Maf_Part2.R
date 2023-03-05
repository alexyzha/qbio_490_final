setwd("/Users/kojiro/Desktop/qbio490/qbio_490_KainingFeng/analysis_data/")
maf_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging
and Masking")
#GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf,
                       clinicalData = clinical,
                       isTCGA = TRUE)
# Start of Data Analysis
# Only take the MAF info of patients who have radiation_therapy information
has_rad_mask <- ifelse(maf_object@clinical.data$radiation_therapy != "",
                       TRUE,
                       FALSE)
has_rad_barcode <- maf_object@clinical.data[has_rad_mask,Tumor_Sample_Barcode]
maf_has_rad <- subsetMaf(maf = maf_object,
                         tsb = has_rad_barcode)
maf_has_rad@clinical.data$radiation_therapy <-
  as.factor(maf_has_rad@clinical.data$radiation_therapy)
# Split the has radiation_therapy info group into 2 categories
rad_yes_mask <- ifelse(maf_has_rad@clinical.data$radiation_therapy == "YES",
                       TRUE,
                       FALSE)
rad_yes_barcode <- maf_has_rad@clinical.data[rad_yes_mask, Tumor_Sample_Barcode]
maf_yes_rad <- subsetMaf(maf = maf_has_rad, tsb = rad_yes_barcode)
rad_no_mask <- ifelse(maf_has_rad@clinical.data$radiation_therapy == "NO",
                      TRUE,
                      FALSE)
rad_no_barcode <- maf_has_rad@clinical.data[rad_no_mask, Tumor_Sample_Barcode]
maf_no_rad <- subsetMaf(maf = maf_has_rad, tsb = rad_no_barcode)
# Determine top 10 most mutated genes in the group
oncoplot(maf = maf_yes_rad, top = 10)
# Creating Oncoplot
coOncoplot(m1 = maf_yes_rad, 
           m2 = maf_no_rad,
           genes = c("PIK3CA", "TP53", "CDH1", "TTN", "RYR2", "KMT2C", "MAP3K1", 
                     "TBX3", "ANKFN1", "SYNE1"),
           m1Name = "Patients that underwent radiation therapy", 
           m2Name = "Patients that did not undergo radiation therapy", 
           borderCol = NA)
# The PIK3CA gene provides instructions for making the p110 alpha (p110α) protein, 
#which is one piece (subunit) of an enzyme called phosphatidylinositol 3-kinase (PI3K).
#The discrepancy may due to the fact that the PIK3CA is an oncogene, which is easy to mutate.

#3 Create a contingency table
barcode_rad_yes <- maf_yes_rad@clinical.data$Tumor_Sample_Barcode
len_rad_yes <- length(barcode_rad_yes)
genePIK_maf <- subsetMaf(maf = maf_has_rad,
                         genes = "PIK3CA")
barcode_PIK_yes <- genePIK_maf@clinical.data$Tumor_Sample_Barcode
len_PIK_yes <- length(barcode_PIK_yes)
both_PIK_and_Rad <- intersect(barcode_rad_yes,barcode_PIK_yes)
len_PIK_and_Rad <- length(both_PIK_and_Rad)
len_rad_only <- len_rad_yes - len_PIK_and_Rad
len_PIK_only <- len_PIK_yes - len_PIK_and_Rad
len_no_PIK_rad <- length(maf_has_rad@clinical.data$Tumor_Sample_Barcode) - len_PIK_and_Rad - len_PIK_only - len_rad_only

contig <- matrix(c(len_PIK_and_Rad, 
                   len_rad_only,
                   len_PIK_only,
                   len_no_PIK_rad), 
                 nrow=2)
contig
mosaicplot(contig)

#Run a Fisher’s Exact Test
fisher_test <- fisher.test(contig)
fisher_test
# The odds ratio = 1.6033, which is larger than 0, meaning that the gene mutation 
# and the radiation therapy are cooccurence. However, the p-value is 0.2612, which
# is much higher than the 0.05 cut-off. This means the data cannot support the dependence 
# of the two variables. 


# co-lollipop
lollipopPlot2(m1 = maf_yes_rad, 
           m2 = maf_no_rad, 
           m1_name = "Patients that underwent radiation therapy", 
           m2_name = "Patients that did not receive radiation therapy", 
           gene = "PIK3CA")
# According to the co-lollipop graph, there's a unit called PI3K_class_1_alpha,
#where the patients didn't receive radiation therapy have a higher mutation rate.

# KM plot
maf_object@clinical.data$Overall_Survival_Status <-
  ifelse(maf_object@clinical.data$vital_status == "Alive",
         maf_object@clinical.data$Overall_Survival_Status <- TRUE,
         maf_object@clinical.data$Overall_Survival_Status <- FALSE)
mafSurvival(maf = maf_object,
            genes = "PIK3CA", ## pick a gene of your choosing
            time = "days_to_last_followup", ## name of the column in maf_object@clinical.data containing survival time
            Status = "Overall_Survival_Status", ## name of the column that contains a boolean value for death events, you may need to recreate this... 
            isTCGA = TRUE)
# The survival probability of the patient mutant gene was lower than the people with 
# wild type gene. But for the most part of the graph, there doesn't seem to have a big difference between them. 
# According to the research, the dependence between the two variables are not significant enough and 
#may need further research to prove their relationships. 

