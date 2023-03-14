if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools")
library(BiocManager)
library(TCGAbiolinks)
library(maftools)

if(!require(survival)){
  install.packages("survival")
}
library(survival)
if(!require(survminer)){
  install.packages("survminer")
}
library(survminer)
if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)

knitr::opts_knit$set(root.dir = normalizePath("/Users/kojiro/Desktop/qbio490/qbio_490_KainingFeng/analysis_data/"))
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")
clinic <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient",
                            directory = "/Users/kojiro/desktop/GDCdata")
clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug", directory = "/Users/kojiro/desktop/GDCdata")
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation",directory= "/Users/kojiro/desktop/GDCdata")

clinic$age_at_initial_pathologic_diagnosis
sum(is.na(clinic$age_at_initial_pathologic_diagnosis))
#I choose the discrete variable "the age of patients at initial pathology diagnosed"

clinical_rad$radiation_dosage
is.na(clinical_rad$radiation_dosage)
#I choose the continuous variable "the dosage of radiation therapy the patients use"

#Hypothesis 1: The younger the patients was initially diagnosed, the lower dosage the patients need for radiation therapy
#Hypothesis 2: The younger the patients was initially diagnosed, the higher survival rate the patients would have 
#Hypothesis 3: The more dosage of radiation therapy the patients use, the higher survival rate the patients would have. 

unique_rad_bcr <- clinical_rad[match(unique(clinical_rad$bcr_patient_barcode), clinical_rad$bcr_patient_barcode), ]
clinic_has_rad_info <- subset(clinic, clinic$has_radiations_information == "YES")
unique_clinic_has_rad <- clinic_has_rad_info[match(unique(clinic_has_rad_info$bcr_patient_barcode), 
                                                   clinic_has_rad_info$bcr_patient_barcode),]
plot(unique_clinic_has_rad$age_at_initial_pathologic_diagnosis, unique_rad_bcr$radiation_dosage)

#remove patients with no age info
age_na_mask <-ifelse(is.na(clinic$age_at_initial_pathologic_diagnosis),F,T)
age_cleaned_clinical <- clinic[age_na_mask,]

# determine what it means to be young, middle and old
young_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis<=50,T,F)
middle_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis>=50 
& age_cleaned_clinical$age_at_initial_pathologic_diagnosis <=70, T,F)
old_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis >=70,T,F)
age_cleaned_clinical$age_status <- ifelse(young_mask,"Young",ifelse(middle_mask,"Middle","Old"))
#making a survival time column for survival plots
age_cleaned_clinical$survival_time <- ifelse(is.na(age_cleaned_clinical$days_to_death),
                                             age_cleaned_clinical$survival_time<- age_cleaned_clinical$days_to_last_followup,
                                             age_cleaned_clinical$survival_time<- age_cleaned_clinical$days_to_death)
#remove -inf values in survival_time
inf_mask <- ifelse(age_cleaned_clinical$survival_time =="-Inf",F,T)
age_cleaned_clinical <- age_cleaned_clinical[inf_mask,]
# adding death event
age_cleaned_clinical$death_event <- ifelse(age_cleaned_clinical$vital_status=="Alive",age_cleaned_clinical$death_event <-FALSE,
                                           age_cleaned_clinical$death_event <- TRUE)
survival_object <- Surv(ifelse(age_cleaned_clinical$vital_status =="Alive", age_cleaned_clinical$days_to_last_followup,age_cleaned_clinical$days_to_death),
                        ifelse(age_cleaned_clinical$vital_status =="Alive", F,T ))
fit_object <- survfit(survival_object ~ age_cleaned_clinical$age_status, data =age_cleaned_clinical)
# Making the Plot
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme =
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend =
                         "right")
KM_plot <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))
KM_plot

# Second KM Plot
# subset out all the blank data
rad_na_mask <- ifelse(unique_rad_bcr$radiation_dosage == "", F, T)
rad_no_na <- unique_rad_bcr[rad_na_mask, ]

# Create Columns for high, middle and low dosage
quantile(as.numeric(rad_no_na$radiation_dosage))
# low dosage <= 78
# medium dosage <= 102
# high dosage > 102
low_mask <- ifelse(as.numeric(rad_no_na$radiation_dosage) <= 78, T, F)
medium_mask <- ifelse(as.numeric(rad_no_na$radiation_dosage) > 78 & 
                        as.numeric(rad_no_na$radiation_dosage) <= 102,
                      T, F)
high_mask <- ifelse(as.numeric(rad_no_na$radiation_dosage) > 102, T, F)
rad_no_na$radiation_status <- ifelse(low_mask, "Low", ifelse(medium_mask, "Medium", "High"))
# merge info from clinic to rad_no_na by barcode
rad_no_na <- merge(rad_no_na, clinic, by.x = "bcr_patient_barcode",
                   by.y = "bcr_patient_barcode")
# making survival time column
rad_no_na$survival_time <- ifelse(is.na(rad_no_na$days_to_death),
                                  rad_no_na$survival_time <- rad_no_na$days_to_last_followup,
                                  rad_no_na$survival_time <- rad_no_na$days_to_death)

# remove all -Inf in survival time
inf_mask <- ifelse(rad_no_na$survival_time == -Inf, F,T)
rad_no_na <- rad_no_na[inf_mask,]

# adding death event
rad_no_na$death_event <- ifelse(rad_no_na$vital_status == "Alive",
                                rad_no_na$death_event <- FALSE,
                                rad_no_na$death_event <- TRUE)
# Making the Plot
survival_object <- Surv(time = rad_no_na$survival_time,
                        event = rad_no_na$death_event)
fit_object <- survfit(survival_object ~ rad_no_na$radiation_status,
                      data = rad_no_na)
survplot <- ggsurvplot(fit_object, pval= TRUE, ggtheme = 
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend =
                         "right")
KM_plot_dosage <- survplot$plot + theme_bw() + theme(axis.title =
                                                       element_text(size=20), axis.text = element_text(size=16),
                                                     legend.title = element_text(size=14), legend.text =
                                                       element_text(size=12))
KM_plot_dosage

#For my first KM plot between "age being diagnosed" and "survival rate," 
#I found that the patients that being diagnosed at an old age have obviously lower survival rate comparing to 
#young and middle aged patient. The patients being diagnosed at a middle age have a slightly lower survival rate than those young age patients. 
# In this case, the research supports the hypothesis that the "The younger the patients was initially diagnosed, the higher survival rate the patients would have"
# Since the p-value is 0.0053, which is much smaller than 0.05, I can conclude that the relationship between these two factors are very significant.

#For my second KM plot between "dosage of radiation therapy" and "survival rate," 
#I found that the patients being treat by low dosage of radiation therapy have obviously lower survival rate comparing to 
#those recieve high and medium dosage. But the patients recieved medium and high dosage of radiation therapy have similar survival rate. 
# Since the p-value is 0.17, which is much larger than 0.05, I cannot conclude that the relationship between these two factors are very significant, and further 
# researches are still needed to find out their relationships.

