#Analysis is PCA of measured proteome by malaria (RDT pos vs. neg) status 

#Clear R environment
rm(list = ls())

#Import the RESERVE-U-1-EBB dataset for all patients 
reserve <- read.csv(file.choose(), header=TRUE)

#PIDs to fixed row identifiers
rownames(reserve) = reserve$pid
reserve$pid = NULL

#Omit those with unknown malaria RDT status
reserve <- reserve[!is.na(reserve$malariardtresult),]

#Restrict to patients with qSOFA>=1
reserve <- subset(reserve, qsofa1p==1)

#Drop PIDs 361 - extreme outlier in protein expression
reserve <- reserve[row.names(reserve) != "361",]

#Select and view the log2-transformed Olink data and covariables 
biomarkers <- reserve[1:238, c(which(colnames(reserve) == 'age'),  which(colnames(reserve) == 'sex'), which(colnames(reserve) == 'illnessduration_enroll'), which(colnames(reserve) == 'hivrdtresult'), which(colnames(reserve) == 'microtbdx'),
                               156, 225:408)]
names(biomarkers)

#Drop biomarkers with <20% of NPX values above each panel’s estimated limits of detection in either of the RESERVE-U cohorts 
library(dplyr)
biomarkers <- dplyr::select(biomarkers, -c("il1alpha",
                                           "il2", 
                                           "il33", 
                                           "il4", 
                                           "il13",
                                           "prcp",
                                           "ltbp2",
                                           "sod1",
                                           "itgam",
                                           "fap",
                                           "mfap5"))  

#Reformat column names
colnames(biomarkers) <- base::toupper(colnames(biomarkers))
colnames(biomarkers)[colnames(biomarkers) == "IL8"] = "IL-8"
colnames(biomarkers)[colnames(biomarkers) == "IL2"] = "IL-2"
colnames(biomarkers)[colnames(biomarkers) == "IL7"] = "IL-7"
colnames(biomarkers)[colnames(biomarkers) == "IL6"] = "IL-6"
colnames(biomarkers)[colnames(biomarkers) == "IL18"] = "IL-18"
colnames(biomarkers)[colnames(biomarkers) == "IL15"] = "IL-15"
colnames(biomarkers)[colnames(biomarkers) == "IL5"] = "IL-5"
colnames(biomarkers)[colnames(biomarkers) == "IL10"] = "IL-10"
colnames(biomarkers)[colnames(biomarkers) == "IFNGAMMA"] = "IFN-\u03B3"
colnames(biomarkers)[colnames(biomarkers) == "IL12RB1"] = "IL-12RB1"
colnames(biomarkers)[colnames(biomarkers) == "IL12"] = "IL-12"
colnames(biomarkers)[colnames(biomarkers) == "IL7R"] = "IL-7R"
colnames(biomarkers)[colnames(biomarkers) == "PDGFSUBUNITB"] = "PDGFB"
colnames(biomarkers)[colnames(biomarkers) == "LAPTGFBETA1"] = "LAPTGFB1"
colnames(biomarkers)[colnames(biomarkers) == "PDL1"] = "PD-L1"
colnames(biomarkers)[colnames(biomarkers) == "PDL2"] = "PD-L2"

#First 2 PCs stratified by malaria status with density plot
library("mixOmics")
library("PLSDAbatch")

pcacols <- c("#d95f02","#7570b3")

#Recode malarial sepsis variable 
biomarkers[biomarkers$MALARIARDTRESULT == 0, "etiology"] <- "Non-malarial"
biomarkers[biomarkers$MALARIARDTRESULT == 1, "etiology"] <- "Malarial"

pcav2 <- pca(biomarkers[1:238, c(7:179)], ncomp = 3, scale = TRUE)
pcadensity <- Scatter_Density(object = pcav2, 
                              batch = biomarkers$etiology,
                              trt = NULL,
                              color.set = pcacols,
                              batch.legend.title = "Sepsis etiology",
                              density.lwd = 0.1,
                              title = NULL,
                              title.cex = 1.5,
                              legend.cex = 1.5,
                              legend.title.cex = 1.5)
pcadensity < pcadensity + theme(text = element_text(size = 18),
                                axis.title = element_text(size = 18),
                                axis.text = element_text(size = 18))
pcadensity