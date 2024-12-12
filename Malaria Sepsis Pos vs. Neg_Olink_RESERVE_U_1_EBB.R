#Analysis is proteomic comparisons of sepsis patients with and without malaria (RDT pos vs. neg)

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

#Differential protein expression
#Select and view the log2-transformed Olink data and covariables
biomarkers <- reserve[1:238, c(which(colnames(reserve) == 'age'),  which(colnames(reserve) == 'sex'), which(colnames(reserve) == 'illnessduration_enroll'), which(colnames(reserve) == 'hivrdtresult'), which(colnames(reserve) == 'microtbdx'),
                               156, 225:408)]
names(biomarkers)

#Omit those with unknown malaria RDT status, illness duration, HIV status, TB status
biomarkers <- biomarkers[!is.na(biomarkers$malariardtresult),]
biomarkers <- biomarkers[!is.na(biomarkers$illnessduration_enroll),]
biomarkers <- biomarkers[!is.na(biomarkers$hivrdtresult),]
biomarkers <- biomarkers[!is.na(biomarkers$microtbdx),]

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

biomarkers[biomarkers$MALARIARDTRESULT==1, "MALARIARDTRESULT_2"] <- "Malaria"
biomarkers[biomarkers$MALARIARDTRESULT==0, "MALARIARDTRESULT_2"] <- "NoMalaria"
biomarkers$MALARIARDTRESULT_2 <- as.factor(biomarkers$MALARIARDTRESULT_2)
biomarkers$SEX <- as.factor(biomarkers$SEX)
biomarkers$HIVRDTRESULT <- as.factor(biomarkers$HIVRDTRESULT)
biomarkers$MICROTBDX <- as.factor(biomarkers$MICROTBDX)

data = biomarkers

#Remove malaria columns and transpose the data so that proteins are rows and samples are columns
expression_data <- biomarkers[1:235, c(7:179)] 
names(expression_data)
expression_data <- t(expression_data)  

group_assignments = factor(biomarkers$MALARIARDTRESULT_2, levels=c("NoMalaria", "Malaria"))
group_assignments

sex_covariate = factor(biomarkers$SEX)
age_covariate = biomarkers$AGE
illnessduration_covariate = biomarkers$ILLNESSDURATION_ENROLL
hiv_covariate = factor(biomarkers$HIVRDTRESULT)
tb_covariate = factor(biomarkers$MICROTBDX)

# Create the design matrix for the linear model
design <- model.matrix(~ 0 + group_assignments +  age_covariate + sex_covariate + illnessduration_covariate + hiv_covariate + tb_covariate)
colnames(design)[1:2] <- levels(group_assignments)
colnames(design)
# run limma
library(limma)
fit <- lmFit(expression_data, design)

# Apply empirical Bayes statistics
fit <- eBayes(fit)

# Create a contrast matrix for the comparison of malaria positive vs malaria negative  
contrast.matrix <- makeContrasts(MalariavsNoMalaria = Malaria - NoMalaria, levels = design)

# Fit the contrast model
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get the table of top differentially expressed proteins
results <- topTable(fit2, coef="MalariavsNoMalaria", adjust.method="BH", sort.by="P", number=Inf)

#Volcano plot
library(EnhancedVolcano)
# Create a data frame for the volcano plot
results$Prot <- toupper(rownames(results))  # Create a column for gene names if not already present
names(results$logFC) <- toupper(results$Prot)
names(results$adj.P.Val) <- toupper(results$Prot)
# 

keyvals <- ifelse(
  results$adj.P.Val < 0.05, 'royalblue',
  ifelse(results$adj.P.Val > 0.05, 'grey',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'royalblue'] <- 'S'

results$selectedLab <- ifelse(results$adj.P.Val < 0.05, as.character(results$Prot), NA)
EnhancedVolcano(
  results,
  lab = NA,
  x = 'logFC',
  y = 'adj.P.Val',
  xlim = c(-1.5,1.5),
  ylim = c(0,3.5),
  title = paste('Volcano Plot of Malarial vs. Non-malarial sepsis', sep = " "),
  pCutoff = 0.05,
  colCustom = keyvals,
  colAlpha = 1,
  cutoffLineType = 'blank',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.8,
  pointSize = 6,
  labSize = 10,
  legendLabSize = 25.0,
  axisLabSize = 25.0,
  max.overlaps = 150,
  drawConnectors = TRUE
)+ geom_text_repel(
  data = subset(results, !is.na(selectedLab)),  # Only use rows where labels are not NA
  aes(label = selectedLab, x = logFC, y = -log10(adj.P.Val)),  # Ensure 'label' aesthetic is correctly mapped
  size = 8,  # Ensure size matches labSize in EnhancedVolcano
  box.padding = unit(0.3, "lines"),  # Adjust padding to avoid overlapping
  point.padding = unit(0.2, "lines"),  # Increase point padding
  max.overlaps = 150,  # Allow more overlaps
  #  segment.color = NA
)

#Discrimination and calibration of protein signature for malarial sepsis
#Model
proteinsignature <- glm(malariardtresult ~ timd4 + 
                           il10 + 
                           lilrb1 + 
                           kir3dl1 +
                          lag3,
                         data = reserve, family = "binomial")
summary(proteinsignature)
exp(proteinsignature$coefficients)
exp(confint.default(proteinsignature))

#ROC
library(pROC)

predicted <- predict(proteinsignature, reserve, type="response")
auc(reserve$malariardtresult, predicted)


rocobj <- plot.roc(reserve$malariardtresult, proteinsignature$fitted.values,
                   percent=FALSE,
                   ci = TRUE,                  # compute AUC (of AUC by default)
                   print.auc = TRUE,
                   print.auc.x = 0.4, 
                   print.auc.y = 0.05,
                   legacy.axes = TRUE,
                   print.auc.pattern = "AUROC %.2f (%.2f-%.2f)")           # print the AUC (will contain the CI)
ciobj <- ci.se(rocobj, boot.n=10000,                         # CI of sensitivity
               specificities = seq(0, 1, 0.05)) # over a select set of specificities
plot(ciobj, type = "shape", col = "#8da0cb")     # plot as a blue shape

#Calibration
library(CalibrationCurves)
reserve$pred_malaria <- predict(proteinsignature, type = 'response')

cal_malaria = val.prob.ci.2(reserve$pred_malaria, reserve$malariardtresult, 
                          smooth = "loess", 
                          connect.smooth = FALSE, 
                          CL.smooth = "fill",
                          cex=1,
                          cex.leg=1,
                          dostats = c("Brier", "Eavg"),
                          cl.level = 0.95,
                          col.ideal = "#0000FF",
                          lwd.ideal = 2)

#Compute probabilities 
library(caret)

#Set CV parameters and seed
train.control <- trainControl(method="repeatedcv", number=10, repeats = 100, 
                              savePredictions = "all", classProbs = TRUE)                

set.seed(12345)

#Recode malaria status to factor
reserve[reserve$malariardtresult==1, "malariardtresult_factor"] <- "yes"
reserve[reserve$malariardtresult==0, "malariardtresult_factor"] <- "no"
reserve$malariardtresult_factor <- as.factor(reserve$malariardtresult_factor)
levels(reserve$malariardtresult_factor)
reserve$malariardtresult_factor <- relevel(reserve$malariardtresult_factor, "yes")
levels(reserve$malariardtresult_factor)

#Logistic model
malaria_model <- caret::train(malariardtresult_factor ~ timd4 + 
                                il10 + 
                                lilrb1 + 
                                kir3dl1 +
                                lag3,
                              data = reserve, method = "glm",
                              family = "binomial", trControl = train.control)

#Extract probabilities 
malaria_model_predict <- predict(malaria_model, type = "prob")

#Merge probabilities 
malaria_probs_reserve <- merge(reserve,malaria_model_predict,all=T,by='row.names')
rownames(malaria_probs_reserve) = malaria_probs_reserve$Row.names
malaria_probs_reserve$Row.names = NULL

library(cutpointr)

cp <- cutpointr(malaria_probs_reserve, yes, malariardtresult_factor, 
                method = maximize_metric, metric = accuracy)
summary(cp)
cp_plot <- plot_sensitivity_specificity(cp) + theme_bw()
cp_plot <- cp_plot + theme(plot.title = element_blank()) + theme(legend.title=element_blank()) + theme(legend.position="top") + theme(text=element_text(size=16))
cp_plot <- cp_plot + scale_color_manual(values=c("#CC6666", "#9999CC")) + geom_line(size = 1)
cp_plot

#Compare expression of proteins in classifier model by malaria status
library(ggpubr)
library(ggprism)
library(rstatix)
#IL-10
il10_plot <- ggstripchart(reserve, x = "malariardtresult_factor", y = "il10", 
                           title = "IL-10",
                           color = "malariardtresult_factor", fill = "malariardtresult_factor", alpha = 0.5, size = 4, palette =  c("#374e55","#df8f44"),
                           order = c("yes", "no"), repel = TRUE) 
il10_plot  <- il10_plot   + stat_summary(fun.y = median, 
                                           fun.ymin = median, 
                                           fun.ymax = median, 
                                           geom = "crossbar", 
                                           width = 0.5)
il10_plot  <- il10_plot  + theme(legend.position = "none")
il10_plot <- il10_plot + theme(plot.title = element_text(size = 12)) + labs(x ="", y = "Protein expression (log2 NPX)") + scale_x_discrete(labels=c("yes" = "Malaria positive", "no" = "Malaria negative"))
il10_plot <- il10_plot + theme_prism(base_size = 16) + theme(legend.position = "none")
il10_plot

il10_df_p_val <- rstatix::wilcox_test(reserve, il10 ~ malariardtresult_factor) %>% 
  rstatix::add_xy_position()

il10_plot <- il10_plot + add_pvalue(il10_df_p_val, 
                                      label = "p = {p}",
                                      label.size = 5,
                                      fontface = "bold",
                                      remove.bracket = TRUE, x = 1.5)
il10_plot

#TIMD4
timd4_plot <- ggstripchart(reserve, x = "malariardtresult_factor", y = "timd4", 
                          title = "TIMD4",
                          color = "malariardtresult_factor", fill = "malariardtresult_factor", alpha = 0.5, size = 4, palette =  c("#374e55","#df8f44"),
                          order = c("yes", "no"), repel = TRUE) 
timd4_plot  <- timd4_plot   + stat_summary(fun.y = median, 
                                         fun.ymin = median, 
                                         fun.ymax = median, 
                                         geom = "crossbar", 
                                         width = 0.5)
timd4_plot  <- timd4_plot  + theme(legend.position = "none")
timd4_plot <- timd4_plot + theme(plot.title = element_text(size = 12)) + labs(x ="", y = "Protein expression (log2 NPX)") + scale_x_discrete(labels=c("yes" = "Malaria positive", "no" = "Malaria negative"))
timd4_plot <- timd4_plot + theme_prism(base_size = 16) + theme(legend.position = "none")
timd4_plot

timd4_df_p_val <- rstatix::wilcox_test(reserve, timd4 ~ malariardtresult_factor) %>% 
  rstatix::add_xy_position()

timd4_plot <- timd4_plot + add_pvalue(timd4_df_p_val, 
              label = "p = {p}",
              label.size = 5,
              fontface = "bold",
              remove.bracket = TRUE, x = 1.5)
timd4_plot

#LILRB1
lilrb1_plot <- ggstripchart(reserve, x = "malariardtresult_factor", y = "lilrb1", 
                             title = "LILRB1",
                             color = "malariardtresult_factor", fill = "malariardtresult_factor", alpha = 0.5, size = 4, palette =  c("#374e55","#df8f44"),
                             order = c("yes", "no"), repel = TRUE) 
lilrb1_plot  <- lilrb1_plot   + stat_summary(fun.y = median, 
                                               fun.ymin = median, 
                                               fun.ymax = median, 
                                               geom = "crossbar", 
                                               width = 0.5)
lilrb1_plot  <- lilrb1_plot  + theme(legend.position = "none")
lilrb1_plot <- lilrb1_plot + theme(plot.title = element_text(size = 12)) + labs(x ="", y = "Protein expression (log2 NPX)") + scale_x_discrete(labels=c("yes" = "Malaria positive", "no" = "Malaria negative"))
lilrb1_plot <- lilrb1_plot + theme_prism(base_size = 16) + theme(legend.position = "none")
lilrb1_plot

lilrb1_df_p_val <- rstatix::wilcox_test(reserve, lilrb1 ~ malariardtresult_factor) %>% 
  rstatix::add_xy_position()

lilrb1_plot <- lilrb1_plot + add_pvalue(lilrb1_df_p_val, 
                                          label = "p = {p}",
                                          label.size = 5,
                                          fontface = "bold",
                                        remove.bracket = TRUE, x = 1.5)
lilrb1_plot

#KIR3DL1
kir3dl1_plot <- ggstripchart(reserve, x = "malariardtresult_factor", y = "kir3dl1", 
                            title = "KIR3DL1",
                            color = "malariardtresult_factor", fill = "malariardtresult_factor", alpha = 0.5, size = 4, palette =  c("#374e55","#df8f44"),
                            order = c("yes", "no"), repel = TRUE) 
kir3dl1_plot  <- kir3dl1_plot   + stat_summary(fun.y = median, 
                                             fun.ymin = median, 
                                             fun.ymax = median, 
                                             geom = "crossbar", 
                                             width = 0.5)
kir3dl1_plot  <- kir3dl1_plot  + theme(legend.position = "none")
kir3dl1_plot <- kir3dl1_plot + theme(plot.title = element_text(size = 12)) + labs(x ="", y = "Protein expression (log2 NPX)") + scale_x_discrete(labels=c("yes" = "Malaria positive", "no" = "Malaria negative"))
kir3dl1_plot <- kir3dl1_plot + theme_prism(base_size = 16) + theme(legend.position = "none")
kir3dl1_plot

kir3dl1_df_p_val <- rstatix::wilcox_test(reserve, kir3dl1 ~ malariardtresult_factor) %>% 
  rstatix::add_xy_position()

kir3dl1_plot <- kir3dl1_plot + add_pvalue(kir3dl1_df_p_val, 
                                        label = "p = {p}",
                                        label.size = 5,
                                        fontface = "bold",
                                        remove.bracket = TRUE, x = 1.5)
kir3dl1_plot

#LAG3
lag3_plot <- ggstripchart(reserve, x = "malariardtresult_factor", y = "lag3", 
                           title = "LAG3",
                           color = "malariardtresult_factor", fill = "malariardtresult_factor", alpha = 0.5, size = 4, palette =  c("#374e55","#df8f44"),
                           order = c("yes", "no"), repel = TRUE) 
lag3_plot  <- lag3_plot   + stat_summary(fun.y = median, 
                                           fun.ymin = median, 
                                           fun.ymax = median, 
                                           geom = "crossbar", 
                                           width = 0.5)
lag3_plot  <- lag3_plot  + theme(legend.position = "none")
lag3_plot <- lag3_plot + theme(plot.title = element_text(size = 12)) + labs(x ="", y = "Protein expression (log2 NPX)") + scale_x_discrete(labels=c("yes" = "Malaria positive", "no" = "Malaria negative"))
lag3_plot <- lag3_plot + theme_prism(base_size = 16) + theme(legend.position = "none")
lag3_plot

lag3_df_p_val <- rstatix::wilcox_test(reserve, lag3 ~ malariardtresult_factor) %>% 
  rstatix::add_xy_position()

lag3_plot <- lag3_plot + add_pvalue(lag3_df_p_val, 
                                      label = "p = {p}",
                                      label.size = 5,
                                      fontface = "bold",
                                    remove.bracket = TRUE, x = 1.5)
lag3_plot

#Combine plots
library(ggpubr)
top_row = ggarrange(il10_plot, timd4_plot, lilrb1_plot, ncol = 3, labels = c("", "", ""))
bottom_row = ggarrange(NULL, kir3dl1_plot, lag3_plot, NULL, ncol = 4, labels = c("", "", ""), widths = c(1,2,2,1))
final_plot = ggarrange(top_row, bottom_row, ncol = 1)
final_plot

