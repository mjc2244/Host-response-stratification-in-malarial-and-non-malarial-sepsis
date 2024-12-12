#Analysis is clinical comparisons of sepsis patients with and without malaria (RDT pos vs. neg)

#Clear R environment
rm(list = ls())

#Import the RESERVE-U-1-EBB dataset for all patients 
reserve <- read.csv(file.choose(), header=TRUE)

#PIDs to fixed row identifiers
rownames(reserve) = reserve$pid
reserve$pid = NULL

#Omit those with unknown malaria RDT status
reserve <- reserve[!is.na(reserve$malariardtresult),]

#Drop PID 361 - extreme outlier in protein expression
reserve <- reserve[row.names(reserve) != "361",]

#Restrict to patients with qSOFA>=1
reserve <- subset(reserve, qsofa1p==1)

#Set malaria and sex as factors
reserve$sex <- as.factor(reserve$sex)
reserve$malariardtresult <- as.factor(reserve$malariardtresult)

#Differences in demographic and clinical characteristics across malaria RDT status 
library(gmodels)
library(dplyr)
library(skimr)
library(gtsummary)

#Re-order malaria RDT result 
reserve$malariardtresult = factor(reserve$malariardtresult, levels=c("1", "0"))

#Create dataframe 
malariaclinical_df <- reserve[, c("sex", 
                                  "age",
                                  "illnessduration_enroll",
                                  "abxmalarialprior", 
                                  "qsofa2p", 
                                  "MEWS_score",
                                  "UVA_score", 
                                  "hivrdtresult", 
                                  "hivstage34", 
                                  "malariardtresult", 
                                  "microtbdx", 
                                  "influenzapcrresult", 
                                  "ivfluid",
                                  "antimalarials",
                                  "quinine",
                                  "artesunate",
                                  "coartem",
                                  "duocortexin",
                                  "antibiotics",
                                  "antitbdrugs",
                                  "hospoutcome",
                                  "death30d")]

theme_gtsummary_compact()

malariaclinical_df %>%
  tbl_summary(
    by = malariardtresult,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"
    ),
    type = list(sex ~ "dichotomous", 
                age ~ "continuous",
                illnessduration_enroll ~ "continuous",
                abxmalarialprior ~ "dichotomous", 
                qsofa2p ~ "dichotomous", 
                MEWS_score ~ "continuous",
                UVA_score ~ "continuous", 
                hivrdtresult ~ "dichotomous", 
                hivstage34 ~ "dichotomous", 
                microtbdx ~ "dichotomous", 
                influenzapcrresult ~ "dichotomous", 
                ivfluid ~ "dichotomous",
                antimalarials ~ "dichotomous",
                quinine ~ "dichotomous",
                artesunate ~ "dichotomous",
                coartem ~ "dichotomous",
                duocortexin ~ "dichotomous",
                antibiotics ~ "dichotomous",
                antitbdrugs ~ "dichotomous",
                hospoutcome ~ "categorical",
                death30d ~ "dichotomous"),
    digits = list(all_continuous() ~ 0, all_dichotomous() ~ c(0,0,0), all_categorical() ~ c(0,0,0)),
    missing_text = "(Missing)",
    label = list(sex ~ "Female sex", 
                 age ~ "Age, years",
                 illnessduration_enroll ~ "Duration of illness prior to enrollment, days", 
                 abxmalarialprior ~ "Reported anti-malarial treatment prior to hospitalization", 
                 qsofa2p ~ "qSOFA score ≥2", 
                 MEWS_score ~ "Modified Early Warning Score",
                 UVA_score ~ "Universla Vital Assessment score", 
                 hivrdtresult ~ "Living with HIV", 
                 hivstage34 ~ "HIV clinical stage 3 or 4", 
                 microtbdx ~ "Microbiological TB positive", 
                 influenzapcrresult ~ "Influenza PCR positive", 
                 ivfluid ~ "Received intravenous fluid",
                 antimalarials ~ "Received antimalarial agent(s)",
                 quinine ~ "Received quinine",
                 artesunate ~ "Received artesunate",
                 coartem ~ "Received artemether-lumefantrine",
                 duocortexin ~ "Received dihydroartemisinin–piperaquine",
                 antibiotics ~ "Received antibacterial agent(s)",
                 antitbdrugs ~ "Received anti-TB agent(s)",
                 hospoutcome ~ "Hospital outcome",
                 death30d ~ "Death at 30-days")) %>% add_n() %>% add_overall() %>% add_p(pvalue_fun = ~style_sigfig(., digits = 3)) %>% as_gt() 



