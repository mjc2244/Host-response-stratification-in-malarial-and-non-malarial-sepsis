#Clear R environment
rm(list = ls())

#Import RESERVE-U-2-TOR dataset 
reserve_tor <- read.csv(file.choose(), header=TRUE)

#PIDs to fixed row identifiers
rownames(reserve_tor) = reserve_tor$pid
reserve_tor$pid = NULL

#Omit those with unknown malaria RDT status
reserve_tor <- reserve_tor[!is.na(reserve_tor$malaria_rdt_result),]

#Restrict to patients with qSOFA score >=1
reserve_tor <- subset(reserve_tor, qsofa_1p_avpu==1)

#Set malaria and sex as factors
reserve_tor$gender <- as.factor(reserve_tor$gender)
reserve_tor$malaria_rdt_result <- factor(reserve_tor$malaria_rdt_result, levels=c('1', '0'))

#Table - Demographic and clinical characteristics by malaria status 
library(gtsummary)
reserve_tor_u_2_table_df <- reserve_tor[, c(  "age",
                                                  "gender",
                                                  "muac_imputed",
                                                  "illnessduration_to_enroll",
                                                  "abx_antimalarial_prior",
                                                  "admit_perf_16_and_above",
                                                  "assist_oob_walk",
                                                  "lactate_result_imputed",
                                                  "lactate4plus",
                                                  "qsofa_2p_avpu", 
                                                  "mews_score",
                                                  "uva_score_avpu", 
                                                  "hiv_rdt_result", 
                                                  "hiv_stage_3_4_all", 
                                                  "cd4_imputed",
                                                  "suppressedviralload_imputed",
                                                  "malaria_rdt_result", 
                                                  "microtbdx", 
                                                  "cryto_ag_result",
                                                  "flu_pcr_uvri_result", 
                                                  "covid_pcr_uvri_result",
                                                  "ivf_hosp",
                                                  "antimalarials_hosp",
                                                  "antimalarials_hosp_med___1",
                                                  "antimalarials_hosp_med___2",
                                                  "antimalarials_hosp_med___3",
                                                  "antimalarials_hosp_med___4",
                                                  "abx_hosp",
                                                  "anti_tb_hosp",                     
                                                  "hosp_outcome", 
                                                  "death30d",
                                                  "death60d")]

theme_gtsummary_compact()

reserve_tor_u_2_table_df %>%
  tbl_summary(
    by = malaria_rdt_result,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"
    ),
    type = list(age ~ "continuous",
                gender ~ "categorical",
                muac_imputed ~ "continuous",
                illnessduration_to_enroll ~ "continuous",
                abx_antimalarial_prior ~ "categorical",
                admit_perf_16_and_above ~ "continuous",
                assist_oob_walk ~ "dichotomous",
                lactate_result_imputed~ "continuous",
                lactate4plus ~ "dichotomous",
                qsofa_2p_avpu ~ "dichotomous", 
                mews_score ~ "continuous",
                uva_score_avpu ~ "continuous", 
                hiv_rdt_result ~ "dichotomous",
                hiv_stage_3_4_all ~ "dichotomous", 
                cd4_imputed ~ "continuous", 
                suppressedviralload_imputed ~ "dichotomous",
                microtbdx ~ "dichotomous",
                cryto_ag_result ~ "dichotomous",
                flu_pcr_uvri_result ~ "dichotomous",
                covid_pcr_uvri_result ~ "dichotomous",
                ivf_hosp ~ "dichotomous",
                antimalarials_hosp ~ "dichotomous",
                antimalarials_hosp_med___1 ~ "dichotomous",
                antimalarials_hosp_med___2 ~ "dichotomous",
                antimalarials_hosp_med___3 ~ "dichotomous",
                antimalarials_hosp_med___4 ~ "dichotomous",
                abx_hosp ~ "dichotomous",
                anti_tb_hosp ~ "dichotomous",                    
                hosp_outcome ~ "categorical",
                death30d ~ "dichotomous",
                death60d ~ "dichotomous"),
    digits = list(all_continuous() ~ 0, all_dichotomous() ~ c(0,0,1), all_categorical() ~ c(0,0,1), lactate_result_imputed ~ c(1, 1)),
    missing_text = "(Missing)",
    label = list(age ~ "Age",
                 gender ~ "Sex",
                 muac_imputed ~ "Mid-upper arm circumference, cm",
                 illnessduration_to_enroll ~ "Illness duration prior to enrollment, days",
                 abx_antimalarial_prior ~ "Received antimalarial treatment prior to hospitalization",
                 admit_perf_16_and_above ~ "KPS at hospital admission",
                 assist_oob_walk ~ "Unable to ambulate without assistance",
                 lactate_result_imputed~ "Whole-blood lactate, mmol/L",
                 lactate4plus ~ "Whole-blood lactate \u2265 4 mmol/L",
                 qsofa_2p_avpu ~ "qSOFA \u2265 2", 
                 mews_score ~ "Modified Early Warning Score",
                 uva_score_avpu ~ "Universal Vital Assessment score", 
                 hiv_rdt_result ~ "Living with HIV",
                 hiv_stage_3_4_all ~ "WHO HIV clinical stage 3 or 4", 
                 cd4_imputed ~ "CD4 count, cells/mm3",
                 suppressedviralload_imputed ~ "HIV-1 viral suppression",
                 microtbdx ~ "Microbiological TB positive",
                 cryto_ag_result ~ "Cryptococcal ag. positive if living with HIV",
                 flu_pcr_uvri_result ~ "Influenza PCR positive",
                 covid_pcr_uvri_result ~ "SARS-CoV-2 PCR positive",
                 ivf_hosp ~ "Received intravenous fluids",
                 antimalarials_hosp ~ "Received anti-malarial agent(s)",
                 antimalarials_hosp_med___1 ~ "Received quinine",
                 antimalarials_hosp_med___2 ~ "Received artesunate",
                 antimalarials_hosp_med___3 ~ "Received artemether-lumefantrine",
                 antimalarials_hosp_med___4 ~ "Received dihydroartemisinin-piperaquine",
                 abx_hosp ~ "Received anti-bacterial agent(s)",
                 anti_tb_hosp ~ "Received anti-TB agent(s)",                    
                 hosp_outcome ~ "Hospital outcome",
                 death30d ~ "Death at 30-days",
                 death60d ~ "Death at 60-days")) %>% add_n() %>% add_overall() %>% add_p(pvalue_fun = ~style_sigfig(., digits = 3)) %>% as_gt() 

