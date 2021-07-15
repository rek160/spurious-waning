library(survival)
library(tidyverse)

setwd("/n/home00/rkahn/VE")
args=(commandArgs(TRUE))
for (i in 1:length(args)) { 
  eval (parse (text = args[[i]] )) 
} 

j <- j
i <- j

length_detectable <- 14 # assumptions
n_match <- 4 # number of controls to match to cases 
  
names <- c("master_vaxreinf","master0.1infected","masterR02","masterR04","master0.65VES","master","masterVEP0","master_aon","master0VES1pS","master0VES") #,"master0.2pS") #change to 0.2 for 51-100

for (name in names[1]){
results_analysis_master <- read.csv(paste0(j,"_results_analysis_",name,".csv"))

analysis_master <- c()
track <- 0
for (R in c(1.25,1.5,2)){ 
  for (VE in c(0.9,0.7,0,0.65)){
      for(pr in c(0,0.1,0.2,0.3)){ #0.20
          for(ri in c(0.05)){ #0.5
            for(relVE in c(1)){ #1.1
              for (sens in c(0.8,0.9,1)){
                for (spec_case in c(0.2,0.6,1)){
                  for (t in c(25,50,75,100,125,150,175,200,225,250)){
                    #for(i in js){
                    track <- track+1
                    
                    #cat(track,"\n")
                    
                    results_analysis_master %>%
                      subset(R0==R & direct_VE==VE & prop_rec==pr & reinf==ri & rel_VE_PI==relVE & j == i) -> results_analysis
                    sum(results_analysis$eventstatus)
                    
                    if (nrow(results_analysis)>0){
                    
                    results_analysis %>%
                      group_by(InfectedNode) %>%
                      arrange(DayInfected) %>%
                      mutate(RI = row_number() - 1,
                             RI = case_when(PI==1 & eventstatus == 1 ~ 1,
                                            TRUE ~ RI)) %>%
                      ungroup() -> results_analysis # note reinfections
                    
                    sum(results_analysis$RI)
                    
                    results_analysis %>%
                      subset(eventstatus==1 & DayInfected<=t) %>%
                      nrow() -> cum_VES
                    
                    results_analysis %>%
                      subset(eventstatus==1 & DayInfected<=t & Symptomatic == 1 & !is.na(Symptomatic)) %>%
                      nrow() -> cum_VESP
                      
                    # Cohort outcome: VES
                    results_analysis %>%  # update results for day t
                      mutate(eventstatus = ifelse(DayInfected > t,0,eventstatus),
                             DayInfected = ifelse(eventstatus==0,t,DayInfected)) -> results_analysis_VES

                    # Cox model VE
                    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis_VES),silent=T)
                    usesurvmod <- !inherits(survmodel, 'try-error')
                    if (usesurvmod==TRUE){
                      vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
                    } else if (results_analysis_VES %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                                 results_analysis_VES %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                      vaccEffEst_total<-c(1,1,1)
                    } else if (results_analysis_VES %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                               results_analysis_VES %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                      vaccEffEst_total<-c(0,0,0)
                    } else {
                      vaccEffEst_total<-c(NA,NA,NA)
                    }

                    VES_cohort_b <- vaccEffEst_total[1]
                    VES_cohort_b
                    
                    # Cox model VE control for risk but don't restrict to past infection
                    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus + Risk,results_analysis_VES),silent=T)
                    usesurvmod <- !inherits(survmodel, 'try-error')
                    if (usesurvmod==TRUE){
                      vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
                    } else if (results_analysis_VES %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                               results_analysis_VES %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                      vaccEffEst_total<-c(1,1,1)
                    } else if (results_analysis_VES %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                               results_analysis_VES %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                      vaccEffEst_total<-c(0,0,0)
                    } else {
                      vaccEffEst_total<-c(NA,NA,NA)
                    }
                    
                    VES_cohort_br <- vaccEffEst_total[1]
                    VES_cohort_br
                    
                    # restrict to no past infection
                    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis_VES %>% subset(PI==0 & RI==0)),silent=T)
                    usesurvmod <- !inherits(survmodel, 'try-error')
                    if (usesurvmod==TRUE){
                      vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
                    } else if (results_analysis_VES %>% subset(PI==0 & RI==0 & TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                               results_analysis_VES %>% subset(PI==0 & RI==0 & TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                      vaccEffEst_total<-c(1,1,1)
                    } else if (results_analysis_VES %>% subset(PI==0 & RI==0 & TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                               results_analysis_VES %>% subset(PI==0 & RI==0 & TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                      vaccEffEst_total<-c(0,0,0)
                    } else {
                      vaccEffEst_total<-c(NA,NA,NA)
                    }
                    
                    VES_cohort_epi <- vaccEffEst_total[1]
                    VES_cohort_epi
                    
                    
                    # restrict to no past infection and control for risk
                    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus+Risk,results_analysis_VES %>% subset(PI==0 & RI==0)),silent=T)
                    usesurvmod <- !inherits(survmodel, 'try-error')
                    if (usesurvmod==TRUE){
                      vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
                    } else if (results_analysis_VES %>% subset(PI==0 & RI==0 & TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                               results_analysis_VES %>% subset(PI==0 & RI==0 & TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                      vaccEffEst_total<-c(1,1,1)
                    } else if (results_analysis_VES %>% subset(PI==0 & RI==0 & TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                               results_analysis_VES %>% subset(PI==0 & RI==0 & TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                      vaccEffEst_total<-c(0,0,0)
                    } else {
                      vaccEffEst_total<-c(NA,NA,NA)
                    }
                    
                    VES_cohort_epir <- vaccEffEst_total[1]
                    VES_cohort_epir
                    
                  
                    # Cohort outcome: VESP
                    results_analysis %>%  # update results for day t
                      mutate(eventstatus = ifelse(DayInfected > t | Symptomatic==0,0,eventstatus),
                             DayInfected = ifelse(eventstatus==0,t,DayInfected)) -> results_analysis_VESP # if not infected make day censored last day
                    
                    # Cox model VE
                    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis_VESP),silent=T)
                    usesurvmod <- !inherits(survmodel, 'try-error')
                    if (usesurvmod==TRUE){
                      vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
                    } else if (results_analysis_VESP %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                               results_analysis_VESP %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                      vaccEffEst_total<-c(1,1,1)
                    } else if (results_analysis_VESP %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                               results_analysis_VESP %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                      vaccEffEst_total<-c(0,0,0)
                    } else {
                      vaccEffEst_total<-c(NA,NA,NA)
                    }
                    
                    VESP_cohort_b <- vaccEffEst_total[1]
                    VESP_cohort_b
                    
                    # Cox model VE control for risk
                    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus + Risk,results_analysis_VESP),silent=T)
                    usesurvmod <- !inherits(survmodel, 'try-error')
                    if (usesurvmod==TRUE){
                      vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
                    } else if (results_analysis_VESP %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                               results_analysis_VESP %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                      vaccEffEst_total<-c(1,1,1)
                    } else if (results_analysis_VESP %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                               results_analysis_VESP %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                      vaccEffEst_total<-c(0,0,0)
                    } else {
                      vaccEffEst_total<-c(NA,NA,NA)
                    }
                    
                    VESP_cohort_br <- vaccEffEst_total[1]
                    VESP_cohort_br
                    
                    # restrict to no past infection
                    # censor people when they get infected so don't update to day t
                    results_analysis %>%  # update results for day t
                      mutate(eventstatus = ifelse(DayInfected > t | Symptomatic==0,0,eventstatus)) -> results_analysis_VESP # if not infected make day censored last day
                    
                    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus,results_analysis_VESP %>% subset(PI==0 & RI==0)),silent=T)
                    usesurvmod <- !inherits(survmodel, 'try-error')
                    if (usesurvmod==TRUE){
                      vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
                    } else if (results_analysis_VESP %>% subset(PI==0 & RI==0 & TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                               results_analysis_VESP %>% subset(PI==0 & RI==0 & TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                      vaccEffEst_total<-c(1,1,1)
                    } else if (results_analysis_VESP %>% subset(PI==0 & RI==0 & TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                               results_analysis_VESP %>% subset(PI==0 & RI==0 & TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                      vaccEffEst_total<-c(0,0,0)
                    } else {
                      vaccEffEst_total<-c(NA,NA,NA)
                    }
                    
                    VESP_cohort_epi <- vaccEffEst_total[1]
                    VESP_cohort_epi
                    
                    
                    # restrict to no past infection and control for risk
                    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus+Risk,results_analysis_VESP %>% subset(PI==0 & RI==0)),silent=T)
                    usesurvmod <- !inherits(survmodel, 'try-error')
                    if (usesurvmod==TRUE){
                      vaccEffEst_total <- 1-exp(survmodel$coefficients + c(0, 1.96, -1.96)*sqrt(survmodel$var))
                    } else if (results_analysis_VESP %>% subset(PI==0 & RI==0 & TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                               results_analysis_VESP %>% subset(PI==0 & RI==0 & TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                      vaccEffEst_total<-c(1,1,1)
                    } else if (results_analysis_VESP %>% subset(PI==0 & RI==0 & TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                               results_analysis_VESP %>% subset(PI==0 & RI==0 & TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                      vaccEffEst_total<-c(0,0,0)
                    } else {
                      vaccEffEst_total<-c(NA,NA,NA)
                    }
                    
                    VESP_cohort_epir <- vaccEffEst_total[1]
                    VESP_cohort_epir
                    
                  
                    

                    
                    
                    
                    # TND outcome: VES
                    results_analysis_infected <- results_analysis %>% # those who will test positive
                      subset(eventstatus == 1 & DayInfected <= t & DayInfected > t-length_detectable) %>%
                      mutate(test = 1)
                    
                    n_VES <- nrow(results_analysis_infected)
                    
                    nm <- ifelse(nrow(results_analysis_infected)>nrow(results_analysis)/n_match,1,n_match)
                    
                    results_analysis_notinfected <- results_analysis %>% # those who would test negative and take a random sample (1-1)
                      subset(DayInfected > t | DayInfected <= t - length_detectable) %>% 
                      sample_n(nrow(results_analysis_infected)*nm) %>%
                      mutate(test = 0)
                    
                    results_analysis_sample <- results_analysis_infected %>% bind_rows(results_analysis_notinfected)
                    
                    if (nrow(results_analysis_sample)>0){
                      
                      # check to make sure no duplciates (should be 0)
                      nrow(results_analysis_sample) - length(unique(results_analysis_sample$InfectedNode))
                      
                      # Calculate OR with everyone
                      VES_TND_b <- 1 - (results_analysis_sample %>% subset(test==1 & TrialStatus==1) %>% nrow() / (results_analysis_sample %>% subset(test==0 & TrialStatus==1) %>% nrow())) / 
                        (results_analysis_sample %>% subset(test==1 & TrialStatus==0) %>% nrow() / (results_analysis_sample %>% subset(test==0 & TrialStatus==0) %>% nrow()))
                      
                      if (is.na(VES_TND_b)){
                        if (results_analysis_sample %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                            results_analysis_sample %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                          VES_TND_b<-1
                        } else if (
                          results_analysis_sample %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                          results_analysis_sample %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                          VES_TND_b<-0
                        } else {
                          VES_TND_b<-NA
                        }
                      }
                      VES_TND_b
                      
                      
                      # baseline + risk
                      model1 <- glm(test ~ TrialStatus + Risk, data=results_analysis_sample,family=binomial)
                      VES_TND_br <- 1-exp(model1$coefficients[2])
                      
                      if (is.na(VES_TND_br)){
                        if (results_analysis_sample %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                            results_analysis_sample %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                          VES_TND_br<-1
                        } else if (
                          results_analysis_sample %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                          results_analysis_sample %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                          VES_TND_br<-0
                        } else {
                          VES_TND_br<-NA
                        }
                      }
                      VES_TND_br
                      
                      # Test for past infection 
                      
                      results_analysis_sample %>%
                        rowwise() %>%
                        mutate(Past_infection = case_when(DayInfected  <= t - length_detectable | (test==1 & RI == 1) | PI==1 ~ rbinom(1,1,sens), # imperfect sensitivity for past infections; includes those infected earlier in simulation, those reinfected on this day or those previously infected
                                                          test == 1 & RI == 0 ~ rbinom(1,1,1-spec_case), # imperfect specificity for past infections, focusing on cases only who have not been previously infected
                                                          TRUE ~ as.integer(0))) -> results_analysis_sample
                      
                      results_analysis_sample %>%
                        subset(Past_infection == 0) -> results_analysis_sample_exclude_pastinfection
                      
                      results_analysis_sample %>%
                        subset(RI == 1) -> results_analysis_sample_pastinfection
                      
                      # Exclude those with past infection and calculate OR = (T+vax+/T-vax+N-)/(T-vax-/T-vax-,N-)
                      VES_TND_epi <- 1 - (results_analysis_sample_exclude_pastinfection %>% subset(test==1 & TrialStatus==1) %>% nrow() / (results_analysis_sample_exclude_pastinfection %>% subset(test==0 & TrialStatus==1) %>% nrow())) / 
                        (results_analysis_sample_exclude_pastinfection %>% subset(test==1 & TrialStatus==0) %>% nrow() / (results_analysis_sample_exclude_pastinfection %>% subset(test==0 & TrialStatus==0) %>% nrow()))
                      
                      if (is.na(VES_TND_epi)){
                        if (results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                            results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                          VES_TND_epi<-1
                        } else if (
                          results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                          results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                          VES_TND_epi<-0
                        } else {
                          VES_TND_epi<-NA
                        }
                      }
                      VES_TND_epi
                      
                      if (nrow(results_analysis_sample_exclude_pastinfection)>0){
                          
                        # control for risk
                        model2 <- glm(test ~ TrialStatus + Risk, data=results_analysis_sample_exclude_pastinfection,family=binomial)
                        VES_TND_epir <- 1-exp(model2$coefficients[2])
                        VES_TND_epir
                      } else{
                        VES_TND_epir <- NA
                      }
                      
                      if (is.na(VES_TND_epir)){
                        if (results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                            results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                          VES_TND_epir<-1
                        } else if (
                          results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                          results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                          VES_TND_epir<-0
                        } else {
                          VES_TND_epir<-NA
                        }
                      }
                      VES_TND_epir
                      
                    
                    } else{
                      VES_TND_b <- NA
                      VES_TND_epi <- NA
                      VES_TND_epir <- NA
                      VES_TND_br <- NA
                    }
                    
                    # TND outcome: VESP
                    results_analysis_infected <- results_analysis %>% # symptomatic
                      subset(Symptomatic == 1) %>%
                      subset(eventstatus == 1 & DayInfected <= t & DayInfected > t-length_detectable) %>%
                      mutate(test = 1)
                    
                    n_VESP <- nrow(results_analysis_infected)
                    
                    nm <- ifelse(nrow(results_analysis_infected)>nrow(results_analysis %>% subset(Symptomatic == 1) )/n_match,1,n_match)
                    nm <- ifelse(nrow(results_analysis_infected)>nrow(results_analysis %>% subset(Symptomatic == 0 | is.na(Symptomatic)) )/n_match,1,n_match)
                    
                    results_analysis_notinfected <- results_analysis %>% # those who would test negative and take a random sample (1-1)
                      subset(Symptomatic == 0 | is.na(Symptomatic) | DayInfected > t | DayInfected <= t - length_detectable) %>%
                      sample_n(nrow(results_analysis_infected)*nm) %>%
                      mutate(test = 0)
                    
                    results_analysis_sample <- results_analysis_infected %>% bind_rows(results_analysis_notinfected)
                    
                    
                    if (nrow(results_analysis_sample)>0){
                      
                      sample_vax <- results_analysis_sample %>% subset(TrialStatus==1) %>% nrow()
                      sample_control <- results_analysis_sample %>% subset(TrialStatus==0) %>% nrow()
                      
                      # check to make sure no duplciates (should be 0)
                      nrow(results_analysis_sample) - length(unique(results_analysis_sample$InfectedNode))
                      
                      # Calculate OR with everyone
                      VESP_TND_b <- 1 - (results_analysis_sample %>% subset(test==1 & TrialStatus==1) %>% nrow() / (results_analysis_sample %>% subset(test==0 & TrialStatus==1) %>% nrow())) / 
                        (results_analysis_sample %>% subset(test==1 & TrialStatus==0) %>% nrow() / (results_analysis_sample %>% subset(test==0 & TrialStatus==0) %>% nrow()))
                      
                      if (is.na(VESP_TND_b)){
                        if (results_analysis_sample %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                            results_analysis_sample %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                          VESP_TND_b<-1
                        } else if (
                          results_analysis_sample %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                          results_analysis_sample %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                          VESP_TND_b<-0
                        } else {
                          VESP_TND_b<-NA
                        }
                      }
                      VESP_TND_b
                      
                      # baseline + risk
                      model1 <- glm(test ~ TrialStatus + Risk, data=results_analysis_sample,family=binomial)
                      VESP_TND_br <- 1-exp(model1$coefficients[2])
                      
                      if (is.na(VESP_TND_br)){
                        if (results_analysis_sample %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                            results_analysis_sample %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                          VESP_TND_br<-1
                        } else if (
                          results_analysis_sample %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                          results_analysis_sample %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                          VESP_TND_br<-0
                        } else {
                          VESP_TND_br<-NA
                        }
                      }
                      VESP_TND_br
                      
                      # Test for past infection 
                      results_analysis_sample %>%
                        rowwise() %>%
                        mutate(Past_infection = case_when(DayInfected  <= t - length_detectable | (test==1 & RI == 1) | PI==1 ~ rbinom(1,1,sens), # imperfect sensitivity for past infections; includes those infected earlier in simulation, those reinfected on this day or those previously infected
                                                          test == 1 & RI == 0 ~ rbinom(1,1,1-spec_case), # imperfect specificity for past infections, focusing on cases only who have not been previously infected
                                                          TRUE ~ as.integer(0))) -> results_analysis_sample
                      
                      results_analysis_sample %>%
                        subset(Past_infection == 0) -> results_analysis_sample_exclude_pastinfection
                      
                      sample_vax_exclude <- results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1) %>% nrow()
                      sample_control_exclude <- results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0) %>% nrow()
                      
                      results_analysis_sample %>%
                        subset(Past_infection == 1) -> results_analysis_sample_pastinfection
                      
                      # Exclude those with past infection and calculate OR = (T+vax+/T-vax+N-)/(T-vax-/T-vax-,N-)
                      VESP_TND_epi <- 1 - (results_analysis_sample_exclude_pastinfection %>% subset(test==1 & TrialStatus==1) %>% nrow() / (results_analysis_sample_exclude_pastinfection %>% subset(test==0 & TrialStatus==1) %>% nrow())) / 
                        (results_analysis_sample_exclude_pastinfection %>% subset(test==1 & TrialStatus==0) %>% nrow() / (results_analysis_sample_exclude_pastinfection %>% subset(test==0 & TrialStatus==0) %>% nrow()))
                      
                      
                      if (is.na(VESP_TND_epi)){
                        if (results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 & 
                            results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                          VESP_TND_epi<-1
                        } else if (
                          results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 & 
                          results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                          VESP_TND_epi<-0
                        } else {
                          VESP_TND_epi<-NA
                        }
                      }
                      VESP_TND_epi
                      
                      if (nrow(results_analysis_sample_exclude_pastinfection)>0){

                        # control for risk
                        model2 <- glm(test ~ TrialStatus + Risk, data=results_analysis_sample_exclude_pastinfection,family=binomial)
                        VESP_TND_epir <- 1-exp(model2$coefficients[2])
                        VESP_TND_epir

                      } else{
                        VE_OR_test <- NA
                        VESP_TND_epir <- NA
                      }

                      if (is.na(VESP_TND_epir)){
                        if (results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() == 0 &
                            results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() > 0){
                          VESP_TND_epir<-1
                        } else if (
                          results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==1 & eventstatus==1) %>% nrow() > 0 &
                          results_analysis_sample_exclude_pastinfection %>% subset(TrialStatus==0 & eventstatus==1) %>% nrow() == 0){
                          VESP_TND_epir<-0
                        } else {
                          VESP_TND_epir<-NA
                        }
                      }
                      VESP_TND_epir
                    
                    } else{
                      VESP_TND_b <- NA
                      VESP_TND_epi <- NA
                      VESP_TND_epir <- NA
                      VESP_TND_br <- NA
                      sample_vax <- NA
                      sample_control <- NA
                      sample_vax_exclude <- NA
                      sample_control_exclude <- NA
                      
                    }
                    
                    # "_b" is baseline, no stratification
                    # "_b" is baseline, controlling for risk group
                    # "_epi" is excluding for past infection
                    # "_epir" is excluding for past infection and controlling for risk group

                    analysis <- bind_cols(VES_cohort_b = VES_cohort_b, 
                                          VES_cohort_epi = VES_cohort_epi, 
                                          VES_cohort_epir = VES_cohort_epir, 
                                          VES_cohort_br = VES_cohort_br, 
                                          VES_TND_b = VES_TND_b, 
                                          VES_TND_epi = VES_TND_epi, 
                                          VES_TND_epir = VES_TND_epir,
                                          VES_TND_br = VES_TND_br,
                                          VESP_cohort_b = VESP_cohort_b,
                                          VESP_cohort_epi = VESP_cohort_epi,
                                          VESP_cohort_epir = VESP_cohort_epir,
                                          VESP_cohort_br = VESP_cohort_br,
                                          VESP_TND_b = VESP_TND_b,
                                          VESP_TND_epi = VESP_TND_epi, 
                                          VESP_TND_epir = VESP_TND_epir,
                                          VESP_TND_br = VESP_TND_br)
                    
                  
                    
                    analysis %>%
                      cbind(numevents = sum(results_analysis$eventstatus),n_VES=n_VES,n_VESP=n_VESP,cum_VES=cum_VES,cum_VESP=cum_VESP,
                            sample_vax = sample_vax, sample_control=sample_control, sample_vax_exclude = sample_vax_exclude, sample_control_exclude = sample_control_exclude,
                            R0=R,direct_VE=VE,reinf=ri,rel_VE_PI=relVE,prop_rec=pr,
                            sens=sens,spec_case=spec_case,t=t,j=i) %>%
                      bind_rows(analysis_master) -> analysis_master
                    
                    analysis_master %>%
                      mutate(VT = case_when(name=="master_aon"~"All-or-nothing",
                                            TRUE ~ "Leaky"),
                             VEP = case_when(name=="master_aon" | name=="masterVEP0" | name=="master0.2pS" ~ 0,
                                             name=="master0VES" | name== "master0VES1pS" ~ 0.85,
                                             name=="master0.65VES" ~ 6/7,
                                             TRUE ~ 0.5),
                             pS = case_when(name=="master0VES1pS"~1,
                                            name=="master0.2pS" ~ 0.2,
                                            TRUE ~ 0.5)) -> analysis_master
                    
                    write.csv(analysis_master,paste0(j,"analysis_",name,".csv"))
                    #}
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}



