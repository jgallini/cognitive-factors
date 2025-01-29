
#This program runs models to estimate a 3-year risk score for
#dementia with cognitive factor scores as the main predictors of interest.
#It then runs a clinical trial simulation to demonstrate power differences
#in detecting a dementia drug effect when selecting participants based on factor 
#scores found to be predictive of AD in 3 years in the models tested compared
#to selecting participants based on a covariate-only model and random selection

#Required data input: training and test data sets

#libraries
library(ggplot2)
library(ggmagnify)
library(scales)
library(dplyr)
library(caret)
library(epiDisplay)
library(purrr)
library(data.table)
library(table1)
library(kableExtra)
library(powerSurvEpi)
here::i_am("Code/RCT_simulation_gh.R")
library(here)

#bringing in data
load(here("Data","train_gh.Rda"))
load(here("Data","test_gh.Rda"))

#renaming data
train<-train_gh
test<-test_gh

#specifying number of iterations of simulation to run
its<-10000


##### PREDICTIVE MODELS #####

#3 year AD conversions logistic regressions using each of the 10 cognitive 
#scores as predictors along with covariates. This will demonstrate the need
#for using the factors as predictors instead of individual test scores. 
#Then running a model with each factor one at a time.
mod_func<-function(cscore){
  fml <- reformulate(c(cscore,"SEX","EDUC","NACCAGEB","RACE","NACCNE4S",
                       "HYPERT","DIABETES","SMOKYRS",
                       "BMI_bin","NACCTBI","DEPever"), response="DemEver")
  mod <- glm(fml,family=binomial(link='logit'),data=train)
  pred<-predict(mod, newdata = 
                  train[,c("Fac1_std","Fac2_std","Fac3_std","Fac4_std",
                          "SEX","EDUC","NACCAGEB","RACE","NACCNE4S",
                          "HYPERT","DIABETES","SMOKYRS","BMI_bin",
                          "NACCTBI","DEPever","LOGIMEM_3","DIGIF_3",
                          "DIGIB_3","ANIMALS_3","VEG_3","TRAILA_3",
                          "TRAILB_3","WAIS_3","MEMUNITS_3","BOSTON_3")])
  prob<-exp(pred)/(1+exp(pred))
  dem<-ifelse(prob>=0.15,1,0)
  df<-as.data.frame(cbind(pred,prob,dem,train$DemEver))
  sens<-sensitivity(as.factor(df$dem),as.factor(df$V4))
  spec<-specificity(as.factor(df$dem),as.factor(df$V4))
  return(c(sens,spec,logistic.display(mod)))
}
#Running the function that creates the 3 year AD model for each of the 
#10 tests and 4 factor scores
mod_func(cscore="LOGIMEM_3"); mod_func(cscore="MEMUNITS_3"); 
mod_func(cscore="DIGIF_3"); mod_func(cscore="DIGIB_3"); 
mod_func(cscore="ANIMALS_3"); mod_func(cscore="VEG_3"); 
mod_func(cscore="BOSTON_3"); mod_func(cscore="TRAILA_3"); 
mod_func(cscore="TRAILB_3"); mod_func(cscore="WAIS_3"); 
mod_func(cscore="Fac1_std"); mod_func(cscore="Fac2_std");
mod_func(cscore="Fac3_std"); mod_func(cscore="Fac4_std");


#Now examining 3 year AD conversion logistic regression model with all 4 factors
sensspec<-train
mallfac <- glm(DemEver ~ Fac1_std+Fac2_std+Fac3_std+Fac4_std+SEX+EDUC+
              NACCAGEB+RACE+NACCNE4S+HYPERT+DIABETES+SMOKYRS+
              BMI_bin+NACCTBI+DEPever,family=binomial(link='logit'),data=train)
sensspec$mallfac_pred<-predict(mallfac, newdata = 
                                 train[,c("Fac1_std","Fac2_std","Fac3_std",
                                          "Fac4_std","SEX","EDUC","NACCAGEB",
                                          "RACE","NACCNE4S","HYPERT",
                                          "DIABETES","SMOKYRS","BMI_bin",
                                          "NACCTBI","DEPever")])
sensspec$mallfac_prob<-exp(sensspec$mallfac_pred)/(1+exp(sensspec$mallfac_pred))
sensspec$mallfac_dem<-ifelse(sensspec$mallfac_prob>=0.15,1,0)
sensitivity(as.factor(sensspec$mallfac_dem),as.factor(sensspec$DemEver))
specificity(as.factor(sensspec$mallfac_dem),as.factor(sensspec$DemEver))
logistic.display(mallfac)
#language and memory factors highly significant here, 
#psychomotor speed also significant, leaving all
#factors in to avoid missing any confounding
mfac_final<-mallfac


#Now examining an only covariates model, no cognitive factors or cognitive tests
#included as predictors
mcov <- glm(DemEver ~ SEX+EDUC+NACCAGEB+RACE+NACCNE4S+HYPERT+DIABETES+
            SMOKYRS+BMI_bin+NACCTBI+DEPever,family=binomial(link='logit'),
            data=train)
sensspec$mcov_pred<-predict(mcov, newdata = 
                              train[,c("Fac1_std","Fac2_std","Fac3_std",
                                       "Fac4_std","SEX","EDUC","NACCAGEB",
                                       "RACE","NACCNE4S","HYPERT","DIABETES",
                                       "SMOKYRS","BMI_bin","NACCTBI",
                                       "DEPever")])
sensspec$mcov_prob<-exp(sensspec$mcov_pred)/(1+exp(sensspec$mcov_pred))
sensspec$mcov_dem<-ifelse(sensspec$mcov_prob>=0.15,1,0)
sensitivity(as.factor(sensspec$mcov_dem),as.factor(sensspec$DemEver))
specificity(as.factor(sensspec$mcov_dem),as.factor(sensspec$DemEver))
logistic.display(mcov)

#renaming test data for use in simulation
d1<-test

#initiating vectors/lists needed for simulation
power_rand2<-vector(); power_pred2<-vector();
power_covpred2<-vector(); power_cogtestpred2<-vector(); powerlist<-list();
n_rand2<-vector(); n_pred2<-vector();
n_covpred2<-vector(); n_cogtestpred2<-vector(); nlist<-list();
hr_rand<-vector(); hr_pred<-vector(); hr_covpred<-vector(); hrlist<-list();
event_rand2<-vector(); event_pred2<-vector(); event_covpred2<-vector();
eventlist<-list();


#### RCT SIMULATION ####

#Runs for a specified number of iterations at the beginning of this program
#for each of the specified treatment effects

for (treat_effect in c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.40,0.45,0.50)){
  
  #first, sampling from test set of patients
  for (i in 1:its){
    #First taking truly random sample
    sids<-unique(d1$NACCID)
    ids_rand<-as.data.frame(sample(sids, 1000, replace = TRUE, prob = NULL))
    colnames(ids_rand)<-"NACCID"
    
    #Now taking sample selected based on factor score model
    #first take factor scores from everyone in test set
    #and predict probability of AD in 3 years using mfac_final 
    #from above code
    d1$pred<-predict(mfac_final, newdata = 
                       d1[,c("Fac1_std","Fac2_std","Fac3_std","Fac4_std",
                            "SEX","EDUC","NACCAGEB","RACE","NACCNE4S",
                            "HYPERT","DIABETES","SMOKYRS","BMI_bin",
                            "NACCTBI","DEPever")])
    d1$prob<-exp(d1$pred)/(1+exp(d1$pred))
    #now sampling from subjects with >15% probability of developing 
    #AD in 3 years
    predids<-d1[d1$prob>0.15,][["NACCID"]]
    ids_pred<-as.data.frame(sample(predids, 1000, replace = TRUE, prob = NULL))
    colnames(ids_pred)<-"NACCID"
    
    #sample selected based only on covariates model, no factor scores
    #following the same selection process as with the factor model except
    #now using mcov model from above instead of mfac_final model
    d1$covpred<-predict(mcov, newdata = 
                          d1[,c("Fac1_std","Fac2_std","Fac3_std","Fac4_std",
                                "SEX","EDUC","NACCAGEB","RACE","NACCNE4S",
                                "HYPERT","DIABETES","SMOKYRS","BMI_bin",
                                "NACCTBI","DEPever")])
    d1$covprob<-exp(d1$covpred)/(1+exp(d1$covpred))
    #again sampling from subjects with >15% probability of 
    #developing AD in 3 years
    covpredids<-d1[d1$covprob>0.15,][["NACCID"]]
    ids_covpred<-as.data.frame(sample(covpredids, 1000, 
                                      replace = TRUE, prob = NULL))
    colnames(ids_covpred)<-"NACCID"
    
    #merging data with sampled ids for all 3 groups
    dat_rand<-merge(ids_rand,d1,by="NACCID")
    dat_pred<-merge(ids_pred,d1,by="NACCID")
    dat_covpred<-merge(ids_covpred,d1,by="NACCID")
    
    #simulating randomization to treatment arms for all cases
    dat_rand$treatment<-rbinom(nrow(dat_rand),1,0.5)
    dat_pred$treatment<-rbinom(nrow(dat_pred),1,0.5)
    dat_covpred$treatment<-rbinom(nrow(dat_covpred),1,0.5)
    
    #updating outcomes for those in the treatment arm for each group
    #basically just taking those who actually did get dementia by the end 
    #of the theoretical study and randomly drawing
    #a yes/no outcome from binomial distribution with p=1-treatment effect
    #also updating follow-up time for would be converters turned non-converters
    #to be three years instead of their original time to dementia
    dat_rand$Status_trial<-
      ifelse(dat_rand$Status==1 & dat_rand$treatment==1,
              rbinom(nrow(dat_rand[dat_rand$Status==1&dat_rand$treatment==1,]),
                     1,1-treat_effect),
             dat_rand$Status)
    dat_rand$Time_trial<-ifelse(dat_rand$Status==1 & dat_rand$Status_trial==0,
                                 365.25*3,dat_rand$Time)
    
    dat_pred$Status_trial<-
      ifelse(dat_pred$Status==1 & dat_pred$treatment==1,
              rbinom(nrow(dat_pred[dat_pred$Status==1&dat_pred$treatment==1,]),
                     1,1-treat_effect),
             dat_pred$Status)
    dat_pred$Time_trial<-ifelse(dat_pred$Status==1 & dat_pred$Status_trial==0,
                                365.25*3,dat_pred$Time)
    
    dat_covpred$Status_trial<-
      ifelse(dat_covpred$Status==1 & dat_covpred$treatment==1,
              rbinom(nrow(dat_covpred[dat_covpred$Status==1&
                          dat_covpred$treatment==1,]),1,1-treat_effect),
             dat_covpred$Status)
    dat_covpred$Time_trial<-ifelse(dat_covpred$Status==1 & 
                                  dat_covpred$Status_trial==0,
                                  365.25*3,dat_covpred$Time)
    
    #modeling treatment effect on time to AD for all 3 types of trials
    #adjusted for only sex, education, and APOE4
    cox_rand <-coxph(Surv(Time_trial, Status_trial)~treatment+SEX+EDUC+NACCNE4S,
                      data = dat_rand)
    cox_pred <-coxph(Surv(Time_trial, Status_trial)~treatment+SEX+EDUC+NACCNE4S,
                      data = dat_pred)
    cox_covpred <-coxph(Surv(Time_trial, Status_trial)~treatment+SEX+EDUC+
                          NACCNE4S,data = dat_covpred)
    
   
    #calculation of power for all three selection methods (Schoenfeld formula)
    
    #power with random selection
    sum_rand<-summary(cox_rand)
    hr_rand[i]<-sum_rand[["coefficients"]][1,2]
    zpower_rand<-
      sqrt(sum_rand[["nevent"]]*(nrow(dat_rand[dat_rand$treatment==1,])/
                                nrow(dat_rand))*
                    (1-(nrow(dat_rand[dat_rand$treatment==1,])/nrow(dat_rand)))*
                    (sum_rand[["coefficients"]][1,1]**2))-1.645
    power_rand<-ifelse(hr_rand[i]>1,0,pnorm(zpower_rand))
    power_rand2[i]<-power_rand
    
    #power with factor model selection
    sum_pred<-summary(cox_pred)
    hr_pred[i]<-sum_pred[["coefficients"]][1,2]
    zpower_pred<-
      sqrt(sum_pred[["nevent"]]*(nrow(dat_pred[dat_pred$treatment==1,])/
                                nrow(dat_pred))*
                    (1-(nrow(dat_pred[dat_pred$treatment==1,])/nrow(dat_pred)))*
                    (sum_pred[["coefficients"]][1,1]**2))-1.645
    power_pred<-ifelse(hr_pred[i]>1,0,pnorm(zpower_pred))
    power_pred2[i]<-power_pred
    
    #power with covariate model selection
    sum_covpred<-summary(cox_covpred)
    hr_covpred[i]<-sum_covpred[["coefficients"]][1,2]
    zpower_covpred<-
     sqrt(sum_covpred[["nevent"]]*(nrow(dat_covpred[dat_covpred$treatment==1,])/
                                     nrow(dat_covpred))*
                      (1-(nrow(dat_covpred[dat_covpred$treatment==1,])/
                            nrow(dat_covpred)))*
                      (sum_covpred[["coefficients"]][1,1]**2))-1.645
    power_covpred<-ifelse(hr_covpred[i]>1,0,pnorm(zpower_covpred))
    power_covpred2[i]<-power_covpred
    
    
    #calculating needed sample size for 0.8 power for each of the 
    #3 selection methods
    
    #using 0.84 as zbeta for 80% power and 1.65 as z alpha for one sided 0.05 
    #test, overwriting in case of hazard ratio over 1
    n_rand<-ifelse(hr_rand[i]>1,Inf,
      ((qnorm(0.8)+1.65)**2)/
      ((nrow(dat_rand[dat_rand$treatment==1,])/nrow(dat_rand))*
         (1-(nrow(dat_rand[dat_rand$treatment==1,])/nrow(dat_rand)))*
         (sum_rand[["coefficients"]][1,1]**2)))
    n_rand2[i]<-n_rand

    n_pred<-ifelse(hr_pred[i]>1,Inf,
      ((qnorm(0.8)+1.65)**2)/
      ((nrow(dat_pred[dat_pred$treatment==1,])/nrow(dat_pred))*
         (1-(nrow(dat_pred[dat_pred$treatment==1,])/nrow(dat_pred)))*
         (sum_pred[["coefficients"]][1,1]**2)))
    n_pred2[i]<-n_pred

    n_covpred<-ifelse(hr_covpred[i]>1,Inf,
      ((qnorm(0.8)+1.65)**2)/
          ((nrow(dat_covpred[dat_covpred$treatment==1,])/nrow(dat_covpred))*
          (1-(nrow(dat_covpred[dat_covpred$treatment==1,])/nrow(dat_covpred)))*
          (sum_covpred[["coefficients"]][1,1]**2)))
    n_covpred2[i]<-n_covpred

    #keeping track of total number of events for each of the 3 selection methods
    #for number needed to recruit calculation later
    event_rand<-sum_rand[["nevent"]]
    event_rand2[i]<-event_rand
    
    event_pred<-sum_pred[["nevent"]]
    event_pred2[i]<-event_pred
    
    event_covpred<-sum_covpred[["nevent"]]
    event_covpred2[i]<-event_covpred
    
    #keeping track of where we are in the simulation as it runs
    if((i %% 10) == 0)print(paste("iteration", i, "treatment effect",
                                  treat_effect))
  }
  
  #storing all of the necessary output into lists
  powerlist[[treat_effect*20]]<-as.data.frame(
    cbind(power_pred2,power_rand2,power_covpred2,treat_effect))
  nlist[[treat_effect*20]]<-as.data.frame(
    cbind(n_pred2,n_rand2,n_covpred2,treat_effect))
  hrlist[[treat_effect*20]]<-as.data.frame(
    cbind(hr_pred,hr_rand,hr_covpred,treat_effect))
  eventlist[[treat_effect*20]]<-as.data.frame(
    cbind(event_pred2,event_rand2,event_covpred2,treat_effect))
}

#saving output
save(powerlist,file=here("Data","powerlist_gh.Rda"))
save(nlist,file=here("Data","nlist_gh.Rda"))
save(hrlist,file=here("Data","hrlist_gh.Rda"))
save(eventlist,file=here("Data","eventlist_gh.Rda")

