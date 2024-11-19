

# This program creates tables and figures for the main results of the factor
# analysis and the randomized controlled trial simulation

# The only table not created by this program is the table summarizing 
# 3 year AD risk model fits. Those table results are taken from the 
# RCT simulation program.

#libraries
library(table1)
library(Hmisc)
library(dplyr)
library(kableExtra)
library(purrr)
library(MASS) 
library(reshape2) 
library(reshape)
library(data.table)
library(ggplot2)
library(ggmagnify)
library(scales)

#bringing in data
load("test_gh.Rda")
load("train_gh.Rda")
load("powerlist_gh.Rda")
load("nlist_gh.Rda")
load("eventlist_gh.Rda")
load("hrlist_gh.Rda")

#renaming data
test<-test_gh
train<-train_gh


###### Table 1 #########

#cleaning and labeling test data to be ready for table 1 function
test$SEX <- factor(test$SEX, levels=c(1,2),labels=c("Male",
                 "Female"))
test$Status<-factor(test$Status,levels=c(0,1),labels=c(
  "AD 3-Year Non-converter","AD 3-Year Converter"))
test$RACE<-factor(test$RACE,levels=c(1,2,3,5),labels=c(
  "White","Black or African-American","American Indian/Alaska Native or
  Native Hawaiian/Other Pacific Islander","Asian"))
test$NACCNE4S<-factor(test$NACCNE4S,levels=c(0,1,2),labels=c(
  "0 e4 alleles","1 e4 allele","2 e4 alleles"))
test$HYPERTEN<-factor(test$HYPERTEN,levels=c(0,1,2),labels=c(
  "Never","Recent/Active","Remote/Inactive"))
test$DIABETES<-factor(test$DIABETES,levels=c(0,1),labels=c(
  "Never/Inactive","Recent/Active"))
test$BMI_bin<-factor(test$BMI_bin,levels=c(0,1),labels=
                       c("Non-Obese","Obese"))
test$NACCTBI<-factor(test$NACCTBI,levels=c(0,1),labels=
                       c("No","Yes"))
test$DEPever<-factor(test$DEPever,levels=c(0,1),labels=
                       c("No","Yes"))
label(test$SEX)<-"Sex at Birth"; label(test$EDUC)<-"Years of Education";
label(test$NACCAGEB)<-"Age at Baseline";label(test$RACE)<-"Race";
label(test$NACCNE4S)<-"Number of e4 alleles"; 
label(test$HYPERTEN)<-"Hypertension";
label(test$DIABETES)<-"Diabetes"; label(test$SMOKYRS)<-"Smoking Years";
label(test$BMI_bin)<-"BMI"; label(test$NACCTBI)<-"TBI Ever";
label(test$DEPever)<-"Depression Ever"
label(test$Time)<-"Follow-Up Time"

#creating table 1
tab1<-table1(~ Time+SEX+EDUC+NACCAGEB+RACE+NACCNE4S+HYPERTEN+DIABETES+
         SMOKYRS+BMI_bin+NACCTBI+DEPever | Status, data=test)

#creating LaTex code to input table 1 into paper
tab1 %>%
  kbl(caption="Summary Statistics of Baseline Covariates by AD Status",
      format="latex",
      #col.names = c("Gender","Education","Count","Mean","Median","SD"),
      align="r") %>%
  kable_paper(full_width = F,  html_font = "Source Sans Pro")



####### Data cleaning for Power and N plots ########

#cleaning data from simulation to be plotted
powerdf<-list_rbind(powerlist)
powerdf2 <- melt(setDT(powerdf), id.vars = c("treat_effect"), 
                 variable.name = "trial_type")
ndf<-list_rbind(nlist)
ndf2 <- melt(setDT(ndf), id.vars = c("treat_effect"), 
             variable.name = "trial_type")
eventdf<-list_rbind(eventlist)
eventdf2<-melt(setDT(eventdf),id.vars = c("treat_effect"), 
               variable.name = "trial_type")
hrf<-list_rbind(hrlist)
hrf2<-melt(setDT(hrf), id.vars = c("treat_effect"), 
           variable.name = "trial_type")


#calculating bootstrapped median by group for graphing (power)
graphdf<- powerdf2 %>% group_by(trial_type,treat_effect) %>% 
  mutate(median=median(value))
graphdf2<-unique(graphdf[,c("treat_effect","trial_type","median")])
graphdf2$trial_type_clean<-ifelse(graphdf2$trial_type=="power_pred2",
                                  "Factor Model",
                                ifelse(graphdf2$trial_type=="power_rand2",
                                       "Random","Covariate Model"))
                                                              
#calculating bootstrapped median by group for graphing (N needed)
eventdf2$needdeaths<-ndf2$value
eventdf2$recruitnum<-eventdf2$needdeaths/(eventdf2$value/1000)
graphevent<- eventdf2 %>% group_by(trial_type,treat_effect) %>% 
  mutate(median=median(recruitnum))
graphevent$trial_type_clean<-
  ifelse(graphevent$trial_type=="event_pred2","Factor Model",
      ifelse(graphevent$trial_type=="event_rand2","Random","Covariate Model"))


##### Power and N needed plots #####

#plotting the difference in median power over treatment effects
#saving as .tiff file
tiff("powerplot.tiff", units="in", width=6, height=5, res=1000)
ggplot(graphdf2[graphdf2$trial_type_clean!="Cognitive Test Model",], 
       aes(x=treat_effect, y=median, group=trial_type_clean, 
           color=trial_type_clean)) + 
  geom_line() + geom_point() + ylab("Median Power (10,000 trials)") +
  xlab("Treatment Effect Size") +
  ggtitle("Figure 3. Median RCT Power Using Different Sampling Methods")+
  theme(legend.position = c(0.8, 0.28))+
  scale_color_discrete(name = "Sampling Type")
dev.off()


#plotting sample sizes required for each trial type (median)
tiff("recruitplot.tiff", units="in", width=6, height=5, res=1000)
ggplot(graphevent, aes(x=treat_effect, y=median, group=trial_type_clean, 
                       color=trial_type_clean)) +
  geom_line() + geom_point() + ylab("Median N(10,000 trials)") +
  xlab("Treatment Effect Size") +
  ggtitle("Figure 4. Median N Required to Recruit for 0.8 Power Using Different 
  Sampling Methods")+
  theme(legend.position = c(0.75, 0.65))+
  scale_color_discrete(name = "Sampling Type")+
  geom_magnify(from=c(0.2,0.5,100,7000),to=c(0.2,0.5,30000,120000),
               axes="xy")+
  scale_y_continuous(labels = label_comma())
dev.off()


###### Observed treatment effect density plot #######

#cleaning hazard ratio data
graphdfhr<- hrf2 %>% group_by(trial_type,treat_effect) %>% 
  mutate(median=median(value),)
graphdfhr$trial_type_clean<-ifelse(graphdfhr$trial_type==
                                     "hr_pred",
                                   "Factor Model",
                                ifelse(graphdfhr$trial_type==
                                         "hr_rand","Random",
                                       "Covariate Model"))

# graph showing the varied distributions of hazard ratios from the 
# three methods
#sving as .tiff file
tiff("dens3plot.tiff", units="in", width=6, height=5, res=1000)
ggplot(graphdfhr[graphdfhr$treat_effect==0.2,], 
                  aes(x = value, fill=trial_type_clean)) + 
  geom_density(alpha = 0.5)+geom_vline(xintercept = 0.8, linetype="dotted", 
                                       color = "black", linewidth=1)+
  ggtitle("Figure 2. Observed Hazard Ratios for a Theoretical Hazard Ratio 
          of 0.8 (Treatment Effect of 0.2)")+
  ylab("Density over 10,000 iterations")+xlab("Empirical Hazard Ratio")+
  theme(legend.position = c(0.75, 0.6))+
  guides(fill=guide_legend(title="Sampling Type"))
dev.off()



