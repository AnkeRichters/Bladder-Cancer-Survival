
library(readr)
library(rstpm2)
library(kableExtra)
library(survminer)
library(ggplot2)
library(data.table)
library(table1)
options(scipen=999)

# 
# ## DATA
# **NKR data and CBS lifetables**  
#   
# Data was used from the Netherlands Cancer Registry (NCR) of all patients:  
# * with bladder cancer;  
# * T1 or higher (all bladder cancers except Ta or Tis);  
# * histologically confirmed;  
# * diagnosed in 2003-2012.   
# 
# 
# Data on mortality rates of the Dutch general population were obtained from CBS (Centraal Bureau voor Statistiek).  
# Data sources are merged together.    
# Lifetables were not available for age>99 and year>2017.    
# Patients who survived past 10 years are censored at 10 years.    

dataset <- read_csv("complete dataset.csv")

# The merged data look like this (n=24,169):
# * Vitstat (1=deceased, 0=censored)
# * Fupyrs (time to vitstat in years)
# * Sex (0=male, 1=female)
# * Age = age in years at diagnosis
# * Tstage (numerical)
# * Nstage (dichotomous, N0 / N+)
# * Mstage (dichotomous, M0 / M+)
# * Morf = histology (0=urothelial carcinoma, 1=non-urothelial carcinoma)
# * Rate = background hazard from life-tables
# * Rate_sxx = background hazard from life-tables adapted based on several sensitivity scenarios

## MODELLING


####  Model summaries
# Basic models: only sex as independent variable, 1df for baseline hazard.
b3<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
             sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=1,
           data=dataset, bhazard=dataset$rate)
#baseline spline
df4<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
              sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=4,                            
            data=dataset, bhazard=dataset$rate, link.type="PH")

#TVCs
tvc3 <- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
                sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=4,  
              tvc=list(sex=1, age=4, tstage=2),
              data=dataset, bhazard=dataset$rate, link.type="PH")


#Basic model (sex + covariates + interaction terms between sex and covariates):
summary(b3)

#Basic model + spline for baseline hazard with 4 df:
summary(df4)

#Basic model with baseline hazard spline with time-covariate interactions to allow for non-proportional hazards over time (final model):
summary(tvc3)


## RESULTS
# **Models and figures used in paper**  
# All the standardized survival functions are corrected for age at diagnosis, T-stage, N-stage, M-stage and histology.

### Table 1 
#**Baseline characteristics of a population-based cohort of 24,169 bladder cancer patients**
table1(~ age + factor(tstage) + factor(nstage) + factor(mstage) + factor(morf) | sex, data=dataset)

#*Models:*
basic  <- stpm2(Surv(fupyrs, vitstat)~sex,
                data=dataset, df=4, 
                bhazard=dataset$rate, 
                tvc=list(sex=1))

basic_rel<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72), 
                  data=dataset, df=4, 
                  bhazard=dataset$rate, 
                  tvc=list(sex=1,age=4), link.type="PH")

full <- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
                sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=4,  
              tvc=list(sex=1, age=4, tstage=2),
              data=dataset, 
              bhazard=dataset$rate, link.type="PH")

### Figure 1
#**Modelled overall and relative survival of male and female bladder cancer patients**
  
#Observed and relative survival model:
obs<- stpm2(Surv(fupyrs, vitstat)~sex, data=dataset, df=4,tvc=list(sex=1))
rel<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72), data=dataset, df=4, bhazard=dataset$rate, tvc=list(sex=1,age=4))

plot(obs,newdata=data.frame(sex=1),type="surv",var="sex",add=F,ci=F,line.col=2, lty=2, 
     xlim=c(0,10), ylim=c(0,1), xlab="Time in years since diagnosis", ylab="Survival",rug=FALSE)
plot(obs,newdata=data.frame(sex=0),type="surv",var="sex",add=T,ci=F,line.col=1, lty=2,rug=FALSE)
plot(rel,newdata=transform(dataset, sex=1), type="meansurv",var="sex",add=T,ci=F,line.col=2, lty=1,rug=FALSE)
plot(rel,newdata=transform(dataset, sex=0),type="meansurv",var="sex",add=T,ci=F,line.col=1, lty=1,rug=FALSE)
legend("topright",c("Women (overall survival)","Men (overall survival)", 
                    "Women (relative survival)","Men (relative survival)"),
       col=c(2,1,2,1), lty=c(2,2,1,1))

### Figure 2A
#**Unadjusted and adjusted relative survival by sex in bladder cancer patients**

#Plot relative survival of men and women, only standardized for age
# Women
plot(basic_rel,newdata=transform(dataset, sex=1),type="meansurv",var="sex",add=F,ci=F,line.col=2, lty=1,
     ylim=c(0,1), xlim=c(0,10), xlab="Time in years since diagnosis", ylab="Survival",rug=FALSE)
# Men
plot(basic_rel,newdata=transform(dataset, sex=0),type="meansurv",var="sex",add=T,ci=F,line.col=1, lty=1,rug=FALSE)

# Plot relative survival of women with men's covariate distribution
subset_men<-subset(dataset, sex==0)      # select only men
subset_men$sex<-1                        # set men's covariate distribution to females
plot(full,newdata=transform(subset_men, sex=1), type="meansurv", add=T, line.col=4, lty=1,rug=FALSE)

legend("topright",c("Women (relative survival)","Men (relative survival)", 
                    "Women (standardized relative survival)"),
       col=c(2,1,4), lty=c(1,1,1)) 


### Figure 2B
#**Unadjusted and adjusted excess hazard by sex in bladder cancer patients**
yticks <- c(0,0.1,0.2,0.3,0.4,0.5,0.6)
ylabels <- c(0,100,200,300,400,500,600)

#Plot hazard functions of men and women, only standardized for age

#Women
plot(basic_rel,newdata=transform(dataset, sex=1),type="meanhaz", var="sex",add=F, ci=F, line.col=2,lty=1,
     ylim=c(0,0.6), xlim=c(0,10), xlab="Time in years since diagnosis", ylab="Sex-specific excess hazard (rate/1000py)", yaxt="n",rug=FALSE)
#Men
plot(basic_rel,newdata=transform(dataset, sex=0),type="meanhaz", var="sex",add=T, ci=F, line.col=1,lty=1,rug=FALSE)
axis(2, at=yticks, labels=ylabels)

# Plot relative survival of women with men's covariate distribution
subset_men<-subset(dataset, sex==0)      # select only men
subset_men$sex<-1                        # change 'man' to 'woman' so now it's a group of women with men's covariates
plot(full,newdata=transform(subset_men, sex=1), type="meanhaz", add=T, line.col=4, lty=1,rug=FALSE)

legend("topright",c("Women (excess hazard)","Men (excess hazard)", 
                    "Women (standardized excess hazard)"),
       col=c(2,1,4), lty=c(1,1,1)) 


### Figure 3A
#**Unadjusted and adjusted excess hazard ratio for sex**

#Plot hazard ratio curve in dashed line
plot(basic_rel,newdata=transform(dataset, sex=0),type="meanhr", var="sex",add=F, ci=F, line.col=1,lty=2,
     ylim=c(0,3), xlim=c(0,10), xlab="Time in years since diagnosis", ylab="Excess hazard ratio",rug=FALSE)
abline(h=1, col="grey80",lty=1)

#Plot standardized hazard ratio curve in normal line
plot(full,newdata=transform(dataset, sex=0), type="meanhr", var="sex", add=T, line.col=1, lty=1,rug=FALSE)
legend("topright",c("Women vs men: unadjusted HR","Women vs men: adjusted HR","No difference in excess hazard"),
       col=c(1,1,"grey80"), lty=c(2,1,1))


### Figure 3B
#**Adjusted excess hazard ratios for age, T-stage, N-stage, M-stage and histology in bladder cancer patients**

full <- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
                sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=4,  
              tvc=list(sex=1, age=4, tstage=2),
              data=dataset, bhazard=dataset$rate, link.type="PH")

#Plot standardized hazard ratio curve of each factor
plot(full,newdata=transform(dataset, tstage=1), type="meanhr", var="tstage", add=F, line.col=2, lty=1,
     ylim=c(0,3), xlim=c(0,10), xlab="Time in years since diagnosis", ylab="Excess hazard ratio",rug=FALSE)
plot(full,newdata=transform(dataset, age=72), type="meanhr", var="age", add=T, line.col=3, lty=1,rug=FALSE)
plot(full,newdata=transform(dataset, nstage=0), type="meanhr", var="nstage", add=T, line.col=4, lty=1,rug=FALSE)
plot(full,newdata=transform(dataset, mstage=0), type="meanhr", var="mstage", add=T, line.col=5, lty=1,rug=FALSE)
plot(full,newdata=transform(dataset, morf=0), type="meanhr", var="morf", add=T, line.col="coral4", lty=1,rug=FALSE)

legend("bottomright",c("HR for age (adjusted for sex, T-stage N-stage, M-stage, histology)",
                       "HR for T-stage (adjusted for age, sex, N-stage, M-stage, histology)",
                       "HR for N-stage (adjusted for age, sex, T-stage, M-stage, histology)",
                       "HR for M-stage (adjusted for age, sex, T-stage, N-stage, histology)",
                       "HR for histology (adjusted for age, sex, T-stage, N-stage, M-stage)",
                       "No difference in excess hazard"),
       col=c(3,2,4,5,"coral4","grey80"), lty=c(1,1,1,1,1,1),
       cex=0.85)
abline(h=1, col="grey80",lty=1)


### Figure 4A
#**Adjusted excess hazard ratios for sex in patients with NMIBC and MIBC**
  
# T-stage subgroups: NMIBC vs MIBC
tstage1 <- subset(dataset, tstage==1)
tstage2 <- subset(dataset, tstage>1)

#All factors associated with T-stage are removed in the first model because every patient has T-stage 1 there.
fullt1<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+nstage+mstage+morf+
                 sex:nsx(age, df=3, centre=72)+sex:nstage+sex:mstage+sex:morf, df=4,  
               tvc=list(sex=1, age=4),
               data=tstage1, bhazard=tstage1$rate, link.type="PH")
fullt2<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
                 sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=4,  
               tvc=list(sex=1, age=4, tstage=2),
               data=tstage2, bhazard=tstage2$rate, link.type="PH")

#Plot standardized hazard ratio curve in normal line
plot(fullt1,newdata=transform(tstage1, sex=0), type="meanhr", var="sex", add=F, line.col=2, lty=1,
     ylim=c(0,3), xlim=c(0,10), xlab="Time in years since diagnosis", ylab="Excess hazard ratio",rug=FALSE)
plot(fullt2,newdata=transform(tstage2, sex=0), type="meanhr", var="sex", add=T, line.col=3, lty=1,rug=FALSE)
legend("topright",c("Adjusted HR - NMIBC","Adjusted HR - MIBC","No difference in excess hazard"),
       col=c(2,3,"grey80"), lty=c(1,1,1))
abline(h=1, col="grey80",lty=1)


### Figure 4B
#**Adjusted excess hazard ratios for sex in patients by age groups**
  
# Age subgroups
age1 <- subset(dataset, age<70)
age2 <- subset(dataset, age>69 & age<80)
age3 <- subset(dataset, age>79)


fulla1<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
                 sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=4,  
               tvc=list(sex=1, age=4, tstage=2),
               data=age1, bhazard=age1$rate, link.type="PH")
fulla2<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
                 sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=4,  
               tvc=list(sex=1, age=4, tstage=2),
               data=age2, bhazard=age2$rate, link.type="PH")
fulla3<- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
                 sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, df=4,  
               tvc=list(sex=1, age=4, tstage=2),
               data=age3, bhazard=age3$rate, link.type="PH")

#Plot standardized hazard ratio curve in normal line
plot(fulla1,newdata=transform(age1, sex=0), type="meanhr", var="sex", add=F, line.col=2, lty=1,
     ylim=c(0,3), xlim=c(0,10), xlab="Time in years since diagnosis", ylab="Excess hazard ratio",rug=FALSE)
plot(fulla2,newdata=transform(age2, sex=0), type="meanhr", var="sex", add=T, line.col=3, lty=1,rug=FALSE)
plot(fulla3,newdata=transform(age3, sex=0), type="meanhr", var="sex", add=T, line.col=4, lty=1,rug=FALSE)
legend("topright",c("Adjusted HR - age 0-69","Adjusted HR - age 70-79", "Adjusted HR - age 80+","No difference in excess hazard"),
       col=c(2,3,4,"grey80"), lty=c(1,1,1,1))
abline(h=1, col="grey80",lty=1)


## APPENDIX 1
#**Sensitivity analysis**  
#*Analysis of scenarios about smoking prevalence in bladder cancer patients and gen.pop.*
  
### Table A1
#**Scenarios regarding smoking in male and female general population and bladder cancer patients**

rows<-c("% Smokers among male pts","% Smokers among female pts", 
        "% Smokers among male gen.pop.", "% Smokers among female gen.pop.",
        "Odds ratio of smoking for all-cause mortality")
X1<-c(86,67,59,49,2.0)
X2<-c(92,75,75,57,2.0)

scenarios <- cbind(rows, X1, X2)

kable(scenarios, caption="Scenarios for sensitivity analysis", digits=3,
      col.names=c("Parameter","Scenario 1","Scenario 2"))%>%
  kable_styling(bootstrap_options="striped", position="left", full_width=F)

#*Background hazard rates were adjusted according to these scenarios in SAS, code can be found here:  
# (Note: here, sex is coded as 1=male 2=female)*
# data sens2;
# set sens;
#  
# *SCENARIO 1;
# if sex=1 then do;
#       %let alpha=0.59;        * Proportion of smokers among general population of this gender;
#       %let or=2;              * Odds ratio of smoking for all cause-mortality ;
#       %let props=0.86;        * Proportion of smokers among bladder cancer patients of this gender;
#       %let propns=0.14;       * Proportion of non-smokers among bladder cancer patients of this gender;
#       prob_nonsmoking_s1=(((1+((prob-&alpha)*(1-&or)))*-1)+SQRT((( 1+((prob-&alpha)*(1-&or)))**2)+((4*prob)*((1-&or)*(&alpha-1)))))/(2*((1-&or)*(&alpha-1)));
#       prob_smoking_s1=(((prob_nonsmoking_s1/(1-prob_nonsmoking_s1))*&or)/(((prob_nonsmoking_s1/(1-prob_nonsmoking_s1))*&or)+1));
#             prob_1=((&props*prob_smoking_s1)+(&propns*prob_nonsmoking_s1));
# end;
# if sex=2 then do;
#       %let alpha=0.49;        * Proportion of smokers among general population of this gender;
#       %let or=2;              * Odds ratio of smoking for all cause-mortality ;
#       %let props=0.67;        * Proportion of smokers among bladder cancer patients of this gender;
#       %let propns=0.33;       * Proportion of non-smokers among bladder cancer patients of this gender;
#       prob_nonsmoking_s1=(((1+((prob-&alpha)*(1-&or)))*-1)+SQRT((( 1+((prob-&alpha)*(1-&or)))**2)+((4*prob)*((1-&or)*(&alpha-1)))))/(2*((1-&or)*(&alpha-1)));
#       prob_smoking_s1=(((prob_nonsmoking_s1/(1-prob_nonsmoking_s1))*&or)/(((prob_nonsmoking_s1/(1-prob_nonsmoking_s1))*&or)+1));
#             prob_1=((&props*prob_smoking_s1)+(&propns*prob_nonsmoking_s1));
# end;
# drop prob_nonsmoking_s1 prob_smoking_s1;
#  
# *SCENARIO 2;
# if sex=1 then do;
#       %let alpha=0.75;              * Proportion of smokers among general population of this gender;
#       %let or=2;                  * Odds ratio of smoking for all cause-mortality ;
#       %let props=0.92;        * Proportion of smokers among bladder cancer patients of this gender;
#       %let propns=0.08;       * Proportion of non-smokers among bladder cancer patients of this gender;
#       prob_nonsmoking_s1=(((1+((prob-&alpha)*(1-&or)))*-1)+SQRT((( 1+((prob-&alpha)*(1-&or)))**2)+((4*prob)*((1-&or)*(&alpha-1)))))/(2*((1-&or)*(&alpha-1)));
#       prob_smoking_s1=(((prob_nonsmoking_s1/(1-prob_nonsmoking_s1))*&or)/(((prob_nonsmoking_s1/(1-prob_nonsmoking_s1))*&or)+1));
#             prob_2=((&props*prob_smoking_s1)+(&propns*prob_nonsmoking_s1));
# end;
# if sex=2 then do;
#       %let alpha=0.57;        * Proportion of smokers among general population of this gender;
#       %let or=2;                  * Odds ratio of smoking for all cause-mortality ;
#       %let props=0.75;        * Proportion of smokers among bladder cancer patients of this gender;
#       %let propns=0.25;       * Proportion of non-smokers among bladder cancer patients of this gender;
#       prob_nonsmoking_s1=(((1+((prob-&alpha)*(1-&or)))*-1)+SQRT((( 1+((prob-&alpha)*(1-&or)))**2)+((4*prob)*((1-&or)*(&alpha-1)))))/(2*((1-&or)*(&alpha-1)));
#       prob_smoking_s1=(((prob_nonsmoking_s1/(1-prob_nonsmoking_s1))*&or)/(((prob_nonsmoking_s1/(1-prob_nonsmoking_s1))*&or)+1));
#             prob_2=((&props*prob_smoking_s1)+(&propns*prob_nonsmoking_s1));
# end;
# drop prob_nonsmoking_s1 prob_smoking_s1;
#  
# X_age=_age;
# X_period=_period;
# drop _age _period;
# run;


### Figure A1
#**Excess hazard ratio for sex in bladder cancer patients adjusted for age, T-stage, N-stage, M-stage and histology based on different scenarios regarding smoking prevalence**

full <- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
                sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, data=dataset, df=4, 
              tvc=list(sex=1, age=4, tstage=2),
              bhazard=dataset$rate, link.type="PH")
s1 <- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
              sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, data=dataset, df=4, 
            tvc=list(sex=1, age=4, tstage=2),
            bhazard=dataset$rate_s1a, link.type="PH")
s2 <- stpm2(Surv(fupyrs, vitstat)~sex+nsx(age, df=3, centre=72)+tstage+nstage+mstage+morf+
              sex:nsx(age, df=3, centre=72)+sex:tstage+sex:nstage+sex:mstage+sex:morf, data=dataset, df=4, 
            tvc=list(sex=1, age=4, tstage=2),
            bhazard=dataset$rate_s1b, link.type="PH")

plot(full,newdata=transform(dataset, sex=0),type="meanhr",var="sex",add=F,ci=F,line.col=1, lty=1,
     xlab="Follow-up time in years",ylab="Excess hazard ratio",xlim=c(0,10), ylim=c(0, 3),rug=FALSE)
plot(s1,newdata=transform(dataset, sex=0),type="meanhr",var="sex",add=T,ci=F,line.col=2, lty=1,rug=FALSE)
plot(s2,newdata=transform(dataset, sex=0),type="meanhr",var="sex",add=T,ci=F,line.col=3, lty=1,rug=FALSE)
legend("topright", c("Unadjusted lifetables","Lifetables adjusted to scenario 1","Lifetables adjusted to scenario 2"),
       col=c(1,2,3), 
       lty=c(1,1,1))


## ADDITIONAL RESULTS
#*Not displayed in paper*  
  
### Illustration of the proportional hazards assumption
#**Hazard function for sex without and with TVC (i.e. proportional and non-proportional  hazards assumed), unadjusted.**
  
#Model **with** proportional hazards assumption:

fit1<- stpm2(Surv(fupyrs, vitstat)~sex, data=dataset, df=4, bhazard=dataset$rate)

plot(fit1,newdata=data.frame(sex=1),type="hazard",var="sex",add=F,ci=T,line.col=2, lty=1, rug=F,
     xlab="Follow-up time in years",ylab="Sex-specific excess hazard",xlim=c(0,10), ylim=c(0, 1.2))
plot(fit1,newdata=data.frame(sex=0),type="hazard",var="sex",add=T,ci=T,line.col=1, lty=1, rug=F)
legend("topright",c("Women","Men"),
       col=c(2,1), lty=c(1,1)) 

#Model **without** proportional hazards assumption:
fit2<- stpm2(Surv(fupyrs, vitstat)~sex, data=dataset, df=4, bhazard=dataset$rate, tvc=list(sex=1))

plot(fit2,newdata=data.frame(sex=1),type="hazard",var="sex",add=F,ci=T,line.col=2, lty=2, rug=F,
     xlab="Follow-up time in years",ylab="Sex-specific excess hazard",xlim=c(0,10), ylim=c(0, 1.2))
plot(fit2,newdata=data.frame(sex=0),type="hazard",var="sex",add=T,ci=T,line.col=1, lty=2, rug=F)
legend("topright",c("Women","Men"),
       col=c(2,1), lty=c(2,2)) 

#Comparing hazard ratio curves when assuming vs. not assuming proportional hazards:
plot(fit1,newdata=data.frame(sex=0),type="hr",var="sex",add=F,ci=T,line.col=1, lty=1,
     xlab="Follow-up time in years",ylab="Sex-specific excess hazard",xlim=c(0,10), ylim=c(0, 4))
plot(fit2,newdata=data.frame(sex=0),type="hr",var="sex",add=T,ci=T,line.col=2, lty=1)
legend("topright",c("HR with assuming proportional hazards","HR without assuming proportional hazards"),
       col=c(1,2), lty=c(1,1)) 

#AIC for model **with** propotional hazards:
print(AIC(fit1))

#AIC for model **without** propotional hazards:
print(AIC(fit2))
