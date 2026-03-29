####### <Replication script for manuscript:"Survey mode and nonresponse bias: 
#a meta-analysis based on the data from ISSP 1996–2018 and ESS rounds 1 to 9" PLOS ONE>####
# Copyright (C) <2022>  <Adam Rybak>

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.
###################################################################

#install and load packages
install.packages("remotes")
remotes::install_github("wviechtb/metafor")
library(metafor)

if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("MathiasHarrer/dmetar")

pacman::p_load(metafor,ggplot2,dplyr,MuMIn,dmetar,knitr,data.table,knitr,xtable,gtsummary,stats,tidyr,expss,readr)

#data preparation
dataset <- read_delim("dataset_Rybak_PLOSONE.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)
dataset <- mutate(dataset, sei = sqrt(vi))
dataset<-mutate(dataset, RR_vi = RR_SE ^ 2)
modelevels<-c("Face-to-face", "Mixed-mode (with f2f)", "Mail mode", "Mixed-mode (no f2f)")
framelevels<-c("Individual frame", "Nonindividual frame")
proglevels<-c("ISSP","ESS")
dataset$mod<- factor(dataset$mod, labels=modelevels)
dataset$fra<- factor(dataset$fra, labels=framelevels)
dataset$int<-interaction(dataset$mod,dataset$fra)
dataset$cntry<-as.factor(dataset$cntry)
dataset$region<-as.factor(dataset$region)
dataset$wave<-as.factor(dataset$wave)
dataset$inc<-as.factor(dataset$`Income WB 2020`)
dataset$prog<-factor(dataset$prog, labels=proglevels)
dataset$gdpgroup<-factor(dataset$gdpgroup, labels=c(1:20))
dataset$fra2<-relevel(dataset$fra,ref = "Nonindividual frame")
wavenames<-dataset$wave
regnames<-dataset$region
intnames<-dataset$int
prognames<-dataset$prog
dataset$rep<-dataset$repetitions
var_lab(dataset$rep) = "Repetitions in survey program"
dataset$drop<-ifelse(dataset$surv=="ESS6_SK",1,ifelse(dataset$surv=="ESS5_SK",1,ifelse(dataset$surv=="ESS4_SK",1,ifelse(dataset$surv=="ESS9_LT",1,0))))
dataset1<-setDT(dataset)
dataset$euro<-ifelse(dataset$region=="Western Europe"|dataset$region=="Eastern Europe"|dataset$region=="Northern Europe"|dataset$region=="Southern Europe", 1,0)

#frame and mode distribution
fra_mod <- table(dataset$fra,dataset$mod)
fra_mod <- addmargins(fra_mod)
kable(fra_mod)

#frequencies of modes used in modes categories
mod_composition<-cross_cases(dataset1,cell_vars=mod,mdset(F2F %to% SELF_COMPLETION))
mod_composition
#country and mode distribution
cntry_mod<-table(dataset$cntry,dataset$mod)
cntry_mod <- cntry_mod[order(-rowSums(cntry_mod)),]
cntry_mod <- addmargins(cntry_mod)
kable(cntry_mod)

#region and mode distribution
reg_mod<-table(dataset$region,dataset$mod)
reg_mod <- reg_mod[order(-rowSums(reg_mod)),]
reg_mod <- addmargins(reg_mod)
kable(reg_mod)

#mean NBI values for regions, modes and frames and their interaction
reg_nbi <- dataset1[ ,list(mean=mean(kohler)), by=region]
kable(reg_nbi)
mod_nbi <- dataset1[ ,list(mean=mean(kohler)), by=mod]
kable(mod_nbi)
fra_nbi <- dataset1[ ,list(mean=mean(kohler)), by=fra]
kable(fra_nbi)
int_nbi <- dataset1[ ,list(mean=mean(kohler)), by=int]
kable(int_nbi)

#forestplots for basic random-effects meta-analyses with subgroups, with Knapp-Hartung modification
modes <- rma(yi = kohler, vi = vi, mods = ~mod -1, test="knha", data=dataset)
modes_rr6 <- rma(yi = RR6, sei = RR_SE, mods = ~mod -1, test="knha", data=dataset)
regions <- rma(yi = kohler, vi = vi, mods = ~region -1, test="knha", data=dataset)
interactions <- rma(yi = kohler, vi = vi, mods = ~int -1, test="knha", data=dataset)
waves <- rma(yi = kohler, vi = vi, mods = ~wave -1, test="knha", data=dataset)
programs <-rma(yi = kohler, vi = vi, mods = ~prog -1, test="knha", data=dataset)
forest(coef(modes), sei=modes$se, refline=NA, slab=c("Face-to-face", "Mixed-mode (with f2f)", "Mail mode", "Mixed-mode (no f2f)"), xlab="Nonresponse bias indicator", header=c("Survey Modes"))
forest(coef(modes_rr6), sei=modes_rr6$se, refline=NA, slab=c("Face-to-face", "Mixed-mode (with f2f)", "Mail mode", "Mixed-mode (no f2f)"), xlab="Response rate RR6", header=c("Survey Modes"))
forest(coef(regions), sei=regions$se, refline=NA, slab=levels(regnames), xlab="Nonresponse bias indicator", header=c("World Region"))
forest(coef(interactions), sei=interactions$se, refline=NA, slab=levels(intnames), xlab="Nonresponse bias indicator", header=c("Mode and frame interactions"))
forest(coef(waves), sei=waves$se, refline=NA, slab=levels(wavenames), xlab="Nonresponse bias indicator", header=c("Survey waves"))
forest(coef(programs), sei=programs$se, refline=NA, slab=levels(prognames), xlab="Nonresponse bias indicator", header=c("Survey programs"))

#mixed-effects meta-regressions with Knapp-Hartung modification, NBI as effect size
nbi<-rma(yi = kohler, vi = vi, test="knha",slab=surv, data=dataset)
summary(nbi)
model1<-rma(yi = kohler, vi = vi, mods=~mod*fra, test="knha",slab=surv, data=dataset)
summary(model1)
model2<-rma(yi = kohler, vi = vi, mods=~prog, test="knha",slab=surv, data=dataset)
summary(model2)
model3<-rma(yi = kohler, vi = vi, mods=~repetitions, test="knha", slab=surv,data=dataset)
summary(model3)
model4<-rma(yi = kohler, vi = vi, mods=~mod*fra+prog, test="knha",slab=surv, data=dataset)
summary(model4)
model5<-rma(yi = kohler, vi = vi, mods=~mod*fra+rep, test="knha",slab=surv, data=dataset)
summary(model5)
model6<-rma(yi = kohler, vi = vi, mods=~mod*fra+prog+rep, test="knha",slab=surv, data=dataset)
summary(model6)
model6a<-rma(yi = kohler, vi = vi, mods=~mod*fra+prog+rep+year, test="knha",slab=surv, data=dataset)
summary(model6a)
anova1<-anova(model1,model4, refit=T)
anova1
anova2<-anova(model1,model5, refit=T)
anova2
anova3<-anova(model4,model6, refit=T)
anova3
anova3a<-anova(model6,model6a, refit=T)
anova3a
anova3b<-anova(model1,model6, refit=T)
anova3b
anova(model6, btt="mod")
anova(model6, btt=2:4)

#permutation test of the model 6 (one need to install dev version of metafor package by remotes::install_github("wviechtb/metafor"))
permu1<-permutest(model6)
print(permu1)

#models on data subsets (only Europe, only ISSP, only high income2020)
table(dataset$euro,dataset$mod)
modeleuro<-rma(yi = kohler, vi = vi, mods=~mod*fra+prog+rep, test="knha",subset = dataset$euro==1, data=dataset)
summary(modeleuro)
modelissp<-rma(yi = kohler, vi = vi, mods=~mod*fra+rep, test="knha",subset = dataset$prog=="ISSP", data=dataset)
summary(modelissp)
modelissp2<-rma(yi = kohler, vi = vi, mods=~mod*fra, test="knha",subset = dataset$prog=="ISSP", data=dataset)
summary(modelissp2)
modelinc<-rma(yi = kohler, vi = vi, mods=~mod*fra+prog+rep, test="knha",subset = dataset$inc=="H", data=dataset)
summary(modelinc)

#outliers identification and exclusion (4 surveys out of 795)
qqnorm(nbi, label=5)
qqnorm(model1, label=5)
qqnorm(model6, label=5)
inf0<-influence(nbi)
plot(inf0, layout=c(8,1))
inf1<-influence(model1)
plot(inf1, layout=c(8,1))
inf6<-influence(model6)
plot(inf6, layout=c(8,1))

modelZ<-rma(yi = kohler, vi = vi, test="knha",slab=surv,subset = drop==0,  data=dataset)
summary(modelZ)
modelA<-rma(yi = kohler, vi = vi, mods=~mod*fra, test="knha",slab=surv,subset = drop==0, data=dataset)
summary(modelA)
#Final set of moderators
modelB<-rma(yi = kohler, vi = vi, mods=~mod*fra+prog+rep, test="knha",slab=surv,subset = drop==0, data=dataset)
summary(modelB)
#
modelBa<-rma(yi = kohler, vi = vi, mods=~mod*fra+prog+rep+RR6, test="knha",slab=surv,subset = drop==0, data=dataset)
summary(modelBa)
anovaX<-anova(modelB, modelBa, refit = T)
anovaX


#permutation test of the final set two-level model
permu2<-permutest(modelB)
print(permu2)

#nested mixed-effects meta-regresions models using t-tests and improved method for approximating the degrees of freedom of the t- and F-distributions, based on dataset without 4 outliers
k1 <- rma.mv(yi = kohler,V = vi, random = ~ 1 | wave/surv, dfs="contain", test="t",subset = drop==0, data=dataset)
summary(k1)
k1a<- rma.mv(yi = kohler,V = vi, random = ~ 1 | wave/surv, sigma2 = c(0,NA), dfs="contain", test="t",subset = drop==0, data=dataset)
anovak1<-anova(k1,k1a)
anovak1
k2 <- rma.mv(yi = kohler,V = vi, random = ~ 1 | cntry/surv, dfs="contain", test="t",subset = drop==0, data=dataset)
summary(k2)
k2a <- rma.mv(yi = kohler,V = vi, random = ~ 1 | cntry/surv,sigma2 = c(0,NA), dfs="contain", test="t",subset = drop==0, data=dataset)
summary(k2a)
anovak2<-anova(k2,k2a)
anovak2
k3<-rma.mv(yi = kohler,V = vi, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = drop==0, data=dataset)
summary(k3)
k3a<-rma.mv(yi = kohler,V = vi, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),sigma2 = c(0,0,NA), dfs="contain", test="t",subset = drop==0, data=dataset)
summary(k3a)
k3b<-rma.mv(yi = kohler,V = vi, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),sigma2 = c(0,NA,NA), dfs="contain", test="t",subset = drop==0, data=dataset)
summary(k3b)
anovak31<-anova(k3,k3a)
anovak31
anovak32<-anova(k3,k3b)
anovak32

model7 <- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep, random = ~ 1 | wave/surv, dfs="contain", test="t",subset = drop==0, data=dataset)
summary(model7)
model7a <- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep, random = ~ 1 | wave/surv,sigma2 = c(0,NA),dfs="contain", test="t",subset = drop==0, data=dataset)
anova4<-anova(model7,model7a)
anova4
model8 <- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep, random = ~ 1 | cntry/surv,dfs="contain", test="t",subset = drop==0, data=dataset)
summary(model8)
model8a <- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep, random = ~ 1 | cntry/surv,sigma2 = c(0,NA),dfs="contain", test="t",subset = drop==0, data=dataset)
anova5<-anova(model8,model8a)
anova5

#Model included in paper as Model I
model9b<- rma.mv(yi = kohler,V = vi, mods = ~mod*fra, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = drop==0, data=dataset)
summary(model9b)
anova7<-anova(model9, model9b, refit=T)
anova7
#

#Model included in paper as Model II 
model9<- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = drop==0, data=dataset)
summary(model9)
#
model9a<- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),sigma2 = c(0,NA,NA),dfs="contain", test="t",subset = drop==0, data=dataset)
summary(model9a)
anova6<-anova(model9,model9a)
anova6
#reversed reference level of survey frame 
model9alt<-rma.mv(yi = kohler,V = vi, mods = ~mod*fra2+prog+rep, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = drop==0, data=dataset)
summary(model9alt)

#Model included in paper as Model III
model10<- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep+RR6, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = drop==0, data=dataset)
summary(model10)
anova10<-anova(model9,model10,refit=T)
anova10
#

#nested models on data subsets (only Europe, only ISSP, only high income2020).
nestedeuro<-rma.mv(yi = kohler, V = vi, mods=~mod*fra+prog+rep,random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = (euro==1 & drop==0), data=dataset)
summary(nestedeuro)
nestedissp<-rma.mv(yi = kohler, V = vi, mods=~mod*fra+rep, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = (prog=="ISSP"), data=dataset)
summary(nestedissp)
nestedinc<-rma.mv(yi = kohler, V = vi, mods=~mod*fra+prog+rep, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = (inc=="H" & drop==0), data=dataset)
summary(nestedinc)

# profiles of the restricted log-likelihood. change number of cores, according to parallel::detectCores()
profilemB<-profile(modelB, parallel="multicore", ncpus=8)
profilem7<-profile(model7, parallel="multicore", ncpus=8)
profilem8<-profile(model8, parallel="multicore", ncpus=8)
profilem9<-profile(model9, parallel="multicore", ncpus=8)
plot(profilemB)
plot(profilem7)
plot(profilem8)
plot(profilem9)

#variance distribution in 3-level models
mlm.variance.distribution(x = k1)
mlm.variance.distribution(x = k2)

###r2 - proportional reduction in the total variance
(sum(k1$sigma2) - sum(model7$sigma2)) / sum(k1$sigma2)
(sum(k2$sigma2) - sum(model8$sigma2)) / sum(k2$sigma2)
(sum(k3$sigma2) - sum(model9$sigma2)) / sum(k3$sigma2)
(sum(k3$sigma2) - sum(model9b$sigma2)) / sum(k3$sigma2)
(sum(k3$sigma2) - sum(model10$sigma2)) / sum(k3$sigma2)

### ICC. Calculation of ICC for cross-nested model is problematic
round(k1$sigma2[1] / sum(k1$sigma2), 3)
round(k2$sigma2[1] / sum(k2$sigma2), 3)

#multimodel inference, AICc criterion
dt2 <- dataset[ which(dataset$drop==0), ]

m.inferenceAICc<-  multimodel.inference(TE = "kohler", 
                                      seTE = "sei",
                                      data = dt2,
                                      predictors = c("mod","fra","rep","prog"),
                                      interaction = T,
                                      eval.criterion='AICc')
m.inferenceAICc

#model with RR6 as an effect size
modelRR<- rma.mv(yi = RR6,V = RR_vi, mods = ~mod*fra+prog+rep, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t", data=dataset)
summary(modelRR)
###########END#########