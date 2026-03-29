####### <Addition for manuscript:"Survey mode and nonresponse bias: 
#a meta-analysis based on the data from ISSP 1996–2018 and ESS rounds 1 to 9" PLOS ONE>####
# Copyright (C) <2024>  <Adam Rybak>

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

install.packages("remotes")
remotes::install_github("wviechtb/metafor")
library(metafor)
pacman::p_load(metafor,dplyr,data.table,expss,readr,xlsx)

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
dataset$prog<-factor(dataset$prog, labels=proglevels)
dataset$fra2<-relevel(dataset$fra,ref = "Nonindividual frame")
wavenames<-dataset$wave
regnames<-dataset$region
intnames<-dataset$int
prognames<-dataset$prog
dataset$rep<-dataset$repetitions
var_lab(dataset$rep) = "Repetitions in survey program"
dataset$drop<-ifelse(dataset$surv=="ESS6_SK",1,ifelse(dataset$surv=="ESS5_SK",1,ifelse(dataset$surv=="ESS4_SK",1,ifelse(dataset$surv=="ESS9_LT",1,0))))

dataset <- mutate(dataset, ID = row_number())
duplicates  <- data.table(dataset)
setkeyv(duplicates, c("kohler","cntry"))
duplicates [,dup := .GRP, by = key(duplicates)]
duplicates <- filter(duplicates,kohler > 0)
dupl <- duplicates  %>% 
  group_by(dup) %>% 
  filter(n() > 1) %>% 
  ungroup()  
del <- dupl$ID
set.seed(1234)
dupl2 <- dupl %>% 
  group_by(dup) %>% 
  sample_n(1) %>% 
  ungroup()  
data_dupl2 <- dupl2[,-31]
data_dupl1 <- dataset[!(dataset$ID %in% del),]
data_corrected <- rbind(data_dupl1,data_dupl2)


model9b<- rma.mv(yi = kohler,V = vi, mods = ~mod*fra, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = drop==0, data=data_corrected)
summary(model9b)
model9<- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = drop==0, data=data_corrected)
summary(model9)
model10<- rma.mv(yi = kohler,V = vi, mods = ~mod*fra+prog+rep+RR6, random = list( ~1 | wave,~1 | cntry,  ~1 | surv),dfs="contain", test="t",subset = drop==0, data=data_corrected)
summary(model10)