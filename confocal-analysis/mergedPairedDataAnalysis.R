library(tidyverse)
library(readxl)

###########data processing###############
# load data and set some as factors
filepath <- "pathTomergedPairedData.csv"
data <- read.csv(file = filepath) 
data <- data %>% mutate(Condition = recode(Condition,"pCaMKII" = "Ex-Ex", "PV" = "Ex-PV"))

cols <- c("Stain_post","Stain_pre","Week","Condition","ImageName")
data[cols] <- lapply(data[cols], factor)
data <- data %>% mutate(Condition = fct_relevel(Condition,"Ex-Ex","Ex-PV"))

# Normalize the intensity measures to the mean of each week/stain combination of Ex->Ex
# split into pre and post stains to make this a little easier
datapost <- data[c(1:10,20:22)]
datapost <- datapost %>% group_by(Week,Stain_post) %>% 
  mutate(normMeanIntensity_post = Mean_post/mean(Mean_post[Condition=="Ex-Ex"]), 
         normRawIntDen_post = RawIntDen_post/mean(RawIntDen_post[Condition == "Ex-Ex"]))

datapre <- data[c(11:22)]
datapre <- datapre %>% group_by(Week,Stain_pre) %>% 
  mutate(normMeanIntensity_pre = Mean_pre/mean(Mean_pre[Condition=="Ex-Ex"]), 
         normRawIntDen_pre = RawIntDen_pre/mean(RawIntDen_pre[Condition == "Ex-Ex"]))

# recombine df
alldata <- cbind(datapost,datapre[c(1:9,13:14)])

# make per cell averages
meansByCell <- alldata %>% 
  group_by(Condition,ImageName,Stain_post,Stain_pre) %>% 
  summarise_at(c("Area_post","normMeanIntensity_post","normRawIntDen_post","Circ__post","Area_pre","normMeanIntensity_pre","normRawIntDen_pre","Circ__pre"),mean)


alldata<-alldata %>% arrange(Condition)
savepathSynapses <- "pathToSaveAlldata.csv"
savepathCells <- "pathToSaveByCellData.csv"
write.csv(alldata,savepathSynapses,row.names = FALSE)
write.csv(meansByCell,savepathCells,row.names = FALSE)
