setwd("/Users/angiepoon528/Desktop/FYP/Experiments/matrix_experiment")
library(ggplot2)
library(dplyr)

data<- read.csv("admin_2022-04-01 13-09-30_CT035437-Angie-2.csv", stringsAsFactors = T)
data<- data[-c(1:19),] #remove rows for details 
colnames(data) <- c("Well", "Fluor", "Target", "Content", "Sample", "Cq", "Start.Quantity", "Mean.Cq", "Dev", "Melting.Temp")

FAM<- subset(data, data$Fluor == "FAM")  #SYB green use FAM fluorophore
FAM <- FAM[, c(1,2,6,8,10)]
FAM[FAM=="NaN"] <- NA #convert "NA" Strings into real NA 
FAM<- na.omit(FAM) 
Treatment <- c(rep(c("No", "IL18", "bIL18", "water"),2))
Gene <- c(rep("IFNg", 4), rep("GAPDH",4))
FAM$Cq<- as.numeric(FAM$Cq)

#Add new column on FAM to specift condition 
FAM<- mutate(FAM,
       Treatment = Treatment,
       Gene = Gene)

# Do this when more datas 
ggplot(FAM, aes(x=factor(Gene), y = Cq, color = factor(Treatment)))+
  geom_boxplot()

#No stacked bars 
FAM%>%
ggplot(aes(x=factor(Gene), y = as.numeric(Mean.Cq), fill= factor(Treatment)))+
  geom_col(width=0.5, position = position_dodge(0.6)) + 
  scale_y_continuous(expand=c(0,0), name= "Ct value")+
  scale_x_discrete("Gene")+
  theme_classic()


