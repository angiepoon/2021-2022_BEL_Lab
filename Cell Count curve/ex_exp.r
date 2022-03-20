library(ggplot2)
library(scales)
setwd("/Users/angiepoon528/Desktop/FYP/Expansion experiment ")
#tutorial: https://stats.oarc.ucla.edu/r/faq/how-can-i-make-individual-growth-curves-in-ggplot2/

#Cell count data 
Day0 <- c(50000, 50000, 425000, 425000)
Day6 <- c(80000, 102500, 575000, 452000)
Day12 <- c(62500, 62500, 395000,410000)

Day6 <- Day6*0.5 #cells are more concentrated coz volume roughly like 500 

# set up data frame 
ID <- c(rep(c(1:4), 3))
Day <- c(rep(0,4), rep(6, 4), rep(12,4))
Count <- c(Day0, Day6, Day12)
Treatment <- c(rep(c("Pos", "Neg"),6))
Batch <- c(rep(c(1,1,2,2),3))
DF <- data.frame(ID, Batch, Treatment, Day, Count)

# change ID, Treatment, Batch to factor variables 
DF <- within(DF, {
  ID <- factor(ID)
  Treatment <- factor(Treatment, levels = c("Pos", "Neg"), labels = c("Pos", "Neg"))
  Batch <- factor(Batch)
  Day <- factor(Day)
})

# Individual plots 
ggplot(DF, aes(x= Day, y= Count))+geom_line()+
  facet_wrap(~ID)

# linear 
ggplot(DF, aes(x= Day, y= Count))+geom_point()+
  stat_smooth(method = "lm", se = F)+facet_wrap(~ID)

# Treatment Nested, + batch 
#plot
p<- ggplot(DF, aes(x= Day, y= Count, group= ID, color = Batch))+geom_line()+
  facet_wrap(~Treatment)
require(scales)
p + scale_y_continuous(labels = comma)+scale_x_discrete(expand = expansion(mult = 0.07))+
  theme_bw()

ggplot(DF, aes(x= Day, y= Count, group= ID, color = Batch))+geom_point()+
  stat_smooth(method = "lm", se = F) + facet_wrap(~Treatment) + 
  scale_y_continuous(labels = comma)+scale_x_discrete(expand = expansion(mult = 0.07))+
  theme_bw()

