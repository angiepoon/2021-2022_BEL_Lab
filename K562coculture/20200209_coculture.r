setwd("/Users/angiepoon528/Desktop/FYP/Coculture /20200209_coculture")
library(ggplot2)

# circular 
# hEPSC-NK
ID <- c(1:4)
Percentage_dying <- c((0.56-0.021),(1.69-0.24), (1.81-0.24), (0.76-0.24))
E_T_ratio <- c(0, 1, 0.33, 0.1)

NK <- data.frame(ID, E_T_ratio ,Percentage_dying)

#line graph 
ggplot(data = NK, aes(x= E_T_ratio, y= Percentage_dying, group=1))+
  geom_line()+
  geom_point()+ 
  xlab("E:T ratio")+ylab("Percentage of dying cells(%)")+
  theme_bw()+
  scale_y_continuous(expand= c(0.1,0.1))+
  theme(text = element_text(size=20), 
        panel.grid.major.x = element_blank())

#NK92mi 
ID_2 <- c(1:4)
Percentage_dying_2 <- c((0.56-0.021),(21.4-0.16), (9.49-0.16), (3.81-0.16))
E_T_ratio_2 <- c(0, 1, 0.33, 0.1)
NK92mi <- data.frame(ID_2, E_T_ratio_2 ,Percentage_dying_2)

#line graph 
ggplot(data = NK92mi, aes(x= E_T_ratio_2, y= Percentage_dying_2, group=1))+
  geom_line()+
  geom_point()+ 
  xlab("E:T ratio")+ylab("Percentage of dying cells(%)")+
  theme_bw()+
  scale_y_continuous(expand= c(0.1,0.1))+
  theme(text = element_text(size=20), 
        panel.grid.major.x = element_blank())

# Total gating 
# hEPSC-NK 
ID_3 <- c(1:4)
Percentage_dying_3 <- c((0.56-0.021), (6.01-0.59), (3.63-0.59),(2.32-0.59) )
E_T_ratio_3 <- c(0, 1, 0.33, 0.1)
NK_3 <- data.frame(ID_3, E_T_ratio_3 ,Percentage_dying_3)
ggplot(data = NK_3, aes(x= E_T_ratio_3, y= Percentage_dying_3, group=1))+
  geom_line()+
  geom_point()+ 
  xlab("E:T ratio")+ylab("Percentage of dying cells(%)")+
  theme_bw()+
  scale_y_continuous(expand= c(0.1,0.1))+
  theme(text = element_text(size=20), 
        panel.grid.major.x = element_blank())

# NK92
ID_4 <- c(1:4)
Percentage_dying_4 <- c((0.56-0.021), (37.6-0.24), (24.7-0.24),(7.41-0.24) )
E_T_ratio_4 <- c(0, 1, 0.33, 0.1)
NK92mi_4 <- data.frame(ID_4, E_T_ratio_4 ,Percentage_dying_4)

ggplot(data = NK92mi_4, aes(x= E_T_ratio_4, y= Percentage_dying_4, group=1))+
  geom_line()+
  geom_point()+ 
  xlab("E:T ratio")+ylab("Percentage of dying cells(%)")+
  theme_bw()+
  scale_y_continuous(expand= c(0.1,0.1))+
  theme(text = element_text(size=20), 
        panel.grid.major.x = element_blank())

## Integrating line graph of hEPSC-NK and NK92mi 
### Total gating ones are used 
ID_5 <-c(1:8)
Percentage_dying_5 <- c(Percentage_dying_3, Percentage_dying_4)
E_T_ratio_5 <- c(E_T_ratio_3, E_T_ratio_4)
Group <- c(rep("hEPSC-NK", 4), rep("NK92mi", 4))

Data_sheet <- data.frame(ID, Group, E_T_ratio_5, Percentage_dying_5)

ggplot(Data_sheet, aes(x= E_T_ratio_5, y= Percentage_dying_5, color = Group))+
  geom_line()+
  geom_point()+
  scale_x_continuous("Effector: Target")+
  ylab("Percentage of dying cells(%)")+
  theme_classic()+
  scale_y_continuous(expand= c(0.1,0.1))+
  theme(text = element_text(size=20), 
        panel.grid = element_blank())


