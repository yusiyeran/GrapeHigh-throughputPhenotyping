wd = "G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\10.图-加性效应"
setwd(wd)

SugarAcid = read_excel("SugarAcid二维分析.xlsx",sheet = "data4R")
SugarAcid = as.data.frame(SugarAcid)
SugarAcid$group_acid = factor(SugarAcid$group_acid,levels =  rev(c("TA1_TA2","TA1","TA2","NA")))
SugarAcid$group_sugar = factor(SugarAcid$group_sugar,levels = rev(c("TS1_TS2_TS3","TS1_TS3","TS1_TS2", "TS2_TS3","TS2","TS1","TS3","NA")))
# ##add label
# label_acid = data.frame(Type = c("TA1_TA2","TA1","TA2","NA"),AcidCol = c("blue4","green4","cyan4","lightcyan"))
# 
# label_sugar = data.frame(Type = c("TS1_TS2_TS3","TS1_TS3","TS1_TS2", "TS2_TS3","TS2","TS1","TS3","NA"),SugarNum = 8:1)
# SugarAcid_plot = merge(SugarAcid,label_acid,by.x ="group_acid",by.y= "Type",all.x = T)
# SugarAcid_plot = merge(SugarAcid_plot,label_sugar,by.x ="group_sugar",by.y= "Type",all.x = T)

SugarAcid$TotalAcid2011_percent = SugarAcid$TotalAcid2011/mean(SugarAcid$TotalAcid2011)*100-100
SugarAcid$sugar_2011_percent = SugarAcid$sugar_2011/mean(SugarAcid$sugar_2011,na.rm = T)*100-100
library(ggrepel)
ggplot(data = SugarAcid, aes(x=sugar_2011,y=TotalAcid2011 ))+
  geom_point(aes(size=group_sugar,color=group_acid,alpha=0.6))+
  xlab("Total sugar 2011")+ ylab("Total acid 2011")+
  geom_label_repel(aes(label = group_sugar ), size = 3)+
  
  theme_bw(base_size =20)

pdf("糖酸加性效应二维_V4.0.pdf",width = 9,height = 6)
p=ggplot(data = SugarAcid, aes(x=sugar_2011_percent,y=TotalAcid2011_percent ))+
  geom_point(aes(size=group_sugar,color=group_acid,alpha=0.6))+
  xlab("Variation in the proportion of total sugar (%)")+ ylab("Variation in the proportion of total acid (%)")+
  geom_label_repel(aes(label = group_sugar ), size = 3)+
  geom_vline(xintercept = 0,color="red",alpha=0.5,linetype = 2)+
  geom_hline(yintercept = 0,color="red",alpha=0.5,linetype = 2)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
  )+
  theme_bw(base_size = 18)
print(p)
dev.off()
