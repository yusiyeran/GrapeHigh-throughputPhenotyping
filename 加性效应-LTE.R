library(rTASSEL)
library(readxl)
library(ggplot2)
library(tidyverse)
library(multcompView)
library(ggsci)
wd = "G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\10.图-加性效应"
setwd(wd)


#######promoter haplotype#####
library(dplyr)
library(stringr)
hapfile = "G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\8.newChrFigure\\20230724-new\\haplotype_V3.0\\v2haplotypes.txt"
all_ano = read.table(hapfile,header = F, #F 为后面字符型做准备
                     sep = "\t")
str(all_ano)

all_ano = all_ano[,-6]
colnames(all_ano) = all_ano[1,]
str(all_ano)
all_ano = all_ano[-1,]
str(all_ano)
all_ano$type = str_split_fixed(all_ano$INFO," :: ",2)[,1]
all_ano$gene = str_split_fixed(all_ano$INFO," :: ",2)[,2]

all_ano_select = all_ano[all_ano$type%in%c("exonic","UTR3","upstream","UTR5"),]
all_ano_select$gene = str_split_fixed(all_ano_select$gene," :: ",2)[,1]
unique(all_ano_select$gene)

markerfile = "G:\\OneDrive\\MS2023_1023\\MS2023-1122\\MS2023-1122-修回\\Additional file_1011.xlsx"
SNPs = read_xlsx(markerfile,sheet = 'Table S14',skip = 3) 
SNPs = as.data.frame(SNPs)
str(SNPs)
unique(SNPs$Trait)
# TotalAcidSNP = SNPs[SNPs$Trait == "total acid",]
# TotalAcidSNP = SNPs[SNPs$Trait == "total sugar",]
TotalAcidSNP = SNPs[SNPs$Trait == "LTE",]
TotalAcidSNP = TotalAcidSNP[TotalAcidSNP$Chromsome%in%c(3,7),]
str(TotalAcidSNP)
 
pheno.tangsuan = read.table("表型文件\\tangsuan.txt",header = T,na.strings = "NA")
pheno.tangsuan = read.table("表型文件\\tangsuan.txt",header = T,na.strings = "NA")
pheno.LTE = read_excel("表型文件\\LTE_new.xlsx",na = "NA")
pheno.LTE = as.data.frame(pheno.LTE)
str(pheno.LTE)
pheno.LTE$'2012' = as.numeric(pheno.LTE$'2012' )
pheno.LTE$'2013' = as.numeric(pheno.LTE$'2013' )
pheno.LTE$'2014' = as.numeric(pheno.LTE$'2014' )
pheno.LTE$'2015' = as.numeric(pheno.LTE$'2015' )
pheno.LTE$'2016' = as.numeric(pheno.LTE$'2016' )
pheno.LTE$'2019' = as.numeric(pheno.LTE$'2019' )
# pheno.tangsuan[is.na(pheno.tangsuan)] =0
# pheno.tangsuan$TotalSugar2011 = pheno.tangsuan$Sucrose2011+pheno.tangsuan$Glucose2011+pheno.tangsuan$Fructose2011
# pheno.tangsuan$TotalSugar2012 = pheno.tangsuan$Sucrose2012+pheno.tangsuan$Glucose2012+pheno.tangsuan$Fructose2012
# pheno.tangsuan$TotalSugar2013 = pheno.tangsuan$Sucrose2013+pheno.tangsuan$Glucose2013+pheno.tangsuan$Fructose2013

# TotalAcidSNP2011 = TotalAcidSNP[TotalAcidSNP$Year == 2011,]
# TotalAcidSNP2012 = TotalAcidSNP[TotalAcidSNP$Year == 2012,]
# TotalAcidSNP2013 = TotalAcidSNP[TotalAcidSNP$Year == 2013,]
# TotalAcidSNP2011 =TotalAcidSNP2013
TotalAcidSNP2011 =TotalAcidSNP
TotalAcidSNP$Marker = paste(TotalAcidSNP$Chr,TotalAcidSNP$Position,sep = "-")
all_ano_select$Marker = paste(all_ano_select$CHROM,all_ano_select$POS,sep = "-")
TotalAcidSNP_group = group_by(TotalAcidSNP,Chromsome) %>% summarise_each(funs(max))
# pheno_TotalAcid = pheno.tangsuan[,c(1,19)]
# pheno_TotalAcid = pheno.tangsuan[,c(1,17:19)] #17:TotalSugar2011 18:TotalSugar2012 19:TotalSugar2013
# pheno = pheno_TotalAcid
pheno = pheno.LTE[,c(1,5)]
pheno = pheno.LTE
# pheno = pheno.shapesize[,1:2]
names(pheno)[1] = "Genotype3"


############################################column plot + significant ###############
pdf("add_trait distribution_LTE_V15.0.pdf",width = 10,height = 6)
for (ii in TotalAcidSNP$Marker[c(93,16)]) { #TotalAcidSNP2011$Marker
  # # ii = "6−5939415"
  # "plotdata\\2014_3-5406949_1.csv"
  # "plotdata\\2014_7-4395972_1.csv"
  
  hap_ii = all_ano_select[all_ano_select$Marker == ii,c(476,477,1:4,6:475)]
  if (dim(hap_ii)[1]>0){
    t_hap_ii =as.data.frame(t(hap_ii))
    for (mu_i in names(t_hap_ii)) {
      t_hap_ii %>% mutate(mu_i = recode(mu_i,'/' = 'NA'))
    }
    
    t_hap_ii = na.omit(t_hap_ii)
    #paste hap
    xx = t_hap_ii[,1]
    if (dim(t_hap_ii)[2] >1 ){
      for (paste_i in 2:dim(t_hap_ii)[2]) {
        xx = paste(xx,t_hap_ii[,paste_i],sep = "_")
        
      }
    }
    t_hap_ii$hap = xx
    t_hap_ii0 = t_hap_ii[-6:-1,]
    
    t_hap_ii0$Genotype =rownames(t_hap_ii0)
    
    t_hap_ii1 = merge(pheno,t_hap_ii0,by.x = "Genotype3",by.y = "Genotype",all.x = T) ###***此处Genotype3 根据你的实际命名更改
    t_hap_ii1 = na.omit(t_hap_ii1)
    hap_num = group_by(t_hap_ii1,hap) %>% summarise_each(funs(length))
    hap_num = as.data.frame(hap_num)
    hap_num1 = hap_num[,1:2]
    names(hap_num1)[2] = "num"
    
    hap_num1$marker_num = dim(t_hap_ii)[2]-1
    hap_num1$gene = ii
    
    
    traits =names(pheno)[-1]
    
    for (trait_i in traits) {
      
      
      pheno_data0 = pheno[,c("Genotype3",trait_i)]
      
      hap_pheno = merge(pheno_data0,t_hap_ii0,by.x = "Genotype3",by.y = "Genotype",all.x = T) ###***此处Genotype3 根据你的实际命名更改
      
      
      
      hap_pheno = na.omit(hap_pheno)
      
      #plot
      
      data_plot = hap_pheno
      
      
      data_plot_summmary1 = group_by(data_plot,hap)%>%summarise_each(funs(mean))
      data_plot_summmary1 = as.data.frame(data_plot_summmary1)
      
      data_plot_summmary2 = group_by(data_plot,hap)%>%summarise_each(funs(length))
      data_plot_summmary2 = as.data.frame(data_plot_summmary2)
      
      names(data_plot_summmary2)[2] = "Genotype"
      
      data_plot_summmary = merge(data_plot_summmary1[,c("hap",trait_i)],data_plot_summmary2[,c("hap","Genotype")],by = "hap")
      
      #####删除个体数量少于5%的单倍型个体
      # data_plot_summmary = data_plot_summmary[data_plot_summmary$Genotype>(dim(data_plot)[1]*0.05),]
      
      data_plot_summmary$name = paste("Hap",1:dim(data_plot_summmary)[1],sep = "")
      data_plot_summmary$name2 = paste(data_plot_summmary$name,"\n (n = ",data_plot_summmary$Genotype,")",sep = "")
      
      level = data_plot_summmary[order(data_plot_summmary[,trait_i]),]
      data_plot$hap = factor(data_plot$hap,levels = level$hap)
      
      
      
      
      ############anova########
      
      if (dim(data_plot_summmary)[1]>1){
        anova <- aov(data_plot[,trait_i]~hap,data=data_plot)
        tukey <- TukeyHSD(anova)
        cld <- multcompLetters4(anova,tukey)
        aa = as.data.frame.list(cld$hap) %>% select(1) %>% rownames_to_column("hap")
        
        data_plot_summmary3 = group_by(data_plot,hap)%>%summarise_each(funs(max))
        data_plot_summmary3 = as.data.frame(data_plot_summmary3)
        
        anova_letter0 = merge(data_plot_summmary3,aa,by = "hap",all.x = T)
        anova_letter = merge(anova_letter0,data_plot_summmary[,c("hap","name2")],by = "hap")
        
        if (length(unique(anova_letter$Letters))>0){
          ################
          data_plot = na.omit(data_plot)
          data_plot0 = merge(data_plot,data_plot_summmary[,c("hap","name2")],by = "hap")
          
          data_plot0$name2 = factor(data_plot0$name2,levels = level$name2 )
          #plot
          p = ggplot( data = data_plot0,aes(x=name2, y=data_plot0[,trait_i], fill=name2)) +
            geom_boxplot() +
            geom_text(data =anova_letter,aes(label = Letters, x = name2,y = anova_letter[,trait_i]), size = 5,vjust = -0.5)+
            
            ylim(min(data_plot0[,trait_i]),max(data_plot0[,trait_i])*1.05)+
            theme_bw(base_size = 25) +
            scale_fill_aaas()+
            theme(
              legend.position="none",
            ) +
            
            ylab(trait_i)+
            ggtitle(ii) +
            xlab("")
          # for (iix in filenames[25:26]) {
          #   data_iix = read.csv(paste("plotdata\\",iix,sep = ""))
          #   data_iix1 = merge(data_plot0,data_iix,by = "Genotype3")
          #   p1 = p+geom_jitter(data = data_iix1,aes(x=name2.x,y=TotalAcid2013.x ,fill =name2.y ,colour = name2.y,size=1,alpha=0.5))+
          #     labs(subtitle =strsplit(iix,"_")[[1]][2])
          # }
          # 
          # print(p1)
          print(p)
        }
        write.csv(data_plot0,paste(trait_i,ii,"1.csv",sep = "_"))
      }
    }
    
  }
  
}

dev.off()



#################################################################
# filenames = list.files("G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\10.图-加性效应\\plotdata-LTE1")
# chr1 = read.csv("plotdata\\TotalSugar2011_1-20723744_1.csv")
# chr11 = read.csv("plotdata\\TotalSugar2011_11-4841924_1.csv")
# chr18 = read.csv("plotdata\\TotalSugar2011_18-7426816_1.csv")
setwd("G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\10.图-加性效应\\plotdata-LTE2")

pdf("add_trait distribution_LTE_V7.0.pdf",width = 8,height = 6)
for (i_year in c(2012:2016,2019)) {
  

chr2 = read.csv(paste(i_year,"_3-5406949_1.csv",sep = ""))
chr7 = read.csv(paste(i_year,"_7-4305464_1.csv",sep = "")) #7-4404321


# chr1_group=group_by(chr1, hap) %>% summarize_each(funs(mean))
# chr11_group=group_by(chr11, hap) %>% summarize_each(funs(mean))
# chr18_group=group_by(chr18, hap) %>% summarize_each(funs(mean))
# 
# chr1_group2=group_by(chr1, hap) %>% summarize_each(funs(length))
# chr11_group2=group_by(chr11, hap) %>% summarize_each(funs(length))
# chr18_group2=group_by(chr18, hap) %>% summarize_each(funs(length))

chr2_group=group_by(chr2, hap) %>% summarize_each(funs(mean))
chr7_group=group_by(chr7, hap) %>% summarize_each(funs(mean))

chr2_group2=group_by(chr2, hap) %>% summarize_each(funs(length))
chr7_group2=group_by(chr7, hap) %>% summarize_each(funs(length))


group = c()
for (iii in 1:dim(chr2)[1]) {
  if ( chr2[iii,"hap"] %in%c("G/G","A/G") ){
    group = c(group,"LTE1")
  }else{
    group = c(group,"NA")
  }
}
chr2$group = group

group = c()
for (iii in 1:dim(chr7)[1]) {
  if ( chr7[iii,"hap"] %in% c("A/G") ){
    group = c(group,"LTE2")
  }else{
    group = c(group,"NA")
  }
}
chr7$group = group


chr2_7 = merge(chr2,chr7,by="Genotype3")

chr2_7[,4] = chr2_7[,4]*-1
# view(chr2_7)
# chr1_11_18[(chr1_11_18)=="NA"] = NA
# aaa=c()
# for (i in 1:dim(chr1_11_18)[1]) {
#   aa = chr1_11_18[i,"group.x"]
#   bb = chr1_11_18[i,"group.y"]
#   cc = chr1_11_18[i,"group.z"]
#   xx = c(aa,bb,cc)
#   xx0 = xx[is.na(xx)==F]
#   if (length(xx0)==3){
#   xxx = paste(xx0[1],xx0[2],xx0[3],sep = "_")
#   aaa = c(aaa,xxx)
#   }
#   if (length(xx0)==2){
#     xxx = paste(xx0[1],xx0[2],sep = "_")
#     aaa = c(aaa,xxx)
#   }
#   if (length(xx0)==1){
#     xxx =xx0[1]
#     aaa = c(aaa,xxx)
#   }
#   if (length(xx0)==0){
#     xxx ="NA"
#     aaa = c(aaa,xxx)
#   }
# }
# # chr1_11_18$group = aaa
# chr1_11_18_group=group_by(chr1_11_18, group) %>% summarize_each(funs(mean))
# chr1_11_18_count=group_by(chr1_11_18, group) %>% summarize_each(funs(length))
# chr1_11_18_group=as.data.frame(chr1_11_18_group)
tem_group = c()
for (i_group in 1:dim(chr2_7)[1]) {
  if (chr2_7[i_group,"group.x"] == "NA"){
    if (chr2_7[i_group,"group.y"] == "NA"){
      tem_group = c(tem_group,"NA")
    }else{
      tem_group = c(tem_group,chr2_7[i_group,"group.y"])
    }
  }else{
    if (chr2_7[i_group,"group.y"] == "NA"){
      tem_group = c(tem_group,chr2_7[i_group,"group.x"])
    }else{
      tem_group = c(tem_group,paste(chr2_7[i_group,"group.x"],chr2_7[i_group,"group.y"],sep = "_"))
    }
  }
}
chr2_7$group = tem_group
# chr2_7$group = paste(chr2_7$group.x,chr2_7$group.y,sep = "_")
chr2_7_group0=group_by(chr2_7, group) %>% summarize_each(funs(mean))
chr2_7_count=group_by(chr2_7, group) %>% summarize_each(funs(length))
chr2_7_max0=group_by(chr2_7, group) %>% summarize_each(funs(max))
chr2_7_group0=as.data.frame(chr2_7_group0)

chr2_7_group = chr2_7_group0[order(chr2_7_group0[,5]),]

chr2_7$group = factor(chr2_7$group,levels = chr2_7_group$group)
chr2_7_max =merge(chr2_7_max0,chr2_7_group,by="group")


############anova########
data_plot = chr2_7

anova <- aov(data_plot[,4]~group,data=data_plot)
tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova,tukey)
anova_letter = as.data.frame(cld$group$Letters)
anova_letter$name = rownames(anova_letter)
anova_letter$letter = anova_letter$`cld$group$Letters`
chr2_7_max=as.data.frame(chr2_7_max)
# chr2_7_group=as.data.frame(chr2_7_group)
anova_letter = merge(anova_letter,chr2_7_max,by.x = "name",by.y = "group")
anova_letter$yy = anova_letter[,7]
anova_letter =anova_letter[order(anova_letter[,26],decreasing = T),]

parent1 = chr2_7[chr2_7$Genotype3 %in%c("bei-hong"),]
parent2 = chr2_7[chr2_7$Genotype3 %in%c("ES7-11-49"),]
p <- ggplot(data = chr2_7,aes(x=group, y=chr2_7[,4]),na.rm = T) + geom_violin()+
  geom_jitter(shape=16,alpha =0.5, position=position_jitter(0.1))+#ylim(0,35)+
  stat_summary(fun= mean, geom="point", size=2, color="red",shape=15,alpha =1,na.rm = T)+
  geom_point(data = parent1,aes(x=group, y=parent1[,4]),na.rm = T,color="purple",shape=15,size=2)+
  geom_point(data = parent2,aes(x=group, y=parent2[,4]),na.rm = T,color="green",shape=15,size = 2)+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 15,colour = "black", family="serif"))+
  geom_text(data =anova_letter,aes(label = letter, x = name,y = yy), size = 5,vjust = -0.5)+
  theme_bw(base_size = 20)+xlab(NULL)+
  theme(axis.text.x = element_text(color = "black"), #angle = 45, hjust = 1,
        axis.text.y = element_text(color = "black")
        )+
  ylab(paste("LTE",i_year,sep = " "))
print(p)
}
dev.off()





write.csv(chr1_11_18,"TotalSugar_chr1_11_18.csv")
write.csv(chr6_17,"TotalAcid_chr6_17.csv")
