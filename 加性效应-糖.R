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
TotalAcidSNP = SNPs[SNPs$Trait == "total acid",]
TotalAcidSNP = SNPs[SNPs$Trait == "total sugar",]
str(TotalAcidSNP)
 
pheno.tangsuan = read.table("表型文件\\tangsuan.txt",header = T,na.strings = "NA")
pheno.tangsuan[is.na(pheno.tangsuan)] =0
pheno.tangsuan$TotalSugar2011 = pheno.tangsuan$Sucrose2011+pheno.tangsuan$Glucose2011+pheno.tangsuan$Fructose2011
pheno.tangsuan$TotalSugar2012 = pheno.tangsuan$Sucrose2012+pheno.tangsuan$Glucose2012+pheno.tangsuan$Fructose2012
pheno.tangsuan$TotalSugar2013 = pheno.tangsuan$Sucrose2013+pheno.tangsuan$Glucose2013+pheno.tangsuan$Fructose2013

# TotalAcidSNP2011 = TotalAcidSNP[TotalAcidSNP$Year == 2011,]
# TotalAcidSNP2012 = TotalAcidSNP[TotalAcidSNP$Year == 2012,]
# TotalAcidSNP2013 = TotalAcidSNP[TotalAcidSNP$Year == 2013,]
TotalAcidSNP2011 =TotalAcidSNP
TotalAcidSNP$Marker = paste(TotalAcidSNP$Chr,TotalAcidSNP$Position,sep = "-")
all_ano_select$Marker = paste(all_ano_select$CHROM,all_ano_select$POS,sep = "-")
TotalAcidSNP_group = group_by(TotalAcidSNP,Chromsome) %>% summarise_each(funs(max))
# pheno_TotalAcid = pheno.tangsuan[,c(1,19)]
pheno_TotalAcid = pheno.tangsuan[,c(1,17:19)] #17:TotalSugar2011 18:TotalSugar2012 19:TotalSugar2013
pheno = pheno_TotalAcid
# pheno = pheno.shapesize[,1:2]
names(pheno)[1] = "Genotype3"


############################################
pdf("add_trait distribution_totalsuagr_V14.0.pdf",width = 10,height = 6)
for (ii in c("1-20723744","11-4841924","18-7426816")) { #TotalAcidSNP2011$Marker #TotalAcidSNP$Marker
  # ii = "6−5939415"
  # 1-20723744
  # 11-4841924
  # 18-7426816
  
  
  
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
      
      
      data_plot_summmary = data_plot_summmary[data_plot_summmary$Genotype>(dim(data_plot)[1]*0.05),]
      
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
filenames = list.files("G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\10.图-加性效应\\plotdata-old3")

setwd("G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\10.图-加性效应\\plotdata-TotalSugar1")

pdf("add_trait distribution_TotalSugar_V9.0.pdf",width = 8,height = 6)
for (i_year in c(2011:2013)) {
  
chr1 = read.csv(paste("TotalSugar",i_year,"_1-20723744_1.csv",sep = ""))
chr11 = read.csv(paste("TotalSugar",i_year,"_11-4841924_1.csv",sep = ""))
chr18 = read.csv(paste("TotalSugar",i_year,"_18-7426816_1.csv",sep = ""))


group = c()
for (iii in 1:dim(chr1)[1]) {
  if ( chr1[iii,"hap"] == "G/G" ){
    group = c(group,"TS1")
  }else{
    group = c(group,"NA")
  }
}
chr1$group = group

group = c()
for (iii in 1:dim(chr11)[1]) {
  if ( chr11[iii,"hap"] == "A/A" ){
    group = c(group,"TS2")
  }else{
    group = c(group,"NA")
  }
}
chr11$group = group

group = c()
for (iii in 1:dim(chr18)[1]) {
  if ( chr18[iii,"hap"] == "A/A" ){
    group = c(group,"TS3")
  }else{
    group = c(group,"NA")
  }
}
chr18$group = group

chr1_group=group_by(chr1, hap) %>% summarize_each(funs(mean))
chr11_group=group_by(chr11, hap) %>% summarize_each(funs(mean))
chr18_group=group_by(chr18, hap) %>% summarize_each(funs(mean))

chr1_group2=group_by(chr1, hap) %>% summarize_each(funs(length))
chr11_group2=group_by(chr11, hap) %>% summarize_each(funs(length))
chr18_group2=group_by(chr18, hap) %>% summarize_each(funs(length))


chr1_11 = merge(chr1,chr11,by="Genotype3")
chr1_11_18 = merge(chr1_11,chr18,by="Genotype3")
names(chr1_11_18)[19] = "group.z"
chr1_11_18[(chr1_11_18)=="NA"] = NA
if (length(which(chr1_11_18[,4]==0))>0){
chr1_11_18 =chr1_11_18[-which(chr1_11_18[,4]==0),]
}
aaa=c()
for (i in 1:dim(chr1_11_18)[1]) {
  aa = chr1_11_18[i,"group.x"]
  bb = chr1_11_18[i,"group.y"]
  cc = chr1_11_18[i,"group.z"]
  xx = c(aa,bb,cc)
  xx0 = xx[is.na(xx)==F]
  if (length(xx0)==3){
  xxx = paste(xx0[1],xx0[2],xx0[3],sep = "_")
  aaa = c(aaa,xxx)
  }
  if (length(xx0)==2){
    xxx = paste(xx0[1],xx0[2],sep = "_")
    aaa = c(aaa,xxx)
  }
  if (length(xx0)==1){
    xxx =xx0[1]
    aaa = c(aaa,xxx)
  }
  if (length(xx0)==0){
    xxx ="NA"
    aaa = c(aaa,xxx)
  }
}
chr1_11_18$group = aaa
chr1_11_18_group0=group_by(chr1_11_18, group) %>% summarize_each(funs(mean))
chr1_11_18_max0=group_by(chr1_11_18, group) %>% summarize_each(funs(max))
chr1_11_18_count=group_by(chr1_11_18, group) %>% summarize_each(funs(length))
chr1_11_18_group1=as.data.frame(chr1_11_18_group0)

chr1_11_18_group = chr1_11_18_group1[order(chr1_11_18_group1[,5],decreasing = T),]

chr1_11_18$group = factor(chr1_11_18$group,levels = chr1_11_18_group$group)
chr1_11_18_max = merge(chr1_11_18_max0,chr1_11_18_group,by="group")

############anova########
data_plot = chr1_11_18

anova <- aov(data_plot[,4]~group,data=data_plot)
tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova,tukey)
anova_letter = as.data.frame(cld$group$Letters)
anova_letter$name = rownames(anova_letter)
anova_letter$letter = anova_letter$`cld$group$Letters`
chr1_11_18_max=as.data.frame(chr1_11_18_max)
# chr2_7_group=as.data.frame(chr2_7_group)
# anova_letter = merge(anova_letter,chr1_11_18_max)
# anova_letter$yy = chr1_11_18_max[order(chr1_11_18_max[,30],decreasing = T),5]

anova_letter = merge(anova_letter,chr1_11_18_max,by.x = "name",by.y = "group")
anova_letter$yy = anova_letter[,7]
#names(anova_letter)
anova_letter =anova_letter[order(anova_letter[,26],decreasing = T),]

parent1 = chr1_11_18[chr1_11_18$Genotype3 %in%c("bei-hong"),]
parent2 = chr1_11_18[chr1_11_18$Genotype3 %in%c("ES7-11-49"),]
p <- ggplot(data = chr1_11_18,aes(x=group, y=chr1_11_18[,4]),na.rm = T) + geom_violin()+
  geom_jitter(shape=16,alpha =0.5, position=position_jitter(0.1))+#ylim(0,35)+
  stat_summary(fun= mean, geom="point", size=2, color="red",shape=15,alpha =1,na.rm = T)+
  geom_point(data = parent1,aes(x=group, y=parent1[,4]),na.rm = T,color="purple",shape=15,size=2)+
  geom_point(data = parent2,aes(x=group, y=parent2[,4]),na.rm = T,color="green",shape=15,size = 2)+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 15,colour = "black", family="serif"))+
  geom_text(data =anova_letter,aes(label = letter, x = name,y = yy), size = 5,vjust = -0.5)+
  theme_bw(base_size = 20)+xlab(NULL)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"), #angle = 45, hjust = 1,
        axis.text.y = element_text(color = "black")
  )+
  ylab(paste("Total sugar",i_year,sep = " "))
print(p)



}
dev.off()

