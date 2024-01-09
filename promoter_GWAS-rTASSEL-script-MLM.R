library(rTASSEL)
library(readxl)
library(ggplot2)
library(tidyverse)
library(MASS)#boxcox变换  正态分布变换
# library(ggforce)
wd = "G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\out_V1.0"
setwd(wd)
# Create rTASSEL data object

options(java.parameters = c("-Xmx15g", "-Xms2g")) ##rtassel memory set

genoPath = "G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\065.missing_maf.vcf"
phenoPath = "G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\2021Out_V2.0"

geno = readGenotypeTableFromPath(genoPath)

file.names = list.files(phenoPath)
pdf("zyy_GWAS-rTASSEL_MLM_V1.0.pdf",width = 10,height = 5.5)
for (i in file.names) {
  

pheno0 = read.csv(paste(phenoPath,"\\",i,sep = ""))
geno_ID = read_excel("G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\06-5-Berrysize-表型文件-02.xlsx",sheet = "GWAS")
pheno1 = merge(geno_ID,pheno0[,-1],by.x = "Genotype2",by.y = "Genotype",all.x = T)
#pheno1 = merge(geno_ID,zyy,by.x = "Genotype2",by.y = "Genotype",all.x = T)
if (i == "P_data_PCs.csv"){
  pheno2 = pheno1[,-1]
}else{
  pheno2 = pheno1[,c(2:4)]
}


for (ii_trait in names(pheno2)[-1]){
  pheno3 = pheno2[,c("Genotype3",ii_trait)]
  pheno33= pheno3

# pcaRes = pca(geno,nComponents = 2)
# pcaRes$PC_Datum
# pheno2 = merge(pheno1,pcaRes$PC_Datum,by.x = "Genotype3",by.y = "Taxa" )
# pheno = readPhenotypeFromDataFrame(pheno2[,c(1,3,4,(dim(pheno2)[2]-1),dim(pheno2)[2])],taxaID = "Genotype3",attributeTypes = c('data','data',"covariate","covariate"))

 
  dis_test = shapiro.test (pheno3[,2])
  if (dis_test$p.value>0.05){
    p0 = ggplot(data = pheno33,aes(x=(pheno33[,2])))+geom_histogram(bins = 20,aes(y=..density..))+xlab(ii_trait)+theme_bw(base_size = 18)
    p01 = ggplot(data = pheno33,aes(sample=(pheno33[,2])))+stat_qq()+stat_qq_line()+theme_bw(base_size = 18)+xlab("theoretical")+ylab ("smple")
  }else{
    p0 = ggplot(data = pheno33,aes(x=(pheno33[,2])))+geom_histogram(bins = 20,aes(y=..density..))+xlab(ii_trait)+theme_bw(base_size = 18)
    p01 = ggplot(data = pheno33,aes(sample=(pheno33[,2])))+stat_qq()+stat_qq_line()+theme_bw(base_size = 18)+xlab("theoretical")+ylab ("smple")
    
    min_phe = min(pheno3[,2],na.rm=T)
    if (min_phe <= 0){
      pheno3[,2] = pheno3[,2]+min_phe*-1.05
      pheno3[which(pheno3[,2]==0),2] = NA
    }

    
    #find optimal lambda for Box-Cox transformation 
    bc <- boxcox(pheno3[,2] ~ c(1:length(pheno3[,2])))
    lambda <- bc$x[which.max(bc$y)]
    cox_phe = (pheno3[,2]^lambda-1)/lambda
    pheno3[,2] = cox_phe
    pheno4 = pheno3
    p00 = ggplot(data = pheno4,aes(x=(pheno4[,2])))+geom_histogram(bins = 20,aes(y=..density..))+xlab(ii_trait)+theme_bw(base_size = 18)
    p001 = ggplot(data = pheno4,aes(sample=(pheno4[,2])))+stat_qq()+stat_qq_line()+theme_bw(base_size = 18)+xlab("theoretical")+ylab ("smple")
    #fit new linear regression model using the Box-Cox transformation
    # new_model <- lm(((lm_y^lambda-1)/lambda) ~lm_x)
  }
  
  pheno = readPhenotypeFromDataFrame(pheno3,taxaID = "Genotype3")


tasGenoPheno <- readGenotypePhenotype(
  genoPathOrObj    = geno,
  phenoPathDFOrObj = pheno
)
tasKin <- kinshipMatrix(tasObj = tasGenoPheno)
tasKinRMat <- as.matrix(tasKin)


# Calculate MLM
tasGenoPheno = filterGenotypeTableSites(tasGenoPheno,
                                        siteMinCount = round(dim(pheno3)[1]*0.9),
                                        siteMinAlleleFreq = 0.05,
                                        # removeMinorSNPStates = T,
                                        # removeSitesWithIndels = T,
                                        siteRangeFilterType = c("none", "sites", "position"),)


tasMLM <- assocModelFitter(
  tasObj = tasGenoPheno,  # <- our prior TASSEL object
  formula = . ~ .,    # <- run only meana.mean2
  fitMarkers = TRUE,      # <- set this to TRUE for GLM
  kinship = tasKin,       # <- our prior kinship object
  # kinship = NULL,
  maxP = 1,
  fastAssociation = FALSE,
  # biallelicOnly = T,
  # appendAddDom = T

)
plot.data0 = tasMLM@results$MLM_Stats
# View visualization
for (j in unique(plot.data0$Trait)) {
  
  
  
  plot.data =plot.data0[plot.data0$Trait == j,]
  plot.data$Pos = as.numeric(plot.data$Pos)
  plot.data$Chr = as.numeric(plot.data$Chr)
  # str(plot.data)
  
  chr.len = plot.data %>% group_by(Chr) %>% summarise(chr.len = max(Pos))
  chr.len = na.omit(chr.len)
  chr.len = as.data.frame(chr.len)
  chr.pos = chr.len %>% mutate(total = cumsum(chr.len)-chr.len) %>% dplyr::select(-chr.len)
  
  snp.pos = chr.pos %>% left_join(plot.data,.,by="Chr")%>%arrange(Chr,Pos)%>% mutate(Poscum = Pos + total)
  
  X.axis = snp.pos %>% group_by(Chr) %>% summarize(center = (max(Poscum)+min(Poscum))/2)
  X.axis = na.omit(X.axis)
  
  data = snp.pos %>% mutate(is_highlight = ifelse(-log10(p)>4.22,"yes","no"))
  
  p = ggplot(snp.pos, aes(x=Poscum,y=-log10(p)))+
    geom_point(aes(color = as.factor(Chr)),alpha=0.8,size=1.3)+
    scale_color_manual(values = rep(c("grey","skyblue"),20))+
    scale_x_continuous(label = X.axis$Chr,breaks = X.axis$center)+
    geom_point(data = subset(data,is_highlight == "yes"),color="orange",size=2)+
    geom_hline(yintercept = 4.22,color= "red",linetype = 2,linewidth = 1)+xlab("Chromosome")+
    # facet_zoom(x=Poscum >= 30000000 & Poscum <= 41000000)+
    # scale_y_continuous(expand = c(0,0))+
    theme_bw(base_size = 18)+#ggtitle(j)+
    theme(
      legend.position = 'none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}
##########raw#########
pheno_raw = readPhenotypeFromDataFrame(pheno33,taxaID = "Genotype3")
tasGenoPheno_raw <- readGenotypePhenotype(
  genoPathOrObj    = geno,
  phenoPathDFOrObj = pheno_raw
)
tasKin_raw <- kinshipMatrix(tasObj = tasGenoPheno_raw)
tasKinRMat_raw <- as.matrix(tasKin_raw)

# Calculate MLM
tasGenoPheno_raw = filterGenotypeTableSites(tasGenoPheno_raw,
                                        siteMinCount = round(dim(pheno33)[1]*0.9),
                                        siteMinAlleleFreq = 0.05,
                                        # removeMinorSNPStates = T,
                                        # removeSitesWithIndels = T,
                                        siteRangeFilterType = c("none", "sites", "position"),)


tasMLM_raw <- assocModelFitter(
  tasObj = tasGenoPheno_raw,  # <- our prior TASSEL object
  formula = . ~ .,    # <- run only meana.mean2
  fitMarkers = TRUE,      # <- set this to TRUE for GLM
  kinship = tasKin_raw,       # <- our prior kinship object
  # kinship = NULL,
  maxP = 1,
  fastAssociation = FALSE,
  # biallelicOnly = T,
  # appendAddDom = T
  
)
# manhattanEH <- manhattanPlot(
#   assocStats = tasMLM$GLM_Stats,
#   trait      = "meana.mean2",
#   threshold  = 5
# )
# 
# manhattanEH
# plot.data0 = tasMLM$MLM_Stats
plot.data0_raw = tasMLM_raw@results$MLM_Stats
# View visualization
for (j in unique(plot.data0_raw$Trait)) {
  


plot.data_raw =plot.data0_raw[plot.data0_raw$Trait == j,]
plot.data_raw$Pos = as.numeric(plot.data_raw$Pos)
plot.data_raw$Chr = as.numeric(plot.data_raw$Chr)
# str(plot.data)

chr.len = plot.data_raw %>% group_by(Chr) %>% summarise(chr.len = max(Pos))
chr.len = na.omit(chr.len)
chr.len = as.data.frame(chr.len)
chr.pos = chr.len %>% mutate(total = cumsum(chr.len)-chr.len) %>% dplyr::select(-chr.len)

snp.pos = chr.pos %>% left_join(plot.data_raw,.,by="Chr")%>%arrange(Chr,Pos)%>% mutate(Poscum = Pos + total)

X.axis = snp.pos %>% group_by(Chr) %>% summarize(center = (max(Poscum)+min(Poscum))/2)
X.axis = na.omit(X.axis)

data_raw = snp.pos %>% mutate(is_highlight = ifelse(-log10(p)>4.22,"yes","no"))

pp = ggplot(snp.pos, aes(x=Poscum,y=-log10(p)))+
  geom_point(aes(color = as.factor(Chr)),alpha=0.8,size=1.3)+
  scale_color_manual(values = rep(c("grey","skyblue"),20))+
  scale_x_continuous(label = X.axis$Chr,breaks = X.axis$center)+
  geom_point(data = subset(data_raw,is_highlight == "yes"),color="orange",size=2)+
  geom_hline(yintercept = 4.22,color= "red",linetype = 2,linewidth = 1)+xlab("Chromosome")+
  # facet_zoom(x=Poscum >= 30000000 & Poscum <= 41000000)+
  # scale_y_continuous(expand = c(0,0))+
  theme_bw(base_size = 18)+#ggtitle(j)+
  theme(
    legend.position = 'none',
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
}
if (dis_test$p.value>0.05){
  p2=cowplot::plot_grid(p0,p01,nrow=1)
  p3 = cowplot::plot_grid(p2,p,nrow=2)
  png(paste(wd,"\\MLM_V3.0\\",i,"_",ii_trait,".png",sep = ""),width = 10,height = 10, units = "in",res = 300)
  # png(paste(wd,"\\","_",j,".png",sep = ""),width = 10,height = 5.5, units = "in",res = 200)
  print(p3)
  
  dev.off()
}else{
  p2=cowplot::plot_grid(p0,p01,nrow=1)
  p22=cowplot::plot_grid(p00,p001,nrow=1)
  p3 = cowplot::plot_grid(p2,pp,p22,p,nrow=4)
  png(paste(wd,"\\MLM_V3.0\\",i,"_",ii_trait,".png",sep = ""),width = 10,height = 18, units = "in",res = 300)
  # png(paste(wd,"\\","_",j,".png",sep = ""),width = 10,height = 5.5, units = "in",res = 200)
  print(p3)
  
  dev.off()
}






don = plot.data0
don = na.omit(don)

don_2 = don[-log10(don$p)>4.22,]

write.csv(don_2,paste(wd,"\\MLM_V3.0\\",i,"_",ii_trait,".csv",sep = ""))

# for (ii in unique(don_2$Chr)) {
#   don_ii = don_2[don_2$Chr == ii,]
#   don_ii$Pos = as.numeric(don_ii$Pos)
#   don_range = range(don_ii[-log10(don_ii$p)>4.22,"Pos"],na.rm = T)
#   
#   tasGenoPhenoFilt <- filterGenotypeTableSites(
#     tasObj              = tasGenoPheno,
#     siteRangeFilterType = "position",
#     startPos            = don_range[1],
#     endPos              = don_range[2],
#     startChr            = as.numeric(ii),
#     endChr              = as.numeric(ii)
#   )
#   out.ii = siteSummary(tasGenoPhenoFilt)
#   
# }

}
}
dev.off()


#### matching gene ######

all_snp = read.csv("C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\out_V1.0\\data_merge.csv")
all_snp = read.csv("C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\out_V1.0\\GLM_V1.0_DATA\\data_merge.csv")
str(all_snp)
all_snp = unique(all_snp)

gene_list = read_excel("C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\参考基因组-AllGenes.xlsx",skip = 1,sheet = "Sheet1")
gene_list= as.data.frame(gene_list)
str(gene_list)


data_all = data.frame()
for (i_snp in 1:dim(all_snp)[1]) {
  
  i_chr = all_snp[i_snp,]$Chr
  i_pos = all_snp[i_snp,]$Pos
  
  gene_chr = all_hap[all_hap$Chromosome == as.character(i_chr),c(1,2,185,186)]
  
  gene_chr = unique(gene_chr)
  
  gene_pos = gene_chr[gene_chr$position== i_pos,"Gene.refGene"]
  
  gene_pos1 =  na.omit(gene_pos)
  if  ( length(gene_pos1)>0 ){
    for (jj in gene_pos1){
  data_temp = all_snp[i_snp,]
  data_temp2 = gene_chr[gene_chr$Gene.refGene == jj,]
  data_all0 = cbind(data_temp,data_temp2)
  data_all = rbind(data_all,data_all0)
    }
  }
}
write.csv(data_all,"PC_data_merge_geneID_V2.0.csv")

data_all = read.csv("C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\out_V1.0\\data_merge_geneID_V2.3_selected.csv")

data_all0 = distinct(data_all,Func.refGene,Gene.refGene,Trait,Marker,.keep_all = T)

data_all = data_all0
#####
#######promoter haplotype#####
library(dplyr)
all_hap = read.csv("C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\原始基因型数据_碱基.csv")
str(all_hap)
all_ano = read.csv("C:\\Users\\wyong\\OneDrive\\RNAseq-xuxiaobo\\基因注释\\PN_V2\\PN06-5.pinotnoir_multianno.csv")

all_hap$merge_name = paste(all_hap$Chromosome,all_hap$position,sep = "-")
all_ano$merge_name = paste(all_ano$Chr,all_ano$Start,sep = "-")
#merge annotation and hap file
all_hap_ano = merge(all_hap,all_ano[all_ano$merge_name%in%all_hap$merge_name,c(6:10,205)],by = "merge_name")
all_hap0 = all_hap
str(all_hap_ano)
dim(all_hap_ano)

all_hap = all_hap_ano[all_hap_ano$Func.refGene%in%c("UTR5","upstream","downstream","exonic","UTR3"),c(3:188)]
dim(all_hap)

genes = unique(data_all$geneID)
hap_num_all = data.frame()
for (ii in genes) {
  
  chr = unique(data_all[data_all$geneID == ii,]$Chr)
  start = unique(data_all[data_all$geneID == ii,]$start)
  end = unique(data_all[data_all$geneID == ii,]$end)
  
  
  hap_ii = all_hap[(all_hap$Chromosome==as.character(chr))&(all_hap$position%in%((start-3000):start)),]
  if (dim(hap_ii)[1]>0){
  t_hap_ii =as.data.frame(t(hap_ii))
  t_hap_ii = na.omit(t_hap_ii)
  #paste hap
  xx = t_hap_ii[,1]
  if (dim(t_hap_ii)[2] >1 ){
  for (paste_i in 2:dim(t_hap_ii)[2]) {
    xx = paste(xx,t_hap_ii[,paste_i],sep = "_")
    
  }
  }
  t_hap_ii$hap = xx
  t_hap_ii0 = t_hap_ii[-5:-1,]
  t_hap_ii0$Genotype =substring(rownames(t_hap_ii0),7)
  pheno2 = na.omit(pheno1)
  t_hap_ii1 = merge(pheno2[,1:2],t_hap_ii0,by.x = "Genotype2",by.y = "Genotype",all.x = T)
  t_hap_ii1 = na.omit(t_hap_ii1)
  hap_num = group_by(t_hap_ii1,hap) %>% summarise_each(funs(length))
  hap_num = as.data.frame(hap_num)
  hap_num1 = hap_num[,1:2]
  names(hap_num1)[2] = "num"
  
  hap_num1$marker_num = dim(t_hap_ii)[2]-1
  hap_num1$gene = ii
  
  
  hap_num_all = rbind(hap_num_all ,hap_num1)
  }
  
}
write.csv(hap_num_all,"hap_num_all_V1.0.csv")
##########################################

#####trait distribution###################
library(tidyverse)
library(hrbrthemes)
library(viridis)
genes = merge_ana_gwas_genes$gene
genes = unique(hap_num_all$gene)
genes = unique(data_all$Gene.refGene)

pdf("Promoter trait distribution_zyy_sbx_PCs_V1.2.pdf",width = 10,height = 6)
#V2:更新注释文件，只要内含子和5‘及3’
for (ii in genes) {
  
  # chr = unique(data_all[data_all$geneID == ii,]$Chr)
  # start = unique(data_all[data_all$geneID == ii,]$start)
  # end = unique(data_all[data_all$geneID == ii,]$end)
  
  
  hap_ii = all_hap[all_hap$Gene.refGene == ii,c(185,1:184)]
  
  t_hap_ii =as.data.frame(t(hap_ii))
  t_hap_ii = na.omit(t_hap_ii)
  #paste hap
  xx = t_hap_ii[,1]
  if (dim(t_hap_ii)[2] >1 ){
    for (paste_i in 2:dim(t_hap_ii)[2]) {
      xx = paste(xx,t_hap_ii[,paste_i],sep = "_")
      
    }
  }
  t_hap_ii$hap = xx
  t_hap_ii0 = t_hap_ii[-5:-1,]
  t_hap_ii0$Genotype =substring(rownames(t_hap_ii0),7)
  # pheno2 = na.omit(pheno1)
  t_hap_ii1 = merge(geno_ID,t_hap_ii0,by.x = "Genotype2",by.y = "Genotype",all.x = T)
  t_hap_ii1 = na.omit(t_hap_ii1)
  hap_num = group_by(t_hap_ii1,hap) %>% summarise_each(funs(length))
  hap_num = as.data.frame(hap_num)
  hap_num1 = hap_num[,1:2]
  names(hap_num1)[2] = "num"
  
  hap_num1$marker_num = dim(t_hap_ii)[2]-1
  hap_num1$gene = ii
  
  traits =unique(data_all[data_all$Gene.refGene == ii,]$Trait)
  # data_path = "C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\2021Out_V2.0_read"
  for (trait_i in traits) {
    # read_trait = strsplit(trait_i,split = "[.]")
    # pheno_data = read.csv(paste(data_path,"\\",read_trait[[1]][1],".csv",sep = ""))
    
    
    # pheno_data0 = pheno_data[,c("Genotype",trait_i)]
    # pheno_data0 = pheno_data[,c("Genotype",trait_i)]
    #pheno_data0 = read.csv("C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\sbx\\2021skin.csv")
    #pheno_data0 = read.csv("C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\data_PCs.csv")
    # m_zyy_sbx = merge(sbx,zyy,by = "Genotype",all.x = T)
    # m_all = merge(m_zyy_sbx,PCs,by = "Genotype",all.x = T)
    
    hap_pheno = merge(t_hap_ii1,pheno_data0[,c("Genotype",trait_i)],by.x = "Genotype2",by.y = "Genotype",all.x = T)
    
    # hap_pheno = merge(t_hap_ii1,zyy[,c("Genotype",trait_i)],by.x = "Genotype2",by.y = "Genotype",all.x = T)
    # hap_pheno = na.omit(hap_pheno)
    #plot
    
    hap_pheno = na.omit(hap_pheno)
    #plot
    
    data_plot = hap_pheno
    
    
    data_plot_summmary1 = group_by(data_plot,hap)%>%summarise_each(funs(mean))
    data_plot_summmary1 = as.data.frame(data_plot_summmary1)
    
    data_plot_summmary2 = group_by(data_plot,hap)%>%summarise_each(funs(length))
    data_plot_summmary2 = as.data.frame(data_plot_summmary2)
    
    names(data_plot_summmary2)[2] = "Genotype"
    
    data_plot_summmary = merge(data_plot_summmary1[,c("hap",trait_i)],data_plot_summmary2[,c("hap","Genotype")],by = "hap")
    
    # data_plot_summmary = data_plot_summmary[order(data_plot_summmary[,trait_i]),]
    
    data_plot_summmary = data_plot_summmary[data_plot_summmary$Genotype>(dim(data_plot)[1]*0.05),]
    
    data_plot_summmary$name = paste("Hap",1:dim(data_plot_summmary)[1],sep = "")
    data_plot_summmary$name2 = paste(data_plot_summmary$name,"\n (n = ",data_plot_summmary$Genotype,")",sep = "")
    
    level = data_plot_summmary[order(data_plot_summmary[,trait_i]),]
    # name_level = data_plot_summmary[order(data_plot_summmary[,trait_i]),]$name
    # num_level = data_plot_summmary[order(data_plot_summmary[,trait_i]),]$Genotype
    data_plot$hap = factor(data_plot$hap,levels = level$hap)
    
    
    
    
    ############anova########
    # library(tidyverse)
    # library(multcompView)
    # library(ggsci)
    if (dim(data_plot_summmary)[1]>1){
    anova <- aov(data_plot[,trait_i]~hap,data=data_plot)
    tukey <- TukeyHSD(anova)
    cld <- multcompLetters4(anova,tukey)
    aa = as.data.frame.list(cld$hap) %>% select(1) %>% rownames_to_column("hap")
    
    data_plot_summmary3 = group_by(data_plot,hap)%>%summarise_each(funs(max))
    data_plot_summmary3 = as.data.frame(data_plot_summmary3)
    
    anova_letter0 = merge(data_plot_summmary3,aa,by = "hap",all.x = T)
    anova_letter = merge(anova_letter0,data_plot_summmary[,c("hap","name2")],by = "hap")
    
    if (length(unique(anova_letter$Letters))>1){
    ################
    data_plot = na.omit(data_plot)
    data_plot0 = merge(data_plot,data_plot_summmary[,c("hap","name2")],by = "hap")
    # x_names = paste(name_level,"\n (n = ",num_level,")",sep = "")
    #factor level
    data_plot0$name2 = factor(data_plot0$name2,levels = level$name2 )
    #plot
    p = ggplot( data = data_plot0,aes(x=name2, y=data_plot0[,trait_i], fill=name2)) +
      geom_boxplot() +
      geom_text(data =anova_letter,aes(label = Letters, x = name2,y = anova_letter[,trait_i]), size = 5,vjust = -0.5)+
      # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      ylim(min(data_plot0[,trait_i]),max(data_plot0[,trait_i])*1.05)+
      theme_bw(base_size = 25) +
      scale_fill_aaas()+
      theme(
        legend.position="none",
      ) +
      # scale_x_discrete(labels =x_names)+
      # ylab(strsplit(trait_i,"2")[[1]][1])+
      ylab(trait_i)+
      ggtitle(ii) +
      xlab("")
    print(p)
    }
    }
  }

  
}

dev.off()


#########merge kenwn genes#############
ana_genes = read.csv("C:\\Users\\wyong\\OneDrive\\C5.23=基因芯片和高通量表型组\\5.GWAS-rTASSEL\\Ana_pathway_id_allgene.csv")
gwas_genes = data.frame(gene = genes,ID = 1)

merge_ana_gwas_genes = merge(gwas_genes,ana_genes,by.x = "gene",by.y = "gene_id")
