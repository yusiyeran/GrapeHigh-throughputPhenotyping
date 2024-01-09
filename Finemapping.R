#fine mapping
library(formattable)

getwd()
wd = "G:\\OneDrive\\C5.23=基因芯片和高通量表型组\\10.图-加性效应\\FineMapping-ShapeIndex"
setwd(wd)
formattable(mtcars,list(mpg=color_bar("pink",proportion)))

formattable(mtcars,list(mpg=color_tile("pink","orange")))

chr4 = read_excel("过滤后基因型数据_碱基.xls")

shapeindex_4 = read_excel("Finemapping-ShapeIndex-Chr4.xlsx",sheet = "Sheet9 (3)order",skip = 2,na="NA")

shapeindex_4 = as.data.frame(shapeindex_4)

shapeindex_phe = read_excel("Finemapping-ShapeIndex-Chr4.xlsx",sheet = "Phenotype",skip = 0,na="NA")
shapeindex_phe = as.data.frame(shapeindex_phe)
shapeindex_phe = na.omit(shapeindex_phe)

shapeindex_phe = shapeindex_phe[order(shapeindex_phe$Phenotype,decreasing = T),]

shapeindex_phe = shapeindex_phe[-which(shapeindex_phe$Genotype%in%c("bei_feng","yan_73")),]
k=0.3
max_shapeindex_phe = shapeindex_phe[1:round(dim(shapeindex_phe)[1]*k),]
min_shapeindex_phe = shapeindex_phe[round(dim(shapeindex_phe)[1]*(1-k)):dim(shapeindex_phe)[1],]
#父母本已筛选
max_shapeindex_4 = shapeindex_4[,max_shapeindex_phe$Genotype]

max_data = data.frame()
for (i in 1:dim(max_shapeindex_4)[1]){
  temp1 = max_shapeindex_4[i,]
  temp1_t = t(temp1)
  temp1_t = na.omit(temp1_t)
  temp1_t = data.frame(Genotype = rownames(temp1_t),group = temp1_t[,1])
  
  temp1_t_group = group_by(temp1_t, group) %>% summarize_each(funs(length))
  temp1_t_group = as.data.frame(temp1_t_group)
  temp1_t_group$fre = temp1_t_group$Genotype/dim(max_shapeindex_4)[2]
  temp1_t_group = temp1_t_group[order(temp1_t_group$fre,decreasing = T),]
  temp1_t_group$chr = shapeindex_4[i,"Chromosome"]
  temp1_t_group$pos = shapeindex_4[i,"position"]
  max_data = rbind(max_data,temp1_t_group)
}

