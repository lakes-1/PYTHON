# GSE68719
# Deseq 
# 
rm(list = ls())
library(limma)
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")
Searchdeg = function(exprSet,design,contrast.matrix){
  ##step1,lmfit起过滤的手续
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  
  ##step3
  # 有了比较矩阵后，coef=1，而number=Inf是把所有结果都打印出来
  tempOutput = topTable(fit2, coef=1, number =Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}
make_limmaContrast=function(data,group_list){
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  head(design)
  exprSet=data
  mode(exprSet)="numeric"
  rownames(design)=colnames(exprSet)
  design
  
  # 比较矩阵
  # 这个矩阵声明，我们要把 npc 组跟 Normal 进行差异分析比较
  contrast.matrix<-makeContrasts("pd-normal",
                                 levels = design)
  contrast.matrix
  
  
  deg = Searchdeg(exprSet,design,contrast.matrix)
  
  return(deg)
}
#dataType1指的是数据的行是独立的，没有重复，其实后面可以把datatype1和datatypeN整合在一块，加一个选项
data_name="GSE68719protein.txt"
label="GSE68719proteinExpadata.txt"
searchDEG_dataType1=function(data_name,label,loc){
data=as.matrix(read.csv(file =data_name, header = T,sep = "\t",fill=T,stringsAsFactors=F,row.names = 1))
label= read.csv(file =label, header = T,sep = "\t",fill=T,stringsAsFactors=F)
data=data[,loc:dim(data)[2]]
mode(data)="numeric"
dat=normalizeBetweenArrays(data)
# limma
library(limma)
# 方法2：制作比较矩阵
# 可以随心所欲的指定任意两组进行比较
# 判断标签是否与样本的列名一致
group_listTemp=label$label
group_list=group_listTemp[group_listTemp%in%c(0,1)]
data=data[,group_listTemp%in%c(0,1)]
group_list[which(group_list==0)]="normal"
group_list[which(group_list==1)]="pd"
degs=make_limmaContrast(data,group_list)

}
brain_protein_68719=searchDEG_dataType1(data_name,label,3)

# dataType2是指第一列的探针或基因名可以重复，但我在内部操作让其均值化的操作
# 其实1，2类可以合并，但更为重要的是输入数据的一致性，还能保证用较少的代码实行比较通泛的功能

data_name="GSE68719_mlpd_PCG_DESeq2_norm_counts.txt"
label="GSE68719rnalabel.txt"
searchDEG_dataType2=function(data_name,label,loc){
  data=as.matrix(read.csv(file =data_name, header = T,sep = "\t",fill=T,stringsAsFactors=F))
  label= read.csv(file =label, header = T,sep = "\t",fill=T,stringsAsFactors=F)
  
  
  ids=data[,1]
  data=data[,loc:dim(data)[2]]
  exp=as.matrix(tapply(as.matrix(as.numeric(data[,1])),ids,simplify = TRUE,mean))
  for(j in 2:ncol(data)){
    b=as.matrix(tapply(as.matrix(as.numeric(data[,j])),ids,simplify = TRUE,mean))
    exp=cbind(exp,b)
  }
  colnames(exp)=colnames(data)
  data=exp
  
  
  data=normalizeBetweenArrays(data)
  # limma
  library(limma)
  group_listTemp=label$label
  group_list=group_listTemp[group_listTemp%in%c(0,1)]
  data=data[,group_listTemp%in%c(0,1)]
  group_list[which(group_list==0)]="normal"
  group_list[which(group_list==1)]="pd"

  degs=make_limmaContrast(data,group_list)
  
}
brain_rna_68719=searchDEG_dataType2(data_name,label,3)

# dataType3，还利用到了平台探针信息
data_name="GSE20292_series_matrix.txt"
label="GSE20292label.txt"
platform="GPL96-57554.txt"
loc=2
searchDEG_dataType3=function(data_name,label,platform,loc){

data=as.matrix(read.csv(file =data_name, header = T,sep = "\t",fill=T,stringsAsFactors=F))
label= read.csv(file =label, header = T,sep = "\t",fill=T,stringsAsFactors=F)
platform= read.csv(file =platform, header = T,sep = "\t",fill=T,stringsAsFactors=F)

# 取共同交叠的探针
ids=data[,1]
row.names(data)=ids
idsInPlatform=intersect(ids,platform$ID)
data=data[idsInPlatform,]
ids_symbol=platform[match(idsInPlatform,platform$ID),]$Gene.Symbol

# 对共同的探针均值化
data=data[,loc:dim(data)[2]]
exp=as.matrix(tapply(as.matrix(as.numeric(data[,1])),ids_symbol,simplify = TRUE,mean))
for(j in 2:ncol(data)){
  b=as.matrix(tapply(as.matrix(as.numeric(data[,j])),ids_symbol,simplify = TRUE,mean))
  exp=cbind(exp,b)
}
colnames(exp)=colnames(data)
data=exp
mode(data)="numeric"
data=normalizeBetweenArrays(data)
# limma
library(limma)
# 方法2：制作比较矩阵
# 可以随心所欲的指定任意两组进行比较
# 判断标签是否与样本的列名一致
group_listTemp=label$label
group_list=group_listTemp[group_listTemp%in%c(0,1)]
data=data[,group_listTemp%in%c(0,1)]
group_list[which(group_list==0)]="normal"
group_list[which(group_list==1)]="pd"

degs=make_limmaContrast(data,group_list)

}
brainSN_rna_20292=searchDEG_dataType3(data_name,label,platform,2)

## GSE20146
data_name="GSE20146_series_matrix.txt"
label="GSE20146label.txt"
platform="GSE20146_platform.txt"
brainGpi_rna_20146=searchDEG_dataType3(data_name,label,platform,2)
## GSE7621
data_name="GSE7621_series_matrix.txt"
label="GSE7621label.txt"
platform="GSE20146_platform.txt"
brainSn_rna_7621=searchDEG_dataType3(data_name,label,platform,2)
## GSE99039
data_name="GSE99039_series_matrix.txt"
label="GSE99039label.txt"
platform="GSE20146_platform.txt"
blood_rna_99039=searchDEG_dataType3(data_name,label,platform,2)


# dataType4,蛋白质GSE62283了
path_GSE62283="D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\ppa_result\\GSE62283_exp.txt"
data_62283=as.matrix(read.csv(file =path_GSE62283, header = T,sep = "\t",fill=T,stringsAsFactors=F,row.names = 1))
gene_GSE62283="D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\ppa_result\\GSE62283_gene.txt"
gene= read.csv(file =gene_GSE62283,header = T,sep = "\t",fill=T,stringsAsFactors=F)
control=which(gene$Description=="Internal Control")
data_62283=data_62283[-control,]
# 标签
load("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\ppa_result\\R_data\\step3_dataNor_label.Rdata")
label_62283=label$class[match(colnames(data_62283),label$dataset)]
group_list=gsub("[1-6.]","",label_62283)
group_list[which(group_list=="PD")]='pd'
group_list[which(group_list=="CT")]='normal'
blood_protein_62283=make_limmaContrast(data_62283,group_list)
# 基因描述转变成基因
gene=gene[-control,]
library(stringr)
gene_name=str_extract_all(gene$Description,'\\(.*\\)')
gene_name1=str_extract_all(gene_name,'\\([A-Z]{1,}[A-Za-z0-9]{1,}\\)')
gene_name3=sapply(gene_name1,function(x){if(length(x)>1){x[-1]}
  else{x[1]}
  
})
gene_name4=gsub(pattern = "[\\(\\)]",replacement = "",x = gene_name3)
gene_name4[which(is.na(gene_name4))]=gene$Description[which(is.na(gene_name4))]
gene_id=apply(as.matrix(gene$ID),1,function(x){strsplit(x,"~")[[1]][2]})
gene_name4=as.matrix(gene_name4)
rownames(gene_name4)=gene_id
# save(blood_protein_62283,blood_rna_99039,brain_protein_68719,brainSN_rna_20292,brainSn_rna_7621,
#      brain_rna_68719,gene_name4,file="R_data/01_degs.Rdata")
save(blood_protein_62283,blood_rna_99039,brain_protein_68719,brainSN_rna_20292,brainSn_rna_7621,
     brain_rna_68719,gene_name4,file="R_data/011_degs.Rdata")