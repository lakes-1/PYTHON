options(stringsAsFactors=F)
rm(list = ls())
library(limma)
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")





data_name="GSE68719protein.txt"
label="GSE68719proteinExpadata.txt"
loc=3
plotDEG_dataType1=function(data_name,label,loc){
  data=as.matrix(read.csv(file =data_name, header = T,sep = "\t",fill=T,stringsAsFactors=F,row.names = 1))
  label= read.csv(file =label, header = T,sep = "\t",fill=T,stringsAsFactors=F)
  data=data[,loc:dim(data)[2]]
  mode(data)="numeric"
  dat=normalizeBetweenArrays(data)
  group_listTemp=label$label
  group_list=group_listTemp[group_listTemp%in%c(0,1)]
  data=data[,group_listTemp%in%c(0,1)]
  group_list=gsub("0","CT",group_list)
  group_list=gsub("1","PD",group_list)
  data2=rbind(data,group_list)
  return(data2)
}
brain_protein_68719=plotDEG_dataType1(data_name,label,3)


data_name="GSE68719_mlpd_PCG_DESeq2_norm_counts.txt"
label="GSE68719rnalabel.txt"
plotDEG_dataType2=function(data_name,label,loc){
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
  group_listTemp=label$label
  group_list=group_listTemp[group_listTemp%in%c(0,1)]
  data=data[,group_listTemp%in%c(0,1)]
  group_list[which(group_list==0)]="CT"
  group_list[which(group_list==1)]="PD"
  data2=rbind(data,group_list)
  return(data2)

  
  
}
brain_rna_68719=plotDEG_dataType2(data_name,label,3)


####
# dataType3，还利用到了平台探针信息
data_name="GSE20292_series_matrix.txt"
label="GSE20292label.txt"
platform="GPL96-57554.txt"
loc=2
plotDEG_dataType3=function(data_name,label,platform,loc){
  
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

  group_listTemp=label$label
  group_list=group_listTemp[group_listTemp%in%c(0,1)]
  data=data[,group_listTemp%in%c(0,1)]
  group_list[which(group_list==0)]="CT"
  group_list[which(group_list==1)]="PD"
  data2=rbind(data,group_list)
  return(data2)
 
  
  
}
brainSN_rna_20292=plotDEG_dataType3(data_name,label,platform,2)


## GSE7621
data_name="GSE7621_series_matrix.txt"
label="GSE7621label.txt"
platform="GSE20146_platform.txt"
brainSn_rna_7621=plotDEG_dataType3(data_name,label,platform,2)
## GSE99039
data_name="GSE99039_series_matrix.txt"
label="GSE99039label.txt"
platform="GSE20146_platform.txt"
blood_rna_99039=plotDEG_dataType3(data_name,label,platform,2)


## GSE62283
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
blood_protein_62283=as.data.frame(cbind(data1,group_list))


save(blood_protein_62283,blood_rna_99039,brain_protein_68719,brainSN_rna_20292,brainSn_rna_7621,
     brain_rna_68719,file="R_data/05_data_exp.Rdata")

