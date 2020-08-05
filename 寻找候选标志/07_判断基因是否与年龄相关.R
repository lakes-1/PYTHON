## 载入数据
rm(list = ls())
options(stringsAsFactors=F)
library(limma)
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")
load("R_data/01_degs.Rdata")
load("R_data/05_data_exp.Rdata")
# 判断基因是否与年龄相关

# 1.载入年龄信息 
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data_duplicate\\result\\01_preprocess")
label = read.csv("df_all_label.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F)
age = read.csv("data_age_sex.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F)

# 2.载入基因表达谱，可以根据05_修改寻找血液一致且fdr显著但rna肉眼可观.R脚本跑一遍获得
gene=c("NDUFB3","OPA1","NDUFA9","NDUFB6","NDUFS4"
       ,"NNT","RIMS3","SRPRB","SYT1","UVRAG")

blood_rna_99039=blood_rna_99039[c(gene,"group_list"),]
brainSN_rna_20292=brainSN_rna_20292[c(gene,"group_list"),]
brainSn_rna_7621=brainSn_rna_7621[c(gene,"group_list"),]
brain_rna_68719=brain_rna_68719[c(gene,"group_list"),]



# 3.对基因表达值与年龄做相关


# 下载临床信息吧
#安装并加载GEOmirror包，去github下载zip，再本地安装
# rm(list = ls())
# if(T){library(devtools)
#   #install_github("jmzeng1314/GEOmirror")
#   #devtools::install('D:\\BaiduNetdiskDownload\\GEOmirror-master.zip\\')
#   #devtools::install_local("D:\\BaiduNetdiskDownload\\GEOmirror-master.zip\\")
#   #devtools::install_local("D:\\BaiduNetdiskDownload\\jmzeng1314-idmap3-59efe3f.tar.gz\\")
#   library(GEOmirror)
#   library(idmap1)
#   library(idmap2)
#   library(idmap3)
#   library("GEOquery")
# }
# 
# 
# 
# gset <- getGEO(a[1],
#                AnnotGPL = F,     ## 注释文件
#                getGPL = F)    
# gset=gset[[1]]
# pd_99039=pData(gset)
# write.table(pd_99039,file="GSE99039pdata.txt",row.names=TRUE,col.names=TRUE,sep="\t")
# pd_20292=pData(gset)
# write.table(pd_20292,file="GSE20292pdata.txt",row.names=TRUE,col.names=TRUE,sep="\t")
# pd_7621=pData(gset)
# write.table(pd_7621,file="GSE7621pdata.txt",row.names=TRUE,col.names=TRUE,sep="\t")
# pd_68719=pData(gset)
# write.table(pd_68719,file="GSE68719pdata.txt",row.names=TRUE,col.names=TRUE,sep="\t")

# #探针注释，提取四个芯片共有的探针数据
# a='GSE99039'
# a='GSE20292'
# a='GSE7621' 
# a="GSE68719"

# 99039这套数据集,762这套数据集没有年龄信息，就只分析其他基因与年龄的关系好了

pd_99039_age=data.frame(pd_99039$geo_accession,pd_99039$`age at 1st symptoms:ch1`,pd_99039$characteristics_ch1.2)
pd_20292_age=cbind(pd_20292$geo_accession,pd_20292$characteristics_ch1.1)
pd_7621_age=cbind(pd_7621$geo_accession,pd_7621$characteristics_ch1.1)
pd_68719_age=cbind(pd_68719$geo_accession,pd_68719$characteristics_ch1.4)



# 99039,提取
n=sapply(pd_99039_age$pd_99039.characteristics_ch1.2,function(x){strsplit(x = x,split = ": ")[[1]][2]})
pd_99039_age=pd_99039_age[n%in%c("CONTROL","IPD"),]
pd_99039_age1=pd_99039_age[!is.na(pd_99039_age$pd_99039..age.at.1st.symptoms.ch1.),]
sample_99039=blood_rna_99039[,pd_99039_age1$pd_99039.geo_accession]
age_information=pd_99039_age1$pd_99039..age.at.1st.symptoms.ch1.
temp=apply(sample_99039[-11,], 1, function(x){cor.test(as.numeric(x),as.numeric(age_information),method = "spearman")})
p_value=sapply(temp,function(x){x$p.value})
cor_value=sapply(temp,function(x){x$estimate})
fdr_value<- p.adjust(p_value,method="fdr")
result=cbind(cor_value,fdr_value)
save(pd_99039_age,pd_20292_age,pd_7621_age,pd_68719_age ,file = "age_data.Rdata")
# 
load("age_data.Rdata")

n=sapply(pd_99039_age$pd_99039.characteristics_ch1.2,function(x){strsplit(x = x,split = ": ")[[1]][2]})
pd_99039_age=pd_99039_age[n%in%c("CONTROL"),]
pd_99039_age1=pd_99039_age[!is.na(pd_99039_age$pd_99039..age.at.1st.symptoms.ch1.),]
sample_99039=blood_rna_99039[,pd_99039_age1$pd_99039.geo_accession]
age_information=pd_99039_age1$pd_99039..age.at.1st.symptoms.ch1.
temp=apply(sample_99039[-11,], 1, function(x){cor.test(as.numeric(x),as.numeric(age_information),method = "spearman")})
p_value=sapply(temp,function(x){x$p.value})
cor_value=sapply(temp,function(x){x$estimate})
fdr_value<- p.adjust(p_value,method="fdr")
result=cbind(cor_value,fdr_value)
