options(stringsAsFactors=F)
rm(list = ls())
library(limma)
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")
load("R_data/03_degs.Rdata")

library(ggplot2)
library(gridExtra)# split the grid 
data_name=c(blood_protein_62283,blood_rna_99039,brain_protein_68719,brainSN_rna_20292,brainSn_rna_7621,
            brain_rna_68719)
blood_rna_99039$group_list[blood_rna_99039$group_list=="normal"]="CT"
blood_rna_99039$group_list[blood_rna_99039$group_list=="pd"]="PD"
brainSN_rna_20292$group_list[brainSN_rna_20292$group_list=="normal"]="CT"
brainSN_rna_20292$group_list[brainSN_rna_20292$group_list=="pd"]="PD"
brainSn_rna_7621$group_list[brainSn_rna_7621$group_list=="normal"]="CT"
brainSn_rna_7621$group_list[brainSn_rna_7621$group_list=="pd"]="PD"
colnames(blood_protein_62283)[1:2]=c("NDUFB3","OPA1")
#¶¨Òå»­Í¼º¯Êý
plot_boxplot=function(data1,group_list,name){

  ggplot(data = data_name[1], aes(x=group_list, y=as.numeric(NDUFB3), color=group_list))+
    geom_boxplot()
  
  name='pic/boxplot1.png'
  ggsave(filename = "temp.jpg")
}


gene_name=c("NDUFB3","NDUFB3")
p1=ggplot(data = blood_rna_99039, aes(x=group_list, y=as.numeric(NDUFB3), color=group_list))+
  geom_boxplot()+labs(y="",title="RNA:GSE99039")+theme(legend.position="none")
p2=ggplot(data = brain_protein_68719, aes(x=group_list, y=as.numeric(NDUFB3), color=group_list))+
  geom_boxplot()+labs(y="",title="PRO:GSE68719")
p3=ggplot(data = brain_rna_68719, aes(x=group_list, y=as.numeric(NDUFB3), color=group_list))+
  geom_boxplot()+labs(y="NDUFB3",title = "RNA:GSE68719")+theme(legend.position="none")
p4=ggplot(data = brainSN_rna_20292, aes(x=group_list, y=as.numeric(NDUFB3), color=group_list))+
  geom_boxplot()+labs(y="",title="RNA:GSE20292")+theme(legend.position="none")
p5=ggplot(data = brainSn_rna_7621, aes(x=group_list, y=as.numeric(NDUFB3), color=group_list))+
  geom_boxplot()+labs(y="",title = "RNA:GSE7621")+theme(legend.position="none")
p6=ggplot(data = blood_protein_62283, aes(x=group_list, y=as.numeric(NDUFB3), color=group_list))+
  geom_boxplot()+labs(y="",title = "PRO:GSE62283")
p6
grid.arrange(p3,p4,p5,p2,p1,p6,nrow=2,ncol=4)
g <- arrangeGrob(p3,p4,p5,p2,p1,p6,nrow=2,ncol=4) #generates g
ggsave(file="pic/NDUFB3.png", g) 






gene_name=c("OPA1","OPA1")
p1=ggplot(data = blood_rna_99039, aes(x=group_list, y=as.numeric(OPA1), color=group_list))+
  geom_boxplot()+labs(y="OPA1",title="RNA:GSE99039")
p2=ggplot(data = brain_protein_68719, aes(x=group_list, y=as.numeric(OPA1), color=group_list))+
  geom_boxplot()+labs(y="OPA1",title="PRO:GSE68719")
p3=ggplot(data = brain_rna_68719, aes(x=group_list, y=as.numeric(OPA1), color=group_list))+
  geom_boxplot()+labs(y="OPA1",title = "RNA:GSE68719")
p4=ggplot(data = brainSN_rna_20292, aes(x=group_list, y=as.numeric(OPA1), color=group_list))+
  geom_boxplot()+labs(y="OPA1",title="RNA:GSE20292")
p5=ggplot(data = brainSn_rna_7621, aes(x=group_list, y=as.numeric(OPA1), color=group_list))+
  geom_boxplot()+labs(y="OPA1",title = "RNA:GSE7621")
p6=ggplot(data = blood_protein_62283, aes(x=group_list, y=as.numeric(OPA1), color=group_list))+
  geom_boxplot()+labs(y="OPA1",title = "PRO:GSE62283")
grid.arrange(p3,p4,p5,p2,p1,p6,nrow=2,ncol=4)
g <- arrangeGrob(p3,p4,p5,p2,p1,p6,nrow=2,ncol=4) #generates g
ggsave(file="pic/OPA1.png", g) 

