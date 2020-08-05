rm(list = ls())
## 1.找到蛋白质部分的基因
options(stringsAsFactors=F)
library(limma)
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")
load("R_data/01_degs.Rdata")
deg1=brain_protein_68719
deg2=blood_protein_62283
fdr_cutoff=0.2
deg1=deg1[deg1$adj.P.Val<fdr_cutoff,]
# caution!!!the result had change the rank of gene.
gene_name4=gene_name4[match(rownames(deg2),rownames(gene_name4)),1]
gene_name4=gene_name4[deg2$adj.P.Val<fdr_cutoff]
deg2=deg2[deg2$adj.P.Val<fdr_cutoff,]
dire_gene=intersect(rownames(deg1),gene_name4)
length(unique(gene_name4[gene_name4%in%dire_gene]))#128,等于dire_gene个数，即交叠的差异基因，刚好都只有一个

# View(gene_name4[gene_name4%in%dire_gene])
a=deg1[dire_gene,]$logFC*deg2[match(dire_gene,gene_name4),]$logFC
b=dire_gene[which(a<0)]
brainVerseBloodProtein=deg1[b,]


length(unique(gene_name4[gene_name4%in%rownames(brainVerseBloodProtein)]))





### 找rna部分
### 就是得注意GSE62283的探针问题，还有就是合并时保留group信息
############-protein consistant with rna show  inverse relation between PD and CT############

deg_62283_loc=gene_name4[match(rownames(brainVerseBloodProtein),gene_name4)]#
# 1.加载工作区间，给他转置
load("R_data/05_data_exp.Rdata")
blood_protein_62283=blood_protein_62283
blood_rna_99039=t(blood_rna_99039)
brain_protein_68719=t(brain_protein_68719)
brainSN_rna_20292=t(brainSN_rna_20292)
brainSn_rna_7621=t(brainSn_rna_7621)
brain_rna_68719=t(brain_rna_68719)

# 2.找在所有数据集都出现过的基因，还得保留group信息
data_all_list=list(blood_rna_99039,brain_protein_68719,brainSN_rna_20292,brainSn_rna_7621,
          brain_rna_68719)
con_gene=colnames(data_all_list[[1]])

# 2.1 提取共有基因，合并表达谱
for (i in (1:length(data_all_list))){
  con_gene=intersect(con_gene,colnames(data_all_list[[i]]))
  print(length(con_gene))
}
con_gene=intersect(con_gene,c(deg_62283_loc,"group_list"))
data_all_list_con=sapply(data_all_list, function(x){x[,con_gene]})
deg_62283_loc=as.matrix(deg_62283_loc)[match(setdiff(con_gene,"group_list"),deg_62283_loc),1]
blood_protein_62283=blood_protein_62283[,c(names(deg_62283_loc),"group_list")]
colnames(blood_protein_62283)=c(con_gene)
blood_rna_99039=as.data.frame(data_all_list_con[[1]])
brain_protein_68719=as.data.frame(data_all_list_con[[2]])
brainSN_rna_20292=as.data.frame(data_all_list_con[[3]])
brainSn_rna_7621=as.data.frame(data_all_list_con[[4]])
brain_rna_68719=as.data.frame(data_all_list_con[[5]])


# 先归一化


blood_protein_62283=cbind(apply(blood_protein_62283[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))}),blood_protein_62283$group_list)
brain_protein_68719=cbind(apply(brain_protein_68719[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))}),brain_protein_68719$group_list)
blood_rna_99039=cbind(apply(blood_rna_99039[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))}),blood_rna_99039$group_list)
brainSN_rna_20292=cbind(apply(brainSN_rna_20292[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))}),brainSN_rna_20292$group_list)
brainSn_rna_7621=cbind(apply(brainSn_rna_7621[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))}),brainSn_rna_7621$group_list)
brain_rna_68719=cbind(apply(brain_rna_68719[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))}),brain_rna_68719$group_list)

blood_protein_62283[,1:73]=apply(blood_protein_62283[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))})
brain_protein_68719[,1:73]=apply(brain_protein_68719[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))})
blood_rna_99039[,1:73]=apply(blood_rna_99039[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))})
brainSN_rna_20292[,1:73]=apply(brainSN_rna_20292[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))})
brainSn_rna_7621[,1:73]=apply(brainSn_rna_7621[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))})
brain_rna_68719[,1:73]=apply(brain_rna_68719[,-74], 2, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))})




colnames(blood_protein_62283)[74]="group_list"
colnames(brain_protein_68719)[74]="group_list"
colnames(blood_rna_99039)[74]="group_list"
colnames(brainSN_rna_20292)[74]="group_list"
colnames(brainSn_rna_7621)[74]="group_list"
colnames(brain_rna_68719)[74]="group_list"


## 批量画图了

i=1
library(ggplot2)
library(gridExtra)# split the grid 
gene_name=con_gene
for (i in 1:(dim(blood_protein_62283)[2]-1)){
  

p1=ggplot(data = as.data.frame(blood_rna_99039), aes(x=group_list, y=as.numeric(blood_rna_99039[,i]), color=group_list))+
  geom_boxplot()+labs(y="",title="RNA:GSE99039")+theme(legend.position="none")
p2=ggplot(data = as.data.frame(brain_protein_68719), aes(x=group_list, y=as.numeric(brain_protein_68719[,i]), color=group_list))+
  geom_boxplot()+labs(y="",title="PRO:GSE68719")
p3=ggplot(data = as.data.frame(brain_rna_68719), aes(x=group_list, y=as.numeric(brain_rna_68719[,i]), color=group_list))+
  geom_boxplot()+labs(y=gene_name[i],title = "RNA:GSE68719")+theme(legend.position="none")
p4=ggplot(data = as.data.frame(brainSN_rna_20292), aes(x=group_list, y=as.numeric(brainSN_rna_20292[,i]), color=group_list))+
  geom_boxplot()+labs(y="",title="RNA:GSE20292")+theme(legend.position="none")
p5=ggplot(data = as.data.frame(brainSn_rna_7621), aes(x=group_list, y=as.numeric(brainSn_rna_7621[,i]), color=group_list))+
  geom_boxplot()+labs(y="",title = "RNA:GSE7621")+theme(legend.position="none")
p6=ggplot(data = as.data.frame(blood_protein_62283), aes(x=group_list, y=as.numeric(blood_protein_62283[,i]), color=group_list))+
  geom_boxplot()+labs(y="",title = "PRO:GSE62283")
p6
#grid.arrange(p3,p4,p5,p2,p1,p6,nrow=2,ncol=4)
g <- arrangeGrob(p3,p4,p5,p2,p1,p6,nrow=2,ncol=4) #generates g
save_file=paste("pic/",gene_name[i],".png",sep = "")
ggsave(file=save_file, g) 

}

# gse99039
group_list1=blood_rna_99039$group_list
which
temp=apply(blood_rna_99039[,-74],2,function(x){median(x[which(group_list1=="PD")])-median(x[which(group_list1=="CT")])})
