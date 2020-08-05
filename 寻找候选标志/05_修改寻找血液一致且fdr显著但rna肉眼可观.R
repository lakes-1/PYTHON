rm(list = ls())
gene=c("NDUFB3","OPA1","NDUFA9","NDUFB6","NDUFS4"
       ,"NNT","RIMS3","SRPRB","SYT1","UVRAG")

## 查看基因在GSE99039的表达方向
gene=c("NDUFB3","OPA1","NDUFA9","NDUFB6","NDUFS4"
       ,"NNT","RIMS3","UVRAG")
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")
load("R_data/05_data_exp.Rdata")
blood_rna_99039=blood_rna_99039[c(gene,"group_list"),]

blood_rna_990391=rbind(t(apply(blood_rna_99039[-9,], 1, function(x){x=as.numeric(x);
(x-min(x))/(max(x)-min(x))})),blood_rna_99039[9,])


apply(blood_rna_99039[-9,],1,function(x){x=as.numeric(x);
median(x[blood_rna_990391[9,]%in%"PD"])-median(x[blood_rna_990391[9,]%in%"CT"])
})



rm(list = ls())
## 1.找到蛋白质部分的基因
options(stringsAsFactors=F)
library(limma)
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")
load("R_data/01_degs.Rdata")




# 查看这些在图中看似差异的基因为什么在表达谱中是非差异

data_all_list=list(brain_protein_68719,brainSN_rna_20292,brainSn_rna_7621,
                   brain_rna_68719)

blood_rna_99039=blood_rna_99039[gene,]
brainSN_rna_20292=brainSN_rna_20292[gene,]
brainSn_rna_7621=brainSn_rna_7621[gene,]
brain_rna_68719=brain_rna_68719[gene,]


# 顺手画个通路富集图
options(stringsAsFactors = F)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

gene_down=gene
gene_down_name <- bitr(gene_down, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                       toType = c("ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                       OrgDb = org.Hs.eg.db)#

gene_down_entrezid=gene_down_name$ENTREZID
# KEGG富集分析得到结果,这边只做了上调
if (T) {
  
  enrichKK <- enrichKEGG(gene         =  gene_down_entrezid,
                         organism     = 'hsa',
                         pvalueCutoff = 0.1,
                         qvalueCutoff =0.2)
  #save(enrichKK,file = "data/enrichkk.rdata")
  enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
}
p_adj=-log10(enrichKK@result$p.adjust)
pathway=enrichKK@result$Description
result=as.data.frame(cbind(p_adj,pathway))


ggplot(result,aes(x=reorder(pathway,as.numeric(p_adj)),round(as.numeric(p_adj),2)))+
  geom_bar(stat = "identity")+coord_flip()

       