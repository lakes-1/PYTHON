options(stringsAsFactors=F)
rm(list = ls())
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")

# the ppi with the direaction  haven't the two gene:OPA1,ndufb3
# so we observe the result of no direation PPI 
data_name="OPA1 protein table.tsv"
data_name="NDUFB3 protein table.tsv"
data=as.matrix(read.csv(file =data_name, header = T,sep = "\t",fill=T,stringsAsFactors=F))
apo1=data[,1]
gene_down=apo1
ene_down_name <- bitr(gene_down, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                      toType = c("ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                      OrgDb = org.Hs.eg.db)#

gene_down_entrezid=ene_down_name[,2]
# KEGG富集分析得到结果,这边只做了上调
if (T) {
  
  enrichKK <- enrichKEGG(gene         =  gene_down_entrezid,
                         organism     = 'hsa',
                         pvalueCutoff = 0.1,
                         qvalueCutoff =0.2)
  #save(enrichKK,file = "data/enrichkk.rdata")
}

head(enrichKK)[,1:6] 
# 打开网页看相关KEGG通路图
browseKEGG(enrichKK, 'hsa00190')

# 将数据中的entrz-id变成symbol
# 更为易读
enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
enrichKK 



## 可视化
#条带图
if (T) {
  # par(mfrow=c(2,1))
  barplot(enrichKK)
  #ggsave("pic/barplot.png")
}

#气泡图
if (T) {
  dotplot(enrichKK)
  ggsave("pic/dotplot.png")
}

#下面的图需要映射颜色，设置和示例数据一样的geneList

# 展示top5通路的共同基因，要放大看。
#Gene-Concept Network
if (T) {
  cnetplot(enrichKK, colorEdge = TRUE, circular = F)
  ggsave("pic/cnetplot.png")
  cnetplot(enrichKK, colorEdge = TRUE, circular = T)
  ggsave("pic/cnetplot_circular.png")
}


#Enrichment Map
if (T) {
  emapplot(enrichKK)
  ggsave("pic/Enrichment_Map.png")
}

#(4)展示通路关系,仅仅是针对于GO数据库结果。
goplot(enrichKK)
#(5)Heatmap-like functional classification
if (T) {
  heatplot(enrichKK)
  ggsave("pic/Enrichment_Heatmap.png")
}