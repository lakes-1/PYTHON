setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\")
gene_down=c("WDR48",'RABL4','AKT1','C19orf18','RALGPS2','GSTO2','DCAMKL2','CARTPT',
'DERL2','MGC24125','C6orf81','PRKACA','MAP3K5')
gene_up=c("")
options(stringsAsFactors = F)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


gene_down_name <- bitr(gene_down, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db)#
gene_down_setdiff=c(11020,166614,728175,221481)
gene_down_search=c(3439,3447,25837,3558,3101,54937, 3112)
gene_down_entrezid=setdiff(c(gene_down_name[,2],gene_down_setdiff,gene_down_search),55103)
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
  browseKEGG(enrichKK, 'hsa05150')
  
  # 将数据中的entrz-id变成symbol
  # 更为易读
  enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  enrichKK 



## 可视化
#条带图
if (T) {
  # par(mfrow=c(2,1))
  barplot(enrichKK,showCategory=20)
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
