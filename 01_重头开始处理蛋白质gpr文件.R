setwd("D:\\搜狗高速下载\\GSE29654_RAW_subset")
dir_names=c("D:\\test\\GSE29654_RAW","D:\\test\\GSE62283_RAW/","D:\\test\\GSE74763_RAW/","D:\\test\\GSE95718_RAW/")
dir="D:\\test\\GSE29654_RAW"
data_gpr=c("GSE29654_gpr.txt","GSE62283_gpr.txt","GSE74763_gpr.txt","GSE95718_gpr.txt")
save_filenames=c("GSE29654_exp.txt","GSE62283_exp.txt","GSE74763_exp.txt","GSE95718_exp.txt")
save_genenames=c("GSE29654_gene.txt","GSE62283_gene.txt","GSE74763_gene.txt","GSE95718_gene.txt")
library(limma)
library(PAA)
i=1
for (i in 1:length(dir_names)){
  setwd(dir_names[i])
  print(getwd())
  filenames=as.matrix(list.files(dir_names[i],pattern = ".gpr.gz"))
  colnames(filenames)=c('FileName')
  write.table(filenames,file=data_gpr[i],row.names=FALSE,col.names=TRUE,sep="\t")
  gpr <- dir_names[i]
  targets=data_gpr[i]
  elist <- loadGPR(gpr.path=gpr, targets.path=targets, array.type="ProtoArray",aggregation = "mean")
  # 标化的时候已经做log2
  normalized.elist <- normalizeArrays(elist=elist, method="rlm",controls="both")
  data_EXP=normalized.elist[["E"]]
  gene=normalized.elist[["genes"]]
  sample_name=normalized.elist[["targets"]]
  write.table(data_EXP,file=save_filenames[i],row.names=FALSE,col.names=TRUE,sep="\t")
  write.table(gene,file=save_genenames[i],row.names=FALSE,col.names=TRUE,sep="\t")
  
} 



log2_signal=log2(head(elist[["E"]])-head(elist[['Eb']]))
log2_e_signal=log2(head(elist[["E"]]))

