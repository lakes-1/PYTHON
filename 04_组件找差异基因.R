rm(list=ls())
pwd="D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\ppa_result"
setwd(pwd)
load(file = 'R_data//step3_dataNor_label.Rdata')
library(limma)

# label process ,else didn't satisfiy the name of contrast matrix 
group_list=label$full_class
group_list=gsub("[1-6.]","",group_list)
group_list[group_list%in%"BC0-"]="BCearly"
group_list[group_list%in%"BC-"]="BClate"


# 1.construct the design and contrast matrix 
# 1.1 design matrix 
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
head(design)
exprSet=as.matrix(data_all_raw_nor)
mode(exprSet)="numeric"
rownames(design)=colnames(exprSet)
head(design)

# 1.2比较矩阵
res=t(combn(unique(group_list),2))
res=res[order(res[,1], res[,2]),]

contrast_class=paste(res[,1],res[,2],sep = "-")
contrast.matrix<-makeContrasts(contrasts = contrast_class,
                               levels = design)


# 2.search the degs 
# 2.1 degs 
  ##step1,lmfit起过滤的手续
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  
data_degs=list()
for (i in 1:66){
  print(i)
  tempOutput1 = topTable(fit2, coef=i, number =Inf)
  data_degs[[i]] = na.omit(tempOutput1) 
}
# 2.2 summarize the deg length
deg_length=sapply(data_degs, function(x){
  sum(x$adj.P.Val<0.05)
})
# 2.3 bind the length and group  
contrast_class_result=cbind(contrast_class,deg_length)
contrast_class_result[,1]=gsub("early","0-2",contrast_class_result[,1])
contrast_class_result[,1]=gsub("late","3-4",contrast_class_result[,1])

write.table(contrast_class_result,file="R_data/contrast_class_result.txt",row.names=FALSE,col.names=TRUE,sep="\t")
