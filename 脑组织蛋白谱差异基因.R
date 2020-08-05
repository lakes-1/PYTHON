##GSE54282考虑去除批次效应
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data/brainTissueProtein/")
data = read.csv("GSE54282_series_matrix.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F,row.names = 1)
platform= read.csv("GSE54282_gpl17047.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F)
label = read.csv("GSE54282label.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F)
# PCA检查
# PCA，图会保存在pic/all_samples_PCA.png
exprSet=data
dat=data
group_list=label$label
group_list=label$label_tissue



# 挑选1000个SD最大的基因画表达量热图
# 结果图存放在pic/heatmap_top1000_sd.png中
if (T) {
    cg=names(tail(sort(apply(dat,1,sd)),1000))# 找到SD最大的1000个基因
    library(pheatmap)
    headmap.dat=dat[cg,]
    # 把每个基因在不同处理和重复中的数据转换为平均值为0，
    # 方差为1的数据，所以在这里也需要先转置(针对基因)。
    headmap.dat=t(scale(t(headmap.dat)))
    headmap.dat[headmap.dat>2]=2 
    headmap.dat[headmap.dat< -2]= -2
    
    # 分组信息设置
    ac=data.frame(group=group_list)
    rownames(ac)=colnames(headmap.dat)
    
    ## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
    pheatmap(headmap.dat,show_colnames =T,show_rownames = F,
             annotation_col=ac,##annotation_col就是分组信息
             filename = 'batch54282heatmap_top1000_sd_group_status.png')
  }

# 利用top500_mad基因画相关性热图热图
# 结果图存放在pic/cor_top500_mad.png中

if (T){
  
  # 利用所有基因画相关性热图
  if (T) {
    ac=data.frame(group=group_list)
    rownames(ac)=colnames(exprSet)
    pheatmap::pheatmap(cor(exprSet),
                       annotation_col = ac,
                       show_rownames = F,
                       filename = 'batch54282cor_all_genes.png')
  }
  
  # 过滤标准
  if (T) {
    dim(exprSet)
    # 过滤标准需要修改,目前是保留每一个基因在>5个人中表达量>1
    exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 5),]
    boxplot(apply(exprSet,1 , mean))
    dim(exprSet)
    # 过滤标准需要修改,目前是去除每一个基因在>5个人中表达量>12
    exprSet=exprSet[!apply(exprSet,1, function(x) sum(x>12) > 5),]
    dim(exprSet)
  }
  
  # 数据normalization处理
  if (T) {
    # 去除文库大小差异normalization，同时利用log做scale
    exprSet=log(edgeR::cpm(exprSet)+1)
    # 取top500 mad 的基因画图
    exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  }
  
  # 利用top500_mad基因画相关性热图热图
  if (T) {
    M=cor(exprSet) 
    pheatmap::pheatmap(M,
                       show_rownames = F,
                       annotation_col = ac,
                       filename = 'batch_cor_top500_mad_status.png')
  }
}

# 来个批次效应搞一搞
# library(BiocManager)
# BiocManager::install("sva")
library(sva)
data <- as.matrix(data)
csif = read.csv("GSE54282batch.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F)
modcombat = model.matrix(~1, data = csif)
batch = csif$batch
combat_edata = ComBat(dat=data, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
#去除批次效应
exprSet=combat_edata
dat=combat_edata


if (T) {
  
  # 过滤标准，这是需要修改的部分
  if (T) {
    dim(exprSet)
    # 过滤标准需要修改,目前是保留每一个基因在>5个人中表达量>1
    exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 5),]
    #boxplot(apply(exprSet,1 , mean))
    dim(exprSet)
    # 过滤标准需要修改,目前是去除每一个基因在>5个人中表达量>12
    exprSet=exprSet[!apply(exprSet,1, function(x) sum(x>12) > 5),]
    dim(exprSet)
  }
  # 去除文库大小差异normalization，同时利用log做scale
  exprSet=log(edgeR::cpm(exprSet)+1)
  dat=exprSet
  dat=as.data.frame(t(dat)) # 画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换。格式要求data.frame
  library("FactoMineR")# 计算PCA
  library("factoextra")# 画图展示
  
  dat.pca <- PCA(dat, graph = F)
  # fviz_pca_ind按样本  fviz_pca_var按基因
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # c("point", "text)2选1
               col.ind = group_list, # color by groups
               # palette = c("#00AFBB", "#E7B800"),# 自定义颜色
               addEllipses = T, # 加圆圈
               legend.title = "Groups"# 图例名称
  )
  #ggsave('all_samples_PCA.png')
}
group_list=as.character(group_list)
