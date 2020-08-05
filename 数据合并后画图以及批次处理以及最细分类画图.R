##导入数据
rm(list=ls())
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data_duplicate\\result/01_preprocess/")
data = read.csv("data_all.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F,row.names = 1)
label = read.csv("df_all_label.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F)
# PCA检查
# PCA，图会保存在pic/all_samples_PCA.png
exprSet=data
dat=t(data)
group_list=as.character(label$class)

# 这个data要求是行为基因，列为样本的,画了全部基因的pca和前1000个sd的热图
# plot_pca_clust 这个是比较简单的就以分组变量为标注，而group——list就直接是分组
# plot_pca_clust_mutipleLabel 就以多个分组作为标注，里面可以直接在ac修改，group——list也不是一维变量
plot_pca_clust=function(data,group_list,name){
  dat=data
  group_list=group_list
  # pca
  if (T) {
    dat=as.data.frame(scale(t(dat)))
    # 画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换。格式要求data.frame
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
    temp=paste(name,'PCA_noNormallization.png')
    ggsave(temp)
  }
  
  # 挑选1000个SD最大的基因画表达量热图
  # 结果图存放在heatmap_top1000_sd.png中
  
  if (T) {
    dat=data
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
    temp=paste(name,'heatmap_top1000_sd.png')
    ## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
    pheatmap(headmap.dat,show_colnames =T,show_rownames = F,
             annotation_col=ac,##annotation_col就是分组信息
             filename = temp)
    
  }
  
  # 利用top500_mad基因画相关性热图热图
  # 结果图存放在cor_top500_mad.png中
  if (T) {
    exprSet=data
    
    # 利用所有基因画相关性热图
    temp=paste(name,"cor_all_genes.png")
    if (T) {
      ac=data.frame(group=group_list)
      rownames(ac)=colnames(exprSet)
      pheatmap::pheatmap(cor(exprSet),
                         annotation_col = ac,
                         show_rownames = F,
                         filename = temp)
    }
    
    # 取top1000 mad 的基因画图
    exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:1000]),]
    temp=paste(name,"cor_top1000_mad.png")
    if (T) {
      M=cor(exprSet) 
      pheatmap::pheatmap(M,
                         show_rownames = F,
                         annotation_col = ac,
                         filename = temp)
    }
  }
  
  
}
plot_pca_clust_mutipleLabel=function(data,group_list,name){
  dat=data
  
  # pca
  if (T) {
    dat=as.data.frame(scale(t(dat)))
    # 画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换。格式要求data.frame
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
    temp=paste(name,'PCA_noNormallization.png')
    ggsave(temp)
  }
  
  # 挑选1000个SD最大的基因画表达量热图
  # 结果图存放在heatmap_top1000_sd.png中
  
  if (T) {
    dat=data
    cg=names(tail(sort(apply(dat,1,sd)),1000))# 找到SD最大的1000个基因
    library(pheatmap)
    
    headmap.dat=dat[cg,]
    # 把每个基因在不同处理和重复中的数据转换为平均值为0，
    # 方差为1的数据，所以在这里也需要先转置(针对基因)。
    headmap.dat=t(scale(t(headmap.dat)))
    headmap.dat[headmap.dat>2]=2 
    headmap.dat[headmap.dat< -2]= -2
    
    
    
    # 分组信息设置
    ac=data.frame(group_six=group_list$class,dataset=group_list$batch,group_fourteen=group_list$full_class)
    rownames(ac)=colnames(headmap.dat)
    temp=paste(name,'heatmap_top1000_sd.png')
    ## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
    pheatmap(headmap.dat,show_colnames =T,show_rownames = F,
             annotation_col=ac,##annotation_col就是分组信息
             filename = temp)
    
  }
  
  # 利用top500_mad基因画相关性热图热图
  # 结果图存放在cor_top500_mad.png中
  if (T) {
    exprSet=data
    
    # 利用所有基因画相关性热图
    temp=paste(name,"cor_all_genes.png")
    if (T) {
      ac=data.frame(group=group_list)
      rownames(ac)=colnames(exprSet)
      pheatmap::pheatmap(cor(exprSet),
                         annotation_col = ac,
                         show_rownames = F,
                         filename = temp)
    }
    
    # 取top1000 mad 的基因画图
    exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:1000]),]
    temp=paste(name,"cor_top1000_mad.png")
    if (T) {
      M=cor(exprSet) 
      pheatmap::pheatmap(M,
                         show_rownames = F,
                         annotation_col = ac,
                         filename = temp)
    }
  }
  
  
}

# 前面两个是属于只能比较死板的用setwd切换目录然后保存图片的方式，
# 后面两个则是可以根据path直接使用，但是也只是相对路径，绝对路径不行好像,我放到后面去了


## 全部数据
name='all_'
plot_pca_clust(dat,group_list,name)


## 画更细致的分类
group_list=as.character(label$full_class)
name='full_class'
plot_pca_clust(dat,group_list,name)


# 画图切目录把，文件夹快爆了
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data_duplicate\\pic\\01_prebatch")
# CT 
# PD样本
# MS样本
# BC样本
# MCI
# AD

class_name=unique(label$class)
for (i in  class_name){
  class_list=label[label$class==i,]
  group_list=class_list$full_class
  data_class=t(data[class_list$dataset,])
  name=i
  plot_pca_clust(data_class,group_list,name)
}




## 全部6个大类的批次效应
library(sva)
data_prebatch <- dat
batch_PD=label$batch
csif = as.data.frame(cbind(label$dataset,batch_PD))
colnames(csif)=c("group_name","batch")
modcombat = model.matrix(~1, data = csif)
batch = csif$batch
# 这边par.prior和priors.plot都是默认参数，其中par.prior代表是否要进行参数调整
# prior.plots is wherther to plot before combatch
combat_edata = ComBat(dat=data_prebatch, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
#去除批次效应
group_list=label$class
name='all_batch'
plot_pca_clust(combat_edata,group_list,name)


## 希望根据六个大类的子类分别做批次效应的去除,然后画图
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data_duplicate\\pic\\02_batch/")
library(sva)
class_name=unique(label$class)
pic_prefix='batch'
# 由于MCI只有一个批次所以报错，可以加个循环
for (i in  class_name){
  class_list=label[label$class==i,]
  group_list=class_list$full_class
  data_class=t(data[class_list$dataset,])
  data_prebatch=data_class
  csif = data.frame(group_name=class_list$dataset,batch=class_list$batch)
  batch = csif$batch
  if (length(unique(batch))>1){
  modcombat = model.matrix(~1, data = csif)
  combat_edata = ComBat(dat=data_prebatch, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
  name=paste(pic_prefix,i,sep = "")
  plot_pca_clust(combat_edata,group_list,name)
  }
  else {next}
}



# 观察批次和多标签的情况
plot_pca_clust=function(data,group_list,name,file_path){
  dat=data
  group_list=group_list
  # pca
  if (T) {
    dat=as.data.frame(scale(t(dat)))
    # 画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换。格式要求data.frame
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
    temp=paste(dir,name,'PCA_noNormallization.png',sep = "")
    ggsave(temp)
  }
  
  # 挑选1000个SD最大的基因画表达量热图
  # 结果图存放在heatmap_top1000_sd.png中
  
  if (T) {
    dat=data
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
    temp=paste(dir,name,'heatmap_top1000_sd.png',sep = "")
    ## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
    pheatmap(headmap.dat,show_colnames =T,show_rownames = F,
             annotation_col=ac,##annotation_col就是分组信息
             filename = temp)
    
  }
  
  # 利用top500_mad基因画相关性热图热图
  # 结果图存放在cor_top500_mad.png中
  if (T) {
    exprSet=data
    
    # 利用所有基因画相关性热图
    temp=paste(dir,name,"cor_all_genes.png",sep = "")
    
    if (T) {
      ac=data.frame(group=group_list)
      rownames(ac)=colnames(exprSet)
      pheatmap::pheatmap(cor(exprSet),
                         annotation_col = ac,
                         show_rownames = F,
                         filename = temp)
    }
    
    # 取top1000 mad 的基因画图
    exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:1000]),]
    temp=paste(dir,name,"cor_top1000_mad.png",sep = "")
    
    if (T) {
      M=cor(exprSet) 
      pheatmap::pheatmap(M,
                         show_rownames = F,
                         annotation_col = ac,
                         filename = temp)
    }
  }
  
  
}
plot_pca_clust_mutipleLabel=function(data,group_list,name,file_path){
  dat=data
  dir=file_path
  # pca
  if (T) {
    dat=as.data.frame(scale(t(dat)))
    # 画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换。格式要求data.frame
    library("FactoMineR")# 计算PCA
    library("factoextra")# 画图展示
    
    dat.pca <- PCA(dat, graph = F)
    # fviz_pca_ind按样本  fviz_pca_var按基因
    fviz_pca_ind(dat.pca,
                 geom.ind = "point", # c("point", "text)2选1
                 col.ind = group_list$full_class, # color by groups
                 # palette = c("#00AFBB", "#E7B800"),# 自定义颜色
                 addEllipses = T, # 加圆圈
                 legend.title = "Groups"# 图例名称
    )
    temp=paste(dir,name,'PCA_noNormallization.png',sep = "")
    ggsave(temp)
  }
  
  # 挑选1000个SD最大的基因画表达量热图
  # 结果图存放在heatmap_top1000_sd.png中
  
  if (T) {
    dat=data
    cg=names(tail(sort(apply(dat,1,sd)),1000))# 找到SD最大的1000个基因
    library(pheatmap)
    
    headmap.dat=dat[cg,]
    # 把每个基因在不同处理和重复中的数据转换为平均值为0，
    # 方差为1的数据，所以在这里也需要先转置(针对基因)。
    headmap.dat=t(scale(t(headmap.dat)))
    headmap.dat[headmap.dat>2]=2 
    headmap.dat[headmap.dat< -2]= -2
    
    
    
    # 分组信息设置
    ac=data.frame(group_six=group_list$class,dataset=as.numeric(group_list$batch),group_fourteen=group_list$full_class)
    rownames(ac)=colnames(headmap.dat)
    temp=paste(dir,name,'heatmap_top1000_sd.png',sep = "")
    ## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
    pheatmap(headmap.dat,show_colnames =T,show_rownames = F,
             annotation_col=ac,##annotation_col就是分组信息
             filename = temp)
    
  }
  
  # 利用top500_mad基因画相关性热图热图
  # 结果图存放在cor_top500_mad.png中
  if (T) {
    exprSet=data
    
    # 利用所有基因画相关性热图
    temp=paste(dir,name,"cor_all_genes.png",sep = "")
    if (T) {
      ac=data.frame(group_six=group_list$class,dataset=group_list$batch,group_fourteen=group_list$full_class)
      rownames(ac)=colnames(exprSet)
      pheatmap::pheatmap(cor(exprSet),
                         annotation_col = ac,
                         show_rownames = F,
                         filename = temp)
    }
    
    # 取top1000 mad 的基因画图
    exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:1000]),]
    temp=paste(dir,name,"cor_top1000_mad.png",sep = "")
    if (T) {
      M=cor(exprSet) 
      pheatmap::pheatmap(M,
                         show_rownames = F,
                         annotation_col = ac,
                         filename = temp)
    }
  }
  
  
}
dir.create("tmp") 
file_path="tmp/"
## 校正批次前
class_name=unique(label$class)
dir.create("prebatch_mutipleLabel")
file_path="prebatch_mutipleLabel/"
for (i in  class_name){
  class_list=label[label$class==i,]
  group_list=class_list
  data_class=t(data[class_list$dataset,])
  name=i
  plot_pca_clust_mutipleLabel(data_class,group_list,name,file_path)
} 


## 校正批次后
dir.create("batch_mutipleLabel") 
file_path="batch_mutipleLabel/"
for (i in  class_name){
  class_list=label[label$class==i,]
  group_list=group_list=class_list
  data_class=t(data[class_list$dataset,])
  data_prebatch=data_class
  csif = data.frame(group_name=class_list$dataset,batch=class_list$batch)
  batch = csif$batch
  if (length(unique(batch))>1){
    modcombat = model.matrix(~1, data = csif)
    combat_edata = ComBat(dat=data_prebatch, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
    name=paste(pic_prefix,i,sep = "")
    plot_pca_clust_mutipleLabel(combat_edata,group_list,name,file_path)
  }
  else {next}
}
