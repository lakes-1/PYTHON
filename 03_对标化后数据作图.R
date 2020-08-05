rm(list=ls())
pwd="D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\ppa_result"
setwd(pwd)
load(file = 'R_data//step1_data_label.Rdata')

# 1.DATA PROCESS 
# 归一化(针对log2后的数据)
data_colsum=colSums(data_all)
data_all_nor=10**6*data_all/data_colsum

# 针对log2前的数据
data_colsum=colSums(data_all_raw)
data_all_raw_nor=10**6*data_all_raw/data_colsum

# 画图
label=label[match(colnames(data_all_raw_nor),label$dataset),]


# 2.LABEL PROCESS 
# 控制顺序，我只能用数字来控制顺序了，啊哈哈哈
data=data_all_raw_nor
if (T){
# pca的顺序
group_list=label$class
group_list[which(group_list=="CT")]="1.CT"
group_list[which(group_list=="MCI")]="2.MCI"
group_list[which(group_list=="PD")]="3.PD"
group_list[which(group_list=="AD")]="4.AD"
group_list[which(group_list=="MS")]="5.MS"
group_list[which(group_list=="BC")]="6.BC"
label$class=group_list
# 热图的顺序
group_list=label$full_class
group_list[group_list%in%c("CT","CY","CO")]="1.CT"
group_list[group_list%in%c("EMCI")]="2.1EMCI"
group_list[group_list%in%c("LMCI")]="2.2LMCI"
group_list[group_list%in%c("ESPD")]="3.1ESPD"
group_list[group_list%in%c("MMPD")]="3.2MMPD"
group_list[group_list%in%c("ESAD")]="4.1EMMAD"
group_list[group_list%in%c("LSAD")]="4.2LMMAD"
group_list[group_list%in%c("MS")]="5.MS"
group_list[group_list%in%c("RRMS")]="5.1RRMS"
group_list[group_list%in%c("SPMS")]="5.2SPMS"
group_list[group_list%in%c("BC0-2")]="6.BC0-2"
group_list[group_list%in%c("BC3-4")]="6.BC3-4"
label$full_class=group_list
}


# 3.PLOT 
# the path of output pictures 
name="pic/all_gene"

# 3.1.PCA图
  
  if (T) {
    dat=data
    dat=as.data.frame(scale(t(dat)))
    # 画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换。格式要求data.frame
    library("FactoMineR")# 计算PCA
    library("factoextra")# 画图展示
    
    dat.pca <- PCA(dat, graph = F)
    # fviz_pca_ind按样本  fviz_pca_var按基因
    fviz_pca_ind(dat.pca,
                 geom.ind = "point", # c("point", "text)2选1
                 col.ind = label$class, # color by groups
                 # palette = c("#00AFBB", "#E7B800"),# 自定义颜色
                 addEllipses = T, # 加圆圈
                 legend.title = "Groups"# 图例名称
    )
    #p + scale_fill_discrete(limits=c("CT", "MCI", "PD","AD","BC","MS"))
    temp=paste(name,'PCA_noNormallization.png')
    ggsave(temp)
  }
  
# 3.2.热图
  
  if (T) {
    dat=data
    #cg=names(tail(sort(apply(dat,1,sd)),1000))# 找到SD最大的1000个基因
    library(pheatmap)
    
    #headmap.dat=dat[cg,]
    headmap.dat=dat
    # 把每个基因在不同处理和重复中的数据转换为平均值为0，
    # 方差为1的数据，所以在这里也需要先转置(针对基因)。
    headmap.dat=t(scale(t(headmap.dat)))
    headmap.dat[headmap.dat>2]=2 
    headmap.dat[headmap.dat< -2]= -2
    
    
    
    # 分组信息设置
    group_list=label
    ac=data.frame(group_six=group_list$class,dataset=group_list$batch,group_fourteen=group_list$full_class)
    rownames(ac)=colnames(headmap.dat)
    temp=paste(name,'heatmap_top1000_sd.png')
    ## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
    pheatmap(headmap.dat,show_colnames =T,show_rownames = F,
             annotation_col=ac,##annotation_col就是分组信息
             filename = temp)
    
  }
  
# 3.3.相关图
  
  if (T) {
    exprSet=data
    
    # 利用所有基因画相关性热图
    temp=paste(name,"cor_all_genes.png")
    if (T) {
      
      
      pheatmap::pheatmap(cor(exprSet),
                         annotation_col = ac,
                         show_rownames = F,
                         filename = temp)
    }
  }

save(data_all_nor,data_all_raw_nor,label,file="./R_data/step3_dataNor_label.Rdata")



# 4.MCI、PD、CT
data=data_all_raw_nor[,label$class%in%c("1.CT","2.MCI","3.PD")]
label=label[label$class%in%c("1.CT","2.MCI","3.PD"),]
name="pic/all_gene_ctMciPD"


# 5.MCI\PD\CT\AD
data=data_all_raw_nor[,label$class%in%c("1.CT","2.MCI","3.PD",'4.AD')]
label=label[label$class%in%c("1.CT","2.MCI","3.PD","4.AD"),]
label=label[order(label$class),]
data1=data[,match(label$dataset,colnames(data))]
name="pic/all_gene_ctMciPdAd"
if (T) {
  dat=data1
  #cg=names(tail(sort(apply(dat,1,sd)),1000))# 找到SD最大的1000个基因
  library(pheatmap)
  
  #headmap.dat=dat[cg,]
  headmap.dat=dat
  # 把每个基因在不同处理和重复中的数据转换为平均值为0，
  # 方差为1的数据，所以在这里也需要先转置(针对基因)。
  headmap.dat=t(scale(t(headmap.dat)))
  headmap.dat[headmap.dat>2]=2 
  headmap.dat[headmap.dat< -2]= -2
  
  
  
  # 分组信息设置
  group_list=label
  ac=data.frame(group_six=group_list$class,dataset=group_list$batch,group_fourteen=group_list$full_class)
  rownames(ac)=colnames(headmap.dat)
  temp=paste(name,'heatmap_top1000_sd.png')
  ## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
  pheatmap(headmap.dat,show_colnames =T,show_rownames = F,cluster_rows = FALSE,
           cluster_cols = FALSE,annotation_col=ac,##annotation_col就是分组信息
           filename = temp)
  
}