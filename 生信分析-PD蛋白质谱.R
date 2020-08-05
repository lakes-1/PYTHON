rm(list = ls())
# 加载所有关于PD,正常数据
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data")


# 读取所有文件
filelist_ct = dir(pattern="CT.txt")
filelist_pd = dir(pattern="PD.txt")
filelist_ndc = dir(pattern="NDC.txt")
##1.读取正常样本
if (T){
file_list <- lapply(filelist_ct,function(x){
  exp<-as.matrix(read.table(x,header=T,sep='\t',stringsAsFactors=F,fill=T,quote = ""))
  rownames(exp)=exp[,1]
  exp=exp[,4:dim(exp)[2]]
  return(exp)
}
)
##读取过程中发现GSE95718这个数据的探针比较少，故而不能合并，需进一步处理
file1_probe=rownames(file_list[[1]])
file4_probe=rownames(file_list[[4]])
con_probe=intersect(file1_probe,file4_probe)
file_ct <- lapply(file_list,function(x){
  x=x[con_probe,]
  return(x)
}
)
file_ct=do.call(cbind, file_ct)

}
##2.读取pd样本
if (T){
  file_list <- lapply(filelist_pd,function(x){
    exp<-as.matrix(read.table(x,header=T,sep='\t',stringsAsFactors=F,fill=T,quote = ""))
    rownames(exp)=exp[,1]
    exp=exp[,4:dim(exp)[2]]
    return(exp)
  }
  )
  ##读取过程中发现GSE95718这个数据的探针比较少，故而不能合并，需进一步处理
  file1_probe=rownames(file_list[[1]])
  file4_probe=rownames(file_list[[7]])
  con_probe=intersect(file1_probe,file4_probe)
  file_pd <- lapply(file_list,function(x){
    x=x[con_probe,]
    return(x)
  }
  )
  file_pd=do.call(cbind, file_pd)
  
}
##3.读取NDC样本
if (T){
  file_list <- lapply(filelist_ndc,function(x){
    exp<-as.matrix(read.table(x,header=T,sep='\t',stringsAsFactors=F,fill=T,quote = ""))
    rownames(exp)=exp[,1]
    exp=exp[,4:dim(exp)[2]]
    return(exp)
  }
  )
  ##读取过程中发现GSE95718这个数据的探针比较少，故而不能合并，需进一步处理
  con_probe=intersect(file1_probe,file4_probe)
  file_ndc <- lapply(file_list,function(x){
    x=x[con_probe,]
    return(x)
  }
  )
  file_ndc=do.call(cbind, file_ndc)
  
}
#合并NDC与正常作为正常
file_allCT=cbind(file_ct,file_ndc)
#na值用整个所在数据的平均值补好了
mode(file_allCT)='numeric'
mode(file_pd)='numeric'
file_allCT[is.na(file_allCT)]=mean(file_allCT,na.rm=T)
file_pd[is.na(file_pd)]=mean(file_pd,na.rm=T)



# 找差异基因，并观察相对富集的通路,正常374个样本，PD274个样本
deg_result=matrix(1,dim(file_pd)[1],1)

for (i in 1:length(deg_result))
          {
 deg_result[i,1]=wilcox.test(as.numeric(file_allCT[i,]), as.numeric(file_pd[i,]), alternative = "two.sided")$p.value
}
pfdr=as.matrix(p.adjust(deg_result[,1],method="fdr",length(deg_result[,1])),ncol=1);
sum(pfdr<0.05)

# 利用4个聚类试下结果
RF_10<-as.matrix(read.table("forest_reg_feature_rank10.txt",header=T,sep='\t',stringsAsFactors=F,fill=T))
PD_CT=cbind(file_allCT,file_pd)
headmap_scale=t(scale(t(PD_CT)))

PD_CT_rank4=headmap_scale[RF_10[1:4,1],]
## 聚类
library(pheatmap)





headmap_log2=log2(headmap.dat)
filter_data_log2=temp

# 分组信息设置
headmap.dat=PD_CT_rank4
group_list=c(rep(0,374),rep(1,274))
ac=data.frame(group=group_list)
colnames(headmap.dat)=c(1:648)
rownames(ac)=colnames(headmap.dat)
dev.new()
pheatmap(headmap.dat,show_colnames =T,show_rownames = F,
         annotation_col=ac)
##annotation_col就是分组信息
         filename = 'heatmap_top1000_sd.png')
pheatmap(headmap.dat,cluster_cols=T,cluster_rows=T)



library(factoextra)

  ####数据是否需要进行处理？都为整数
  ####1.不进行处理
  if(T){
    
    filter_data=as.matrix(headmap.dat)
    d <- dist(t(filter_data))
    hc1 <- hclust(d,"single")
    plot(hc1)
    #plot(hc1,hang=-1,type="tirangle")     
    rect.hclust(hc1,k=2)
    
  }
