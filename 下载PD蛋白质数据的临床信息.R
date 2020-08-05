#安装并加载GEOmirror包，去github下载zip，再本地安装
rm(list = ls())
if(T){library(devtools)
  #install_github("jmzeng1314/GEOmirror")
  #devtools::install('D:\\BaiduNetdiskDownload\\GEOmirror-master.zip\\')
  #devtools::install_local("D:\\BaiduNetdiskDownload\\GEOmirror-master.zip\\")
  #devtools::install_local("D:\\BaiduNetdiskDownload\\jmzeng1314-idmap3-59efe3f.tar.gz\\")
  library(GEOmirror)
  library(idmap1)
  library(idmap2)
  library(idmap3)
  library("GEOquery")
}



gset <- getGEO(a[1],
               AnnotGPL = F,     ## 注释文件
               getGPL = F)    
gset=gset[[1]]
pd_62283=pData(gset)
write.table(pd_62283,file="GSE62283pdata.txt",row.names=TRUE,col.names=TRUE,sep="\t")
pd_74763=pData(gset)
write.table(pd_74763,file="GSE74763pdata.txt",row.names=TRUE,col.names=TRUE,sep="\t")
pd_95718=pData(gset)
write.table(pd_95718,file="GSE95718pdata.txt",row.names=TRUE,col.names=TRUE,sep="\t")

#探针注释，提取四个芯片共有的探针数据
a='GSE62283'
a='GSE74763'
a='GSE95718'
