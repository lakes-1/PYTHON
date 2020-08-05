# 定义目录
dir="D://BaiduNetdiskDownload//bioinformatics_data//PD/data_duplicate/"
setwd(dir)
# 拿到当前文件夹的所有txt文件
fileNames=list.files(dir,pattern = ".txt")

# 定义函数，返回每个数据每个样本的基因表达的均值、方差、最大、最小，小于0的基因个数
summary_statatic=function(fileNames,i){
  data = as.matrix(read.csv(file = fileNames[i],header = T,sep = "\t",fill=T,stringsAsFactors=F))
  data_noGeneInfo=data[,4:dim(data)[2]]
  mode(data_noGeneInfo)='numeric'
  temp=apply(data_noGeneInfo,2,function(x){
   
    c(N_length=dim(data_noGeneInfo)[1],value_mean=mean(x),value_sd=sd(x),value_median=median(x),
      value_max=max(x),value_min=min(x),value_lesss0=sum(x<0),value_more0_min=min(x[x>0]))
      })

  return(temp)
}
result_summary=matrix(,8,1)
for (i in 1:length(fileNames)){
  result_subset=summary_statatic(fileNames,i)
  result_summary=cbind(result_summary,result_subset)
} 
result_summary=result_summary[,-1]



# plot the difference disturbution between min and max
dataset <- data.frame(value =c(result_summary[5,],result_summary[2,],abs(result_summary[6,])), 
                      group = rep(c("max","mean","min"),each=dim(result_summary)[2])
)
boxplot( value ~ group, dataset)


dataset <- data.frame(value =c(result_summary[2,],abs(result_summary[6,])), 
                      group = rep(c("mean","min"),each=dim(result_summary)[2])
)
boxplot( value ~ group, dataset)

sum(result_summary[1,])