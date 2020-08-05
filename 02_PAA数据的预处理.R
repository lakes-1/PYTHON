rm(list=ls())
pwd="D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\ppa_result"
setwd(pwd)
datanames=c("GSE29654_exp.txt","GSE62283_exp.txt","GSE74763_exp.txt","GSE95718_exp.txt")
genenames=c("GSE29654_gene.txt","GSE62283_gene.txt","GSE74763_gene.txt","GSE95718_gene.txt")
i=datanames[1]
# 

data_all_list=list()


# 1.合并表达谱
# 注意，这边合并错误同样的是因为GSE95718这一套探针比较少
# 1.1读取所有文件表达谱
for (i in (1:length(datanames))){
  data = read.csv(file =datanames[i],header = T,sep = "\t",fill=T,stringsAsFactors=F,row.names = 1)
  gene = read.csv(file =genenames[i],header = T,sep = "\t",fill=T,stringsAsFactors=F)
  control=which(gene$Description=="Internal Control")
  data=data[-control,]
  data_all_list[[i]]=data
}
# 1.2 提取共有基因，合并表达谱
con_gene=rownames(data_all_list[[1]])
for (i in (1:length(data_all_list))){
  con_gene=intersect(con_gene,rownames(data_all_list[[i]]))
  print(length(con_gene))
}
head(data_all_list[[1]][con_gene,])
data_all_list_con=sapply(data_all_list, function(x){x[con_gene,]})
data_all=do.call(cbind,data_all_list_con)



# 2.修改label以吻合数据
label_file="df_all_label.txt"
label= read.csv(file =label_file,header = T,sep = "\t",fill=T,stringsAsFactors=F,row.names = 1)
con_sample=intersect(label$dataset,colnames(data_all))
setdiff(label$dataset,con_sample)
setdiff(colnames(data_all),con_sample)
# 2.2得到原因是，ADNI部分的空格
label$dataset=gsub(" ","_",label$dataset)
# 2.3 AST的大小写问题
label$dataset=gsub("AST","Ast",label$dataset)
# 2.4 bios的.问题
label$dataset=gsub("Bios.","Bios_",label$dataset)
# 2.1 得到原因是因为表读入进来后-会变成.号，对label做部分变换即可
label$dataset=gsub("-",".",label$dataset)
# 2.5 z
label$dataset[label$dataset=="Bios_CT13"]="Bios.CT13"

# 3.fit in 2^^x
data_all_raw=2**data_all
save(data_all,data_all_raw,label,file='R_data/step1_data_label.Rdata')
