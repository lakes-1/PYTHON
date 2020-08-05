rm(list=ls())
pwd="D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\ppa_result"
setwd(pwd)
load(file = 'R_data//step3_dataNor_label.Rdata')
library(limma)

# 1.subset the CT,MCI,PD and the correspodant label 
group_list=label$class
group_list=gsub("[1-6.]","",group_list)
data_subset=data_all_raw_nor[,which(group_list%in%c("CT","MCI","PD"))]
group_list=group_list[group_list%in%c("CT","MCI","PD")]


# 2.search degs 
contrast_class=c("PD-MCI","PD-CT","CT-MCI")
search_limma_degs=function(data_subset,group_list,contrast_class){
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  head(design)
  exprSet=as.matrix(data_subset)
  mode(exprSet)="numeric"
  rownames(design)=colnames(exprSet)
  head(design)
  
  # 2.2 contrast.matrix
  
  contrast.matrix<-makeContrasts(contrasts = contrast_class,
                                 levels = design)
  
  # 3.search the degs 
  # 3.1 degs 
  ##step1,lmfit起过滤的手续
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  
  data_degs=list()
  for (i in 1:length(contrast_class)){
    print(i)
    tempOutput1 = topTable(fit2, coef=i, number =Inf)
    data_degs[[i]] = na.omit(tempOutput1) 
  }
  return(data_degs)
}
data_degs=search_limma_degs(data_subset = data_subset,group_list,contrast_class)
deg_pd_mci=data_degs[[1]]
deg_pd_ct=data_degs[[2]]
deg_ct_mci=data_degs[[3]]

# 3.the union of three class degs 
deg1=deg_ct_mci
deg2=deg_pd_ct
deg3=deg_pd_mci
fdr_cutoff=0.05
fc_cutoff=1.2
# 3.1 这个函数,enenen,没啥软子用
twodegs_union=function(deg1,deg2,fdr_cutoff=0.2,fc_cutoff=1.2){

  # the differentiate express genes
  deg1_cutoff=deg1[deg1$adj.P.Val<fdr_cutoff&deg1$logFC>fc_cutoff,]
  deg2_cutoff=deg2[deg2$adj.P.Val<fdr_cutoff&deg2$logFC>fc_cutoff,]
  dire_gene=union(rownames(deg1_cutoff),rownames(deg2_cutoff))
  gene_direction=deg1[dire_gene,]
  return(gene_direction)
}
con_1=twodegs_union(deg_ct_mci,deg_pd_ct,fdr_cutoff = fdr_cutoff)
deg1_cutoff=deg1[deg1$adj.P.Val<fdr_cutoff&deg1$logFC>fc_cutoff,]
deg2_cutoff=deg2[deg2$adj.P.Val<fdr_cutoff&deg2$logFC>fc_cutoff,]
deg3_cutoff=deg3[deg3$adj.P.Val<fdr_cutoff&deg3$logFC>fc_cutoff,]
union_gene_three=unique(c(rownames(deg1_cutoff),rownames(deg2_cutoff),rownames(deg3_cutoff)))




#union_gene_three=twodegs_union(con_1,deg_pd_mci,fdr_cutoff = fdr_cutoff)

# 4.write the table 
write.table(union_gene_three,file="R_data/ct_pd_mci_uniondegs.txt",row.names=TRUE,col.names=TRUE,sep="\t")
write.table(data_subset,file="R_data/ct_pd_mci_exp.txt",row.names=TRUE,col.names=TRUE,sep="\t")
write.table(group_list,file="R_data/ct_pd_mci_label.txt",row.names=FALSE,col.names=TRUE,sep="\t")


data_noGeneInfo=data_all_raw_nor
temp=apply(data_noGeneInfo,2,function(x){
  
  c(N_length=dim(data_noGeneInfo)[1],value_mean=mean(x),value_sd=sd(x),value_median=median(x),
    value_max=max(x),value_min=min(x),value_lesss0=sum(x<0),value_more0_min=min(x[x>0]))
})
