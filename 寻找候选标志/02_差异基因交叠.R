rm(list = ls())
library(limma)
setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data\\brainTissueProtein")
load("R_data/01_degs.Rdata")
load("R_data/011_degs.Rdata")
# 免疫基因的载入
immune_file="immune_genes.csv"
immune= read.csv(file =immune_file, header = T,sep = ",",fill=T,stringsAsFactors=F)
immune_gene=immune[,1]


############-RNA is inverse between PD and CT############
# all blood RNA degs overlap 
# all brain's sn rna with blood and in the immune_gene,fdr=0.2
fdr_cutoff=0.2
# fdr_cutoff=0.05
twodegs_direction=function(deg1,deg2,dire_con=TRUE,fdr_cutoff=0.2){
  
  # the differentiate express genes 
  deg1=deg1[deg1$adj.P.Val<fdr_cutoff,]
  deg2=deg2[deg2$adj.P.Val<fdr_cutoff,]
  dire_gene=intersect(rownames(deg1),rownames(deg2))
  a=deg1[dire_gene,]$logFC*deg2[dire_gene,]$logFC
  # overlalp and the direation are same or verse 
  if(dire_con)
  {
    b=dire_gene[which(a>0)]
    gene_direction=deg1[b,]
  }
  else{
    
    b=dire_gene[which(a<0)]
    gene_direction=deg1[b,]
  }
  return(gene_direction)
}
brainSN_rna_con=twodegs_direction(brainSN_rna_20292,brainSn_rna_7621,fdr_cutoff = fdr_cutoff)
# brain sn + ba9 verse blood and in the immune_gene,fdr=0.2
brainSN_rna_con=brain_rna_68719
brain_con=twodegs_direction(brainSN_rna_con,brain_rna_68719,fdr_cutoff = fdr_cutoff)
brainVerseBloodRNA=twodegs_direction(brain_con,blood_rna_99039,dire_con = FALSE,fdr_cutoff = fdr_cutoff)



############-protein is inverse between PD and CT############
deg1=brain_protein_68719
deg2=blood_protein_62283
deg1=deg1[deg1$adj.P.Val<fdr_cutoff,]
# caution!!!the result had change the rank of gene.
gene_name4=gene_name4[match(rownames(deg2),rownames(gene_name4)),1]
gene_name4=gene_name4[deg2$adj.P.Val<fdr_cutoff]
deg2=deg2[deg2$adj.P.Val<fdr_cutoff,]
dire_gene=intersect(rownames(deg1),gene_name4)
length(unique(gene_name4[gene_name4%in%dire_gene]))#128,等于dire_gene个数，即交叠的差异基因，刚好都只有一个

# View(gene_name4[gene_name4%in%dire_gene])
a=deg1[dire_gene,]$logFC*deg2[match(dire_gene,gene_name4),]$logFC
b=dire_gene[which(a<0)]
brainVerseBloodProtein=deg1[b,]


############-protein consistant with rna show  inverse relation between PD and CT############
brainVerseBlood=twodegs_direction(brainVerseBloodProtein,brainVerseBloodRNA)
brainVerseBlood
deg_62283_loc=gene_name4[match(rownames(brainVerseBlood),gene_name4)]#
# deg_62283_loc=c("B19R15C07","B02R07C07")
save(brainVerseBlood,blood_protein_62283,deg_62283_loc,file = "R_data/step02_intersect_degs.Rdata")


brain_rna_con_immune=intersect(rownames(brainVerseBlood),immune_gene)
brain_rna_con_immune