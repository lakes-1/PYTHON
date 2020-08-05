setwd("D:\\BaiduNetdiskDownload\\bioinformatics_data\\PD\\data_duplicate\\result\\01_preprocess")
label = read.csv("df_all_label.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F)
age = read.csv("data_age_sex.txt",header = T,sep = "\t",fill=T,stringsAsFactors=F)
label1=age[match(label$dataset,age$Sample.name),]
write.table(label1,file="label12.txt",row.names=FALSE,col.names=TRUE,sep="\t")
## 后在excel中把label和label1合并在一块


## 根据患者计算年龄分布
install.packages('plyr')
library('plyr')
label_konw=label[label$age!='UNK',]
label_konw$age=as.numeric(label_konw$age)
data_class<-ddply(label_konw,c('class'),summarise,
             N=length(age),
             age_mean=mean(age),
             age_sd=sd(age),
             age_median=median(age),
             age_max=max(age),
             age_min=min(age)
             )

data_fullclass<-ddply(label_konw,c('full_class'),summarise,
                  N=length(age),
                  age_mean=mean(age),
                  age_sd=sd(age),
                  age_median=median(age),
                  age_max=max(age),
                  age_min=min(age)
)

data_ct=label_konw[label_konw$class=='CT',]
ct_co=data_ct[data_ct$full_class%in%c('CT','CO'),]
data_ctco<-ddply(ct_co,c('class'),summarise,
                      N=length(age),
                      age_mean=mean(age),
                      age_sd=sd(age),
                      age_median=median(age),
                      age_max=max(age),
                      age_min=min(age)
)

ct_cy=data_ct[data_ct$full_class%in%c('CT','CY'),]
data_ctcy<-ddply(ct_cy,c('class'),summarise,
                 N=length(age),
                 age_mean=mean(age),
                 age_sd=sd(age),
                 age_median=median(age),
                 age_max=max(age),
                 age_min=min(age)
)
