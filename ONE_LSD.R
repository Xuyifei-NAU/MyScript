# library(dplyr)
# 
# list.files()
# dat <- read.delim("data_new.txt",header = T,row.names = 1,stringsAsFactors = F)
# design <- read.delim("metadata.txt",header = T)
# names(dat) <- c ("SampleId","16SrRNA" ,"aadA","tetG","tetM","sul1","intI2","bla.oxa")
# dat$SampleId <- factor(dat$SampleId,levels= unique(dat$SampleId))
# design$SampleId <- factor(design$SampleId,levels= unique(design$SampleId))
# design$group <- factor(design$group,levels= unique(design$group))
# 
# 
# ## 转换数据格式
# dat_m <- reshape2::melt(dat,id.vars=c("SampleId"),variable.name = "gene",value.name = "value")
# ## 合并处理与设计文件
# qpcrdat <- merge(dat_m,design,by="SampleId")
# ## 删除值为0的样品
# qpcrdat <- qpcrdat[qpcrdat$value != 0,]
# argsdat <- subset(qpcrdat,gene!="16SrRNA")

# ---------------------------------------------------------------------
## p1

# 
# 
# a <- data.frame(stringsAsFactors = F)
# type <- unique(argsdat$gene)
# for (i in type)
# {
#   sub_dat <- subset(argsdat,gene == i)
#   library(agricolae)
#   fit <- aov(value~SampleId,sub_dat)
#   out <- LSD.test(fit,'SampleId')
#   
#   out$groups$type <- i
#   out$groups$SampleId <- rownames(out$groups)
#   
#   a <- rbind(a,merge(out$means[,1:2], out$groups,by='value'))
# }
# 

# -------------------------------------------------------------------------
# data <- argsdat
# group="gene"
# compare <- "location"
# value <- "value"
ONE_LSD <- function(data,group,compare,value){
  library(agricolae)

a <- data.frame(stringsAsFactors = F)
type <- unique(data[,group])
for (i in type)
{
 # sub_dat <- subset(data,group == i)
  sub_dat <- data[data[,group]==i,]
 # fit <- aov(value ~ compare,sub_dat)
  fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
  out <- LSD.test(fit,'sub_dat[, compare]',p.adj='BH')
  
  out$groups$type <- i
  out$groups$compare <- rownames(out$groups)
  
  a <- rbind(a,merge(out$means[,1:2], out$groups,by='sub_dat[, value]'))
  print(a)
}
names(a) <- c('mean','std','lsd',group,compare)
return(a)
}
#res <- ONE_LSD(data=argsdat,group='gene',compare='trt',value='value')


