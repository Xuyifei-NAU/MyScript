#!/usr/bin/env Rscript

# Author Xuyifei

# 1.2 解析命令行
# 设置清华源加速下载
#site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse")
  require("optparse",character.only=T) 
}

# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list <- c("vegan")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 解析参数-h显示帮助信息
option_list <- list(
  make_option(c("-i", "--otutab"), type="character", default="",
              help="Input reads count file; such as OTU table [default %default]"),
  make_option(c("-t", "--taxtab"), type="character", default="",
              help="Input taxonomy table [default %default]"),
  make_option(c("-d", "--design"), type="character", default="",
              help="Input design table [default %default]"),
  make_option(c("-p", "--position"), type="numeric", default="",
              help="Which column is the treatment in design table [default %default]"),
  make_option(c("-T", "--threshold"), type="numeric", default=1,
              help="Threshold of abundance percentage, such as 0.1% [default %default]"),
  make_option(c("-o", "--output"), type="character", default="",
              help="Output directory and filename [default %default]"),
  make_option(c("-r", "--relativeab"), type="character", default="",
              help="Output directory and filename [default %default]"),
  make_option(c("-a", "--outtax"), type="character", default=1,
              help="Output directory and filename [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The input feature table is ", opts$input,  sep = ""))
print(paste("Threshold of abundance percentage: ", opts$threshold,  sep = ""))
print(paste("Output directory and filename: ", opts$output,  sep = ""))

# 读取文件 --------------------------------------------------------------------
otutab = read.table(opts$otutab, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
design <- read.table(opts$design,header = T,row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
tax <- read.table(opts$taxtab,header = T,row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
#otutab <- read.table('~/Test_Data/18S/otu_table.txt',header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)

# 筛选统计 ----------------------------------------------------------------------


otutab1 <- t(otutab)
otutab1 <- otutab1/rowSums(otutab1)
otutab1 <- t(otutab1)
filter_result <- otutab1[which(rowSums(otutab1)/length(colnames(otutab))>=(opts$threshold)/100 ),]*100

# OTU丰度筛选阈值，默认0.1%，主要为圈图展示合适数据的OTU

## 提取otu的id
id <- rownames(filter_result)
len <- length(id)
filter_tax <- tax[rownames(tax)%in%id,]
filter_otu <- otutab[rownames(otutab)%in%id,]
#filter_result$sum 
#design

#t(filter_result)

## 合并设计表与筛选后的丰度表
dat_m <- merge(t(filter_result),design,by='row.names')[,-1]
sum <- aggregate(dat_m[,c(1:len)],by= list(dat_m[,opts$position+len-1]),FUN=mean)
#mean <- aggregate(dat_m[,c(1:len)],by= list(dat_m[,1+len]),FUN=mean)
rownames(sum) <- sum[,1]
sum <- sum[,-1]
## 合并丰度以及分类
filter_tax <- merge(filter_tax,t(sum),by='row.names')
rownames(filter_tax) <- filter_tax[,1]
filter_tax <- filter_tax[,-1]
## 添加一个otu在所有样本中的总丰度
#filter_tax$mean <- rowSums(filter_result)/len
filter_tax <- merge(filter_tax, data.frame(mean=rowSums(filter_result)/(2*len)),by='row.names')
rownames(filter_tax) <- filter_tax[,1]
filter_tax <- filter_tax[,-1]


# 输出文件 --------------------------------------------------------------------

#write.table(filter_result, file=paste(opts$output,sep=""),append = F, quote = F, eol = "", row.names = T, col.names = T)
write.table(data.frame ("OUTID"= rownames(filter_otu), filter_otu),file=paste(opts$output,sep="") ,append = F, sep = '\t',quote = F,  row.names = F, col.names = T)
write.table(data.frame ("OUTID"= rownames(filter_tax), filter_tax),file=paste(opts$outtax,sep="") ,append = F, sep = '\t',quote = F,  row.names = F, col.names = T)
write.table(data.frame ("OUTID"= rownames(filter_result), filter_result),file=paste(opts$relativeab,sep="") ,append = F, sep = '\t',quote = F,  row.names = F, col.names = T)



# 分割线 ---------------------------------------------------------------------

#
# #
# setwd("~/Test_data/18S")
# design <- read.table('metadata.txt',header = T,row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# otutab = read.table('otu_table.txt', header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# tax <- read.table('taxonomy.txt',header = T,row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# 
# 
# otutab <- t(otutab)
# otutab <- otutab/rowSums(otutab)
# otutab <- t(otutab)
# filter_result <- otutab[which(rowSums(otutab)/length(colnames(otutab))>=(opts$threshold)/100 ),]
# 
# ## 提取otu的id
# id <- rownames(filter_result)
# len <- length(id)
# filter_tax <- tax[rownames(tax)%in%id,]
# 
# ##filter_result$sum
# ##design
# 
# ##t(filter_result)
# 
# ## 合并设计表与筛选后的丰度表
# dat_m <- merge(t(filter_result),design,by='row.names')[,-1]
# sum <- aggregate(dat_m[,c(1:len)],by= list(dat_m[,1+len]),FUN=sum)
# ##mean <- aggregate(dat_m[,c(1:len)],by= list(dat_m[,1+len]),FUN=mean)
# rownames(sum) <- sum[,1]
# sum <- sum[,-1]
# ## 合并丰度以及分类
# filter_tax <- merge(filter_tax,t(sum),by='row.names')
# rownames(filter_tax) <- filter_tax[,1]
# filter_tax <- filter_tax[,-1]
# ## 添加一个otu在所有样本中的总丰度
# #filter_tax$mean <- rowSums(filter_result)/len
# filter_tax <- merge(filter_tax, data.frame(mean=rowSums(filter_result)/len),by='row.names')
# rownames(filter_tax) <- filter_tax[,1]
# filter_tax <- filter_tax[,-1]


# -------------------------------------------------------------------------


# library(vegan)
# test1 <- data.frame(
#   S1=c(1,2,3),
#   S2=c(4,5,6),
#   S3=c(7,8,9)
# )
# rownames(test1) <- c('O1','O2','O3')
# 
# norm = t(t(test1)/colSums(test1,na=T)) * 100
# 
# (test1 <- t(test1))
# (test1 <- test1/rowSums(test1))


