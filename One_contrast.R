
# DESeq2 ------------------------------------------------------------------
One_DESeq2 <- function(dds,design,treatment,sig,padj)
{
# setwd('~/Test_Data/18S')
# tax_sum <- read.delim('taxonomy.txt',header = T,stringsAsFactors = F,row.names = 1)
# tax_filter <- read.delim('without_unknow_taxonomy.txt',header = T,stringsAsFactors = F,row.names = 1)
# otu <- read.delim('otu_table.txt',header = T,stringsAsFactors = F,row.names = 1)
# design <- read.delim('metadata.txt',header = T,stringsAsFactors = F,row.names = 1)
# ptree <- read.tree('tree.nwk')

  # head(design) # 实验设计表格
  # head(otus) # 最最原始的otu表格
  # # 判别otu表和设计表样本是否对应
  # idx <- rownames(design)%in%colnames(otus)
  # subdesign <- design[idx,][,c(1,2,4)]
  # count <- otus[,idx]
  # ## 将筛选好的otu表和设计表转化为矩阵形式（或者不转也行？）然后构建DESeq对象
  # count <- as.matrix(count)
  # subdesign <- as.matrix(subdesign)
  # dds_prec <- DESeqDataSetFromMatrix(countData = count,colData = subdesign,design = ~ precipitation)
  # ## DESeq函数进行分析
  # dds_prec2 <- DESeq(dds_prec)


# dds=dds_prec2
# design=design
# treatment='precipitation'
# sig=0.05
# padj='BH'

library(DESeq2)
trts =unique(design[,which(colnames(design)==treatment)])

for (i in 1:length(trts))
{
  for (j in 1:length(trts))
  {
    if (j != i)
    {
      result <- results(dds,contrast = c(treatment,trts[i],trts[j]),alpha = sig,pAdjustMethod=padj)
      result = result[order(result$padj),] # 按照p值排序
      write.table(data.frame ("OUTID"= rownames(result), result),file=paste0('DESeq2_',trts[i],'-',trts[j],".txt"),sep="\t",quote=F,row.names=F)
      if (i==1&j==2)
      {
        result$contrast <- paste0(trts[i],'-',trts[j])
        aa <- as.data.frame(result)
      }
      else 
      {
        result$contrast <- paste0(trts[i],'-',trts[j])
        aa <- rbind(aa,as.data.frame(result))
      }
    }
  }
}
return(aa)
}




#  edgeR ------------------------------------------------------------------


One_edgeR <- function(fit,design_matrix,padj,sig)
{
  
  
  library(edgeR)
  library(dplyr)
  # head(design) # 实验设计表格
  # head(otus) # 最最原始的otu表格
  # 
  # # 判别otu表和设计表样本是否对应
  # idx <- rownames(design)%in%colnames(otus)
  # subdesign <- design[idx,][,c(1,2,4)]
  # count <- otus[,idx]
  # 
  # # 整理DGEList 
  # d_prec <- DGEList(counts = count,group = subdesign$precipitation) %>% calcNormFactors()
  # # 整理设计分组矩阵
  # design.prec <- model.matrix(~0+d_prec$samples$group)
  # colnames(design.prec) <-levels(d_prec$samples$group)
  # 
  # ## 进行分析
  # fit_prec <- estimateGLMCommonDisp(d_prec,design.prec) %>% estimateTagwiseDisp() %>% glmFit(design.prec)
  # 
  
  
  
# fit = fit_prec
# design_matrix=design.prec
# padj='BH'
# sig=0.05

trts <- colnames(design_matrix)

for (i in 1:length(trts))
{
  for (j in 1:length(trts))
  {
    if (j != i)
    {
      contrast_list <- makeContrasts(contrasts = paste0(trts[i],'-',trts[j]),levels = design_matrix)
      LRD <- glmLRT(fit,contrast =contrast_list )
      res <- decideTestsDGE(LRD,p.value = sig,adjust.method = padj,lfc = 0)
      
      res_final <- as.data.frame(LRD$table)
      res_final <- merge(res_final,as.data.frame(res),by='row.names')
      res_final <- cbind(res_final,padj=p.adjust(res_final$PValue,method = padj))
      res_final$level <- ifelse(res_final[,6]==1,'enriched',ifelse(res_final[,6]==-1,'depleted','nosig')) # 报错
      write.table(data.frame (res_final),file= paste0('edgeR_',trts[i],'-',trts[j],".txt"),sep="\t",quote=F,row.names=F)
      
      if (i==1&j==2)
      {
        res_final$contrast <- paste0(trts[i],'-',trts[j])
        bb <- as.data.frame(res_final[,-6])
      } 
      else 
        {
        res_final$contrast <- paste0(trts[i],'-',trts[j])
        bb <- rbind(bb,as.data.frame(res_final[,-6]))
        }
    }
  }
}
return(bb)
}



