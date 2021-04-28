One_wilcox <- function(data,group,compare,value){
  # data=rep_nw_uT
  # group='G'
  # compare='group'
  # value='posi_pro'
  
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (ele in type){
    # sub_dat <- subset(data,group == i)
    sub_dat <- data[data[,group]==ele,]
   
    ## 两两比较
    GP <- unique(sub_dat[[compare]])
    group1 <- NULL
    group2 <- NULL
    p <- NULL
    for (i in 1:(length(GP) - 1)) {
      for (j in (i + 1):length(GP)) {
        group1 <- c(group1, GP[i])
        group2 <- c(group2, GP[j])
        group_ij <- subset(sub_dat, sub_dat[[compare]] %in% c(GP[i], GP[j]))
        pos <- which(colnames(group_ij) %in% compare)
        group_ij[,pos] <- factor(group_ij[,pos], levels = c(GP[i], GP[j]))
        
        wilcox_test <- wilcox.test(group_ij[[value]]~group_ij[[compare]],  alternative = 'two.sided', conf.level = 0.95)
        p <- c(p, wilcox_test$p.value)
      }
    }
    ## 整理结果
    result <- data.frame(group1, group2,  p)
    result$padj <- p.adjust(result$p, method = 'BH')   #推荐加上 p 值校正，这里使用 Benjamini 方法校正 p 值
    result$sig <- ifelse(result$p <= 0.01,'**',ifelse(result$p <= 0.05,'*', 'no sig'))
    result$sig.adj <- ifelse(result$padj <= 0.01,'**',ifelse(result$padj <= 0.05,'*', 'no sig'))

    result$type <- ele
    ## 合并结果
    a <- rbind(a,result)
    #print(a)
  }
  names(a) <- c('compare1','compare2','p','padj','sig','sigadj',group)
  return(a)
}

