
One_tax_abundance <- function(tax=tax,otu=otu,tar='',all=F){
  # tax<- read.delim('~/Test_Data/18S/taxonomy.txt',header = T,stringsAsFactors = F,row.names = 1)
  # otu <- read.delim('~/Test_Data/18S/otu_table.txt',header = T,stringsAsFactors = F,row.names = 1)
  # tax=tax
  # otu=otu
  # tar='genus'
  all_tab <- merge(tax,otu,by='row.names') ## 合并otu和tax
  start <- dim(tax)[2]+2 # otu表开始位置
  end <- dim(all_tab)[2] # otu表结束位置
  
  if (all==FALSE){
    if (tar=='') tar <- names(tax)[2]
    
    select_tab <- data.frame(all_tab[[tar]],all_tab[,start:end]) ## 合并数据集
    na.omit(select_tab) # 删除NA值所在的行
    select_tab[,1] <- factor(select_tab[,1]) ## 转化为因子型
    res <- aggregate(select_tab[,-1],by=list(select_tab[,1]),sum)
    names(res)[1] <- tar # 重命名
    
    rownames(res) <- res[,1]
    res <- as.data.frame(res)[,-1]
    
    return(res)
  } else {
    class <- names(tax)
    res_list <- list()
    for (tar in class){
      select_tab <- data.frame(all_tab[[tar]],all_tab[,start:end]) ## 合并数据集
      na.omit(select_tab) # 删除NA值所在的行
      select_tab[,1] <- factor(select_tab[,1]) ## 转化为因子型
      res <- aggregate(select_tab[,-1],by=list(select_tab[,1]),sum)
      names(res)[1] <- tar # 重命名
      
      rownames(res) <- res[,1]
      res <- as.data.frame(res)[,-1]
      
      res_list[[tar]] <- res
    }
    return(res_list)
  }

}
# res1 <- One_tax_abundance(tax=tax,otu=otu,tar='phylum')
# res2 <- One_tax_abundance(tax=tax,otu=otu,all=T)
# res2[['phylum']]


One_filter_abundance <- function(otu=otu,mt=0.01,num=""){
  # otu <- read.delim('~/Test_Data/18S/otu_table.txt',header = T,stringsAsFactors = F,row.names = 1)
  # otu=otu
  # mt=0.01
  # num=""
  
  ## 计算相对丰度
  re_ab <-  t(t(otu)/colSums(otu)) 
  ab_all <- data.frame(abundance=rowMeans(re_ab),ID=rownames(re_ab))
  
  ## 默认按照丰度阈值进行提取
  if (num ==""){
  ## 按照丰度进行提取
  filter_ID <- ab_all[ab_all$abundance>=mt,][,2]
  } else {
  ## 按照数目进行提取
  filter_ID <-ab_all [order(ab_all$abundance,decreasing = T),][1:num,2]
  }
  
  
  ## 根据id提取表格
  ID <- rownames(otu) %in% filter_ID
  final_table <- re_ab[ID,] %>% as.data.frame()
  return(final_table)
}

# res1 <- One_filter_abundance(otu,mt=0.001)
# res2 <- One_filter_abundance(otu,num=200)


One_filter_abundance_profile <- function(otu=otu,div=1000,th_num=500){
  ## otu <- read.delim('~/Test_Data/18S/otu_table.txt',header = T,stringsAsFactors = F,row.names = 1)
  ## 计算相对丰度
  re_ab <-  t(t(otu)/colSums(otu)) 
  ab_all <- data.frame(abundance=rowMeans(re_ab),ID=rownames(re_ab))
  
  m1 <- max(ab_all$abundance)
  m2 <- min(ab_all$abundance)
  
  mt_l <- round(seq(m2,m1,(m1-m2)/div),5)
  for (mt in mt_l ) {
    
    filter_ID <- ab_all[ab_all$abundance>=mt,][,2]
    len <- length(filter_ID)
    if (mt==mt_l[1]){
      res_list <- data.frame('Abundance threshold'=mt,'Numbers'=len)
    } else{
      res_list <- rbind(res_list,data.frame('Abundance threshold'=mt,'Numbers'=len))
    }
  }
  
  ## 去掉 number 小于20的
  ## th_num =200
  res_list <- subset(res_list,Numbers>th_num)
  p <- ggplot(data=res_list,aes(`Abundance.threshold`,Numbers))+geom_bar(stat = 'identity')
  
  return(list(res_list,p))
}

# res1 <- One_filter_abundance_profile(otu,div=1000,th_num=200)
# res1[[1]]
# res1[[2]]