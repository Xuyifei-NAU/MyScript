
One_lefse_pack <- function(taxonomy,otutab,tar,path,group,compare,method='LSD',write=F){
  library(reshape2)
  library(do)
  # ## 变量区：
  # otutab = read.table('~/Test_Data/18S/lefse/otu_table.txt', header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
  # # 2. 读取物种注释
  # taxonomy = read.table("~/Test_Data/18S/lefse/taxonomy.txt", header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F) 
  # # 3. 读取样本分组信息
  # metadata = read.table("~/Test_Data/18S/lefse/metadata.txt", header=T, row.names=1, sep="\t", comment.char="")
  # 
  # taxonomy = taxonomy
  # otutab = otutab
  # tar = 'Order'
  # path = "~/Test_Data/18S/lefse/temp/input.res"
  # group = 'temperature'
  # compare = 'precipitation'
  # method = 'LSD'
  # write = F
  
  
  taxonomy <- Replace(taxonomy,from = '-','_') ## 因为lefse的结果中，-都会变成_
  source('~/MyScript/One_tax_abundance.R')
  otu_tax_family <- One_tax_abundance(taxonomy,otutab,tar = tar) ## 按照tar水平整合counts
  norm_family <- as.data.frame((t(t(otu_tax_family)/colSums(otu_tax_family)))*100) ## 计算相对丰度
  
  
  ## 2.1 读入lefse结果 ----------
  res_temp <- read.delim(path, header=FALSE, row.names=NULL) %>% 
    subset(V3!='') %>% Replace(from = paste0(tolower(strsplit(tar,"")[[1]][1]),'_'),to='')
  
  
  ## 合并丰度表与lefse结果，随后melt数据框,再与metadata合并
  res_table <- merge(norm_family,res_temp[,c(1,3)],by.x='row.names',by.y='V1') %>% 
    melt(id.vars=c('Row.names','V3')) %>% merge(metadata,by.x='variable',by.y='row.names')
  
  ## 差异检验
  res_list <- list()
  
  gp_id <- which(names(res_table)==group)
  gp_ele <- unique(res_table[['V3']])
  
  source('~/MyScript/One_LSD.R')
  for (gp in gp_ele){
    core <- subset(res_table,V3==gp)
    if (method=='LSD') {
      res <- ONE_LSD(core,'Row.names',compare,'value')
    } else if(method=='KW'){
      res <- One_wilcox(core,'Row.names',compare,'value')
    }
    res_list[[gp]] <- res
  }
  
  
  if (write == T){
    name_temp <- merge(taxonomy,res_temp[,c(1,3)],by.y = 'V1',by.x='Family')
    dir.create('lefse_result')
    openxlsx::write.xlsx(name_temp,paste0('lefse_result/',group,'_',tar,'.xlsx'))
  }
  res_list[['table']] <- res_table
  return(res_list)
}



