One_format2lefse <- function (otutab, taxonomy, metadata, thre = 0.01, groupID = "Group", sub_group='None',rank_start='None',rank='None',
          output = "LEfSe.txt") 
{
  idx = rownames(otutab) %in% rownames(taxonomy)
  otutab = otutab[idx, ]
  tax = taxonomy[rownames(otutab), ]
  
  
  p_list = c("dplyr",'tidyverse','stringr')
  for (p in p_list) {
    if (!requireNamespace(p)) {
      install.packages(p)
    }
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
  norm = t(t(otutab)/colSums(otutab, na = T)) * 100
  idx = rowMeans(norm) > thre
  HA = norm[idx, ]
  tax = tax[rownames(HA), ]
  tax$Phylum = paste(tax$Kingdom, tax$Phylum, sep = "|")
  tax$Class = paste(tax$Phylum, tax$Class, sep = "|")
  tax$Order = paste(tax$Class, tax$Order, sep = "|")
  tax$Family = paste(tax$Order, tax$Family, sep = "|")
  tax$Genus = paste(tax$Family, tax$Genus, sep = "|")
  tax$Species = paste(tax$Genus, tax$Species, sep = "|")
  grp <- tax[rownames(tax), "Kingdom", drop = F]
  merge = cbind(HA, grp)
  HA_Kingdom = merge %>% group_by(Kingdom) %>% summarise_all(sum)
  colnames(HA_Kingdom)[1] = "Class"
  grp <- tax[rownames(tax), "Phylum", drop = F]
  merge = cbind(HA, grp)
  HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
  colnames(HA_Phylum)[1] = "Class"
  grp <- tax[rownames(tax), "Class", drop = F]
  merge = cbind(HA, grp)
  HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
  colnames(HA_Class)[1] = "Class"
  grp <- tax[rownames(tax), "Order", drop = F]
  merge = cbind(HA, grp)
  HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
  colnames(HA_Order)[1] = "Class"
  grp <- tax[rownames(tax), "Family", drop = F]
  merge = cbind(HA, grp)
  HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
  colnames(HA_Family)[1] = "Class"
  grp <- tax[rownames(tax), "Genus", drop = F]
  merge = cbind(HA, grp)
  HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
  colnames(HA_Genus)[1] = "Class"
  all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family, 
              HA_Genus)
  metadata$group = metadata[, groupID]
  
  ## 针对 rank 进行数据整理
  if (rank != 'None' & (rank %in% names(taxonomy)) & rank_start != 'None' & (rank_start %in% names(taxonomy))){
    class_id <- which(names(taxonomy)==rank) ## 获取rank的位置
    star_pos <- which(names(taxonomy)==rank_start) ## 获取rank起始位置
    
    all <- all[str_count(all$Class,"\\|") ==(class_id-1),] ## 按照rank位置拆分分类表
    
    all_split <- separate(all,col = 'Class',sep='\\|',into = names(taxonomy))[,star_pos:class_id] # 按照 ｜ 拆分列
    
    filter_id <- all_split[[rank]] != "" ## 获取最细分类级水平下，注释到的行的索引
    
    all_unite <- unite(all_split[filter_id,], "Class", sep = "|", remove = T)
    
    all <- cbind(all_unite,all[,-1][filter_id, ])
    
  }
  
  if (rank != 'None' & (rank %in% names(taxonomy)) & rank_start == 'None'){
    class_id <- which(names(taxonomy)==rank) ## 获取rank的位置

    all <- all[str_count(all$Class,"\\|") ==(class_id-1),] ## 按照rank位置拆分分类表

    all_split <- separate(all,col = 'Class',sep='\\|',into = names(taxonomy))[,class_id] # 按照 ｜ 拆分列

    all <- cbind(all_split,all[,-1])[all_split != "", ]

    names(all)[1] <- 'Class'
  }

  
  
  if(sub_group == 'None'){
    colnames(all)[2:dim(all)[2]] = as.character(metadata[colnames(all)[2:dim(all)[2]], 
    ]$group)
  }
  

  if(sub_group!='None' & (sub_group %in% names(metadata)) ){
  colnames(all)[2:dim(all)[2]] = as.character(metadata[colnames(all)[2:dim(all)[2]], 
    ]$group)
  all <- rbind(sub_class=c('sub_class',as.character(metadata[[sub_group]])),all)
  }
  
  write.table(all, file = paste(output, sep = ""), append = FALSE,
              sep = "\t", quote = F, row.names = F, col.names = T)
  
  return(all)
}


