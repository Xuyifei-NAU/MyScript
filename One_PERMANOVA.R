
## 改编自 小白鱼的生统笔记
One_PERMANOVA <- function(otu,group,dis_method,perm,p.adj)
{
  ##dis_method可选 "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" ,"mahalanobis".
  ##p.adj可选'bonferroni','holm','hochberg','hommel','fdr'或'BH','BY','none'  

library(vegan)
otu <- data.frame(t(otu))

## PERMANOVA 分析（所有分组间比较，即整体差异）

#根据 group$site 这一列分组进行 PERMANOVA 分析检验组间差异，基于 999 次置换，详情 ?adonis

#（1）直接输入 OTU 丰度表，在参数中指定距离类型
#使用dis_method距离测度
adonis_result <- adonis(otu~group[,2], group, distance = dis_method, permutations =perm)


#查看结果，上述 3 条命令所计算的内容一致，以其中一个为例

adonis_result$aov.tab 

#对于各项统计内容：Df，自由度，其值=所比较的分组数量-1；SumsOfSqs，即Sums of squares，总方差，又称离差平方和；MeanSqs，即Mean squares，均方（差）；F.Model，F检验值；R2，即Variation (R2)，方差贡献，表示不同分组对样品差异的解释度，即分组方差与总方差的比值，R2越大表示分组对差异的解释度越高；Pr (>F)，显著性p值，默认p<0.05即存在显著差异。


## 可选输出

#可选输出
otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'PERMANOVA.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')


## PERMANOVA 分析(使用循环处理，进行小分组间比较，如两组间)
##PERMANOVA 分析（使用循环处理，进行小分组间比较，如两组间）
#推荐使用 OTU 丰度表作为输入数据，每次筛选分组后重新计算样本距离，避免由于样本数减少可能导致的距离变动而造成误差
group_name <-  unique(group[,2])

adonis_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(group, group[,2] %in% c(group_name[i], group_name[j]))
    otu_ij <- otu[group_ij[,1], ]
    adonis_result_otu_ij <- adonis(otu_ij~group_ij[,2], group_ij, permutations = perm, distance = dis_method)	#指定距离测度，基于 999 次置换
    adonis_result_two <- rbind(adonis_result_two, c(paste(group_name[i], group_name[j], sep = '/'), dis_method, unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
  }
}
adonis_result_two <- data.frame(adonis_result_two, stringsAsFactors = FALSE)
names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
subset(group, group$Type =='BS')
#可选添加 p 值校正，例如 Benjamini 校正
adonis_result_two$'Pr (>F)' <- as.numeric(adonis_result_two$'Pr (>F)')
adonis_result_two <- cbind(adonis_result_two,p.adjust(adonis_result_two$'Pr (>F)', method = p.adj))
names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)',paste0('P_adj_',p.adj))

#输出
write.table(adonis_result_two, 'PERMANOVA.result_two.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

PERMANOVA.result_all <- otuput
PERMANOVA.result_two <- adonis_result_two
result <- list(PERMANOVA.result_all,PERMANOVA.result_two)
return(result)
}




