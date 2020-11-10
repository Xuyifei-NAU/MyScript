
##改编自小白鱼的生统笔记代码

One_ANOSIM <- function(otu,group,dis_method,p.adj,plot)
{
  
  ##dis_method可选 "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" ,"mahalanobis".
  ##p.adj可选'bonferroni','holm','hochberg','hommel','fdr'或'BH','BY','none'  
  
dir.create('anosim_result', recursive = TRUE)
library(vegan)
otu <- data.frame(t(otu))

#（1）直接输入 OTU 丰度表，在参数中指定距离类型
#使用 Bray-Curtis 距离测度
anosim_result <- anosim(otu, group$site, distance = 'bray', permutations = 999)

output <- NULL
output <- cbind(anosim_result$permutations,dis_method,anosim_result$signif,anosim_result$statistic)
anosim_result$signif	#p 值
# R>0，说明组间差异大于组内差异，即组间差异显著；R<0，说明组内差异大于组间差异；R值的绝对值越大表明相对差异越大
anosim_result$statistic	#R 值
colnames(output) <- c('permutations','distance','P','R')

write.table(output, file = 'anosim_result/anosim.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')


if (plot==TRUE){
#生成颜色
count <-length(unique(group[,2]))
library(randomcoloR)
set.seed(1120)
palette <- distinctColorPalette(count)
#作图展示
#pdf(paste('anosim.all.pdf', sep = ''), width = 10, height = 5)
png(paste('anosim_result/anosim.all.png', sep = ''), width = 800, height = 400)
plot(anosim_result, col = c('gray', palette))
dev.off()
}


##ANOSIM 分析（使用循环处理，进行小分组间比较，如两组间）
#推荐使用 OTU 丰度表作为输入数据，每次筛选分组后重新计算样本距离，避免由于样本数减少可能导致的距离变动而造成误差
group_name <- unique(group[,2])

anosim_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(group, site %in% c(group_name[i], group_name[j]))
    otu_ij <- otu[group_ij$names, ]
    anosim_result_otu_ij <- anosim(otu_ij, group_ij$site, permutations = 999, distance = 'bray')	#Bray-Curtis 距离测度，基于 999 次置换
    anosim_result_two <- rbind(anosim_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', anosim_result_otu_ij$statistic, anosim_result_otu_ij$signif))
    
    if (plot==TRUE){
    #每次循环输出图片
    #pdf(paste('anosim_two/anosim.', group_name[i], '_', group_name[j], '.pdf', sep = ''), width = 7, height = 5)
    png(paste('anosim_result/anosim.', group_name[i], '_', group_name[j], '.png', sep = ''), width = 600, height = 400)
    plot(anosim_result_otu_ij, col = c('gray', palette[3],palette[4])) # 3和4的颜色好看一些
    dev.off()
    }
  }
}

#带 R 值和 p 值的表格
anosim_result_two <- data.frame(anosim_result_two, stringsAsFactors = FALSE)
names(anosim_result_two) <- c('group', 'distance', 'R', 'P_value')



#可选添加 p 值校正过程，例如 Benjamini 校正
anosim_result_two$P_value <- as.numeric(anosim_result_two$P_value)
anosim_result_two <- cbind(anosim_result_two,p.adjust(anosim_result_two$P_value, method = p.adj))
names(anosim_result_two) <- c('group', 'distance', 'R', 'P_value',paste0('P_adj_',p.adj))

write.table(anosim_result_two, 'anosim_result/ANOSIM.result_two.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

result <- list(anosim_result,anosim_result_two)
return(result)
}

