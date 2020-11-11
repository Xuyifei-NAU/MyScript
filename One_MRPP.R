## 改编自 小白鱼的生统笔记



One_MRPP <- function(otu,group,dis_method,perm,p.adj)
{

library(vegan)
otu <- data.frame(t(otu))

# 进行总体的mrpp分析
mrpp_result <- mrpp(otu, group[,2], distance = dis_method, permutations =perm)

# 查看结果
mrpp_result


#生成表格
write.table(data.frame(Group = 'all', Distance = dis_method, A = mrpp_result$A, 
                       Observe_delta = mrpp_result$delta, Expect_delta = mrpp_result$E.delta, P_value = mrpp_result$Pvalue), 
            file = 'MRPP.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE)

##MRPP ??????ʹ??ѭ?????���????С???????Ƚϣ???��???䣩
#?Ƽ?ʹ?? OTU ???ȱ???Ϊ???????ݣ?ÿ??ɸѡ?????????¼??????????룬?????????????????ٿ??ܵ??µľ????䶯??????????
group_name <- unique(group[,2])

mrpp_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(group, group[,2] %in% c(group_name[i], group_name[j]))
    otu_ij <- otu[group_ij[,1], ]
    mrpp_result_otu_ij <- mrpp(otu_ij, group_ij[,2], permutations = perm, distance = dis_method)	#Bray-Curtis ???????ȣ????? 999 ???û?
    mrpp_result_two <- rbind(mrpp_result_two, c(paste(group_name[i], group_name[j], sep = '/'), dis_method, mrpp_result_otu_ij$A, mrpp_result_otu_ij$delta, mrpp_result_otu_ij$E.delta, mrpp_result_otu_ij$Pvalue))
  }
}
mrpp_result_two <- data.frame(mrpp_result_two, stringsAsFactors = FALSE)
names(mrpp_result_two) <- c('group', 'distance', 'A', 'Observe_delta', 'Expect_delta', 'P_value')

# 加上校正的p值
mrpp_result_two$P_value <- as.numeric(mrpp_result_two$P_value)
mrpp_result_two <- cbind(mrpp_result_two,p.adjust(mrpp_result_two$'P_value', method = p.adj))
names(mrpp_result_two) <- c('group', 'distance', 'A', 'Observe_delta', 'Expect_delta', 'P_value',paste0('P_adj_',p.adj))

result <- list(mrpp_result,mrpp_result_two)

#????
write.table(mrpp_result_two, 'MRPP.result_two.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
return(result)


}



