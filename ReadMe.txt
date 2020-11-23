2020-11-23 更新

#z_score() 
##对数据进行z_score处理
##用法：z_score(your_data)


#wdnm.sh 
##解决salmon quantmerge后文件中出现NA值（已解决）的问题
##详见：https://www.bilibili.com/read/cv7561884
##wdnm.sh

##Dunn_out()
##用来整理非参数检验(posthoc.kruskal.nemenyi.test和posthoc.kruskal.dunn.test的结果)的多重比较p值
##用法：Dunn_out(out)

##One_PERMANOVA
##进行多元置换方差分析
##otu为otu表，行名为OTUID，列名为sampleID
##group为分组信息，第一列为与otu列名对应的sampleID，第二列为分组分析
##dis_method可选 "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" ,"mahalanobis".
##perm 置换检验次数
##p.adj可选'bonferroni','holm','hochberg','hommel','fdr'或'BH','BY','none'
##输出结果为两个表格，一个是总体的差异，一个是两两比较的差异
Out <- One_PERMANOVA(otu,group,dis_method,perm,p.adj)
Out[[1]] #查看总体的差异
Out[[2]] #查看两两比较的差异



##One_ANOSIM
##进行ANOSIM分析
##otu为otu表，行名为OTUID，列名为sampleID
##group为分组信息，第一列为与otu列名对应的sampleID，第二列为分组分析
##dis_method可选 "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" ,"mahalanobis".
##perm 置换检验次数
##p.adj可选'bonferroni','holm','hochberg','hommel','fdr'或'BH','BY','none'
##plot 为TRUE则输出相对应的箱线图
##输出结果为一个文件夹，两个表格，一个是总体的差异，一个是两两比较的差异，以及总的结果图和两两比较图
Out <- One_ANOSIM(otu,group,dis_method,perm,p.adj,plot)
Out[[1]] #查看总体的差异
Out[[2]] #查看两两比较的差异


##One_MRPP
##进行MRPP分析
##otu为otu表，行名为OTUID，列名为sampleID
##group为分组信息，第一列为与otu列名对应的sampleID，第二列为分组分析
##dis_method可选 "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" ,"mahalanobis".
##perm 置换检验次数
##p.adj可选'bonferroni','holm','hochberg','hommel','fdr'或'BH','BY','none'
##输出结果为一个文件夹，两个表格，一个是总体的差异，一个是两两比较的差异，以及总的结果图和两两比较图
Out <- One_MRPP(otu,group,dis_method,perm,p.adj)
Out[[1]] #查看总体的差异
Out[[2]] #查看两两比较的差异


##downarticle_fromsogou
###爬取R教程的推送
###searchkey:搜索的关键词；output：输出的格式，默认是xlsx，可以改为csv
source('/Users/Xuyifei/MyScript/download_article.R',local = T)
downarticle_fromsogou(searchkey = 'R语言+PCA')
downarticle_fromsogou(searchkey = 'R语言+PCA',output = 'csv')
