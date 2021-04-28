One_network <- function(normtab,tax,rank='OTU',abundance = 0.005,appear='none',num='none',cor_type = 'spearman',
                        P.adj = 'BH',rt=0.7,pt=0.05){
### 准备工作
  # setwd('~/Test_Data/18S')
  # library(ggplot2)
  # library(ggtreeExtra)
  # library(ggtree)
  # library(phyloseq)
  # library(reshape2)
  # library(dplyr)
  # library(ape)
  # library(tidyr)
  # library(ggprism)
  # tax_sum <- read.delim('taxonomy.txt',header = T,stringsAsFactors = F,row.names = 1)
  # tax_filter <- read.delim('without_unknow_taxonomy.txt',header = T,stringsAsFactors = F,row.names = 1)
  # otu <- read.delim('otu_table.txt',header = T,stringsAsFactors = F,row.names = 1)
  # design <- read.delim('metadata.txt',header = T,stringsAsFactors = F,row.names = 1)
  # ptree <- read.tree('tree.nwk')
  # # 设置中文字体
  # par(family=MyFonts[1])
  # # 构建phyloseq对象
  # ps <- phyloseq(
  #   OTU=otu_table(otu,taxa_are_rows = TRUE),
  #   TAX <- tax_table(as.matrix(tax_sum)),
  #   sampledata <- sample_data(design),
  #   ptree
  # )
  # 
  # # 2 抽平 
  # library(vegan)
  # ## 先对数据进行转制
  # otu_t <- t(otu)
  # head(otu_t)
  # ## 计算抽平的阈值
  # raremin <- min(rowSums(otu_t))
  # ## 抽平
  # set.seed(888)
  # otu_rare <- rrarefy(otu_t,raremin)
  # ## 查看抽平后sum是否一致
  # rowSums(otu_rare)
  # ## 稀释曲线
  # rarecurve(otu_rare,step=1000,sample=raremin,xlab = "Reads", ylab = "Species",label = F,main="抽平后") # 抽平之后的情况
  # rarecurve(otu_t,step=10000,xlab = "Reads", ylab = "Species",label = F,col = "blue", cex = 0.6,main="抽平前") # 抽平之前的情况
  # 
  # ## 转制回去
  # otu_dat <- t(otu_rare)
  
  #可选事先过滤一些低丰度或低频的类群
  # normtab = norm_otu
  # tax = tax_sum
  # rank='OTU'
  # abundance = 0.005 # 筛选丰度 常见的还有0.0001（万分之一）
  # appear='none' # 样本出现次数阈值
  # num='none' # otu或其他分类级的总数
  # cor_type = 'spearman'
  # P.adj = 'BH'
  # rt=0.7
  # pt=0.05
  
### 正文
library(Hmisc)
library(psych)
## abundance的筛选
if (abundance != 'none'){
# filtered_otu <- normtab[which(rowSums(normtab) >= abundance), ]    #例如只保留相对丰度总和高于 0.005 的属
filtered_otu <- normtab[which(rowMeans(normtab) >= abundance), ]    #按照平均相对丰度算
}

## appear的筛选
if (appear != 'none'){
filtered_otu1 <- filtered_otu
filtered_otu1[filtered_otu1>0] <- 1
filtered_otu <- filtered_otu[which(rowSums(filtered_otu1) >= appear), ]    #例如只保留在 5 个及以上样本中出现的属
}

## num的筛选
if (num != 'none') {
filtered_otu$mean <- rowMeans(filtered_otu)
filtered_otu <- filtered_otu[order(filtered_otu$mean,decreasing = T),]
filtered_otu <- filtered_otu[1:num,!names(filtered_otu) %in% c("mean")]
filtered_otu <- filtered_otu[,!names(filtered_otu) %in% c("mean")] # 把生成的mean删除
 }

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
  if(ncol(filtered_otu)>4) {
    otu_corr <- rcorr(t(filtered_otu), type = cor_type)
    
    #阈值筛选
    #将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
    r <- otu_corr$r
    # r[is.na(r)] <- 0
    r[abs(r) < rt] <- 0  
    
    #选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
    p <- otu_corr$P
    # p[is.na(p)] <- 0
    p <- p.adjust(p, method = P.adj)    #可选 p 值校正，这里使用 BH 法校正 p 值
    p[p>=pt] <- -1
    p[p<pt & p>=0] <- 1
    p[p==-1] <- 0
  } else{
    otu_corr <- corr.test(t(filtered_otu),method = cor_type,use='complete',adjust =P.adj,ci=F)
    #阈值筛选
    #将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
    r <- otu_corr$r
    r[is.na(r)] <- 0
    r[abs(r) < rt] <- 0 
    #选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
    p <- otu_corr$p
    p[is.na(p)] <- 0
    p[p>=pt] <- -1
    p[p<pt & p>=0] <- 1
    p[p==-1] <- 0
  }

  #根据上述筛选的 r 值和 p 值保留数据
  z <- r * p
  diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
  head(z)[1:6,1:6]
  z[is.na(z)] <- 0
  # 如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
  # write.table(data.frame(z, check.names = FALSE), 'network/otu_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
  

# 15.2 igraph构建网络 --------------------------------------------------------------
## 获得网络
if(sum(z)!=0){
library(igraph)
library(dplyr)
# 将邻接矩阵转化为 igraph 网络的邻接列表
# 构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数
# z <- read.delim('network/otu_corr.matrix.txt', row.names = 1, sep = '\t', check.names = FALSE)
g <- graph.adjacency(as.matrix(z), weighted = TRUE, mode = 'undirected')
g

#自相关也可以通过该式去除
g <- simplify(g)
#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]) )

# 15.3 添加其他信息 -------------------------------------------------------------
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

# 直接将正负相关剥离开
E(g)$state <-  E(g)$correlation 
E(g)$state[E(g)$state<0] <- -1
E(g)$state[E(g)$state>0] <- 1

#为节点（otu）添加属性信息（界门纲目科属水平注释）
if (rank=='OTU'){
  filtered_tax <- tax[as.character(V(g)$name), ]
} else if(rank %in% colnames(tax)){
  pos_id <- which(colnames(tax)==rank)
  delt <- ncol(tax_sum)-pos_id
  filtered_tax <- tax[,1:pos_id] %>% dplyr::distinct() # 对列去冗余
  ## 补充分类级低于rank的列为UN
  for (i in 1:delt){
    filtered_tax[[colnames(tax)[pos_id+i]]] <- 'UN'
  }
  rownames(filtered_tax) <- filtered_tax[[rank]] # 重新上行名
  filtered_tax <- filtered_tax[as.character(V(g)$name), ] ## 排序
}

V(g)$kingdom <- filtered_tax$kingdom
V(g)$phylum <- filtered_tax$phylum
V(g)$class <- filtered_tax$class
V(g)$order <- filtered_tax$order
V(g)$family <- filtered_tax$family
V(g)$genus <- filtered_tax$genus


#节点度（Degree）
#由于本示例是个无向网络，故无出度和入度之分
V(g)$degree <- degree(g)
V(g)$degree

#加权度（Weighted degree）
V(g)$weight_degree <- strength(g)
V(g)$weight_degree

#接近中心性（Closeness centrality）
V(g)$closeness_centrality <- closeness(g)
V(g)$closeness_centrality

#介数中心性（Betweenness centrality）
V(g)$betweenness_centrality <- betweenness(g)
V(g)$betweenness_centrality

#特征向量中心性（Eigenvector centrality）
V(g)$eigenvector_centrality <- evcent(g)$vector
V(g)$eigenvector_centrality

#边介数中心性（Edge betweenness centrality）
E(g)$betweenness_centrality <- edge.betweenness(g)
E(g)$betweenness_centrality

#模块性（modularity），详见 ?cluster_fast_greedy，?modularity，有多种模型
fc <- cluster_fast_greedy(g)
modularity <- modularity(g, membership(fc))


#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(g)$modularity <- membership(cluster_fast_greedy(g))

# 计算模块内连通度（Zi）和模块间连通度（Pi）
source('~/MyScript/zi_pi.r')
#节点属性列表，包含节点所划分的模块
nodes_list <- data.frame(
  nodes_id = V(g)$name, 
  degree = V(g)$degree, 
  modularity = V(g)$modularity
)
rownames(nodes_list) <-  nodes_list$nodes_id
#两个文件的节点顺序要一致
z2 <- z[V(g)$name,V(g)$name]
nodes_list <- nodes_list[rownames(z2), ]
#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, z2, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

V(g)$within_module_connectivities <- zi_pi$within_module_connectivities
V(g)$among_module_connectivities <- zi_pi$among_module_connectivities

## 汇总节点表和边的表格
nodes <- data.frame(
  nodes_id = V(g)$name, 
  degree = V(g)$degree, 
  modularity = V(g)$modularity,
  weight_degree = V(g)$weight_degree,
  closeness_centrality = V(g)$closeness_centrality,
  betweenness_centrality = V(g)$betweenness_centrality,
  eigenvector_centrality = V(g)$eigenvector_centrality,
  within_module_connectivities = V(g)$within_module_connectivities,
  among_module_connectivities = V(g)$among_module_connectivities,
  kingdom = V(g)$kingdom,
  phylum = V(g)$phylum,
  class = V(g)$class,
  order = V(g)$order,
  family = V(g)$family,
  genus = V(g)$genus
)

edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表
edges <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation,
  state = E(g)$state,
  betweenness_centrality = E(g)$betweenness_centrality
)


# 网络特征计算 ------------------------------------------------------------------

#节点数量（number of nodes）和边数量（number of edges）
nodes_num <- length(V(g))
nodes_num

edges_num <- length(E(g))
edges_num

#平均度（average degree）
average_degree <- mean(degree(g))
average_degree

#平均加权度（average weighted degree），仅适用于含权网络
#average_weight_degree <- mean(strength(igraph))

#节点和边连通度（connectivity）
nodes_connectivity <- vertex.connectivity(g)
nodes_connectivity

edges_connectivity <- edge.connectivity(g)
edges_connectivity

#平均路径长度（average path length）
average_path_length <- average.path.length(g, directed = FALSE)
average_path_length

#网络直径（diameter）
graph_diameter <- diameter(g, directed = FALSE)
graph_diameter

#图密度（density）
graph_density <- graph.density(g)
graph_density

#聚类系数（clustering coefficient）
clustering_coefficient <- transitivity(g)
clustering_coefficient

#介数中心性（betweenness centralization)
betweenness_centralization <- centralization.betweenness(g)$centralization
betweenness_centralization

#度中心性（degree centralization）
degree_centralization <- centralization.degree(g)$centralization
degree_centralization

#模块性（modularity），详见 ?cluster_fast_greedy，?modularity，有多种模型
modularity

#同配混合（assortative mixing），例如
# otu_class <- read.delim('node_attribute.txt', row.names = 1, stringsAsFactors = FALSE)
# V(g)$group <- otu_class[V(g)$name,'group']
# assortativity.nominal(g, (V(g)$group == 'class2')+1, directed = FALSE)

#互惠性（reciprocity），仅适用于有向网络
#reciprocity(g, mode = 'default')
#reciprocity(g, mode = 'ratio')

#选择部分做个汇总输出
network_character <- data.frame(
  nodes_num,    #节点数量（number of nodes）
  edges_num,    #边数量（number of edges）
  average_degree,    #平均度（average degree)
  nodes_connectivity,    #节点连通度（vertex connectivity）
  edges_connectivity,    #边连通度（edges connectivity）
  average_path_length,    #平均路径长度（average path length）
  graph_diameter,    #网络直径（diameter）
  graph_density,    #图密度（density）
  clustering_coefficient,    #聚类系数（clustering coefficient）
  betweenness_centralization,    #介数中心性（betweenness centralization)
  degree_centralization,    #度中心性
  modularity    #模块性（modularity）
)
network_character

# -------------------------------------------------------------------------




res <- list(g=g,nodes_list=nodes,edges_list=edges,network_character=network_character,z=z)
} else res = 'Unable to form a network'

return(res)

}