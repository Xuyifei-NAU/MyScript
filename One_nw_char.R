

# boot 节点的特征 --------------------------------------------------------------


node_char <- function(g) {
  
  #计算节点拓扑特征
  all_degree <- degree(g, mode = 'all')  #节点度
  all_betweenness <- betweenness(g, normalized = TRUE)  #节点介数中心性
  all_closeness <- closeness(g, normalized = TRUE)  #节点接近中心性
  all_transitivity <- transitivity(g, 'local', vids = V(g))  #节点连通性
  all_transitivity[is.na(all_transitivity)] <- 0
  names(all_transitivity)<- V(g)$name
  
  #Bootstrapping，10000 次自举重抽样，估计节点特征的理论分布
  set.seed(123)
  boot_degree <- replicate(10000, sample(all_degree, 1, replace = TRUE))
  boot_betweenness <- replicate(10000, sample(all_betweenness, 1, replace = TRUE))
  boot_closeness <- replicate(10000, sample(all_closeness, 1, replace = TRUE))
  boot_transitivity <- replicate(10000, sample(all_transitivity, 1, replace = TRUE))
  
  #合并数据框
  node_characterstics <- data.frame(cbind(boot_degree, boot_betweenness, boot_closeness, boot_transitivity))
  node_characterstics
}



# boot 边的特征 ---------------------------------------------------------------


edge_char <- function(g) {
  
  #计算节点拓扑特征
  all_state <- E(g)$state  #正边还是负边
  all_correlation <- E(g)$correlation  #相关性（无正负）
  all_betweenness <- E(g)$betweenness_centrality
  
  #Bootstrapping，10000 次自举重抽样，估计节点特征的理论分布
  set.seed(123)
  boot_state <- replicate(10000, sample(all_state, 1, replace = TRUE))
  boot_correlation <- replicate(10000, sample(all_correlation, 1, replace = TRUE))
  boot_betweenness <- replicate(10000, sample(all_betweenness, 1, replace = TRUE))
  
  #合并数据框
  edge_characterstics <- data.frame(cbind(boot_state, boot_correlation, boot_betweenness))
  edge_characterstics
}


# 计算自然连接度 -----------------------------------------------------------------


nc <- function(adj_matrix) {
  #获取 0-1 矩阵，1 表示节点间存在边，0 表示不存在边
  adj_matrix <- as.matrix(adj_matrix)
  adj_matrix[abs(adj_matrix) != 0] <- 1
  
  #矩阵的特征分解，获取特征值 λ
  lambda <- eigen(adj_matrix, only.values = TRUE)$values
  lambda <- sort(lambda, decreasing = TRUE)
  
  #计算“平均特征根”，获得自然连通度
  lambda_sum <- 0
  N = length(lambda)
  for (i in 1:N) lambda_sum = lambda_sum + exp(lambda[i])
  lambda_average <- log(lambda_sum/N, base = exp(1))
  lambda_average
}



# 模拟移除网络中的节点 --------------------------------------------------------------

One_boot_remove <- function(adj_mat,num=200,group_name){
  # adj_mat=nw_uT_CP
  # num=200
  # group_name='ut-cp'
  #读取网络邻接矩阵，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
  adj1=adj_mat
  
  #计算自然连通度
  natural_connectivity1 <- nc(adj1)
  
  
  #转化为 igraph 邻接列表，计算节点平均度
  g1 <- graph_from_adjacency_matrix(as.matrix(adj1), mode = 'undirected', diag = FALSE)
  average_degree1 <- mean(degree(g1))
  
  
  #随机移除 200 个节点，并计算上述 3 种网络特征
  for (i in 1:num) {
    ## 构建网络
    g1 <- graph_from_adjacency_matrix(as.matrix(adj1_remove), mode = 'undirected', diag = FALSE)
    
    #自相关也可以通过该式去除
    g1 <- simplify(g1)
    #孤立节点的删除（删除度为 0 的节点）
    g1 <- delete.vertices(g1, names(degree(g1)[degree(g1) == 0]) )
    source <- data.frame(as_edgelist(g1))[[1]]
    target <- data.frame(as_edgelist(g1))[[2]]
    nodelist <- data.frame(nodelist=as.character(c(source,target))) %>% dplyr::distinct()
    
    #在邻接矩阵中随机移除 i 个节点
    remove_node <- nodelist[sample(1:nrow(nodelist), i),]
    remove_node <- which(rownames(adj1)%in%remove_node)
    adj1_remove <- adj1[-remove_node,-remove_node]
    
    #计算自然连通度
    natural_connectivity1 <- c(natural_connectivity1, nc(adj1_remove))
    
    #计算节点平均度
    average_degree1 <- c(average_degree1, mean(degree(g1)))
  }
  
  #整理数据
  dat <- data.frame(remove_node = rep(0:num, 2),
                    variable = c(c(rep('natural_connectivity', num+1), rep('average_degree', num+1))),
                    values = c(natural_connectivity1, average_degree1),
                    network = c(rep(group_name, (num+1)*2)))
  #write.csv(dat, 'dat.csv', row.names = FALSE, quote = FALSE)
  return(dat)
}



# 随机网络构建 ------------------------------------------------------------------
One_random_network <- function(g){
  # g = nw_CT_CP$g
  degree_dist <- table(degree(g))# 不统计degree值的分布，比如4有几个，5有几个
  degree_num <- as.numeric(names(degree_dist)) ## 节点的度有哪些取值
  degree_count <- as.numeric(degree_dist) ## 节点度取值的count
  names(degree_count) <- degree_num
  degs <- rep(degree_num, degree_count) #把 degree_num 重复  degree_count次
  
  
  #获得广义随机图（构建 3 个），并转换为邻接矩阵
  # set.seed(123)
  g_rand1 <- degree.sequence.game(degs, method = 'simple')
  adj_matrix_rand <- as.matrix(get.adjacency(g_rand1))
  return(adj_matrix_rand)
}


# 随机移除部分节点计算平均自然连通度 -------------------------------------------------------

One_prop_remove <- function(adj1,prop=0.5,num=99){
  # adj1=nw_CT_CP$z
  # prop=0.5
  # num=9
  
  natural_connectivity <- ""
  
  for (i in 1:num){
    
    ## 构建网络
    g1 <- graph_from_adjacency_matrix(as.matrix(adj1_remove), mode = 'undirected', diag = FALSE)
    
    #自相关也可以通过该式去除
    g1 <- simplify(g1)
    #孤立节点的删除（删除度为 0 的节点）
    g1 <- delete.vertices(g1, names(degree(g1)[degree(g1) == 0]) )
    source <- data.frame(as_edgelist(g1))[[1]]
    target <- data.frame(as_edgelist(g1))[[2]]
    nodelist <- data.frame(nodelist=as.character(c(source,target))) %>% dplyr::distinct()
    
    #在邻接矩阵中随机移除 i 个节点
    remove_node <- nodelist[sample(1:nrow(nodelist), round(nrow(nodelist)*prop)),]
    remove_node <- which(rownames(adj1)%in%remove_node)
    adj1_remove <- adj1[-remove_node,-remove_node]
    
    #计算自然连通度
    natural_connectivity <- c(natural_connectivity, nc(adj1_remove))
  }
  return(natural_connectivity)
}


