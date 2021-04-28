##定义函数
zi.pi<-function(nodes_bulk, z.bulk, modularity_class, degree){

z.bulk[abs(z.bulk)>0]<-1
module<-which(colnames(nodes_bulk)==modularity_class)
module.max<-max(nodes_bulk[,module])
degree<-which(colnames(nodes_bulk)==degree)

#按照模块将相关矩阵分割
bulk.module<-list(NA)
length(bulk.module)<-module.max

for(i in 1:max(nodes_bulk[,module])){
	bulk.module[[i]]<-z.bulk[which(nodes_bulk[,module]==i),which(nodes_bulk[,module]==i)]
	bulk.module[[i]]<-as.data.frame(bulk.module[[i]])
	rownames(bulk.module[[i]])<-rownames(z.bulk)[which(nodes_bulk[,module]==i)]
	colnames(bulk.module[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
}

# within-module degree z
z_bulk<-list(NA)
length(z_bulk)<-module.max

for(i in 1:length(z_bulk)){
	z_bulk[[i]]<-bulk.module[[i]][,1]
	z_bulk[[i]]<-as.data.frame(z_bulk[[i]])
	colnames(z_bulk[[i]])<-"z"
	rownames(z_bulk[[i]])<-rownames(bulk.module[[i]])
}

#计算z值
for(i in 1:max(nodes_bulk[,module])){
	if(length(bulk.module[[i]])==1){
		z_bulk[[i]][,1]<-0
	}else if(sum(bulk.module[[i]])==0){
		z_bulk[[i]][,1]<-0
	}else{
		k<-rowSums(bulk.module[[i]])
		mean<-mean(k)
		sd<-sd(k)
		if (sd==0){
			z_bulk[[i]][,1]<-0
		}else{
			z_bulk[[i]][,1]<-(k-mean)/sd
		}
	}
}

#z值合并
for(i in 2:max(nodes_bulk[,module])) {
	z_bulk[[i]]<-rbind(z_bulk[[i-1]],z_bulk[[i]])
}
z_bulk<-z_bulk[[module.max]]

#按照模块将相关矩阵列分割
bulk.module1<-list(NA)
length(bulk.module1)<-module.max

for(i in 1:max(nodes_bulk[,module])){
	bulk.module1[[i]]<-z.bulk[,which(nodes_bulk[,module]==i)]
	bulk.module1[[i]]<-as.data.frame(bulk.module1[[i]])
	rownames(bulk.module1[[i]])<-rownames(z.bulk)
	colnames(bulk.module1[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
}

#among-module connectivity c
c_bulk<-list(NA)
length(c_bulk)<-module.max

for(i in 1:length(c_bulk)){
	c_bulk[[i]]<-z.bulk[,1]
	c_bulk[[i]]<-as.matrix(c_bulk[[i]])
	colnames(c_bulk[[i]])<-"c"
	rownames(c_bulk[[i]])<-rownames(z.bulk)
	c_bulk[[i]][,1]<-NA
}

#每个节点各模块连接数平方
for(i in 1:max(nodes_bulk[,module])){
	c_bulk[[i]]<-rowSums(bulk.module1[[i]])
	c_bulk[[i]]<-as.matrix(c_bulk[[i]])
	c_bulk[[i]]<-c_bulk[[i]]*c_bulk[[i]]
	colnames(c_bulk[[i]])<-"c"
	rownames(c_bulk[[i]])<-rownames(z.bulk)
}

#平方和
for(i in 2:max(nodes_bulk[,module])){
	c_bulk[[i]]<-c_bulk[[i]]+c_bulk[[i-1]]
}
c_bulk<-c_bulk[[module.max]]

c_bulk1<-1-(c_bulk/(nodes_bulk[,degree]*nodes_bulk[,degree]))
colnames(c_bulk1)<-"c"

#z,c整合
z_c_bulk<-c_bulk1
z_c_bulk<-as.data.frame(z_c_bulk)
z_c_bulk$z<-z_bulk[match(rownames(c_bulk1),rownames(z_bulk)),]
z_c_bulk<-z_c_bulk[,c(2,1)]
names(z_c_bulk)[1:2]<-c('within_module_connectivities','among_module_connectivities')

z_c_bulk$nodes_id<-rownames(z_c_bulk)
nodes_bulk$nodes_id<-rownames(nodes_bulk)
z_c_bulk<-merge(z_c_bulk,nodes_bulk,by='nodes_id')
z_c_bulk

}
