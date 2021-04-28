# dat <- read.table('~/Desktop/otherdata/vegan 1.txt',row.names = 1,header = T,stringsAsFactors = F)
# design <- read.table('~/Desktop/otherdata/metadata.tsv',row.names = 1,header = T,stringsAsFactors = F)
# dat <- merge(dat,design,by='row.names')2
# library(reshape2)
# dat <- melt(dat,id.vars = -c(2:7),variable.name = 'alpha')
# 
# dat$alpha <- as.factor(dat$alpha)
# names(dat)[7] <- 'v'
# 
# 



# 1 -----------------------------------------------------------------------


ONE_Tukey_HSD1 <- function(data,group,compare,value){
  library(multcomp)
  
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    #fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
    ## 重命名方便后面使用
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    
    fit <- aov(value ~ g1,data = sub_dat )
    Tukey_HSD = TukeyHSD(fit, ordered = TRUE, conf.level = 0.95)
    options(warn = -1)
    tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(g1 = 'Tukey')), decreasing = TRUE)
    Tukey.labels <- data.frame(Letters=tuk$mcletters$Letters, stringsAsFactors = FALSE)
    ## 提取字母分组行名为group组名
    Tukey.labels$compare = rownames(Tukey.labels)
    Tukey.labels$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
    aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    
    a <- rbind(a,merge(mean_sd,Tukey.labels,by='compare'))
  }
 
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}


# 2 -----------------------------------------------------------------------

ONE_Tukey_HSD2 <- function(data,group,compare,value){
  library(multcompView)
  
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    #fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
    ## 重命名方便后面使用
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    
    fit <- aov(value ~ g1,data = sub_dat )
    Tukey_HSD = TukeyHSD(fit, ordered = TRUE, conf.level = 0.95)
    options(warn = -1)
    tuk <- multcompLetters2(value ~ g1, Tukey_HSD$g1[,"p adj"], sub_dat)

    
    #tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(g1 = 'Tukey')), decreasing = TRUE)
    Tukey.labels <- data.frame(tuk['Letters'], stringsAsFactors = FALSE)
    ## 提取字母分组行名为group组名
    Tukey.labels$compare = rownames(Tukey.labels)
    Tukey.labels$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    
    a <- rbind(a,merge(mean_sd,Tukey.labels,by='compare'))
  }
  
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}

# 
# df <- ONE_Tukey_HSD1(data=dat,group='alpha',compare='Group',value='v')
# df <- ONE_Tukey_HSD1(dat,'alpha','Group','v')
# ggplot(df,aes(x=reorder(Group,mean),y=mean))+geom_bar(stat = 'identity')+facet_wrap(.~alpha,scales = "free_y")+
#   geom_text(data=df,aes(x=Group,y=mean+std,label=Letters))
# ggplot(dat)+geom_boxplot(aes(x=Group,y=v,fill=Group))+geom_text(data=df,aes(x=Group,y=mean+1.3*std,label=Letters))+
#   facet_wrap(.~alpha,scales = "free_y")+ labs(x='Group',y='AlphaDiv')+
#   ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
# 
# df <- ONE_Tukey_HSD2(data=dat,group='alpha',compare='Group',value='v')
# df <- ONE_Tukey_HSD2(dat,'alpha','Group','v')
# ggplot(df,aes(x=reorder(Group,mean),y=mean))+geom_bar(stat = 'identity')+facet_wrap(.~alpha,scales = "free_y")+
#   geom_text(data=df,aes(x=Group,y=mean+std,label=Letters))
# ggplot(dat)+geom_boxplot(aes(x=Group,y=v,fill=Group))+geom_text(data=df,aes(x=Group,y=mean+1.3*std,label=Letters))+
#   facet_wrap(.~alpha,scales = "free_y")+ labs(x='Group',y='AlphaDiv')+
#   ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
# 


