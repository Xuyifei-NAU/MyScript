


downarticle_fromsogou <- function(searchkey, maxpage = 2,output = 'xlsx') {
  
  
if (!requireNamespace("tidyverse", quietly=TRUE))
    install.packages("tidyverse")
if (!requireNamespace("rvest", quietly=TRUE))
    install.packages("rvest")
if (!requireNamespace("openxlsx", quietly=TRUE))
    install.packages("openxlsx")
  library(tidyverse)
  library(rvest)
  library(openxlsx)
  
  allurl <- c()
  alltitle <- c()
  allsource <- c()
  
  for (temppage in c(1:maxpage)) {
    cat(temppage, '\r')
    
    url1 <- paste0("https://weixin.sogou.com/weixin?type=2&query=",URLencode(searchkey),"&page=", temppage, "&ie=utf8")
    web <- read_html(url1)
    articlelink <- web %>% html_nodes("ul[class='news-list']") %>% 
      html_nodes("li") %>% html_nodes("h3") %>% 
      html_nodes("a[target='_blank']") %>% html_attr('href')
    articlelink <- sapply(X = articlelink, function(x){paste0('https://weixin.sogou.com', x)})
    names(articlelink) <- NULL
    articletitle <- web %>% html_nodes("ul[class='news-list']") %>% 
      html_nodes("li") %>% html_nodes("h3") %>% html_nodes("a[target='_blank']") %>% html_text()
    
    sourcetitle <- web %>% html_nodes("ul[class='news-list']") %>% 
      html_nodes("li") %>% html_nodes("a[class='account']") %>% html_text()
    
    allurl <- c(allurl, articlelink)
    alltitle <- c(alltitle, articletitle)
    allsource <- c(allsource, sourcetitle)
    
    Sys.sleep(3)
  }
  
  result <- data.frame("title" = alltitle, "url" = allurl, "source" = allsource)
  
  if(output=='csv'){
  
  write.csv(result,paste0(searchkey,'教程',Sys.Date(),'.csv'))
  }else{
  write.xlsx(result,paste0(searchkey,'教程',Sys.Date(),'.xlsx'))
  }
  return(as_tibble(result))

}









