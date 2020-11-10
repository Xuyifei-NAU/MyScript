Dunn_out <- function(out)
{
  b <- NULL
  for (i in 1:length(colnames(out$p.value))){
    for (j in i:length(rownames(out$p.value))){
      a <- cbind(colnames(out$p.value)[i],rownames(out$p.value)[j],out$p.value[j,i])
      b <- rbind(b,a)
    }
  }
  b <- data.frame(b)
  colnames(b) <- c('source','target','p.adj')
  
  #b$p.adj <- round(as.numeric(b$p.adj),digits=5)
  b$p.adj <- as.numeric(b$p.adj)
  sig <- NULL
  for (i in 1:length(b$p.adj)){
    if (0.01 <= b[i,3] & b[i,3] < 0.05){
      sig <- rbind(sig,'*')
    }else if(0.001 <= b[i,3] & b[i,3] < 0.01){
      sig <- rbind(sig,'**')
    }else if(b[i,3] < 0.001){
      sig <- rbind(sig,'***')
    }else if(b[i,3]>=0.05){
      sig <- rbind(sig,'NG')
    }
  }
  ## 合并返回结果
  result <- cbind(b,sig)
  return(result)
}