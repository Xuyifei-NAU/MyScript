## Z-Score
####---------------来源代码---------------------------------------------
# data_A <- rnorm(100, 80, 10)     # randomly create population dataset
# data_B <- rnorm(100, 400, 100)   # randomly create population dataset
# 
# hist(data_A) #histogram
# hist(data_B) #histogram
# 
# #Calculate population mean and standard deviation
# A_data_std <- sd(data_A)*sqrt((length(data_A)-1)/(length(data_A)))
# A_data_mean <- mean(data_A)
# 
# B_data_std <- sd(data_B)*sqrt((length(data_B)-1)/(length(data_B)))
# B_data_mean <- mean(data_B)
# 
# # Provided that A got 92 and B got 610
# A_obs <- 92
# B_obs <- 610
# 
# A_Z_score <-  (A_obs - A_data_mean) / A_data_std
# B_Z_score <-  (B_obs - B_data_mean) / B_data_std
# 
# sd(1:3) ^ 2

###--------------------------------------------------------------------------------------
z_score=function(data)
  {
  data_std=sd(data)*sqrt((length(data)-1)/(length(data)))    
  data_mean <- mean(data)  
  data_z_score <- (data - data_mean) / data_std
  return(data_z_score)
}

