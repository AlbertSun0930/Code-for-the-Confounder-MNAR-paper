rm(list = ls())
library(BB)
library(ggplot2)
library(mice)
library(parallel)
library(doParallel)

####################################################################################################
### For the results in the boxplots: 

generate_data_OCPC <- function(sampleN) {
  Tgamma = c(-0.5,1,1)
  Tbeta = c(1,3,1,-1)
  Talpha = c(1,-2,1,3)
  N = sampleN ## Sample size 500. We can set N=500 when evaluating sample size 500.
  C1 = rnorm(N,0,1)
  C2 = rbinom(N,1,0.5)
  P_A = exp(Tgamma[1]+Tgamma[2]*C1+Tgamma[3]*C2)/(1+exp(Tgamma[1]+Tgamma[2]*C1+Tgamma[3]*C2))
  A = ifelse(P_A>runif(N,0,1),1,0)
  Y = rnorm(N, Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1+Tbeta[4]*C2, 1)
  P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y))
  RC = ifelse(P_RC>runif(N,0,1),1,0)
  RC[is.na(RC)]=1
  O_C1 = ifelse(RC==1,C1,NA)
  list(C1 = C1, C2 = C2, A = A, Y = Y, RC = RC, O_C1 = O_C1)
}

generate_data_OMPC <- function(sampleN) {
  Tgamma = c(-0.5,1,1)
  Tbeta = c(-1,3,2,-3)
  Talpha = c(1,-2,1,3)
  N = sampleN ## Sample size 500. We can set N=500 when evaluating sample size 500.
  C1 = rnorm(N,0,1)
  C2 = rbinom(N,1,0.5)
  P_A = exp(Tgamma[1]+Tgamma[2]*C1+Tgamma[3]*C2)/(1+exp(Tgamma[1]+Tgamma[2]*C1+Tgamma[3]*C2))
  A = ifelse(P_A>runif(N,0,1),1,0)
  Y = rnorm(N, Tbeta[1]+Tbeta[2]*A+Tbeta[3]*exp(C1*C2)+Tbeta[4]*C1*C2, 1)
  P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y))
  RC = ifelse(P_RC>runif(N,0,1),1,0)
  RC[is.na(RC)]=1
  O_C1 = ifelse(RC==1,C1,NA)
  list(C1 = C1, C2 = C2, A = A, Y = Y, RC = RC, O_C1 = O_C1)
}

generate_data_OCPM <- function(sampleN) {
  Tgamma = c(-3,3,3)
  Tbeta = c(1,3,1,-1)
  Talpha = c(1,-2,1,3)
  N = sampleN ## Sample size 500. We can set N=500 when evaluating sample size 500.
  C1 = rnorm(N,0,1)
  C2 = rbinom(N,1,0.5)
  P_A = exp(Tgamma[1]+Tgamma[2]*exp(C1*C2)+Tgamma[3]*C1*C2)/(1+exp(Tgamma[1]+Tgamma[2]*exp(C1*C2)+Tgamma[3]*C1*C2))
  A = ifelse(P_A>runif(N,0,1),1,0)
  Y = rnorm(N, Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1+Tbeta[4]*C2, 1)
  P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y))
  RC = ifelse(P_RC>runif(N,0,1),1,0)
  RC[is.na(RC)]=1
  O_C1 = ifelse(RC==1,C1,NA)
  list(C1 = C1, C2 = C2, A = A, Y = Y, RC = RC, O_C1 = O_C1)
}

generate_data_OMPM <- function(sampleN) {
  Tgamma = c(-3,3,3)
  Tbeta = c(-1,3,2,-3)
  Talpha = c(1,-2,1,3)
  N = sampleN ## Sample size 500. We can set N=500 when evaluating sample size 500.
  C1 = rnorm(N,0,1)
  C2 = rbinom(N,1,0.5)
  P_A = exp(Tgamma[1]+Tgamma[2]*exp(C1*C2)+Tgamma[3]*C1*C2)/(1+exp(Tgamma[1]+Tgamma[2]*exp(C1*C2)+Tgamma[3]*C1*C2))
  A = ifelse(P_A>runif(N,0,1),1,0)
  Y = rnorm(N, Tbeta[1]+Tbeta[2]*A+Tbeta[3]*exp(C1*C2)+Tbeta[4]*C1*C2, 1)
  P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y))
  RC = ifelse(P_RC>runif(N,0,1),1,0)
  RC[is.na(RC)]=1
  O_C1 = ifelse(RC==1,C1,NA)
  list(C1 = C1, C2 = C2, A = A, Y = Y, RC = RC, O_C1 = O_C1)
}

process_data <- function(data, times, verbose = TRUE, file = "output.txt"){
  C1 <- data$C1
  C2 <- data$C2
  A <- data$A
  Y <- data$Y
  RC <- data$RC
  O_C1 <- data$O_C1
  
  fun <- function(alpha) { 
    f <- numeric(length(alpha)) 					
    f[1] <- sum(exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2 - alpha[4]*Y)*ifelse(RC==1,1,0)) - sum(ifelse(RC==1,0,1))
    f[2] <- sum(exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2 - alpha[4]*Y)*ifelse(RC==1,1,0)*A) - sum(ifelse(RC==1,0,1)*A)
    f[3] <- sum(exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2 - alpha[4]*Y)*ifelse(RC==1,1,0)*Y) - sum(ifelse(RC==1,0,1)*Y)
    f[4] <- sum(exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2 - alpha[4]*Y)*ifelse(RC==1,1,0)*A*Y) - sum(ifelse(RC==1,0,1)*A*Y)
    f 
  } 
  startalpha <- c(0,-2,2,3)
  result = dfsane(startalpha,fun,quiet = TRUE)
  alpha = result$par
  
  mw = ifelse(RC==1,1,0)*(1+exp(-alpha[1]-alpha[2]*C1 - alpha[3]*C2-alpha[4]*Y))
  OM <- lm(Y~cbind(A,C1,C2), weights = mw)
  OMATEWEE <- OM$coefficients[2]
  
  PWM <- glm(A~cbind(C1,C2), weights = mw, family = 'binomial')
  ps <- predict(PWM,type = 'response')
  pw <- ifelse(A==1, 1/ps, 1/(1-ps))
  
  IPWM <- lm(Y~A, weights = mw*pw)
  IPWATEWEE = IPWM$coefficients[2]
  
  DRM = lm(Y~cbind(A,C1,C2), weights =  mw*pw)
  DRATEWEE <- DRM$coefficients[2]
  
  #complete case analysis
  OMCC <- lm(Y[RC==1]~cbind(A,C1,C2)[RC==1,])
  OMATECC <- OMCC$coefficients[2]
  PWMCC <- glm(A[RC==1]~cbind(C1,C2)[RC==1,], family = 'binomial')
  psCC <- predict(PWMCC,type = 'response')
  pwCC <- ifelse(A[RC==1]==1, 1/psCC, 1/(1-psCC))
  
  IPWMCC <- lm(Y[RC==1]~A[RC==1], weights = pwCC)
  IPWATECC = IPWMCC$coefficients[2]
  
  DRMCC = lm(Y[RC==1]~cbind(A,C1,C2)[RC==1,], weights = pwCC)
  DRATECC <- DRMCC$coefficients[2]
  
  
  #multiple imputation
  imp = mice(cbind(O_C1,C2,A,Y),m=5,method='pmm' ,print=FALSE)
  OMMI = pool(with(imp,lm(Y~cbind(A,O_C1,C2))))
  OMATEMI = OMMI$pooled[,3][2]
  
  IPW_results = numeric(5)
  for (i in 1:5) {
    data_imputed = complete(imp, i)
    PWM = glm(A ~ O_C1 + C2, data=data_imputed, family='binomial')
    ps = predict(PWM, type='response')
    pw = ifelse(data_imputed$A == 1, 1/ps, 1/(1 - ps))
    
    IPW_model = lm(Y ~ A, data=data_imputed, weights=pw)
    IPW_results[i] = summary(IPW_model)$coefficients[2, 1]
  }
  IPWATEMI = mean(IPW_results)
  
  DR_results = numeric(5)
  for (i in 1:5) {
    data_imputed = complete(imp, i)
    PWM = glm(A ~ O_C1 + C2, data=data_imputed, family='binomial')
    ps = predict(PWM, type='response')
    pw = ifelse(data_imputed$A == 1, 1/ps, 1/(1 - ps))
    
    DR_model = lm(Y ~ A + O_C1 + C2, data=data_imputed, weights=pw)
    DR_results[i] = summary(DR_model)$coefficients[2, 1]
  }
  DRATEMI = mean(DR_results)
  
  cat("Completed Monte Carlo", times, "\n", file = file, append = TRUE)
  list( as.numeric(c(OMATEWEE,IPWATEWEE,DRATEWEE,OMATECC,IPWATECC,DRATECC,OMATEMI,IPWATEMI,DRATEMI)))
}

sampleN=500
MCN = 1000
num_cores <- 70
data_info <- list(
  list(func = generate_data_OCPC, name = "results_df_OCPC_500"),
  list(func = generate_data_OCPM, name = "results_df_OCPM_500"),
  list(func = generate_data_OMPC, name = "results_df_OMPC_500"),
  list(func = generate_data_OMPM, name = "results_df_OMPM_500")
)
for (info in data_info) {
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 1)
  results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
    # 生成唯一的文件名
    output_file <- paste0("resultsrecord.txt")  # 这里可以区分文件名
    data <- info$func(sampleN)  # 调用相应的函数
    process_data(data, times, verbose = TRUE, file = output_file)
  }
  stopCluster(cl)
  results_df <- do.call(rbind, lapply(results, function(x) data.frame(t(x))))
  # 使用 eval 和 paste0 来创建变量
  eval(parse(text = paste0(info$name, " <- results_df")))
  print(info$name)
}


sampleN=2000
MCN = 1000
num_cores <- 70
data_info <- list(
  list(func = generate_data_OCPC, name = "results_df_OCPC_2000"),
  list(func = generate_data_OCPM, name = "results_df_OCPM_2000"),
  list(func = generate_data_OMPC, name = "results_df_OMPC_2000"),
  list(func = generate_data_OMPM, name = "results_df_OMPM_2000")
)
for (info in data_info) {
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 1)
  results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
    # 生成唯一的文件名
    output_file <- paste0("resultsrecord.txt")  # 这里可以区分文件名
    data <- info$func(sampleN)  # 调用相应的函数
    process_data(data, times, verbose = TRUE, file = output_file)
  }
  stopCluster(cl)
  results_df <- do.call(rbind, lapply(results, function(x) data.frame(t(x))))
  # 使用 eval 和 paste0 来创建变量
  eval(parse(text = paste0(info$name, " <- results_df")))
  print(info$name)
}



####################################################################################################
### For the results of the coverage probabilities of the confidence intervals: 

process_data <- function(data, B2, times, verbose = TRUE, file = "output.txt"){
  C1ori <- data$C1
  C2ori <- data$C2
  Aori <- data$A
  Yori <- data$Y
  RCori <- data$RC
  O_C1ori <- data$O_C1
  
  BSrec <- matrix(NA, nrow = B2, ncol=3)
  for(i in 1:B2){
    BSnum <- sample(1:500, 500, replace = TRUE)
    C1 <- C1ori[BSnum]
    C2 <- C2ori[BSnum]
    A <- Aori[BSnum]
    Y <- Yori[BSnum]
    RC <- RCori[BSnum]
    O_C1 <- O_C1ori[BSnum]
    
    fun <- function(alpha) { 
      f <- numeric(length(alpha)) 					
      f[1] <- sum(exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2 - alpha[4]*Y)*ifelse(RC==1,1,0)) - sum(ifelse(RC==1,0,1))
      f[2] <- sum(exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2 - alpha[4]*Y)*ifelse(RC==1,1,0)*A) - sum(ifelse(RC==1,0,1)*A)
      f[3] <- sum(exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2 - alpha[4]*Y)*ifelse(RC==1,1,0)*Y) - sum(ifelse(RC==1,0,1)*Y)
      f[4] <- sum(exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2 - alpha[4]*Y)*ifelse(RC==1,1,0)*A*Y) - sum(ifelse(RC==1,0,1)*A*Y)
      f 
    } 
    startalpha <- c(0,-2,2,3)+rep(rnorm(1,0,0.3),4)
    result = dfsane(startalpha,fun,quiet = TRUE)
    alpha = result$par
    
    mw = ifelse(RC==1,1,0)*(1+exp(-alpha[1]-alpha[2]*C1 - alpha[3]*C2-alpha[4]*Y))
    OM <- lm(Y~cbind(A,C1,C2), weights = mw)
    OMATEWEE <- OM$coefficients[2]
    
    PWM <- glm(A~cbind(C1,C2), weights = mw, family = 'binomial')
    ps <- predict(PWM,type = 'response')
    pw <- ifelse(A==1, 1/ps, 1/(1-ps))
    
    IPWM <- lm(Y~A, weights = mw*pw)
    IPWATEWEE = IPWM$coefficients[2]
    
    DRM = lm(Y~cbind(A,C1,C2), weights =  mw*pw)
    DRATEWEE <- DRM$coefficients[2]
    
    if (i %% 10 == 0) {
      cat("Completed bootstrap sample", i, "\n", file = file, append = TRUE)
    }
    BSrec[i,] <- c(OMATEWEE, IPWATEWEE, DRATEWEE)
  }
  
  # CI025 <- numeric(3)
  # CI975 <- numeric(3)
  # for(k in 1:3){
  #   CI025[k] <- quantile(BSrec[, k], 0.025, na.rm = TRUE)
  #   CI975[k] <- quantile(BSrec[, k], 0.975, na.rm = TRUE)
  # }
  
  CI025 <- apply(BSrec, 2, mean) - 1.96*apply(BSrec, 2, sd)
  CI975 <- apply(BSrec, 2, mean) + 1.96*apply(BSrec, 2, sd)
  
  list(CI025 = CI025, CI975 = CI975)
}

###############################################################
## OCPC
MCN = 1000
num_cores <- 70
B2 <- 500

#########################################
## sample size 500
sampleN=500
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
    # 生成唯一的文件名
    output_file <- paste0("REC_", times, ".txt") # 这里可以区分文件名
    data <- generate_data_OCPC(sampleN)  # 调用相应的函数
    process_data(data,B2, times, verbose = TRUE, file = output_file)
}
stopCluster(cl)

CI500025list = matrix(nrow = MCN, ncol = 3)
CI500975list = matrix(nrow = MCN, ncol = 3)
for (i in 1:MCN) {
  CI500025list[i,] = as.numeric(unlist(results[2*i-1]))
  CI500975list[i,] = as.numeric(unlist(results[2*i]))
}

c(mean((CI500025list[,1]<= 3)& (CI500975list[,1]>=3)), 
  mean((CI500025list[,2]<= 3)& (CI500975list[,2]>=3)),
  mean((CI500025list[,3]<= 3)& (CI500975list[,3]>=3)))

#########################################
## sample size 2000
sampleN=2000
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
  # 生成唯一的文件名
  output_file <- paste0("REC_", times, ".txt") # 这里可以区分文件名
  data <- generate_data_OCPC(sampleN)  # 调用相应的函数
  process_data(data,B2, times, verbose = TRUE, file = output_file)
}
stopCluster(cl)

CI2000025list = matrix(nrow = MCN, ncol = 3)
CI2000975list = matrix(nrow = MCN, ncol = 3)
for (i in 1:MCN) {
  CI2000025list[i,] = as.numeric(unlist(results[2*i-1]))
  CI2000975list[i,] = as.numeric(unlist(results[2*i]))
}

c(mean((CI2000025list[,1]<= 3)& (CI2000975list[,1]>=3)), 
  mean((CI2000025list[,2]<= 3)& (CI2000975list[,2]>=3)),
  mean((CI2000025list[,3]<= 3)& (CI2000975list[,3]>=3)))


###############################################################
## OCPM
MCN = 1000
num_cores <- 70
B2 <- 500

#########################################
## sample size 500
sampleN=500
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
  # 生成唯一的文件名
  output_file <- paste0("REC_", times, ".txt") # 这里可以区分文件名
  data <- generate_data_OCPM(sampleN)  # 调用相应的函数
  process_data(data,B2, times, verbose = TRUE, file = output_file)
}
stopCluster(cl)

CI500025list = matrix(nrow = MCN, ncol = 3)
CI500975list = matrix(nrow = MCN, ncol = 3)
for (i in 1:MCN) {
  CI500025list[i,] = as.numeric(unlist(results[2*i-1]))
  CI500975list[i,] = as.numeric(unlist(results[2*i]))
}

c(mean((CI500025list[,1]<= 3)& (CI500975list[,1]>=3)), 
  mean((CI500025list[,2]<= 3)& (CI500975list[,2]>=3)),
  mean((CI500025list[,3]<= 3)& (CI500975list[,3]>=3)))

#########################################
## sample size 2000
sampleN=2000
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
  # 生成唯一的文件名
  output_file <- paste0("REC_", times, ".txt") # 这里可以区分文件名
  data <- generate_data_OCPM(sampleN)  # 调用相应的函数
  process_data(data,B2, times, verbose = TRUE, file = output_file)
}
stopCluster(cl)

CI2000025list = matrix(nrow = MCN, ncol = 3)
CI2000975list = matrix(nrow = MCN, ncol = 3)
for (i in 1:MCN) {
  CI2000025list[i,] = as.numeric(unlist(results[2*i-1]))
  CI2000975list[i,] = as.numeric(unlist(results[2*i]))
}

c(mean((CI2000025list[,1]<= 3)& (CI2000975list[,1]>=3)), 
  mean((CI2000025list[,2]<= 3)& (CI2000975list[,2]>=3)),
  mean((CI2000025list[,3]<= 3)& (CI2000975list[,3]>=3)))



###############################################################
## OMPC
MCN = 1000
num_cores <- 70
B2 <- 500

#########################################
## sample size 500
sampleN=500
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
  # 生成唯一的文件名
  output_file <- paste0("REC_", times, ".txt") # 这里可以区分文件名
  data <- generate_data_OMPC(sampleN)  # 调用相应的函数
  process_data(data,B2, times, verbose = TRUE, file = output_file)
}
stopCluster(cl)

CI500025list = matrix(nrow = MCN, ncol = 3)
CI500975list = matrix(nrow = MCN, ncol = 3)
for (i in 1:MCN) {
  CI500025list[i,] = as.numeric(unlist(results[2*i-1]))
  CI500975list[i,] = as.numeric(unlist(results[2*i]))
}

c(mean((CI500025list[,1]<= 3)& (CI500975list[,1]>=3)), 
  mean((CI500025list[,2]<= 3)& (CI500975list[,2]>=3)),
  mean((CI500025list[,3]<= 3)& (CI500975list[,3]>=3)))

#########################################
## sample size 2000
sampleN=2000
cl <- makeCluster(num_cores)
registerDoParallel(cl)
results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
  # 生成唯一的文件名
  output_file <- paste0("REC_", times, ".txt") # 这里可以区分文件名
  data <- generate_data_OMPC(sampleN)  # 调用相应的函数
  process_data(data,B2, times, verbose = TRUE, file = output_file)
}
stopCluster(cl)

CI2000025list = matrix(nrow = MCN, ncol = 3)
CI2000975list = matrix(nrow = MCN, ncol = 3)
for (i in 1:MCN) {
  CI2000025list[i,] = as.numeric(unlist(results[2*i-1]))
  CI2000975list[i,] = as.numeric(unlist(results[2*i]))
}

c(mean((CI2000025list[,1]<= 3)& (CI2000975list[,1]>=3)), 
  mean((CI2000025list[,2]<= 3)& (CI2000975list[,2]>=3)),
  mean((CI2000025list[,3]<= 3)& (CI2000975list[,3]>=3)))