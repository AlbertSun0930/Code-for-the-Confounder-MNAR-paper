rm(list = ls())
library(BB)
library(ggplot2)
library(mice)
library(parallel)
library(doParallel)

generate_data_OCPC_sen <- function(sampleN,alpha_A) {
  Tgamma = c(-0.5,1,1)
  Tbeta = c(1,3,1,-1)
  Talpha = c(1,-2,1,3)
  N = sampleN ## Sample size 500. We can set N=500 when evaluating sample size 500.
  C1 = rnorm(N,0,1)
  C2 = rbinom(N,1,0.5)
  P_A = exp(Tgamma[1]+Tgamma[2]*C1+Tgamma[3]*C2)/(1+exp(Tgamma[1]+Tgamma[2]*C1+Tgamma[3]*C2))
  A = ifelse(P_A>runif(N,0,1),1,0)
  Y = rnorm(N, Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1+Tbeta[4]*C2, 1)
  P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y+alpha_A*A)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*C2+Talpha[4]*Y+alpha_A*A))
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

MCN = 1000
num_cores <- 70
sampleN=500
alpha_A_list = c(-0.4,-0.2,0.2,0.4)
res <- c()
for (alpha_A in alpha_A_list) {
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 1)
  results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS', 'BB', 'mice')) %dopar% {
    # 生成唯一的文件名
    output_file <- paste0("resultsrecord.txt")  # 这里可以区分文件名
    data <- generate_data_OCPC_sen(sampleN,alpha_A)  # 调用相应的函数
    process_data(data, times, verbose = TRUE, file = output_file)
  }
  stopCluster(cl)
  results_df <- do.call(rbind, lapply(results, function(x) data.frame(t(x))))
  res <- rbind(res,results_df)
  print(alpha_A)
}

resbiasstdest <- round(cbind(
apply(res[1:1000,],2,mean)-3,
apply(res[1:1000,],2,sd),
apply(res[1001:2000,],2,mean)-3,
apply(res[1001:2000,],2,sd),
apply(res[2001:3000,],2,mean)-3,
apply(res[2001:3000,],2,sd),
apply(res[3001:4000,],2,mean)-3,
apply(res[3001:4000,],2,sd)),3)

resbiasstdmcse <- round(cbind(
  apply(res[1:1000,],2,sd)/sqrt(1000),
  apply(res[1:1000,],2,sd)/sqrt(2*999),
  apply(res[1001:2000,],2,sd)/sqrt(1000),
  apply(res[1001:2000,],2,sd)/sqrt(2*999),
  apply(res[2001:3000,],2,sd)/sqrt(1000),
  apply(res[2001:3000,],2,sd)/sqrt(2*999),
  apply(res[3001:4000,],2,sd)/sqrt(1000),
  apply(res[3001:4000,],2,sd)/sqrt(2*999)),3)



finalres <- cbind(as.numeric(unlist(t(resbiasstdest[,1:2]))),
                  as.numeric(unlist(t(resbiasstdmcse[,1:2]))),
                  as.numeric(unlist(t(resbiasstdest[,3:4]))),
                  as.numeric(unlist(t(resbiasstdmcse[,3:4]))),
                  as.numeric(unlist(t(resbiasstdest[,5:6]))),
                  as.numeric(unlist(t(resbiasstdmcse[,5:6]))),
                  as.numeric(unlist(t(resbiasstdest[,7:8]))),
                  as.numeric(unlist(t(resbiasstdmcse[,7:8]))))

finalres

