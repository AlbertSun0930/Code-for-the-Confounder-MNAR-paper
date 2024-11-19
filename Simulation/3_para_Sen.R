rm(list = ls())
library(BB)
library(ggplot2)
library(mice)
library(parallel)
library(doParallel)

################################################################################

### Binary outcome

Talpha = c(0.5,-1,2)
Tbeta = c(0.5,1.5,-0.5)
Tgamma = c(0.5,0.5)

generate_data_sen_BY <- function(sampleN,alpha_A) {
  N = sampleN ## Sample size 500. We can set N=500 when evaluating sample size 500.
  C1 = rnorm(N,0,1)-0.5
  P_A = exp(Tgamma[1]+Tgamma[2]*C1)/(1+exp(Tgamma[1]+Tgamma[2]*C1))
  A = ifelse(P_A>runif(N,0,1),1,0)
  P_Y = exp(Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1+alpha_A*A)/(1+exp(Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1+alpha_A*A))
  Y = ifelse(P_Y>runif(N,0,1),1,0)
  P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*Y)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*Y))
  RC = ifelse(P_RC>runif(N,0,1),1,0)
  O_C1 = ifelse(RC==1,C1,NA)
  list(C1 = C1, A = A, Y = Y, RC = RC, O_C1 = O_C1)
}

process_data_BY <- function(data, times, verbose = TRUE, file = "output.txt") {
  C1 <- data$C1
  A <- data$A
  Y <- data$Y
  RC <- data$RC
  O_C1 <- data$O_C1
  
  fun <- function(alpha) { 
    f <- numeric(length(alpha))                   
    f[1] <- sum(exp(-alpha[1]-alpha[2]*C1 - alpha[3]*Y)*ifelse(RC==1,1,0)) - sum(ifelse(RC==1,0,1))
    f[2] <- sum(exp(-alpha[1]-alpha[2]*C1 - alpha[3]*Y)*ifelse(RC==1,1,0)*A) - sum(ifelse(RC==1,0,1)*A)
    f[3] <- sum(exp(-alpha[1]-alpha[2]*C1 - alpha[3]*Y)*ifelse(RC==1,1,0)*Y) - sum(ifelse(RC==1,0,1)*Y)
    f 
  } 
  startalpha <- Talpha 
  result = dfsane(startalpha,fun,quiet = TRUE)
  WEEalpha = result$par
  
  mw = ifelse(RC==1,1,0)*(1+exp(-WEEalpha[1]-WEEalpha[2]*C1 - WEEalpha[3]*Y))
  PWM <- glm(A~cbind(C1), weights = mw, family = 'binomial')
  WEEgamma <- PWM$coefficients
  
  OM <- glm(Y~cbind(A,C1), weights = mw, family = 'binomial')
  WEEbeta <- OM$coefficients
  
  #complete case analysis
  cc <- glm(Y[RC==1]~cbind(A,C1)[RC==1,],family = 'binomial')
  CCbeta <- cc$coefficients
  
  ccPWM <- glm(A[RC==1]~cbind(C1)[RC==1,], family = 'binomial')
  CCgamma <- ccPWM$coefficients
  
  #multiple imputation
  imp = mice(cbind(O_C1,A,Y),m=3,method='midastouch',print=FALSE)
  mi = pool(with(imp,glm(Y~A+O_C1,family = 'binomial')))
  MIbeta = mi$pooled[,3]
  
  mi2 = pool(with(imp,glm(A~O_C1,family = 'binomial')))
  MIgamma = mi2$pooled[,3]
  
  cat("Completed Monte Carlo", times, "\n", file = file, append = TRUE)
  list( as.numeric(c(WEEalpha,WEEgamma,WEEbeta,CCgamma,CCbeta,MIgamma,MIbeta )))
}


alpha_A_list = c(-0.2,-0.1,0.1,0.2)
sampleN=500
MCN = 1000
num_cores <- 70

BYres500 <- c()
BYresmcse500 <- c()
for (i in 1:4) {
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 63)
  results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS','BB','mice')) %dopar% {
    # Generate a unique file name for each parallel worker
    output_file <- paste0("resultsrecord.txt")
    data <- generate_data_sen_BY(sampleN,alpha_A_list[i])
    process_data_BY(data,times, verbose = TRUE, file = output_file)
  }
  stopCluster(cl)
  results_df <- do.call(rbind, lapply(results, function(x) data.frame(t(x))))
  
  WEEbias500 <- apply(results_df[,4:8],2,mean) - c(Tgamma,Tbeta)
  WEEstd500 <- apply(results_df[,4:8],2,sd)
  WEEbias_mcse500 <- WEEstd500/sqrt(MCN)
  WEEstd_mcse500 <- WEEstd500/sqrt(2*(MCN-1))
  
  CCbias500 <- apply(results_df[,9:13],2,mean) - c(Tgamma,Tbeta)
  CCstd500 <- apply(results_df[,9:13],2,sd)
  CCbias_mcse500 <- CCstd500/sqrt(MCN)
  CCstd_mcse500 <- CCstd500/sqrt(2*(MCN-1))
  
  MIbias500 <- apply(results_df[,14:18],2,mean) - c(Tgamma,Tbeta)
  MIstd500 <- apply(results_df[,14:18],2,sd)
  MIbias_mcse500 <- MIstd500/sqrt(MCN)
  MIstd_mcse500 <- MIstd500/sqrt(2*(MCN-1))
  
  BYres500 <- rbind(BYres500, rbind(WEEbias500,
                                    WEEstd500,
                                    CCbias500,
                                    CCstd500,
                                    MIbias500,
                                    MIstd500))
  BYresmcse500 <- rbind(BYresmcse500, rbind(WEEbias_mcse500,
                                            WEEstd_mcse500,
                                            CCbias_mcse500,
                                            CCstd_mcse500,
                                            MIbias_mcse500,
                                            MIstd_mcse500))
  print(i)
}

round(cbind(BYres500[,1],BYresmcse500[,1],
      BYres500[,2],BYresmcse500[,2],
      BYres500[,3],BYresmcse500[,3],
      BYres500[,4],BYresmcse500[,4],
      BYres500[,5],BYresmcse500[,5]),3)


################################################################################

### Continuous outcome

Talpha = c(-1,1,1)
Tbeta = c(0.5,1.5,-0.5)
Tgamma = c(0.5,0.5)

generate_data_sen_CY <- function(sampleN,alpha_A) {
  N = sampleN ## Sample size 500. We can set N=500 when evaluating sample size 500.
  C1 = rnorm(N)-0.5
  P_A = exp(Tgamma[1]+Tgamma[2]*C1)/(1+exp(Tgamma[1]+Tgamma[2]*C1))
  A = ifelse(P_A>runif(N,0,1),1,0)
  Y = rnorm(N,Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1,1)
  P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*Y)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*Y))
  RC = ifelse(P_RC>runif(N,0,1),1,0)
  O_C1 = ifelse(RC==1,C1,NA)
  list(C1 = C1, A = A, Y = Y, RC = RC, O_C1 = O_C1)
}

process_data_CY <- function(data, times, verbose = TRUE, file = "output.txt") {
  C1 <- data$C1
  A <- data$A
  Y <- data$Y
  RC <- data$RC
  O_C1 <- data$O_C1
  
  fun <- function(alpha) { 
    f <- numeric(length(alpha))                   
    f[1] <- sum(exp(-alpha[1]-alpha[2]*C1 - alpha[3]*Y)*ifelse(RC==1,1,0)) - sum(ifelse(RC==1,0,1))
    f[2] <- sum(exp(-alpha[1]-alpha[2]*C1 - alpha[3]*Y)*ifelse(RC==1,1,0)*A) - sum(ifelse(RC==1,0,1)*A)
    f[3] <- sum(exp(-alpha[1]-alpha[2]*C1 - alpha[3]*Y)*ifelse(RC==1,1,0)*Y) - sum(ifelse(RC==1,0,1)*Y)
    f 
  } 
  startalpha <- Talpha + rnorm(1,0,0.6)
  result = dfsane(startalpha,fun,quiet = TRUE)
  WEEalpha = result$par
  
  mw = ifelse(RC==1,1,0)*(1+exp(-WEEalpha[1]-WEEalpha[2]*C1 - WEEalpha[3]*Y))
  PWM <- glm(A~cbind(C1), weights = mw, family = 'binomial')
  WEEgamma <- PWM$coefficients
  
  OM <- lm(Y~cbind(A,C1), weights = mw)
  WEEbeta <- OM$coefficients
  
  #complete case analysis
  cc <- lm(Y[RC==1]~cbind(A,C1)[RC==1,])
  CCbeta <- cc$coefficients
  
  ccPWM <- glm(A[RC==1]~cbind(C1)[RC==1,], family = 'binomial')
  CCgamma <- ccPWM$coefficients
  
  #multiple imputation
  imp = mice(cbind(O_C1,A,Y),m=3,method='midastouch',print=FALSE)
  mi = pool(with(imp,lm(Y~A+O_C1)))
  MIbeta = mi$pooled[,3]
  
  mi2 = pool(with(imp,glm(A~O_C1,family = 'binomial')))
  MIgamma = mi2$pooled[,3]
  
  cat("Completed Monte Carlo", times, "\n", file = file, append = TRUE)
  list( as.numeric(c(WEEalpha,WEEgamma,WEEbeta,CCgamma,CCbeta,MIgamma,MIbeta)))
}

alpha_A_list = c(-0.2,-0.1,0.1,0.2)
sampleN=500
MCN = 1000
num_cores <- 70

CYres500 <- c()
CYresmcse500 <- c()
for (i in 1:4) {
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 63)
  results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS','BB','mice')) %dopar% {
    # Generate a unique file name for each parallel worker
    output_file <- paste0("resultsrecord.txt")
    data <- generate_data_sen_CY(sampleN,alpha_A_list[i])
    process_data_CY(data,times, verbose = TRUE, file = output_file)
  }
  stopCluster(cl)
  results_df <- do.call(rbind, lapply(results, function(x) data.frame(t(x))))
  
  WEEbias500 <- apply(results_df[,4:8],2,mean) - c(Tgamma,Tbeta)
  WEEstd500 <- apply(results_df[,4:8],2,sd)
  WEEbias_mcse500 <- WEEstd500/sqrt(MCN)
  WEEstd_mcse500 <- WEEstd500/sqrt(2*(MCN-1))
  
  CCbias500 <- apply(results_df[,9:13],2,mean) - c(Tgamma,Tbeta)
  CCstd500 <- apply(results_df[,9:13],2,sd)
  CCbias_mcse500 <- CCstd500/sqrt(MCN)
  CCstd_mcse500 <- CCstd500/sqrt(2*(MCN-1))
  
  MIbias500 <- apply(results_df[,14:18],2,mean) - c(Tgamma,Tbeta)
  MIstd500 <- apply(results_df[,14:18],2,sd)
  MIbias_mcse500 <- MIstd500/sqrt(MCN)
  MIstd_mcse500 <- MIstd500/sqrt(2*(MCN-1))
  
  CYres500 <- rbind(CYres500, rbind(WEEbias500,
                                    WEEstd500,
                                    CCbias500,
                                    CCstd500,
                                    MIbias500,
                                    MIstd500))
  CYresmcse500 <- rbind(CYresmcse500, rbind(WEEbias_mcse500,
                                            WEEstd_mcse500,
                                            CCbias_mcse500,
                                            CCstd_mcse500,
                                            MIbias_mcse500,
                                            MIstd_mcse500))
  print(i)
}

round(cbind(CYres500[,1],CYresmcse500[,1],
            CYres500[,2],CYresmcse500[,2],
            CYres500[,3],CYresmcse500[,3],
            CYres500[,4],CYresmcse500[,4],
            CYres500[,5],CYresmcse500[,5]),3)














