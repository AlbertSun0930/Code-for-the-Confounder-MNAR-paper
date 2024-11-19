rm(list = ls())
library(BB)
library(ggplot2)
library(mice)
library(parallel)
library(doParallel)

### True values for parameters
Talpha = c(0.5,-1,2)
Tbeta = c(0.5,1.5,-0.5)
Tgamma = c(0.5,0.5)

#################################################################################
# Asymptotic variance

VarN=100000
set.seed(53)
connum = 2
internum = 1

C1 = rnorm(VarN,-0.5,1)
P_A = exp(Tgamma[1]+Tgamma[2]*C1)/(1+exp(Tgamma[1]+Tgamma[2]*C1))
A = ifelse(P_A>runif(VarN,0,1),1,0)
P_Y = exp(Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1)/(1+exp(Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1))
Y = ifelse(P_Y>runif(VarN,0,1),1,0)
P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*Y)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*Y))
RC = ifelse(P_RC>runif(VarN,0,1),1,0)

## Variance for alpha

Mmatrix = matrix(nrow = VarN, ncol = 9)
for (i in 1:VarN) {
  Ti = -RC[i]*exp(-Talpha[1]-Talpha[2]*C1[i]-Talpha[3]*Y[i])
  Ui = matrix(c(1,A[i],Y[i]),ncol=1)%*%(c(1,C1[i],Y[i]))*Ti
  Mmatrix[i,] = as.vector(Ui)
}

M = matrix(apply(Mmatrix,2,mean),nrow =3)
Phimatrix = matrix(nrow = 3, ncol = VarN)
invM = solve(M)
for (i in 1:VarN) {
  mi = matrix(c(1,A[i],Y[i]),ncol = 1)%*%(RC[i]/P_RC[i]-1)
  Tpi = -invM %*% mi
  Phimatrix[,i] = Tpi
}

Vpmatrix = matrix(nrow = VarN, ncol = 9)
for (i in 1:VarN) {
  Vpi = Phimatrix[,i] %*% t(Phimatrix[,i])
  Vpmatrix[i,] = as.vector(Vpi)
}

Valpha = matrix(apply(Vpmatrix,2,mean),ncol=3)

##################################################

## Variance for gamma

qmatrix = matrix(nrow = VarN,ncol = 2)
Qamatrix = matrix(nrow = VarN,ncol = 6)
Qgmatrix = matrix(nrow = VarN,ncol = 4)
Vgtmatrix = matrix(nrow = VarN,ncol = 4)

for (i in 1:VarN) {
  Qgi = - (RC[i]/P_RC[i])*P_A[i]*(1-P_A[i])*matrix(c(1,C1[i],C1[i],C1[i]^2),ncol=2)
  Qgmatrix[i,] = as.vector(Qgi)
}
Qg = matrix(apply(Qgmatrix,2,mean),ncol=2)
invQg = solve(Qg)

for (i in 1:VarN) {
  Qai = -RC[i]*exp(-Talpha[1]-Talpha[2]*C1[i]-Talpha[3]*Y[i])*(A[i]-P_A[i])*matrix(c(1,C1[i]),ncol=1)%*%matrix(c(1,C1[i],Y[i]),ncol=3)
  Qamatrix[i,] = as.vector(Qai)
}
Qa = matrix(apply(Qamatrix,2,mean),ncol=3)

for (i in 1:VarN) {
  qi = matrix((RC[i]/P_RC[i])*(A[i]-P_A[i])*c(1,C1[i]),ncol = 1)
  qmatrix[i,] = as.vector(qi)
  Vgti = qi+Qa %*% matrix(Phimatrix[,i],ncol=1)
  Vgtmatrix[i,] = as.vector(Vgti%*%t(Vgti))
}
Vgt = matrix(apply(Vgtmatrix,2,mean),ncol=2)
Vg = solve(Qg)%*%Vgt%*%t(solve(Qg))

####################################################
## variance for beta

smatrix = matrix(nrow = VarN,ncol = 3)
Sbmatrix = matrix(nrow = VarN,ncol = 9)
Samatrix = matrix(nrow = VarN,ncol = 9)
Vbtmatrix = matrix(nrow = VarN,ncol = 9)

for (i in 1:VarN) {
  Sbi = - (RC[i]/P_RC[i])*P_Y[i]*(1-P_Y[i])*matrix(c(1,A[i],C1[i]),ncol=1)%*%matrix(c(1,A[i],C1[i]),ncol=3)
  Sbmatrix[i,] = as.vector(Sbi)
}
Sb = matrix(apply(Sbmatrix,2,mean),ncol=3)
invSb = solve(Sb)

for (i in 1:VarN) {
  Sai = -RC[i]*exp(-Talpha[1]-Talpha[2]*C1[i]-Talpha[3]*Y[i])*(Y[i]-P_Y[i])*matrix(c(1,A[i],C1[i]),ncol=1)%*%matrix(c(1,C1[i],Y[i]),ncol=3)
  Samatrix[i,] = as.vector(Sai)
}
Sa = matrix(apply(Samatrix,2,mean),ncol=3)

for (i in 1:VarN) {
  si = matrix((RC[i]/P_RC[i])*(Y[i]-P_Y[i])*c(1,A[i],C1[i]),ncol = 1)
  smatrix[i,] = as.vector(si)
  Vbti = si+Sa %*% matrix(Phimatrix[,i],ncol=1)
  Vbtmatrix[i,] = as.vector(Vbti%*%t(Vbti))
}
Vbt = matrix(apply(Vbtmatrix,2,mean),ncol=3)
Vb = solve(Sb)%*%Vbt%*%t(solve(Sb))

sqrt(c(diag(Vg),diag(Vb))/500)
sqrt(c(diag(Vg),diag(Vb))/2000)


# CCasyOM <- glm(Y[RC==1]~cbind(A,C1)[RC==1,],family = 'binomial')
# CCVb <- (summary(CCasyOM)$coefficients[,2])^2*VarN
# CCasyPWM <- glm(A[RC==1]~cbind(C1)[RC==1,],family = 'binomial')
# CCVg <- (summary(CCasyPWM)$coefficients[,2])^2*VarN
# 
# O_C1 = ifelse(RC==1,C1,NA)
# imp = mice(cbind(O_C1,A,Y),m=3,method='pmm',print=FALSE)
# m=3
# mibetares <- with(imp,glm(Y~A+O_C1,family = 'binomial'))
# MIbetaVar1 <- apply(matrix(summary(mibetares)$std.error,ncol=m)^2,1,mean)*VarN   # vary with sample size,
# estbetares <- matrix(as.numeric(summary(mibetares)$estimate),ncol=m )
# MIbetaVar2 <- apply((estbetares - apply(estbetares,1,mean))^2,1,sum)*(m+1)/(m*(m-1))  # constant when sample sizes change, vary with imputation times.
# 
# migammares <- with(imp,glm(A~O_C1,family = 'binomial'))
# MIgammaVar1 <- apply(matrix(summary(migammares)$std.error,ncol=m)^2,1,mean)*VarN   # vary with sample size,
# estgammares <- matrix(as.numeric(summary(migammares)$estimate),ncol=m )
# MIgammaVar2 <- apply((estgammares - apply(estgammares,1,mean))^2,1,sum)*(m+1)/(m*(m-1))  # constant when sample sizes change, vary with imputation times.


#################################################################################

generate_data <- function(sampleN) {
  N = sampleN ## Sample size 500. We can set N=500 when evaluating sample size 500.
  C1 = rnorm(N,0,1)-0.5
  P_A = exp(Tgamma[1]+Tgamma[2]*C1)/(1+exp(Tgamma[1]+Tgamma[2]*C1))
  A = ifelse(P_A>runif(N,0,1),1,0)
  P_Y = exp(Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1)/(1+exp(Tbeta[1]+Tbeta[2]*A+Tbeta[3]*C1))
  Y = ifelse(P_Y>runif(N,0,1),1,0)
  P_RC = exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*Y)/(1+exp(Talpha[1]+Talpha[2]*C1+Talpha[3]*Y))
  RC = ifelse(P_RC>runif(N,0,1),1,0)
  O_C1 = ifelse(RC==1,C1,NA)
  list(C1 = C1, A = A, Y = Y, RC = RC, O_C1 = O_C1)
}

process_data <- function(data, times, verbose = TRUE, file = "output.txt") {
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
  CCbetavar <- (summary(cc)$coefficients[,2])^2
  
  ccPWM <- glm(A[RC==1]~cbind(C1)[RC==1,], family = 'binomial')
  CCgamma <- ccPWM$coefficients
  CCgammavar <- (summary(ccPWM)$coefficients[,2])^2
  
  #multiple imputation
  imp = mice(cbind(O_C1,A,Y),m=3,method='pmm', donors = 50 ,print=FALSE)
  mi = pool(with(imp,glm(Y~A+O_C1,family = 'binomial')))
  MIbeta = mi$pooled[,3]
  MIbetavar = mi$pooled[,6]
  
  mi2 = pool(with(imp,glm(A~O_C1,family = 'binomial')))
  MIgamma = mi2$pooled[,3]
  MIgammavar = mi2$pooled[,6]
  
  cat("Completed Monte Carlo", times, "\n", file = file, append = TRUE)
  list( as.numeric(c(WEEalpha,WEEgamma,WEEbeta,CCgamma,CCbeta,MIgamma,MIbeta,CCgammavar,CCbetavar,MIgammavar,MIbetavar )))
}

##############################

sampleN=500
MCN = 1000
num_cores <- 70
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 87)
results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS','BB','mice')) %dopar% {
  # Generate a unique file name for each parallel worker
  output_file <- paste0("resultsrecord.txt")
  data <- generate_data(sampleN)
  process_data(data,times, verbose = TRUE, file = output_file)
}
stopCluster(cl)
results_df <- do.call(rbind, lapply(results, function(x) data.frame(t(x))))

WEEbias500 <- apply(results_df[,4:8],2,mean) - c(Tgamma,Tbeta)
WEEhstd500 <- c(sqrt(diag(Vg/sampleN)),sqrt(diag(Vb/sampleN)))
WEEstd500 <- apply(results_df[,4:8],2,sd)
WEEcicr500 <- apply(((t(results_df[,4:8]) + 1.96*WEEhstd500)>c(Tgamma,Tbeta))&((t(results_df[,4:8]) - 1.96*WEEhstd500)<c(Tgamma,Tbeta)),1, mean)
WEEbias_mcse500 <- WEEstd500/sqrt(MCN)
WEEstd_mcse500 <- WEEstd500/sqrt(2*(MCN-1))
WEEcicr_mcse500 <- sqrt(WEEcicr500*(1-WEEcicr500)/MCN)

CCbias500 <- apply(results_df[,9:13],2,mean) - c(Tgamma,Tbeta)
CChstd500 <- sqrt(apply(results_df[,19:23],2,mean))
CCstd500 <- apply(results_df[,9:13],2,sd)
CCcicr500 <- apply(((t(results_df[,9:13]) + 1.96*CChstd500)>c(Tgamma,Tbeta))&((t(results_df[,9:13]) - 1.96*CChstd500)<c(Tgamma,Tbeta)),1, mean)
CCbias_mcse500 <- CCstd500/sqrt(MCN)
CChstd_mcse500 <- sqrt(apply(results_df[,19:23],2,var)/(4*MCN*CChstd500^2))
CCstd_mcse500 <- CCstd500/sqrt(2*(MCN-1))
CCcicr_mcse500 <- sqrt(CCcicr500*(1-CCcicr500)/MCN)

MIbias500 <- apply(results_df[,14:18],2,mean) - c(Tgamma,Tbeta)
MIhstd500 <- sqrt(apply(results_df[,24:28],2,mean))
MIstd500 <- apply(results_df[,14:18],2,sd)
MIcicr500 <- apply(((t(results_df[,14:18]) + 1.96*MIhstd500)>c(Tgamma,Tbeta))&((t(results_df[,14:18]) - 1.96*MIhstd500)<c(Tgamma,Tbeta)),1, mean)
MIbias_mcse500 <- MIstd500/sqrt(MCN)
MIhstd_mcse500 <- sqrt(apply(results_df[,19:23],2,var)/(4*MCN*MIhstd500^2))
MIstd_mcse500 <- MIstd500/sqrt(2*(MCN-1))
MIcicr_mcse500 <- sqrt(MIcicr500*(1-MIcicr500)/MCN)

##############################

sampleN=2000
MCN = 1000
num_cores <- 70
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1)
results <- foreach(times = 1:MCN, .combine = 'c', .packages = c('MASS','BB','mice')) %dopar% {
  # Generate a unique file name for each parallel worker
  output_file <- paste0("resultsrecord.txt")
  data <- generate_data(sampleN)
  process_data(data,times, verbose = TRUE, file = output_file)
}
stopCluster(cl)
results_df <- do.call(rbind, lapply(results, function(x) data.frame(t(x))))

WEEbias2000 <- apply(results_df[,4:8],2,mean) - c(Tgamma,Tbeta)
WEEhstd2000 <- c(sqrt(diag(Vg/sampleN)),sqrt(diag(Vb/sampleN)))
WEEstd2000 <- apply(results_df[,4:8],2,sd)
WEEcicr2000 <- apply(((t(results_df[,4:8]) + 1.96*WEEhstd2000)>c(Tgamma,Tbeta))&((t(results_df[,4:8]) - 1.96*WEEhstd2000)<c(Tgamma,Tbeta)),1, mean)
WEEbias_mcse2000 <- WEEstd2000/sqrt(MCN)
WEEstd_mcse2000 <- WEEstd2000/sqrt(2*(MCN-1))
WEEcicr_mcse2000 <- sqrt(WEEcicr2000*(1-WEEcicr2000)/MCN)

CCbias2000 <- apply(results_df[,9:13],2,mean) - c(Tgamma,Tbeta)
CChstd2000 <- sqrt(apply(results_df[,19:23],2,mean))
CCstd2000 <- apply(results_df[,9:13],2,sd)
CCcicr2000 <- apply(((t(results_df[,9:13]) + 1.96*CChstd2000)>c(Tgamma,Tbeta))&((t(results_df[,9:13]) - 1.96*CChstd2000)<c(Tgamma,Tbeta)),1, mean)
CCbias_mcse2000 <- CCstd2000/sqrt(MCN)
CChstd_mcse2000 <- sqrt(apply(results_df[,19:23],2,var)/(4*MCN*CChstd2000^2))
CCstd_mcse2000 <- CCstd2000/sqrt(2*(MCN-1))
CCcicr_mcse2000 <- sqrt(CCcicr2000*(1-CCcicr2000)/MCN)

MIbias2000 <- apply(results_df[,14:18],2,mean) - c(Tgamma,Tbeta)
MIhstd2000 <- sqrt(apply(results_df[,24:28],2,mean))
MIstd2000 <- apply(results_df[,14:18],2,sd)
MIcicr2000 <- apply(((t(results_df[,14:18]) + 1.96*MIhstd2000)>c(Tgamma,Tbeta))&((t(results_df[,14:18]) - 1.96*MIhstd2000)<c(Tgamma,Tbeta)),1, mean)
MIbias_mcse2000 <- MIstd2000/sqrt(MCN)
MIhstd_mcse2000 <- sqrt(apply(results_df[,19:23],2,var)/(4*MCN*MIhstd2000^2))
MIstd_mcse2000 <- MIstd2000/sqrt(2*(MCN-1))
MIcicr_mcse2000 <- sqrt(MIcicr2000*(1-MIcicr2000)/MCN)




t5001 <- rbind(WEEbias500,WEEhstd500,WEEstd500,WEEcicr500,CCbias500,CChstd500,CCstd500,CCcicr500,MIbias500,MIhstd500,MIstd500,MIcicr500)
t5002 <- rbind(WEEbias_mcse500, rep(NA,5), WEEstd_mcse500,WEEcicr_mcse500,
               CCbias_mcse500,CChstd_mcse500,CCstd_mcse500,CCcicr_mcse500,
               MIbias_mcse500,MIhstd_mcse500,MIstd_mcse500,MIcicr_mcse500)
t500 <- data.frame(t5001[, 1], t5002[, 1], t5001[, 2], t5002[, 2],t5001[, 3], t5002[, 3],
                   t5001[, 4], t5002[, 4],t5001[, 5], t5002[, 5])
as.matrix( round(t500,3))


t20001 <- rbind(WEEbias2000,WEEhstd2000,WEEstd2000,WEEcicr2000,CCbias2000,CChstd2000,CCstd2000,CCcicr2000,MIbias2000,MIhstd2000,MIstd2000,MIcicr2000)
t20002 <- rbind(WEEbias_mcse2000, rep(NA,5),WEEstd_mcse2000,WEEcicr_mcse2000,
                CCbias_mcse2000,CChstd_mcse2000,CCstd_mcse2000,CCcicr_mcse2000,
                MIbias_mcse2000,MIhstd_mcse2000,MIstd_mcse2000,MIcicr_mcse2000)
t2000 <- data.frame(t20001[, 1], t20002[, 1], t20001[, 2], t20002[, 2],t20001[, 3], t20002[, 3],
                    t20001[, 4], t20002[, 4],t20001[, 5], t20002[, 5])
as.matrix( round(t2000,3))
