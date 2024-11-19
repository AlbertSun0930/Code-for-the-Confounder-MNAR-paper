library(BB)
library(ggplot2)
library(mice)
library(foreign)

mergeddata <- read.csv('mergeddata.csv')

OMPtable <- matrix(nrow =datasets, ncol = 1)
IPWtable <- matrix(nrow =datasets, ncol = 1)
DRtable <- matrix(nrow =datasets, ncol = 1)
CCtable <- matrix(nrow =datasets, ncol = 1)
MItable <- matrix(nrow =datasets, ncol = 1)
alphatable <- matrix(nrow =datasets, ncol = 2)

ccipw = matrix(nrow = datasets, ncol = 1)
ccaipw = matrix(nrow = datasets, ncol = 1)
miipw = matrix(nrow = datasets, ncol = 1)
miaipw = matrix(nrow = datasets, ncol = 1)


# Bootstrap standard errors
set.seed(3)
for (i in 1:datasets) {
  samplenum <- sample(dim(mergeddata2)[1],replace = TRUE)
  A <- mergeddata2$DMDMARTL[samplenum]
  C1 <- mergeddata2$RIDAGEYR[samplenum]
  C2 <- mergeddata2$RIAGENDR[samplenum]
  C3 <- mergeddata2$INDFMMPI[samplenum]
  C4 <- mergeddata2$RIAGENDR[samplenum]
  Y <- mergeddata2$GM[samplenum]
  RC <- ifelse(is.na(C3),0,1)
  C3[is.na(C3)]=6
  
  fun <- function(alpha) { 
    f <- numeric(length(alpha)) 					
    f[1] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)) - sum(ifelse(RC==1,0,1))
    f[2] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*A) - sum(ifelse(RC==1,0,1)*A)
    f[3] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*Y) - sum(ifelse(RC==1,0,1)*Y)		
    f[4] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*Y^2) - sum(ifelse(RC==1,0,1)*Y^2)
    f[5] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*A*Y) - sum(ifelse(RC==1,0,1)*A*Y)
    f[6] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*A*Y^2) - sum(ifelse(RC==1,0,1)*A*Y^2)	
    f 
  } 
  startalpha <- c(0,0,0,-2,0,0)
  result = dfsane(startalpha,fun,quiet = TRUE)
  alpha = result$par
  alphatable[i,1] <- alpha[4]
  alphatable[i,2] <- alpha[6]
  
  mw = ifelse(RC==1,1,0)*(1+exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y))
  OM <- lm(Y~cbind(A,C1,C2,C3,C4),weights = mw)
  beta <- OM$coefficients
  
  PWM <- glm(A~cbind(C1,C2,C3,C4),weights = mw,family = 'quasibinomial')
  gamma <- PWM$coefficients
  
  OMPtable[i] <- beta[2]
  
  ps <- predict(PWM,type = 'response')
  pw <- ifelse(A==1, 1/ps, 1/(1-ps))
  
  IPWM <- lm(Y~A, weights = mw*pw)
  IPWtable[i] = IPWM$coefficients[2]
  
  DRM = lm(Y~cbind(A,C1,C2,C3,C4), weights =  mw*pw)
  DRtable[i] <- DRM$coefficients[2]
  
  CCtable[i] <- lm(Y[RC==1]~cbind(A,C1,C2,C3,C4)[RC==1,])$coefficients[2]
  
  ccPWM <- glm(A[RC==1]~cbind(C1,C2,C3,C4)[RC==1,], family = 'binomial')
  ccps <- predict(ccPWM,type = 'response')
  ccpw <- ifelse(A[RC==1]==1, 1/ccps, 1/(1-ccps))
  
  ccipw[i] <- lm(Y[RC==1]~A[RC==1], weights = ccpw)$coefficients[2]
  ccaipw[i] <- lm(Y[RC==1]~cbind(A,C1,C2,C3,C4)[RC==1,], weights = ccpw)$coefficients[2]
  
  O_C <- C3
  O_C[O_C==6] = NA
  imp = mice(cbind(O_C,A,C1,C2,C4,Y),m=3,method='pmm',print=FALSE)
  mi = pool(with(imp,lm(Y~A+O_C+C1+C2+C4)))
  MItable[i] = mi$pooled[2,3]
  
  com <- complete(imp,action = 'long')
  comd1 <- com[com$.imp==1,]
  comd2 <- com[com$.imp==2,]
  comd3 <- com[com$.imp==3,]
  
  miPWM1 <- glm(A~cbind(C1,C2,O_C,C4), data  = comd1, family = 'binomial')
  mips1 <- predict(miPWM1,type = 'response')
  mipw1 <- ifelse(A==1, 1/mips1, 1/(1-mips1))
  miipw1 <- lm(Y~A, data  = comd1, weights = mipw1)$coefficients[2]
  midr1 <- lm(Y~cbind(A,C1,C2,O_C,C4), data  = comd1, weights = mipw1)$coefficients[2]
  
  miPWM2 <- glm(A~cbind(C1,C2,O_C,C4), data  = comd2, family = 'binomial')
  mips2 <- predict(miPWM2,type = 'response')
  mipw2 <- ifelse(A==1, 1/mips2, 1/(1-mips2))
  miipw2 <- lm(Y~A, data  = comd2, weights = mipw2)$coefficients[2]
  midr2 <- lm(Y~cbind(A,C1,C2,O_C,C4), data  = comd2, weights = mipw2)$coefficients[2]
  
  miPWM3 <- glm(A~cbind(C1,C2,O_C,C4), data  = comd3, family = 'binomial')
  mips3 <- predict(miPWM3,type = 'response')
  mipw3 <- ifelse(A==1, 1/mips3, 1/(1-mips3))
  miipw3 <- lm(Y~A, data  = comd3, weights = mipw3)$coefficients[2]
  midr3 <- lm(Y~cbind(A,C1,C2,O_C,C4), data  = comd3, weights = mipw3)$coefficients[2]  
  
  miipw[i] <- mean(c(miipw1,miipw2,miipw3))
  miaipw[i] <- mean(c(midr1,midr2,midr3))
}


resu <- data.frame(cbind(alphatable,DRtable, IPWtable, OMPtable,ccipw,ccaipw,CCtable,miipw,miaipw,MItable))
colnames(resu) = c('alpha1','alpha2','WEE-DR','WEE-IPW','WEE-OMP','CC-AIPW','CC-IPW','CC-OMP','MI-AIPW','MI-IPW','MI-OMP')
BSstd <- apply(resu, 2, sd)[-2]



# Point estimate from original data
set.seed(5)
peod  = rep(0,10)
A <- mergeddata2$DMDMARTL
C1 <- mergeddata2$RIDAGEYR
C2 <- mergeddata2$DMDEDUC2
C3 <- mergeddata2$INDFMMPI
C4 <- mergeddata2$RIAGENDR
Y <- mergeddata2$GM
RC <- ifelse(is.na(C3),0,1)
C3[is.na(C3)]=6

fun <- function(alpha) { 
  f <- numeric(length(alpha)) 					
  f[1] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)) - sum(ifelse(RC==1,0,1))
  f[2] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*A) - sum(ifelse(RC==1,0,1)*A)
  f[3] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*Y) - sum(ifelse(RC==1,0,1)*Y)		
  f[4] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*Y^2) - sum(ifelse(RC==1,0,1)*Y^2)
  f[5] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*A*Y) - sum(ifelse(RC==1,0,1)*A*Y)
  f[6] <-  sum(ifelse(RC==1,1,0)*exp(-alpha[1]-alpha[2]*C1-alpha[3]*C2-alpha[4]*C3-alpha[5]*C4 - alpha[6]*Y)*A*Y^2) - sum(ifelse(RC==1,0,1)*A*Y^2)	
  f 
} 
startalpha <- c(0,0,0,-2,0,0)
result = dfsane(startalpha,fun,quiet = TRUE)
alpha2 = result$par
alphatable[i,1] <- alpha2[4]
alphatable[i,2] <- alpha2[6]

peod[1] <- alpha2[4]

mw = ifelse(RC==1,1,0)*(1+exp(-alpha2[1]-alpha2[2]*C1-alpha2[3]*C2-alpha2[4]*C3-alpha2[5]*C4 - alpha2[6]*Y))
OM <- lm(Y~cbind(A,C1,C2,C3,C4),weights = mw)
beta <- OM$coefficients
peod[4] <- beta[2]

PWM <- glm(A~cbind(C1,C2,C3,C4),weights = mw,family = 'quasibinomial')
ps <- predict(PWM,type = 'response')
pw <- ifelse(A==1, 1/ps, 1/(1-ps))

IPWM <- lm(Y~A, weights = mw*pw)
peod[3] = IPWM$coefficients[2]

DRM = lm(Y~cbind(A,C1,C2,C3,C4), weights =  mw*pw)
peod[2] <- DRM$coefficients[2]

peod[7] <- lm(Y[RC==1]~cbind(A,C1,C2,C3,C4)[RC==1,])$coefficients[2]

ccPWM <- glm(A[RC==1]~cbind(C1,C2,C3,C4)[RC==1,], family = 'binomial')
ccps <- predict(ccPWM,type = 'response')
ccpw <- ifelse(A[RC==1]==1, 1/ccps, 1/(1-ccps))

peod[6] <- lm(Y[RC==1]~A[RC==1], weights = ccpw)$coefficients[2]
peod[5] <- lm(Y[RC==1]~cbind(A,C1,C2,C3,C4)[RC==1,], weights = ccpw)$coefficients[2]

O_C <- C3
O_C[O_C==6] = NA
imp = mice(cbind(O_C,A,C1,C2,C4,Y),m=3,method='pmm',print=FALSE)
mi = pool(with(imp,lm(Y~A+O_C+C1+C2+C4)))
peod[10] = mi$pooled[2,3]

com <- complete(imp,action = 'long')
comd1 <- com[com$.imp==1,]
comd2 <- com[com$.imp==2,]
comd3 <- com[com$.imp==3,]

miPWM1 <- glm(A~cbind(C1,C2,O_C,C4), data  = comd1, family = 'binomial')
mips1 <- predict(miPWM1,type = 'response')
mipw1 <- ifelse(A==1, 1/mips1, 1/(1-mips1))
miipw1 <- lm(Y~A, data  = comd1, weights = mipw1)$coefficients[2]
midr1 <- lm(Y~cbind(A,C1,C2,O_C,C4), data  = comd1, weights = mipw1)$coefficients[2]

miPWM2 <- glm(A~cbind(C1,C2,O_C,C4), data  = comd2, family = 'binomial')
mips2 <- predict(miPWM2,type = 'response')
mipw2 <- ifelse(A==1, 1/mips2, 1/(1-mips2))
miipw2 <- lm(Y~A, data  = comd2, weights = mipw2)$coefficients[2]
midr2 <- lm(Y~cbind(A,C1,C2,O_C,C4), data  = comd2, weights = mipw2)$coefficients[2]

miPWM3 <- glm(A~cbind(C1,C2,O_C,C4), data  = comd3, family = 'binomial')
mips3 <- predict(miPWM3,type = 'response')
mipw3 <- ifelse(A==1, 1/mips3, 1/(1-mips3))
miipw3 <- lm(Y~A, data  = comd3, weights = mipw3)$coefficients[2]
midr3 <- lm(Y~cbind(A,C1,C2,O_C,C4), data  = comd3, weights = mipw3)$coefficients[2]  


peod[9] <- mean(c(miipw1,miipw2,miipw3))
peod[8] <- mean(c(midr1,midr2,midr3))


# Wald-type confidence intervals

CI025 <- peod - 1.96*BSstd
CI975 <- peod + 1.96*BSstd

# Final results

round(cbind(peod, BSstd, CI025, CI975),3)

