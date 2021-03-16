###Pre-Processing Functions###
library("MASS", lib.loc="C:/Program Files/R/R-3.5.1/library")
library("pracma", lib.loc="~/R/win-library/3.5")

#1. Centering
centering_icaR<-function(X){
  m_mean=matrix(data=1, nrow=nrow(X))%*%colMeans(X)
  output=X-m_mean
  return(output)}

#2. Covariance Matrix Funtion
cov_matrix<-function(X){
  m=ncol(X) #Number of columns of input matrix
  n=nrow(X) #Number of rows of input matrix
  
  #Compute X-E(X) for each element of the matrix
  m_mean=matrix(data=1, nrow=nrow(X))%*%colMeans(X)
  D=X-m_mean
  
  #Computing each element of the matrix
  C<-((n-1)^-1)*t(D)%*%D #Equivalent to computing [X-E(X)]^2
  return(C)
}

#3. Whitening Functions
###3.1 PCA Whitening
whiten_PCA<-function(X){
  #Centering
  X_c=centering_icaR(X)
  
  #Compute Covariance Matrix
  S=cov_matrix(X_c)
  
  #Eigen Value Decomposition
  D=eigen(S)$values
  for (i in 1:length(D)){
    if(D[i]<0){D[i]=0}
  }
  U=eigen(S)$vectors
  #Multiply by the sign of the eigen values
  U=sweep(U, 2, sign(diag(U)), "*") 
  
  #Regularization step
  epsilon = 1e-10 #To prevent NaNs
  
  #PCA Whitening Matrix
  W=(diag((D+epsilon)^(-1/2)))%*%t(U)
  
  #Xpca
  Xtilde=X%*%t(W)
  output=list(Xtilde,W)
  return (output)
  }

##3.2 Whiten ZCA
whiten_ZCA<-function(X){
  #Centering
  X_c=centering_icaR(X)
  
  #Compute Covariance Matrix
  S=cov_matrix(X_c)
  
  #Eigen Value Decomposition
  D=eigen(S)$values
  for (i in 1:length(D)){
    if(D[i]<0){D[i]=0}
  }
  U=eigen(S)$vectors
  #Multiply by the sign of the eigen values
  U=sweep(U, 2, sign(diag(U)), "*") 
  
  #Regularization step
  epsilon = 1e-10 #To prevent NaNs
  
  #ZCA Whitening Matrix
  W=U%*%(diag((D+epsilon)^(-1/2)))%*%t(U)
  
  #Xpzca
  Xtilde=X%*%t(W)
  output=list(Xtilde,W)
  return (output)
}

###Non-quadratic sub-gaussian or super-gaussian functions###
#fsup and its derivatives
fsup<-function(x,a){
  output=((1/a)*log(cosh(a*x)))-((a/2)*x^2)
  return(output)
}

fsup_1<-function(x,a){
  output=tanh(a*x)-(a*x)
  return(output)
}

fsup_2<-function(x,a){
  output=(a*((sech(a*x))^2)-1)
  return(output)
}

#fsub and its derivatives
fsub<-function(x,b){
  output=(b/4)*(x^4)
  return(output)
}

fsub_1<-function(x,b){
  output=b*x^3
  return(output)
}

fsub_2<-function(x,b){
  output=3*b*x^2
  return(output)
}

###Closeness Metrics###
#g(y,ref) = MSE
g.MSE<-function(y,ref){
  n=length(y)
  output=((y-ref)^2)
  return(output)
}
g.MSE_1<-function(y,ref){
  n=length(y)
  output=2*(y-ref)
  return(output)
}
g.MSE_2<-function(y,ref){
  n=length(y)
  output=2*(ones(n,1))
  return(output)
}
#g(y,ref) = cor
g.cor<-function(y,ref){
  n=length(y)
  output=(1/n)*((y*ref)^2)
  return(output)
}
g.cor_1<-function(y,ref){
  n=length(y)
  output=-2*(((1/n)*ref)/((1/n)*(y*ref)^3))
  return(output)
}
g.cor_2<-function(y,ref){
  n=length(y)
  output=6*((1/n)*(ref)^2)/((1/n)*(y*ref)^4)
  return(output)
}

######4. Lu & Rajapakse ICA with reference ICA-R Prototype######
ICA_R<-function(X, ref, threshold, learningRate, mu, lambda, gamma, a, b, ro, maxIter, OverValue, print = TRUE){
###i. Define dimensions of y (ICs)
ICnum=ncol(X)
IClen=nrow(X)
###ii.a Initialize parameters
mu=mu0
lambda0=lambda
gamma=1
a=a
b=b
ro=ro
ref=ref
v_gauss=rnorm(IClen, mean=0, sd=1)
###ii.a Initialize Weights
w0=as.matrix(runif(ICnum))
w=w0/norm(w0, type=c("F"))
oldw=w
###ii.b Standarize Reference Signal
ref.stand=(ref-mean(ref))/sd(ref)
###iii. Center and Whiten data matrix X
X.center=centering_icaR(X)
X.whiten=whiten_ZCA(X.center)[[1]]
#a. Compute Moore-Penrose generalized inverse of X
pInvRxx=round(ginv(cov_matrix(X.whiten)),14)
#Define flow control parameters
flag=1
loop=1

while(flag==1){
###Compute initial y
y=X.whiten%*%w
y.stand=(y-mean(y))/sd(y)

  
###iv. Compute ro.hat
ro.hat=2*ro*(mean(fsup(a,y)-fsup(a,v_gauss)))

###v. Compute d using g'' as the correlation
d=mean(mu*g.MSE_2(y.stand,ref.stand))+(8*lambda)-mean(ro.hat*fsup_2(a,y))

###vi. Update mu/lambda multipliers
mu=max(-mu, gamma*(g.MSE(y.stand,ref.stand)-(threshold*(1-exp(-loop)))))
lambda=lambda+(gamma*(mean(y^2)-1)^2)

###vii. Compute phi and psi 
phi=(1/IClen)*(t((-ro.hat*fsup_1(a,y)))%*%X.whiten)*(1/d)
psi=(1/IClen)*(t(mu*g.MSE_1(y.stand,ref.stand))%*%X.whiten)*(1/d)
omega=(4*lambda)*(mean(y^2)-1)*((1/IClen)*(t(y)%*%X.whiten))*(1/d)

###viii. Update/Normalize Weights
w=w-t((learningRate*(phi+psi+omega))%*%pInvRxx)
w=w/norm(w, type=c("F"))

###ix. Convergence criteria
wchange=1-as.numeric(abs(t(oldw)%*%w))
if(wchange<OverValue){
  flag=0
  print(sprintf("Converged after %d iteration",loop))
}

if(loop>=maxIter){
  flag=0
  print(sprintf("After %d iteration, still cannot converge",loop))

}
print(wchange)
oldw=w
loop=loop+1

}

y=X.whiten%*%w
print("End of cICA algorithm")

output=list(y,w)
return(output) 

}
#########################################################################################################
######5. Zhillin Zhang ICA with Reference ######
Zhang_cICA<-function(X, ref, threshold, learningRate, mu, lambda0, a, b, gamma, maxIter, OverValue, print=TRUE){
###i. Define dimensions of y (ICs)
ICnum=ncol(X)
IClen=nrow(X)
###ii.a Initialize parameters
mu=mu0
lambda=lambda0
gamma=1
rou=1
a=a
b=b
ref=ref
###ii.a Initialize Weights
w0=matrix(runif(ICnum))
w=w0/norm(w0, type=c("F"))
oldw=w
###ii.b Standarize Reference Signal
ref.stand=(ref-mean(ref))/sd(ref)
###iii. Center and Whiten data matrix X
X.center=centering_icaR(X)  
X.whiten=whiten_ZCA(X.center)[[1]]
#a. Compute Moore-Penrose generalized inverse of X
pInvRxx=round(ginv(cov_matrix(X.whiten)),14)
#Define flow control parameters
flag=1
loop=1
  
#While loop
while (flag==1) {
#1.Output using current iteration and standarize y
y=X.whiten%*%w
y.stand=(y-mean(y))/sd(y)
#2. Generate Gaussian Random Matrix and compute rou.hat
v_gauss=rnorm(IClen, mean=0, 1) 
rou.hat=rou*sign((mean(fsup(a,y)-fsup(a,v_gauss))))
    
#3. Gamma1 and Gamma2
Gamma1=((rou.hat*(t(fsup_1(a,y))%*%X.whiten))/IClen) - ((mu/2)*(t(g.MSE_1(y.stand, ref.stand))%*%X.whiten)/(IClen)) - ((lambda*t(y)%*%X.whiten)/(IClen))
Gamma2=((rou.hat*mean(fsup_2(a,y)))) - ((mu/2)*(mean(g.MSE_2(y.stand, ref.stand)))) - lambda

#4. Update Equations for mu and Lambda
mu=max(0, mu + (gamma * (g.MSE(y.stand,ref.stand)-(threshold*(1-exp(-loop))))))
lambda =lambda+(gamma*(mean(y.stand^2)-1)^2)

#5. Update/Normalize Weight Vector
w=w-(t(learningRate*(Gamma1*(1/Gamma2))%*%pInvRxx))
w=w/norm(w, type=c("F"))

#6.Check for Convergence
wchange=1-as.numeric(abs(t(oldw)%*%w))
if(wchange<OverValue){
  flag=0
  print(sprintf("Converged after %d iteration",loop))}

if(loop>=maxIter){
  flag=0
  print(sprintf("After %d iteration, still cannot converge",loop))
  
}
print(wchange)
oldw=w
loop=loop+1

}

y=X.whiten%*%w
print("End of cICA algorithm")

output=list(y,w)
return(output) 
}
#########################################################################################################

######6. Lin fastICA with Reference ######
Lin_fastICAR<-function(X, ref, threshold, learningRate, mu, a, b, ro, gamma, maxIter, OverValue, print=TRUE){
###i. Define dimensions of y (ICs)
ICnum=ncol(X)
IClen=nrow(X)
###ii.a Initialize parameters
mu=mu0
lambda=lambda0
gamma=gamma
rou=ro
a=a
b=b
ref=ref
###ii.a Initialize Weights
pInvX=round(ginv(X),14)
w0=(pInvX%*%ref)
#w0=as.matrix(runif(ICnum))
w=w0/norm(w0, type=c("F"))
oldw=w
###ii.b Standarize Reference Signal
ref.stand=(ref-mean(ref))/sd(ref)
###iii. Center and Whiten data matrix X
X.center=centering_icaR(X)  
X.whiten=whiten_ZCA(X.center)[[1]]
#a. Compute Moore-Penrose generalized inverse of X
pInvRxx=round(ginv(cov_matrix(X.whiten)),14)
#Define flow control parameters
flag=1
loop=1
  
#While loop
while (flag==1) {
    
#1.Output using current iteration and standarize y
y=X.whiten%*%w
y.stand=(y-mean(y))/sd(y)
#2. Generate Gaussian Random Matrix and compute rou.hat
v_gauss=rnorm(IClen, mean=0, 1) 
rou.hat=(mean(fsup(a,y)-fsup(a,v_gauss)))
    
#3. L and delta
L=rou.hat*((1/IClen)*(t(fsup_1(a,y))%*%X.whiten)) - ((mu/2)*(t(g.MSE_1(y.stand, ref.stand))%*%X.whiten)/(IClen))
delta=((rou.hat*mean(fsup_2(a,y)))) - ((mu/2)*(mean(g.MSE_2(y.stand, ref.stand))))

#4. Update Equations for mu 
mu=max(0, mu + (gamma * (g.MSE(y.stand,ref.stand)-(threshold*(1-exp(-loop))))))

#5. Update/Normalize Weight Vector
w=w-(t(learningRate*L*(1/delta)))
w=w/norm(w, type=c("F"))
    
#6.Check for Convergence
wchange=1-as.numeric(abs(t(oldw)%*%w))
  if(wchange<OverValue){
    flag=0
    print(sprintf("Converged after %d iteration",loop))
    }
    
  if(loop>=maxIter){
    flag=0
    print(sprintf("After %d iteration, still cannot converge",loop))
  }

print(wchange)
oldw=w
loop=loop+1

}
  
y=X.whiten%*%w
print("End of cICA algorithm")
  
output=list(y,w)
return(output) 
}
#########################################################################################################



####7. Signal To Noise Ratio###
SNR<-function(s,c){
  Var.S=var(s)
  MSE=(sum(g.MSE(c,s)))/length(S[,1])
  SNR=10*(log10(Var.S/MSE))
  return(SNR)
}
################################

####8. Individual Performance Index ####
IPI<-function(w,A){
  #i. Compute Permutation Matrix
  Perm.Vec=t(w)%*%A
  
  #ii. Define dimensions and maximum permutation vector
  M=length(Perm.Vec)
  p.k=max(abs(Perm.Vec))
  
  #iii. Compute IPI
  s.IPI=0
  for (j in 1:M){
    p.j=abs(Perm.Vec[,j])
    aux=p.j/p.k
    s.IPI=s.IPI+aux
    #print(aux)
    #print(s.IPI)
  }
  IPI=s.IPI-1
  #i.v Return IPI
  return(IPI)
}

#####################################TEST 1#####################################T
#Test matrixes Simulation 2: Random Frequencies#
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
#SIMULATION HYPERPARAMETERS
N=300
k=seq(1:N)
ts = 1e-4 # sampling period
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1
#Empty Data Frames
aux.SNR.c1=data.frame()
aux.SNR.c2=data.frame()
aux.SNR.c3=data.frame()
aux.SNR.fastICA=data.frame()
aux.SNR.PCA=data.frame()
SNR.df=data.frame()
IPI.df=data.frame()
aux.IPI.c1=data.frame()
aux.IPI.c2=data.frame()
aux.IPI.c3=data.frame()
aux.IPI.fastICA=data.frame()
Total.MSE=data.frame()

###SIMULATION###
for (j in 1:100){
#SIMULATION PARAMETERS#  
ts = 1e-4 # sampling period
f1 = runif(N, min=0.001, max=0.1)/ts    # true frequency
f2 = runif(N, min=0.001, max=0.1)/ts    
f3 = runif(N, min=0.001, max=0.1)/ts
S=matrix(0, N, 5)
  
S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k))
S[,2] = cos(2*pi*f2*ts*k)
S[,3] = cos(2*pi*f3*ts*k + 2)
S[,4] = rnorm(N)
S[,5] = rnorm(N)
  
A = matrix(runif(25), 5,5)
X = S%*%A + matrix(runif(1500), 300,5)
L=3 #Three relevant signals
#cICA parameters
ref1=S[,1]
ref2=S[,2]
ref3=S[,3]
ref=data.frame(ref1, ref2, ref3)
#Run fastICA
#Run fastICA
ic=fastICA(X, L, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
Xhat.ica=t(t(ic$A)%*%t(ic$S))

#PCA
Xpca=prcomp(X)
nComp=L
Xhat.pca.aux=Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
Xhat.pca=scale(Xhat.pca.aux, center = -colMeans(X), scale = FALSE)

#cICA LOOP
#Y Empty Data Frames
c1.y.df=data.frame(seq(1:N))
c2.y.df=data.frame(seq(1:N))
c3.y.df=data.frame(seq(1:N))
#W Empty Data Frames
c1.w.df=data.frame(seq(1:5))
c2.w.df=data.frame(seq(1:5))
c3.w.df=data.frame(seq(1:5))
for (i in 1:L){
  #Run cICA Algorithms
  c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
  c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
  c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
  #Extract ICs
  c1.y.df=cbind.data.frame(c1.y.df, as.numeric(c1[[1]]))
  c2.y.df=cbind.data.frame(c2.y.df, as.numeric(c2[[1]]))
  c3.y.df=cbind.data.frame(c3.y.df, as.numeric(c3[[1]]))
  #Extract W (only for cICA implementations)
  c1.w.df=cbind.data.frame(c1.w.df, c1[[2]])
  c2.w.df=cbind.data.frame(c2.w.df, c2[[2]])
  c3.w.df=cbind.data.frame(c3.w.df, c3[[2]])
  #Compute SNR
  SNR.df[i,1]=as.numeric(SNR(c1[[1]], S[,i]))
  SNR.df[i,2]=as.numeric(SNR(c2[[1]], S[,i]))
  SNR.df[i,3]=as.numeric(SNR(c3[[1]], S[,i]))
  SNR.df[i,4]=as.numeric(SNR(ic$S[,i],S[,i]))
  SNR.df[i,5]=as.numeric(SNR(as.matrix(Xpca$x[,i]), S[,i]))
  
  #Compute IPI (Ad-hoc measurement for ICA framework)
  IPI.df[i,1]=IPI(c1[[2]], A)
  IPI.df[i,2]=IPI(c2[[2]], A)
  IPI.df[i,3]=IPI(c3[[2]], A)
  IPI.df[i,4]=IPI(as.matrix(ic$A[i,]),A)

}
colnames(SNR.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA", "PCA")
colnames(IPI.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA")
#Compute Reconstruction Xhat
WhitenMatrix=whiten_ZCA(X)[[2]]
#c1: Lu-Rajapakse
c1.y=as.matrix(c1.y.df[,-1])
c1.w=as.matrix(c1.w.df[,-1])
Xhat.c1.aux=t(ginv(t(c1.w))%*%t(c1.y))
Xhat.c1=t(inv(WhitenMatrix)%*%t(Xhat.c1.aux))

#c2: Zilling-Zhang
c2.y=as.matrix(c2.y.df[,-1])
c2.w=as.matrix(c2.w.df[,-1])
Xhat.c2.aux=t(ginv(t(c2.w))%*%t(c2.y))
Xhat.c2=t(inv(WhitenMatrix)%*%t(Xhat.c2.aux))

#c3: Lin's
c3.y=as.matrix(c3.y.df[,-1])
c3.w=as.matrix(c3.w.df[,-1])
Xhat.c3.aux=t(ginv(t(c3.w))%*%t(c3.y))
Xhat.c3=t(inv(WhitenMatrix)%*%t(Xhat.c3.aux))

#Store Simulation Results
#SNR
aux.SNR.c1=rbind.data.frame(aux.SNR.c1, SNR.df[,1])
aux.SNR.c2=rbind.data.frame(aux.SNR.c2, SNR.df[,2])
aux.SNR.c3=rbind.data.frame(aux.SNR.c3, SNR.df[,3])
aux.SNR.fastICA=rbind.data.frame(aux.SNR.fastICA, SNR.df[,4])
aux.SNR.PCA=rbind.data.frame(aux.SNR.PCA, SNR.df[,5])
#IPI
aux.IPI.c1=rbind.data.frame(aux.IPI.c1, IPI.df[,1])
aux.IPI.c2=rbind.data.frame(aux.IPI.c2, IPI.df[,2])
aux.IPI.c3=rbind.data.frame(aux.IPI.c3, IPI.df[,3])
aux.IPI.fastICA=rbind.data.frame(aux.IPI.fastICA, IPI.df[,4])
#Reconstruction MSE/Amari Distance
MSE.df=data.frame(round(mse(Xhat.pca,X),5), round(mse(Xhat.ica,X),5), round(mse(Xhat.c1,X),5), round(mse(Xhat.c2,X),5), round(mse(Xhat.c3,X),5))
colnames(MSE.df)=c("PCA", "ICA", "Lu-Rajapakse", "Zhang", "Lin's")
Total.MSE=rbind.data.frame(Total.MSE, colMeans(MSE.df))
}

#Rename dataframes
#SNR
colnames(aux.SNR.c1)=c("S1", "S2", "S3")
colnames(aux.SNR.c2)=c("S1", "S2", "S3")
colnames(aux.SNR.c3)=c("S1", "S2", "S3")
colnames(aux.SNR.fastICA)=c("S1", "S2", "S3")
colnames(aux.SNR.PCA)=c("S1", "S2", "S3")
#IPI
colnames(aux.IPI.c1)=c("S1", "S2", "S3")
colnames(aux.IPI.c2)=c("S1", "S2", "S3")
colnames(aux.IPI.c3)=c("S1", "S2", "S3")
colnames(aux.IPI.fastICA)=c("S1", "S2", "S3")
#MSE
colnames(Total.MSE)=c("pca", "ica", "Lu-Rajapakse", "Zhang", "Lin's")


####Plotting Results####
library("boot", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
#SNR side-to-side comparison
for (i in 1:L){
a1=data.frame("LR",aux.SNR.c1[,i])
colnames(a1)=c("Algorithm", "SNR")
a2=data.frame("Zhang",aux.SNR.c2[,i])
colnames(a2)=c("Algorithm", "SNR")
a3=data.frame("Lin",aux.SNR.c3[,i])
colnames(a3)=c("Algorithm", "SNR")
a4=data.frame("fastICA",aux.SNR.fastICA[,i])
colnames(a4)=c("Algorithm", "SNR")
a5=data.frame("PCA",aux.SNR.PCA[,i])
colnames(a5)=c("Algorithm", "SNR")
a=rbind(a1, a2, a3, a4, a5)
p=ggplot(a, aes(x=SNR, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
  xlab("Signal-to-Noise Ratio") + ylab("Frequency Count")+labs(title = colnames(aux.SNR.c1)[i])
file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/UPDATED Simulation Results Plots/Sim 1 Noiseless/SNR/SNR_benchmark_", 
                  colnames(aux.SNR.c1)[i],".PNG", sep="")
ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
PCAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.SNR.c1[,i]
  bootCI.SNR.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=1000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.SNR.c1$basic[4], median(x1), bootCI.SNR.c1$basic[5]))
  #Zhang
  x2=aux.SNR.c2[,i]
  bootCI.SNR.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=1000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.SNR.c2$basic[4], median(x2), bootCI.SNR.c2$basic[5]))
  #Lin
  x3=aux.SNR.c3[,i]
  bootCI.SNR.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=1000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.SNR.c3$basic[4], median(x3), bootCI.SNR.c3$basic[5]))
  #fastICA
  x4=aux.SNR.fastICA[,i]
  bootCI.SNR.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=1000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.SNR.fastICA$basic[4], median(x4), bootCI.SNR.fastICA$basic[5]))
  #kPCA
  x5=aux.SNR.PCA[,i]
  bootCI.SNR.PCA=boot.ci(boot(x5, function(x,j) median(x[j]), R=1000))
  PCAResults=rbind.data.frame(PCAResults,data.frame(bootCI.SNR.PCA$basic[4], median(x5), bootCI.SNR.PCA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.SNR", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.SNR", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.SNR", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.SNR", "fastICA.upperCI")
colnames(PCAResults)=c("PCA.lowerCI", "PCA.Median.SNR", "PCA.upperCI")

SNR.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults, PCAResults)
write.csv(SNR.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/UPDATES/SNR_Results_Sim1.csv", row.names = FALSE)

#IPI side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.IPI.c1[,i])
  colnames(a1)=c("Algorithm", "IPI")
  a2=data.frame("Zhang",aux.IPI.c2[,i])
  colnames(a2)=c("Algorithm", "IPI")
  a3=data.frame("Lin",aux.IPI.c3[,i])
  colnames(a3)=c("Algorithm", "IPI")
  a4=data.frame("fastICA",aux.IPI.fastICA[,i])
  colnames(a4)=c("Algorithm", "IPI")
  a=rbind(a1, a2, a3, a4)
  p=ggplot(a, aes(x=IPI, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Individual Performance Index") + ylab("Frequency Count")+labs(title = colnames(aux.IPI.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/UPDATED Simulation Results Plots/Sim 1 Noiseless/IPI/IPI_benchmark_", 
                    colnames(aux.IPI.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.IPI.c1[,i]
  bootCI.IPI.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.IPI.c1$basic[4], median(x1), bootCI.IPI.c1$basic[5]))
  #Zhang
  x2=aux.IPI.c2[,i]
  bootCI.IPI.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.IPI.c2$basic[4], median(x2), bootCI.IPI.c2$basic[5]))
  #Lin
  x3=aux.IPI.c3[,i]
  bootCI.IPI.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.IPI.c3$basic[4], median(x3), bootCI.IPI.c3$basic[5]))
  #fastICA
  x4=aux.IPI.fastICA[,i]
  bootCI.IPI.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.IPI.fastICA$basic[4], median(x4), bootCI.IPI.fastICA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.IPI", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.IPI", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.IPI", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.IPI", "fastICA.upperCI")
IPI.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults)
write.csv(IPI.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/UPDATES/IPI_Results_Sim1.csv", row.names = FALSE)

#####################################TEST 2#####################################T
#Test matrixes Simulation 2: Random Frequencies#
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
#SIMULATION HYPERPARAMETERS
N=300
k=seq(1:N)
ts = 1e-4 # sampling period
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1
#Empty Data Frames
aux.SNR.c1=data.frame()
aux.SNR.c2=data.frame()
aux.SNR.c3=data.frame()
aux.SNR.fastICA=data.frame()
aux.SNR.PCA=data.frame()
SNR.df=data.frame()
IPI.df=data.frame()
aux.IPI.c1=data.frame()
aux.IPI.c2=data.frame()
aux.IPI.c3=data.frame()
aux.IPI.fastICA=data.frame()
Total.MSE=data.frame()

###SIMULATION###
for (j in 1:100){
  #SIMULATION PARAMETERS#  
  ts = 1e-4 # sampling period
  f1 = runif(N, min=0.001, max=0.1)/ts    # true frequency
  f2 = runif(N, min=0.001, max=0.1)/ts    
  f3 = runif(N, min=0.001, max=0.1)/ts
  S=matrix(0, N, 5)
  
  S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k)) 
  S[,2] = cos(2*pi*f2*ts*k)
  S[,3] = cos(2*pi*f3*ts*k + 2)
  S[,4] = rnorm(N)
  S[,5] = rnorm(N)
  
  A = matrix(runif(25), 5,5)
  X = S%*%A
  L=3 #Three relevant signals
  
  #cICA parameters
  ref1=S[,1] + rnorm(N, 0, 1)
  ref2=S[,2] + rnorm(N, 0, 1)
  ref3=S[,3] + rnorm(N, 0, 1)
  ref=data.frame(ref1, ref2, ref3)

  #Run fastICA
  #Run fastICA
  ic=fastICA(X, L, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
  Xhat.ica=t(t(ic$A)%*%t(ic$S))
  
  #PCA
  Xpca=prcomp(X)
  nComp=L
  Xhat.pca.aux=Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
  Xhat.pca=scale(Xhat.pca.aux, center = -colMeans(X), scale = FALSE)
  
  #cICA LOOP
  #Y Empty Data Frames
  c1.y.df=data.frame(seq(1:N))
  c2.y.df=data.frame(seq(1:N))
  c3.y.df=data.frame(seq(1:N))
  #W Empty Data Frames
  c1.w.df=data.frame(seq(1:5))
  c2.w.df=data.frame(seq(1:5))
  c3.w.df=data.frame(seq(1:5))
  for (i in 1:L){
    #Run cICA Algorithms
    c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    #Extract ICs
    c1.y.df=cbind.data.frame(c1.y.df, as.numeric(c1[[1]]))
    c2.y.df=cbind.data.frame(c2.y.df, as.numeric(c2[[1]]))
    c3.y.df=cbind.data.frame(c3.y.df, as.numeric(c3[[1]]))
    #Extract W (only for cICA implementations)
    c1.w.df=cbind.data.frame(c1.w.df, c1[[2]])
    c2.w.df=cbind.data.frame(c2.w.df, c2[[2]])
    c3.w.df=cbind.data.frame(c3.w.df, c3[[2]])
    #Compute SNR
    SNR.df[i,1]=as.numeric(SNR(c1[[1]], S[,i]))
    SNR.df[i,2]=as.numeric(SNR(c2[[1]], S[,i]))
    SNR.df[i,3]=as.numeric(SNR(c3[[1]], S[,i]))
    SNR.df[i,4]=as.numeric(SNR(ic$S[,i],S[,i]))
    SNR.df[i,5]=as.numeric(SNR(as.matrix(Xpca$x[,i]), S[,i]))
    
    #Compute IPI (Ad-hoc measurement for ICA framework)
    IPI.df[i,1]=IPI(c1[[2]], A)
    IPI.df[i,2]=IPI(c2[[2]], A)
    IPI.df[i,3]=IPI(c3[[2]], A)
    IPI.df[i,4]=IPI(as.matrix(ic$A[i,]),A)
    
  }
  colnames(SNR.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA", "PCA")
  colnames(IPI.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA")
  #Compute Reconstruction Xhat
  WhitenMatrix=whiten_ZCA(X)[[2]]
  #c1: Lu-Rajapakse
  c1.y=as.matrix(c1.y.df[,-1])
  c1.w=as.matrix(c1.w.df[,-1])
  Xhat.c1.aux=t(ginv(t(c1.w))%*%t(c1.y))
  Xhat.c1=t(inv(WhitenMatrix)%*%t(Xhat.c1.aux))
  
  #c2: Zilling-Zhang
  c2.y=as.matrix(c2.y.df[,-1])
  c2.w=as.matrix(c2.w.df[,-1])
  Xhat.c2.aux=t(ginv(t(c2.w))%*%t(c2.y))
  Xhat.c2=t(inv(WhitenMatrix)%*%t(Xhat.c2.aux))
  
  #c3: Lin's
  c3.y=as.matrix(c3.y.df[,-1])
  c3.w=as.matrix(c3.w.df[,-1])
  Xhat.c3.aux=t(ginv(t(c3.w))%*%t(c3.y))
  Xhat.c3=t(inv(WhitenMatrix)%*%t(Xhat.c3.aux))
  
  #Store Simulation Results
  #SNR
  aux.SNR.c1=rbind.data.frame(aux.SNR.c1, SNR.df[,1])
  aux.SNR.c2=rbind.data.frame(aux.SNR.c2, SNR.df[,2])
  aux.SNR.c3=rbind.data.frame(aux.SNR.c3, SNR.df[,3])
  aux.SNR.fastICA=rbind.data.frame(aux.SNR.fastICA, SNR.df[,4])
  aux.SNR.PCA=rbind.data.frame(aux.SNR.PCA, SNR.df[,5])
  #IPI
  aux.IPI.c1=rbind.data.frame(aux.IPI.c1, IPI.df[,1])
  aux.IPI.c2=rbind.data.frame(aux.IPI.c2, IPI.df[,2])
  aux.IPI.c3=rbind.data.frame(aux.IPI.c3, IPI.df[,3])
  aux.IPI.fastICA=rbind.data.frame(aux.IPI.fastICA, IPI.df[,4])
  #Reconstruction MSE
  MSE.df=data.frame(round(mse(Xhat.pca,X),5), round(mse(Xhat.ica,X),5), round(mse(Xhat.c1,X),5), round(mse(Xhat.c2,X),5), round(mse(Xhat.c3,X),5))
  colnames(MSE.df)=c("PCA", "ICA", "Lu-Rajapakse", "Zhang", "Lin's")
  Total.MSE=rbind.data.frame(Total.MSE, colMeans(MSE.df))
  
}

#Rename dataframes
#SNR
colnames(aux.SNR.c1)=c("S1", "S2", "S3")
colnames(aux.SNR.c2)=c("S1", "S2", "S3")
colnames(aux.SNR.c3)=c("S1", "S2", "S3")
colnames(aux.SNR.fastICA)=c("S1", "S2", "S3")
colnames(aux.SNR.PCA)=c("S1", "S2", "S3")
#IPI
colnames(aux.IPI.c1)=c("S1", "S2", "S3")
colnames(aux.IPI.c2)=c("S1", "S2", "S3")
colnames(aux.IPI.c3)=c("S1", "S2", "S3")
colnames(aux.IPI.fastICA)=c("S1", "S2", "S3")
#MSE
colnames(Total.MSE)=c("pca", "ica", "Lu-Rajapakse", "Zhang", "Lin's")


####Plotting Results####
library("boot", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
#SNR side-to-side comparison
for (i in 1:L){
  a1=data.frame("LR",aux.SNR.c1[,i])
  colnames(a1)=c("Algorithm", "SNR")
  a2=data.frame("Zhang",aux.SNR.c2[,i])
  colnames(a2)=c("Algorithm", "SNR")
  a3=data.frame("Lin",aux.SNR.c3[,i])
  colnames(a3)=c("Algorithm", "SNR")
  a4=data.frame("fastICA",aux.SNR.fastICA[,i])
  colnames(a4)=c("Algorithm", "SNR")
  a5=data.frame("PCA",aux.SNR.PCA[,i])
  colnames(a5)=c("Algorithm", "SNR")
  a=rbind(a1, a2, a3, a4, a5)
  p=ggplot(a, aes(x=SNR, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Signal-to-Noise Ratio") + ylab("Frequency Count")+labs(title = colnames(aux.SNR.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/UPDATED Simulation Results Plots/Sim 2 White Noise/SNR/SNR_benchmark_", 
                    colnames(aux.SNR.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
PCAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.SNR.c1[,i]
  bootCI.SNR.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=1000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.SNR.c1$basic[4], median(x1), bootCI.SNR.c1$basic[5]))
  #Zhang
  x2=aux.SNR.c2[,i]
  bootCI.SNR.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=1000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.SNR.c2$basic[4], median(x2), bootCI.SNR.c2$basic[5]))
  #Lin
  x3=aux.SNR.c3[,i]
  bootCI.SNR.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=1000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.SNR.c3$basic[4], median(x3), bootCI.SNR.c3$basic[5]))
  #fastICA
  x4=aux.SNR.fastICA[,i]
  bootCI.SNR.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=1000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.SNR.fastICA$basic[4], median(x4), bootCI.SNR.fastICA$basic[5]))
  #kPCA
  x5=aux.SNR.PCA[,i]
  bootCI.SNR.PCA=boot.ci(boot(x5, function(x,j) median(x[j]), R=1000))
  PCAResults=rbind.data.frame(PCAResults,data.frame(bootCI.SNR.PCA$basic[4], median(x5), bootCI.SNR.PCA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.SNR", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.SNR", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.SNR", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.SNR", "fastICA.upperCI")
colnames(PCAResults)=c("PCA.lowerCI", "PCA.Median.SNR", "PCA.upperCI")

SNR.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults, PCAResults)
write.csv(SNR.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/UPDATES/SNR_Results_Sim2.csv", row.names = FALSE)

#IPI side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.IPI.c1[,i])
  colnames(a1)=c("Algorithm", "IPI")
  a2=data.frame("Zhang",aux.IPI.c2[,i])
  colnames(a2)=c("Algorithm", "IPI")
  a3=data.frame("Lin",aux.IPI.c3[,i])
  colnames(a3)=c("Algorithm", "IPI")
  a4=data.frame("fastICA",aux.IPI.fastICA[,i])
  colnames(a4)=c("Algorithm", "IPI")
  a=rbind(a1, a2, a3, a4)
  p=ggplot(a, aes(x=IPI, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Individual Performance Index") + ylab("Frequency Count")+labs(title = colnames(aux.IPI.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/UPDATED Simulation Results Plots/Sim 2 White Noise/IPI/IPI_benchmark_", 
                    colnames(aux.IPI.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.IPI.c1[,i]
  bootCI.IPI.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.IPI.c1$basic[4], median(x1), bootCI.IPI.c1$basic[5]))
  #Zhang
  x2=aux.IPI.c2[,i]
  bootCI.IPI.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.IPI.c2$basic[4], median(x2), bootCI.IPI.c2$basic[5]))
  #Lin
  x3=aux.IPI.c3[,i]
  bootCI.IPI.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.IPI.c3$basic[4], median(x3), bootCI.IPI.c3$basic[5]))
  #fastICA
  x4=aux.IPI.fastICA[,i]
  bootCI.IPI.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.IPI.fastICA$basic[4], median(x4), bootCI.IPI.fastICA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.IPI", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.IPI", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.IPI", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.IPI", "fastICA.upperCI")
IPI.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults)
write.csv(IPI.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/UPDATES/IPI_Results_Sim2.csv", row.names = FALSE)

#####################################TEST 3#####################################T
#Test matrixes Simulation 2: Random Frequencies#
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
#SIMULATION HYPERPARAMETERS
N=300
k=seq(1:N)
ts = 1e-4 # sampling period
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1
#Empty Data Frames
aux.SNR.c1=data.frame()
aux.SNR.c2=data.frame()
aux.SNR.c3=data.frame()
aux.SNR.fastICA=data.frame()
aux.SNR.PCA=data.frame()
SNR.df=data.frame()
IPI.df=data.frame()
aux.IPI.c1=data.frame()
aux.IPI.c2=data.frame()
aux.IPI.c3=data.frame()
aux.IPI.fastICA=data.frame()
Total.MSE=data.frame()

###SIMULATION###
for (j in 1:100){
  #SIMULATION PARAMETERS#  
  ts = 1e-4 # sampling period
  f1 = runif(N, min=0.001, max=0.1)/ts    # true frequency
  f2 = runif(N, min=0.001, max=0.1)/ts    
  f3 = runif(N, min=0.001, max=0.1)/ts
  S=matrix(0, N, 5)
  
  S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k)) 
  S[,2] = cos(2*pi*f2*ts*k)
  S[,3] = cos(2*pi*f3*ts*k + 2)
  S[,4] = rnorm(N)
  S[,5] = rnorm(N)
  
  A = matrix(runif(25), 5,5)
  X = S%*%A
  L=3 #Three relevant signals
  #cICA parameters
  ref1=S[,1] + rnorm(N, mean(S[,1]), sd(S[,1]))
  ref2=S[,2] + rnorm(N, mean(S[,2]), sd(S[,2]))
  ref3=S[,3] + rnorm(N, mean(S[,3]), sd(S[,3]))
  ref=data.frame(ref1, ref2, ref3)
  #Run fastICA
  #Run fastICA
  ic=fastICA(X, L, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
  Xhat.ica=t(t(ic$A)%*%t(ic$S))
  
  #PCA
  Xpca=prcomp(X)
  nComp=L
  Xhat.pca.aux=Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
  Xhat.pca=scale(Xhat.pca.aux, center = -colMeans(X), scale = FALSE)
  
  #cICA LOOP
  #Y Empty Data Frames
  c1.y.df=data.frame(seq(1:N))
  c2.y.df=data.frame(seq(1:N))
  c3.y.df=data.frame(seq(1:N))
  #W Empty Data Frames
  c1.w.df=data.frame(seq(1:5))
  c2.w.df=data.frame(seq(1:5))
  c3.w.df=data.frame(seq(1:5))
  for (i in 1:L){
    #Run cICA Algorithms
    c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    #Extract ICs
    c1.y.df=cbind.data.frame(c1.y.df, as.numeric(c1[[1]]))
    c2.y.df=cbind.data.frame(c2.y.df, as.numeric(c2[[1]]))
    c3.y.df=cbind.data.frame(c3.y.df, as.numeric(c3[[1]]))
    #Extract W (only for cICA implementations)
    c1.w.df=cbind.data.frame(c1.w.df, c1[[2]])
    c2.w.df=cbind.data.frame(c2.w.df, c2[[2]])
    c3.w.df=cbind.data.frame(c3.w.df, c3[[2]])
    #Compute SNR
    SNR.df[i,1]=as.numeric(SNR(c1[[1]], S[,i]))
    SNR.df[i,2]=as.numeric(SNR(c2[[1]], S[,i]))
    SNR.df[i,3]=as.numeric(SNR(c3[[1]], S[,i]))
    SNR.df[i,4]=as.numeric(SNR(ic$S[,i],S[,i]))
    SNR.df[i,5]=as.numeric(SNR(as.matrix(Xpca$x[,i]), S[,i]))
    
    #Compute IPI (Ad-hoc measurement for ICA framework)
    IPI.df[i,1]=IPI(c1[[2]], A)
    IPI.df[i,2]=IPI(c2[[2]], A)
    IPI.df[i,3]=IPI(c3[[2]], A)
    IPI.df[i,4]=IPI(as.matrix(ic$A[i,]),A)
    
  }
  colnames(SNR.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA", "PCA")
  colnames(IPI.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA")
  #Compute Reconstruction Xhat
  WhitenMatrix=whiten_ZCA(X)[[2]]
  #c1: Lu-Rajapakse
  c1.y=as.matrix(c1.y.df[,-1])
  c1.w=as.matrix(c1.w.df[,-1])
  Xhat.c1.aux=t(ginv(t(c1.w))%*%t(c1.y))
  Xhat.c1=t(inv(WhitenMatrix)%*%t(Xhat.c1.aux))
  
  #c2: Zilling-Zhang
  c2.y=as.matrix(c2.y.df[,-1])
  c2.w=as.matrix(c2.w.df[,-1])
  Xhat.c2.aux=t(ginv(t(c2.w))%*%t(c2.y))
  Xhat.c2=t(inv(WhitenMatrix)%*%t(Xhat.c2.aux))
  
  #c3: Lin's
  c3.y=as.matrix(c3.y.df[,-1])
  c3.w=as.matrix(c3.w.df[,-1])
  Xhat.c3.aux=t(ginv(t(c3.w))%*%t(c3.y))
  Xhat.c3=t(inv(WhitenMatrix)%*%t(Xhat.c3.aux))
  
  #Store Simulation Results
  #SNR
  aux.SNR.c1=rbind.data.frame(aux.SNR.c1, SNR.df[,1])
  aux.SNR.c2=rbind.data.frame(aux.SNR.c2, SNR.df[,2])
  aux.SNR.c3=rbind.data.frame(aux.SNR.c3, SNR.df[,3])
  aux.SNR.fastICA=rbind.data.frame(aux.SNR.fastICA, SNR.df[,4])
  aux.SNR.PCA=rbind.data.frame(aux.SNR.PCA, SNR.df[,5])
  #IPI
  aux.IPI.c1=rbind.data.frame(aux.IPI.c1, IPI.df[,1])
  aux.IPI.c2=rbind.data.frame(aux.IPI.c2, IPI.df[,2])
  aux.IPI.c3=rbind.data.frame(aux.IPI.c3, IPI.df[,3])
  aux.IPI.fastICA=rbind.data.frame(aux.IPI.fastICA, IPI.df[,4])
  #Reconstruction MSE
  MSE.df=data.frame(round(mse(Xhat.pca,X),5), round(mse(Xhat.ica,X),5), round(mse(Xhat.c1,X),5), round(mse(Xhat.c2,X),5), round(mse(Xhat.c3,X),5))
  colnames(MSE.df)=c("PCA", "ICA", "Lu-Rajapakse", "Zhang", "Lin's")
  Total.MSE=rbind.data.frame(Total.MSE, colMeans(MSE.df))
  
}

#Rename dataframes
#SNR
colnames(aux.SNR.c1)=c("S1", "S2", "S3")
colnames(aux.SNR.c2)=c("S1", "S2", "S3")
colnames(aux.SNR.c3)=c("S1", "S2", "S3")
colnames(aux.SNR.fastICA)=c("S1", "S2", "S3")
colnames(aux.SNR.PCA)=c("S1", "S2", "S3")
#IPI
colnames(aux.IPI.c1)=c("S1", "S2", "S3")
colnames(aux.IPI.c2)=c("S1", "S2", "S3")
colnames(aux.IPI.c3)=c("S1", "S2", "S3")
colnames(aux.IPI.fastICA)=c("S1", "S2", "S3")
#MSE
colnames(Total.MSE)=c("pca", "ica", "Lu-Rajapakse", "Zhang", "Lin's")


####Plotting Results####
library("boot", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
#SNR side-to-side comparison
for (i in 1:L){
  a1=data.frame("LR",aux.SNR.c1[,i])
  colnames(a1)=c("Algorithm", "SNR")
  a2=data.frame("Zhang",aux.SNR.c2[,i])
  colnames(a2)=c("Algorithm", "SNR")
  a3=data.frame("Lin",aux.SNR.c3[,i])
  colnames(a3)=c("Algorithm", "SNR")
  a4=data.frame("fastICA",aux.SNR.fastICA[,i])
  colnames(a4)=c("Algorithm", "SNR")
  a5=data.frame("PCA",aux.SNR.PCA[,i])
  colnames(a5)=c("Algorithm", "SNR")
  a=rbind(a1, a2, a3, a4, a5)
  p=ggplot(a, aes(x=SNR, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Signal-to-Noise Ratio") + ylab("Frequency Count")+labs(title = colnames(aux.SNR.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/UPDATED Simulation Results Plots/Sim 3 Gaussian Noise/SNR/SNR_benchmark_", 
                    colnames(aux.SNR.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
PCAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.SNR.c1[,i]
  bootCI.SNR.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=1000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.SNR.c1$basic[4], median(x1), bootCI.SNR.c1$basic[5]))
  #Zhang
  x2=aux.SNR.c2[,i]
  bootCI.SNR.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=1000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.SNR.c2$basic[4], median(x2), bootCI.SNR.c2$basic[5]))
  #Lin
  x3=aux.SNR.c3[,i]
  bootCI.SNR.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=1000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.SNR.c3$basic[4], median(x3), bootCI.SNR.c3$basic[5]))
  #fastICA
  x4=aux.SNR.fastICA[,i]
  bootCI.SNR.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=1000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.SNR.fastICA$basic[4], median(x4), bootCI.SNR.fastICA$basic[5]))
  #kPCA
  x5=aux.SNR.PCA[,i]
  bootCI.SNR.PCA=boot.ci(boot(x5, function(x,j) median(x[j]), R=1000))
  PCAResults=rbind.data.frame(PCAResults,data.frame(bootCI.SNR.PCA$basic[4], median(x5), bootCI.SNR.PCA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.SNR", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.SNR", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.SNR", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.SNR", "fastICA.upperCI")
colnames(PCAResults)=c("PCA.lowerCI", "PCA.Median.SNR", "PCA.upperCI")

SNR.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults, PCAResults)
write.csv(SNR.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/UPDATES/SNR_Results_Sim3.csv", row.names = FALSE)

#IPI side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.IPI.c1[,i])
  colnames(a1)=c("Algorithm", "IPI")
  a2=data.frame("Zhang",aux.IPI.c2[,i])
  colnames(a2)=c("Algorithm", "IPI")
  a3=data.frame("Lin",aux.IPI.c3[,i])
  colnames(a3)=c("Algorithm", "IPI")
  a4=data.frame("fastICA",aux.IPI.fastICA[,i])
  colnames(a4)=c("Algorithm", "IPI")
  a=rbind(a1, a2, a3, a4)
  p=ggplot(a, aes(x=IPI, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Individual Performance Index") + ylab("Frequency Count")+labs(title = colnames(aux.IPI.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/UPDATED Simulation Results Plots/Sim 3 Gaussian Noise/IPI/IPI_benchmark_", 
                    colnames(aux.IPI.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.IPI.c1[,i]
  bootCI.IPI.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.IPI.c1$basic[4], median(x1), bootCI.IPI.c1$basic[5]))
  #Zhang
  x2=aux.IPI.c2[,i]
  bootCI.IPI.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.IPI.c2$basic[4], median(x2), bootCI.IPI.c2$basic[5]))
  #Lin
  x3=aux.IPI.c3[,i]
  bootCI.IPI.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.IPI.c3$basic[4], median(x3), bootCI.IPI.c3$basic[5]))
  #fastICA
  x4=aux.IPI.fastICA[,i]
  bootCI.IPI.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.IPI.fastICA$basic[4], median(x4), bootCI.IPI.fastICA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.IPI", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.IPI", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.IPI", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.IPI", "fastICA.upperCI")
IPI.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults)
write.csv(IPI.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/UPDATES/IPI_Results_Sim3.csv", row.names = FALSE)

#####################################TEST 4#####################################T
#Test matrixes Simulation 2: Random Frequencies#
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
#SIMULATION HYPERPARAMETERS
N=300
k=seq(1:N)
ts = 1e-4 # sampling period
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1
#Empty Data Frames
aux.SNR.c1=data.frame()
aux.SNR.c2=data.frame()
aux.SNR.c3=data.frame()
aux.SNR.fastICA=data.frame()
aux.SNR.PCA=data.frame()
SNR.df=data.frame()
IPI.df=data.frame()
aux.IPI.c1=data.frame()
aux.IPI.c2=data.frame()
aux.IPI.c3=data.frame()
aux.IPI.fastICA=data.frame()
Total.MSE=data.frame()

###SIMULATION###
for (j in 1:100){
  #SIMULATION PARAMETERS#  
  ts = 1e-4 # sampling period
  f1 = runif(N, min=0.001, max=0.1)/ts    # true frequency
  f2 = runif(N, min=0.001, max=0.1)/ts    
  f3 = runif(N, min=0.001, max=0.1)/ts
  S=matrix(0, N, 5)
  
  S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k)) 
  S[,2] = cos(2*pi*f2*ts*k)
  S[,3] = cos(2*pi*f3*ts*k + 2)
  S[,4] = rnorm(N)
  S[,5] = rnorm(N)
  
  A = matrix(runif(25), 5,5)
  X = S%*%A
  L=3 #Three relevant signals
  
  #cICA parameters
  p1=spec.pgram(X[,1], plot=FALSE, taper=0.1, demean=TRUE)
  freq1=p1$freq[which(p1$spec>=sort(p1$spec, decreasing = TRUE)[3], arr.ind = TRUE)]
  ref1=sin(freq1[1]*k)
  p2=spec.pgram(X[,2], plot=FALSE, taper=0.1, demean=TRUE)
  freq2=p2$freq[which(p2$spec>=sort(p2$spec, decreasing = TRUE)[3], arr.ind = TRUE)]
  ref2=cos(freq2[1]*k)
  p3=spec.pgram(X[,3], plot=FALSE, taper=0.1, demean=TRUE)
  freq3=p3$freq[which(p3$spec>=sort(p3$spec, decreasing = TRUE)[3], arr.ind = TRUE)]
  ref3=cos(freq3[1]*k)
  ref=data.frame(ref1, ref2, ref3)
  #Run fastICA
  #Run fastICA
  ic=fastICA(X, L, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
  Xhat.ica=t(t(ic$A)%*%t(ic$S))
  
  #PCA
  Xpca=prcomp(X)
  nComp=L
  Xhat.pca.aux=Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
  Xhat.pca=scale(Xhat.pca.aux, center = -colMeans(X), scale = FALSE)
  
  #cICA LOOP
  #Y Empty Data Frames
  c1.y.df=data.frame(seq(1:N))
  c2.y.df=data.frame(seq(1:N))
  c3.y.df=data.frame(seq(1:N))
  #W Empty Data Frames
  c1.w.df=data.frame(seq(1:5))
  c2.w.df=data.frame(seq(1:5))
  c3.w.df=data.frame(seq(1:5))
  for (i in 1:L){
    #Run cICA Algorithms
    c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    #Extract ICs
    c1.y.df=cbind.data.frame(c1.y.df, as.numeric(c1[[1]]))
    c2.y.df=cbind.data.frame(c2.y.df, as.numeric(c2[[1]]))
    c3.y.df=cbind.data.frame(c3.y.df, as.numeric(c3[[1]]))
    #Extract W (only for cICA implementations)
    c1.w.df=cbind.data.frame(c1.w.df, c1[[2]])
    c2.w.df=cbind.data.frame(c2.w.df, c2[[2]])
    c3.w.df=cbind.data.frame(c3.w.df, c3[[2]])
    #Compute SNR
    SNR.df[i,1]=as.numeric(SNR(c1[[1]], S[,i]))
    SNR.df[i,2]=as.numeric(SNR(c2[[1]], S[,i]))
    SNR.df[i,3]=as.numeric(SNR(c3[[1]], S[,i]))
    SNR.df[i,4]=as.numeric(SNR(ic$S[,i],S[,i]))
    SNR.df[i,5]=as.numeric(SNR(as.matrix(Xpca$x[,i]), S[,i]))
    
    #Compute IPI (Ad-hoc measurement for ICA framework)
    IPI.df[i,1]=IPI(c1[[2]], A)
    IPI.df[i,2]=IPI(c2[[2]], A)
    IPI.df[i,3]=IPI(c3[[2]], A)
    IPI.df[i,4]=IPI(as.matrix(ic$A[i,]),A)
    
  }
  colnames(SNR.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA", "PCA")
  colnames(IPI.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA")
  #Compute Reconstruction Xhat
  WhitenMatrix=whiten_ZCA(X)[[2]]
  #c1: Lu-Rajapakse
  c1.y=as.matrix(c1.y.df[,-1])
  c1.w=as.matrix(c1.w.df[,-1])
  Xhat.c1.aux=t(ginv(t(c1.w))%*%t(c1.y))
  Xhat.c1=t(inv(WhitenMatrix)%*%t(Xhat.c1.aux))
  
  #c2: Zilling-Zhang
  c2.y=as.matrix(c2.y.df[,-1])
  c2.w=as.matrix(c2.w.df[,-1])
  Xhat.c2.aux=t(ginv(t(c2.w))%*%t(c2.y))
  Xhat.c2=t(inv(WhitenMatrix)%*%t(Xhat.c2.aux))
  
  #c3: Lin's
  c3.y=as.matrix(c3.y.df[,-1])
  c3.w=as.matrix(c3.w.df[,-1])
  Xhat.c3.aux=t(ginv(t(c3.w))%*%t(c3.y))
  Xhat.c3=t(inv(WhitenMatrix)%*%t(Xhat.c3.aux))
  
  #Store Simulation Results
  #SNR
  aux.SNR.c1=rbind.data.frame(aux.SNR.c1, SNR.df[,1])
  aux.SNR.c2=rbind.data.frame(aux.SNR.c2, SNR.df[,2])
  aux.SNR.c3=rbind.data.frame(aux.SNR.c3, SNR.df[,3])
  aux.SNR.fastICA=rbind.data.frame(aux.SNR.fastICA, SNR.df[,4])
  aux.SNR.PCA=rbind.data.frame(aux.SNR.PCA, SNR.df[,5])
  #IPI
  aux.IPI.c1=rbind.data.frame(aux.IPI.c1, IPI.df[,1])
  aux.IPI.c2=rbind.data.frame(aux.IPI.c2, IPI.df[,2])
  aux.IPI.c3=rbind.data.frame(aux.IPI.c3, IPI.df[,3])
  aux.IPI.fastICA=rbind.data.frame(aux.IPI.fastICA, IPI.df[,4])
  #Reconstruction MSE
  MSE.df=data.frame(round(mse(Xhat.pca,X),5), round(mse(Xhat.ica,X),5), round(mse(Xhat.c1,X),5), round(mse(Xhat.c2,X),5), round(mse(Xhat.c3,X),5))
  colnames(MSE.df)=c("PCA", "ICA", "Lu-Rajapakse", "Zhang", "Lin's")
  Total.MSE=rbind.data.frame(Total.MSE, colMeans(MSE.df))
  
}

#Rename dataframes
#SNR
colnames(aux.SNR.c1)=c("S1", "S2", "S3")
colnames(aux.SNR.c2)=c("S1", "S2", "S3")
colnames(aux.SNR.c3)=c("S1", "S2", "S3")
colnames(aux.SNR.fastICA)=c("S1", "S2", "S3")
colnames(aux.SNR.PCA)=c("S1", "S2", "S3")
#IPI
colnames(aux.IPI.c1)=c("S1", "S2", "S3")
colnames(aux.IPI.c2)=c("S1", "S2", "S3")
colnames(aux.IPI.c3)=c("S1", "S2", "S3")
colnames(aux.IPI.fastICA)=c("S1", "S2", "S3")
#MSE
colnames(Total.MSE)=c("pca", "ica", "Lu-Rajapakse", "Zhang", "Lin's")


####Plotting Results####
library("boot", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
#SNR side-to-side comparison
for (i in 1:L){
  a1=data.frame("LR",aux.SNR.c1[,i])
  colnames(a1)=c("Algorithm", "SNR")
  a2=data.frame("Zhang",aux.SNR.c2[,i])
  colnames(a2)=c("Algorithm", "SNR")
  a3=data.frame("Lin",aux.SNR.c3[,i])
  colnames(a3)=c("Algorithm", "SNR")
  a4=data.frame("fastICA",aux.SNR.fastICA[,i])
  colnames(a4)=c("Algorithm", "SNR")
  a5=data.frame("PCA",aux.SNR.PCA[,i])
  colnames(a5)=c("Algorithm", "SNR")
  a=rbind(a1, a2, a3, a4, a5)
  p=ggplot(a, aes(x=SNR, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Signal-to-Noise Ratio") + ylab("Frequency Count")+labs(title = colnames(aux.SNR.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/UPDATED Simulation Results Plots/Sim 4 Rough Ref Signal/SNR/SNR_benchmark_", 
                    colnames(aux.SNR.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
PCAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.SNR.c1[,i]
  bootCI.SNR.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=1000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.SNR.c1$basic[4], median(x1), bootCI.SNR.c1$basic[5]))
  #Zhang
  x2=aux.SNR.c2[,i]
  bootCI.SNR.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=1000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.SNR.c2$basic[4], median(x2), bootCI.SNR.c2$basic[5]))
  #Lin
  x3=aux.SNR.c3[,i]
  bootCI.SNR.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=1000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.SNR.c3$basic[4], median(x3), bootCI.SNR.c3$basic[5]))
  #fastICA
  x4=aux.SNR.fastICA[,i]
  bootCI.SNR.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=1000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.SNR.fastICA$basic[4], median(x4), bootCI.SNR.fastICA$basic[5]))
  #kPCA
  x5=aux.SNR.PCA[,i]
  bootCI.SNR.PCA=boot.ci(boot(x5, function(x,j) median(x[j]), R=1000))
  PCAResults=rbind.data.frame(PCAResults,data.frame(bootCI.SNR.PCA$basic[4], median(x5), bootCI.SNR.PCA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.SNR", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.SNR", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.SNR", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.SNR", "fastICA.upperCI")
colnames(PCAResults)=c("PCA.lowerCI", "PCA.Median.SNR", "PCA.upperCI")

SNR.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults, PCAResults)
write.csv(SNR.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/UPDATES/SNR_Results_Sim4.csv", row.names = FALSE)

#IPI side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.IPI.c1[,i])
  colnames(a1)=c("Algorithm", "IPI")
  a2=data.frame("Zhang",aux.IPI.c2[,i])
  colnames(a2)=c("Algorithm", "IPI")
  a3=data.frame("Lin",aux.IPI.c3[,i])
  colnames(a3)=c("Algorithm", "IPI")
  a4=data.frame("fastICA",aux.IPI.fastICA[,i])
  colnames(a4)=c("Algorithm", "IPI")
  a=rbind(a1, a2, a3, a4)
  p=ggplot(a, aes(x=IPI, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Individual Performance Index") + ylab("Frequency Count")+labs(title = colnames(aux.IPI.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/UPDATED Simulation Results Plots/Sim 4 Rough Ref Signal/IPI/IPI_benchmark_", 
                    colnames(aux.IPI.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.IPI.c1[,i]
  bootCI.IPI.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.IPI.c1$basic[4], median(x1), bootCI.IPI.c1$basic[5]))
  #Zhang
  x2=aux.IPI.c2[,i]
  bootCI.IPI.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.IPI.c2$basic[4], median(x2), bootCI.IPI.c2$basic[5]))
  #Lin
  x3=aux.IPI.c3[,i]
  bootCI.IPI.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.IPI.c3$basic[4], median(x3), bootCI.IPI.c3$basic[5]))
  #fastICA
  x4=aux.IPI.fastICA[,i]
  bootCI.IPI.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.IPI.fastICA$basic[4], median(x4), bootCI.IPI.fastICA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.IPI", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.IPI", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.IPI", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.IPI", "fastICA.upperCI")
IPI.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults)
write.csv(IPI.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/UPDATES/IPI_Results_Sim4.csv", row.names = FALSE)

#####################################TEST 6 HIERARCHICAL DATA SET FULL DATA#####################################
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
library("aTSA", lib.loc="~/R/win-library/3.5")
library("tseries", lib.loc="~/R/win-library/3.5")
library("TSA", lib.loc="~/R/win-library/3.5")
library("forecast", lib.loc="~/R/win-library/3.5")

#Load Simulated Data Sets
Sim.Feeders=data.frame(read.csv("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulated_Feeder.csv", header = TRUE, sep = ","))
Sim.SSEE=data.frame(read.csv("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulated_SSEE.csv", header = TRUE, sep = ","))
Sim.Feeders=Sim.Feeders[,2:dim(Sim.Feeders)[2]]
Sim.SSEE=Sim.SSEE[,2:dim(Sim.SSEE)[2]]

#Decompose SSEE
#Noise
Noise.SSEE=data.frame(seasadj(stl(msts(Sim.SSEE[,2], seasonal.periods=c(24,7*24)), "periodic")))
for (i in 3:dim(Sim.SSEE)[2]){
  Noise.SSEE=data.frame(Noise.SSEE, seasadj(stl(msts(Sim.SSEE[,i], seasonal.periods=c(24,7*24)), "periodic")))
}
Noise.SSEE=data.frame(Sim.SSEE[,1], Noise.SSEE)
names(Noise.SSEE)=c("DateTime", paste0("Noise_",colnames(Sim.SSEE[,2:dim(Sim.SSEE)[2]])))


#3.2 Seasonal Pattern (Substracting the noice from the obtained TS)
Season.Pattern=data.frame(Sim.SSEE[,2]-seasadj(stl(msts(Sim.SSEE[,2], seasonal.periods=c(24,7*24)), "periodic")))
for (i in 3:dim(Sim.SSEE)[2]){
  Season.Pattern=data.frame(Season.Pattern, Sim.SSEE[,i]-seasadj(stl(msts(Sim.SSEE[,i], seasonal.periods=c(24,7*24)), "periodic")))
}
Season.Pattern=data.frame(Sim.SSEE[,1], Season.Pattern)
names(Season.Pattern)=c("DateTime", paste0("Season_",colnames(Sim.SSEE[,2:dim(Sim.SSEE)[2]])))

Transformed.SSEE=data.frame(Noise.SSEE, Season.Pattern[,2:dim(Season.Pattern)[2]])

###1. REFERENCE SIGNAL EXTRACTION FROM TRANSFORMED SSEE
library("dplyr", lib.loc="~/R/win-library/3.5")
library("tsoutliers", lib.loc="~/R/win-library/3.5")
library("TSclust", lib.loc="~/R/win-library/3.5")
library("factoextra", lib.loc="~/R/win-library/3.5")
library("dendextend", lib.loc="~/R/win-library/3.5")
library("cluster", lib.loc="C:/Program Files/R/R-3.5.1/library")
library("NbClust", lib.loc="~/R/win-library/3.5")
library("dtwclust", lib.loc="~/R/win-library/3.5")
library("hydroGOF", lib.loc="~/R/win-library/3.5")
#1.1 SELECT NUMBER OF CLUSTERS
#DoParallel routine
require("doParallel")
# Create parallel workers
workers <- makeCluster(6L)
# Preload dtwclust in each worker; not necessary but useful
invisible(clusterEvalQ(workers, library("factoextra")))
# Register the backend; this step MUST be done
registerDoParallel(workers)
set.seed(12345)
AvgSilhouette=fviz_nbclust(Transformed.SSEE[,2:dim(Transformed.SSEE)[2]], FUN = hcut, hc_func=c("hclust"), hc_method=c("complete"), method = "silhouette")
#gap_stat=clusGap(Final.Data.Set, FUN = hcut, hc_func=c("hclust"), hc_method=c("complete"), hc_metric = "pearson",
#                 K.max = round(0.05*dim(Output.Period.Outliers)[2]), B = 50, method="globalmax")
# Stop parallel workers
stopCluster(workers)
# Go back to sequential computation
registerDoSEQ()
#3. Extract Maximum Silhouette
Cluster.Sil=AvgSilhouette$data
CLUS.NUM=as.numeric(Cluster.Sil[which.max(Cluster.Sil$y),1])

#1.2 RUN AGGLOMERATIVE HCLUST COMPLETE LINKAGE
dist_ts=TSclust::diss(SERIES = t(Transformed.SSEE[,2:dim(Transformed.SSEE)[2]]), METHOD = "INT.PER") # note the dataframe must be transposed
hcPER.complete=hclust(dist_ts, method="complete")
plot(hcPER.complete)
rect.hclust(hcPER.complete, k=CLUS.NUM, border=2:12) #Change Clus.Num
sub_grp.complete=cutree(hcPER.complete, k = CLUS.NUM) #Change Clus.Num
t=data.frame(sub_grp.complete)
Output.HClust=data.frame(row.names(t), t, row.names=NULL)
names(Output.HClust)=c("Alim_ID", "Cluster")

#1.3 Compute Cluster Centroids
#Extract Memberships & Centroids
Clust.Center=list()
for (i in (1:CLUS.NUM)){
  aux1=as.character(t(Output.HClust%>%dplyr::filter(Cluster==i)%>%dplyr::select(Alim_ID)))
  idx=match(aux1, names(Transformed.SSEE))
  aux2=data.frame(Transformed.SSEE[,2],Transformed.SSEE[, idx])
  Clust.Center[[i]]=base::rowMeans(aux2[,2:dim(aux2)[2]])
}

#1.4 Reference Signals
ref=data.frame(Clust.Center[[1]])
for (i in 2:CLUS.NUM){
  ref=cbind.data.frame(ref,Clust.Center[[i]])
}

#2. BUILD BLIND SOURCE SETTING AND RESULTS DATA FRAMES
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1.5
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1
a0=1
b0=1
ro0=1

X=as.matrix(Sim.Feeders[,2:dim(Sim.Feeders)[2]])
mu_X=colMeans(X)
centered_X=centering_icaR(X)
Total.MSE=data.frame()


for (j in 1:100){
#2.1 Empty Data Frames
#Y Empty Data Frames
c1.y.df=data.frame(Sim.SSEE[,1])
c2.y.df=data.frame(Sim.SSEE[,1])
c3.y.df=data.frame(Sim.SSEE[,1])

#W Empty Data Frames
c1.w.df=data.frame(colnames(X))
c2.w.df=data.frame(colnames(X))
c3.w.df=data.frame(colnames(X))

#3. RUN ICA ALGORITMS!
#Benchmarks
#Run fastICA
ic=fastICA(X, CLUS.NUM, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
Xhat.ica=t(t(ic$A)%*%t(ic$S))

#PCA
Xpca=prcomp(X)
nComp=CLUS.NUM
Xhat.pca.aux=Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
Xhat.pca=scale(Xhat.pca.aux, center = -mu_X, scale = FALSE)

for (i in 1:CLUS.NUM){
  #Run cICA Algorithms
  c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
  c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
  c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
  #Extract ICs
  c1.y.df=cbind.data.frame(c1.y.df, as.numeric(c1[[1]]))
  c2.y.df=cbind.data.frame(c2.y.df, as.numeric(c2[[1]]))
  c3.y.df=cbind.data.frame(c3.y.df, as.numeric(c3[[1]]))
  #Extract W (only for cICA implementations)
  c1.w.df=cbind.data.frame(c1.w.df, c1[[2]])
  c2.w.df=cbind.data.frame(c2.w.df, c2[[2]])
  c3.w.df=cbind.data.frame(c3.w.df, c3[[2]])
  
    }
#Compute Reconstruction Xhat
WhitenMatrix=whiten_ZCA(X)[[2]]
#c1: Lu-Rajapakse
c1.y=as.matrix(c1.y.df[,-1])
c1.w=as.matrix(c1.w.df[,-1])
Xhat.c1.aux=t(ginv(t(c1.w))%*%t(c1.y))
Xhat.c1=t(inv(WhitenMatrix)%*%t(Xhat.c1.aux))

#c2: Zilling-Zhang
c2.y=as.matrix(c2.y.df[,-1])
c2.w=as.matrix(c2.w.df[,-1])
Xhat.c2.aux=t(ginv(t(c2.w))%*%t(c2.y))
Xhat.c2=t(inv(WhitenMatrix)%*%t(Xhat.c2.aux))

#c3: Lin's
c3.y=as.matrix(c3.y.df[,-1])
c3.w=as.matrix(c3.w.df[,-1])
Xhat.c3.aux=t(ginv(t(c3.w))%*%t(c3.y))
Xhat.c3=t(inv(WhitenMatrix)%*%t(Xhat.c3.aux))

####GENERATE MSE dataFrames####
#Reconstruction MSE
MSE.df=data.frame(round(mse(Xhat.pca, X),5), round(mse(Xhat.ica, whiten_ZCA(X)[[1]]),5), round(mse(Xhat.c1, whiten_ZCA(X)[[1]]),5), round(mse(Xhat.c2, whiten_ZCA(X)[[1]]),5), 
                  round(mse(Xhat.c3, whiten_ZCA(X)[[1]]),5))
colnames(MSE.df)=c("PCA", "ICA", "Lu-Rajapakse", "Zhang", "Lin's")
Total.MSE=rbind.data.frame(Total.MSE, colMeans(MSE.df))
colnames(Total.MSE)=c("PCA", "ICA", "Lu-Rajapakse", "Zhang", "Lin's")
}

####FIT ARIMA TO FEEDERS AND EVALUATE RMSE, MAE####
Result.M=data.frame()
for (j in 1:dim(X)[2]){
train.index=365-(31+30+31)
x.aux=ts(X[,j], start=1, end=365, frequency=24)
x.train=window(x.aux, start=1, end=train.index)
x.test=window(x.aux, start=train.index+1)
fit=auto.arima(x.train)
fcast=forecast::forecast(fit, h=(365-273))
Result.M=rbind.data.frame(Result.M,accuracy(fcast, x.test)[2,2:5])
}
colnames(Result.M)=c("RMSE", "MAE", "MPE", "MAPE")

#1. FIT ARIMA TO PC's
#Run Forecast
Forecast.M.pca=data.frame(seq(1:((31+30+31))*24))
train.index=(365-(31+30+31))*24
for (k in 1:CLUS.NUM){
x.aux=msts(as.ts(Xpca$x[,k]), seasonal.periods=c(12,24,7*24))
x.train=forecast:::subset.msts(x.aux, start=1, end=train.index)
x.test=forecast:::subset.msts(x.aux, start=train.index+1)
fit=auto.arima(x.train)
fcast=forecast::forecast(fit, h=((365-273)*24))
Forecast.M.pca=data.frame(Forecast.M.pca, fcast$mean)
}
#Reconstruct Xhat
Forecast.X.pca=as.matrix(data.frame(as.numeric(Forecast.M.pca[,2]), as.numeric(Forecast.M.pca[,3]), as.numeric(Forecast.M.pca[,4])))
Forecast.pca.aux=Forecast.X.pca %*% t(Xpca$rotation[,1:CLUS.NUM])
Forecast.Xhat.pca=scale(Forecast.pca.aux, center = colMeans(Forecast.pca.aux), scale = FALSE)
#Compute Performance Metrics
X.test=X[(train.index+1):dim(X)[1],]
Result.M.pca=data.frame()
for (k in 1:dim(X.test)[2]){
Result.M.pca=rbind.data.frame(Result.M.pca, forecast::accuracy(Forecast.Xhat.pca[,k], X.test[,k]))
}

#2. FIT ARIMA TO IC's from FASTICA
#Run Forecast
Forecast.M.fastica=data.frame(seq(1:((31+30+31))*24))
train.index=(365-(31+30+31))*24
for (k in 1:CLUS.NUM){
  x.aux=msts(as.ts(ic$S[,k]), seasonal.periods=c(12,24,7*24))
  x.train=forecast:::subset.msts(x.aux, start=1, end=train.index)
  x.test=forecast:::subset.msts(x.aux, start=train.index+1)
  fit=auto.arima(x.train)
  fcast=forecast::forecast(fit, h=((365-273)*24))
  Forecast.M.fastica=data.frame(Forecast.M.fastica, fcast$mean)
}
#Reconstruct Xhat
Forecast.X.fastica=as.matrix(data.frame(as.numeric(Forecast.M.fastica[,2]), as.numeric(Forecast.M.fastica[,3]), as.numeric(Forecast.M.fastica[,4])))
Forecast.Xhat.fastica=t(t(ic$A)%*%t(Forecast.X.fastica))
#Compute Performance Metrics
X.test=X[(train.index+1):dim(X)[1],]
Result.M.fastica=data.frame()
for (k in 1:dim(X.test)[2]){
  Result.M.fastica=rbind.data.frame(Result.M.fastica, forecast::accuracy(Forecast.Xhat.fastica[,k], X.test[,k]))
}

#3. FIT ARIMA TO IC's from Lu-Rajapakse's
#Run Forecast
Forecast.M.c1=data.frame(seq(1:((31+30+31))*24))
train.index=(365-(31+30+31))*24
for (k in 1:CLUS.NUM){
  x.aux=msts(as.ts(c1.y[,k]), seasonal.periods=c(12,24,7*24))
  x.train=forecast:::subset.msts(x.aux, start=1, end=train.index)
  x.test=forecast:::subset.msts(x.aux, start=train.index+1)
  fit=auto.arima(x.train)
  fcast=forecast::forecast(fit, h=((365-273)*24))
  Forecast.M.c1=data.frame(Forecast.M.c1, fcast$mean)
}
#Reconstruct Xhat
Forecast.X.c1=as.matrix(data.frame(as.numeric(Forecast.M.c1[,2]), as.numeric(Forecast.M.c1[,3]), as.numeric(Forecast.M.c1[,4])))
Forecast.Xhat.c1=t(ginv(t(c1.w))%*%t(Forecast.X.c1))
#Compute Performance Metrics
X.whiten=whiten_ZCA(X)[[1]]
X.test=X.whiten[(train.index+1):dim(X)[1],]
Result.M.c1=data.frame()
for (k in 1:dim(X.test)[2]){
  Result.M.c1=rbind.data.frame(Result.M.c1, forecast::accuracy(Forecast.Xhat.c1[,k], X.test[,k]))
}

#4. FIT ARIMA TO IC's from Zhang
#Run Forecast
Forecast.M.c2=data.frame(seq(1:((31+30+31))*24))
train.index=(365-(31+30+31))*24
for (k in 1:CLUS.NUM){
  x.aux=msts(as.ts(c2.y[,k]), seasonal.periods=c(12,24,7*24))
  x.train=forecast:::subset.msts(x.aux, start=1, end=train.index)
  x.test=forecast:::subset.msts(x.aux, start=train.index+1)
  fit=auto.arima(x.train)
  fcast=forecast::forecast(fit, h=((365-273)*24))
  Forecast.M.c2=data.frame(Forecast.M.c2, fcast$mean)
}
#Reconstruct Xhat
Forecast.X.c2=as.matrix(data.frame(as.numeric(Forecast.M.c2[,2]), as.numeric(Forecast.M.c2[,3]), as.numeric(Forecast.M.c2[,4])))
Forecast.Xhat.c2=t(ginv(t(c2.w))%*%t(Forecast.X.c2))
#Compute Performance Metrics
X.whiten=whiten_ZCA(X)[[1]]
X.test=X.whiten[(train.index+1):dim(X)[1],]
Result.M.c2=data.frame()
for (k in 1:dim(X.test)[2]){
  Result.M.c2=rbind.data.frame(Result.M.c2, forecast::accuracy(Forecast.Xhat.c2[,k], X.test[,k]))
}

#5. FIT ARIMA TO IC's from Lin's
#Run Forecast
Forecast.M.c3=data.frame(seq(1:((31+30+31))*24))
train.index=(365-(31+30+31))*24
for (k in 1:CLUS.NUM){
  x.aux=msts(as.ts(c3.y[,k]), seasonal.periods=c(12,24,7*24))
  x.train=forecast:::subset.msts(x.aux, start=1, end=train.index)
  x.test=forecast:::subset.msts(x.aux, start=train.index+1)
  fit=auto.arima(x.train)
  fcast=forecast::forecast(fit, h=((365-273)*24))
  Forecast.M.c3=data.frame(Forecast.M.c3, fcast$mean)
}
#Reconstruct Xhat
Forecast.X.c3=as.matrix(data.frame(as.numeric(Forecast.M.c3[,2]), as.numeric(Forecast.M.c3[,3]), as.numeric(Forecast.M.c3[,4])))
Forecast.Xhat.c3=t(ginv(t(c3.w))%*%t(Forecast.X.c3))
#Compute Performance Metrics
X.whiten=whiten_ZCA(X)[[1]]
X.test=X.whiten[(train.index+1):dim(X)[1],]
Result.M.c3=data.frame()
for (k in 1:dim(X.test)[2]){
  Result.M.c3=rbind.data.frame(Result.M.c3, forecast::accuracy(Forecast.Xhat.c3[,k], X.test[,k]))
}

#Relative RMSE and MAE#
#RMSE RESULT TABLE
# define a function to remove outliers
FindOutliers <- function(data) {
  lowerq = quantile(data)[2]
  upperq = quantile(data)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  # we identify extreme outliers
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  result <- which(data > extreme.threshold.upper | data < extreme.threshold.lower)
}


Rel.RMSE=data.frame((Result.M.pca%>%select(RMSE)/Result.M%>%select(RMSE)), (Result.M.fastica%>%select(RMSE)/Result.M%>%select(RMSE)),
                    (Result.M.c1%>%select(RMSE)/Result.M%>%select(RMSE)), (Result.M.c2%>%select(RMSE)/Result.M%>%select(RMSE)),
                    (Result.M.c3%>%select(RMSE)/Result.M%>%select(RMSE)))
colnames(Rel.RMSE)=c("pca", "fastICA", "Lu-Rajapakse", "Zhang", "Lin's")
#Deal with infinites and outliers
Rel.RMSE.aux=Rel.RMSE[!is.infinite(rowSums(Rel.RMSE)),]
temp=FindOutliers(Rel.RMSE.aux[,1])
Rel.RMSE.final=Rel.RMSE.aux[-temp,]

#PCA vs the Rest
t.test(Rel.RMSE.final[,1], Rel.RMSE.final[,2], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.RMSE.final[,1], Rel.RMSE.final[,3], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.RMSE.final[,1], Rel.RMSE.final[,4], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.RMSE.final[,1], Rel.RMSE.final[,5], alternative=c("two.sided"), conf.level=0.95)

#fastICA vs Lu-Rajapakse, Zhang & Lin
t.test(Rel.RMSE.final[,2], Rel.RMSE.final[,3], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.RMSE.final[,2], Rel.RMSE.final[,4], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.RMSE.final[,2], Rel.RMSE.final[,5], alternative=c("two.sided"), conf.level=0.95)

#Lu-Rajapakse vs Zhang, Lin
t.test(Rel.RMSE.final[,3], Rel.RMSE.final[,4], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.RMSE.final[,3], Rel.RMSE.final[,5], alternative=c("two.sided"), conf.level=0.95)

#Zhang vs Lin
t.test(Rel.RMSE.final[,4], Rel.RMSE.final[,5], alternative=c("two.sided"), conf.level=0.95)

#ColMeans of the RMSE
colMeans(Rel.RMSE.final)

#MAE
Rel.MAE=data.frame((Result.M.pca%>%select(MAE)/Result.M%>%select(MAE)), (Result.M.fastica%>%select(MAE)/Result.M%>%select(MAE)),
                    (Result.M.c1%>%select(MAE)/Result.M%>%select(MAE)), (Result.M.c2%>%select(MAE)/Result.M%>%select(MAE)),
                    (Result.M.c3%>%select(MAE)/Result.M%>%select(MAE)))
colnames(Rel.MAE)=c("pca", "fastICA", "Lu-Rajapakse", "Zhang", "Lin's")
#Deal with infinites and outliers
Rel.MAE.aux=Rel.MAE[!is.infinite(rowSums(Rel.MAE)),]
temp=FindOutliers(Rel.MAE.aux[,1])
Rel.MAE.final=Rel.MAE.aux[-temp,]

#PCA vs the Rest
t.test(Rel.MAE.final[,1], Rel.MAE.final[,2], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.MAE.final[,1], Rel.MAE.final[,3], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.MAE.final[,1], Rel.MAE.final[,4], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.MAE.final[,1], Rel.MAE.final[,5], alternative=c("two.sided"), conf.level=0.95)

#fastICA vs Lu-Rajapakse, Zhang & Lin
t.test(Rel.MAE.final[,2], Rel.MAE.final[,3], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.MAE.final[,2], Rel.MAE.final[,4], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.MAE.final[,2], Rel.MAE.final[,5], alternative=c("two.sided"), conf.level=0.95)

#Lu-Rajapakse vs Zhang, Lin
t.test(Rel.MAE.final[,3], Rel.MAE.final[,4], alternative=c("two.sided"), conf.level=0.95)
t.test(Rel.MAE.final[,3], Rel.MAE.final[,5], alternative=c("two.sided"), conf.level=0.95)

#Zhang vs Lin
t.test(Rel.MAE.final[,4], Rel.MAE.final[,5], alternative=c("two.sided"), conf.level=0.95)

#T-Test for RMSE and MAE



#Diebold-Mariano Test for RMSE and MAE
#baseline
rmse.base=Result.M%>%select(RMSE)
mae.base=Result.M%>%select(MAE)
#pca
rmse.pca=Result.M.pca%>%select(RMSE)
mae.pca=Result.M.pca%>%select(MAE)
#fastica
rmse.fastica=Result.M.fastica%>%select(RMSE)
mae.fastica=Result.M.fastica%>%select(MAE)
#Lu-Rajapakse cica
rmse.c1=Result.M.c1%>%select(RMSE)
mae.c1=Result.M.c1%>%select(MAE)
#Zhang
rmse.c2=Result.M.c2%>%select(RMSE)
mae.c2=Result.M.c2%>%select(MAE)
#Lin
rmse.c3=Result.M.c3%>%select(RMSE)
mae.c3=Result.M.c3%>%select(MAE)

#####RMSE#####
#c1 vs baseline
round(mean(rmse.c1[,1]-rmse.base[,1]),3)
round(sd(rmse.c1[,1]-rmse.base[,1]),3)
dm.test(rmse.c1[,1], rmse.base[,1], alternative = c("two.sided"), h = 1, power = 2)
#c1 vs pca
round(mean(rmse.c1[,1]-rmse.pca[,1]),3)
round(sd(rmse.c1[,1]-rmse.pca[,1]),3)
dm.test(rmse.c1[,1], rmse.pca[,1], alternative = c("two.sided"), h = 1, power = 2)
#c1 vs fastICA
round(mean(rmse.c1[,1]-rmse.fastica[,1]),3)
round(sd(rmse.c1[,1]-rmse.fastica[,1]),3)
dm.test(rmse.c1[,1], rmse.fastica[,1], alternative = c("two.sided"), h = 1, power = 2)

#c2 vs baseline
round(mean(rmse.c2[,1]-rmse.base[,1]),3)
round(sd(rmse.c2[,1]-rmse.base[,1]),3)
dm.test(rmse.c2[,1], rmse.base[,1], alternative = c("two.sided"), h = 1, power = 2)
#c2 vs pca
round(mean(rmse.c2[,1]-rmse.pca[,1]),3)
round(sd(rmse.c2[,1]-rmse.pca[,1]),3)
dm.test(rmse.c2[,1], rmse.pca[,1], alternative = c("two.sided"), h = 1, power = 2)
#c2 vs fastICA
round(mean(rmse.c2[,1]-rmse.fastica[,1]),3)
round(sd(rmse.c2[,1]-rmse.fastica[,1]),3)
dm.test(rmse.c2[,1], rmse.fastica[,1], alternative = c("two.sided"), h = 1, power = 2)


#c3 vs baseline
round(mean(rmse.c3[,1]-rmse.base[,1]),3)
round(sd(rmse.c3[,1]-rmse.base[,1]),3)
dm.test(rmse.c3[,1], rmse.base[,1], alternative = c("two.sided"), h = 1, power = 2)
#c3 vs pca
round(mean(rmse.c3[,1]-rmse.pca[,1]),3)
round(sd(rmse.c3[,1]-rmse.pca[,1]),3)
dm.test(rmse.c3[,1], rmse.pca[,1], alternative = c("two.sided"), h = 1, power = 2)
#c3 vs fastICA
round(mean(rmse.c3[,1]-rmse.fastica[,1]),3)
round(sd(rmse.c3[,1]-rmse.fastica[,1]),3)
dm.test(rmse.c3[,1], rmse.fastica[,1], alternative = c("two.sided"), h = 1, power = 2)

#c1 vs c2
round(mean(rmse.c1[,1]-rmse.c2[,1]),3)
round(sd(rmse.c1[,1]-rmse.c2[,1]),3)
dm.test(rmse.c1[,1], rmse.c2[,1], alternative = c("two.sided"), h = 1, power = 2)
#c1 vs c3
round(mean(rmse.c1[,1]-rmse.c3[,1]),3)
round(sd(rmse.c1[,1]-rmse.c3[,1]),3)
dm.test(rmse.c1[,1], rmse.c3[,1], alternative = c("two.sided"), h = 1, power = 2)
#c3 vs fastICA
round(mean(rmse.c2[,1]-rmse.c3[,1]),3)
round(sd(rmse.c2[,1]-rmse.c3[,1]),3)
dm.test(rmse.c2[,1], rmse.c3[,1], alternative = c("two.sided"), h = 1, power = 2)


#####MAE#####
#c1 vs baseline
round(mean(mae.c1[,1]-mae.base[,1]),3)
round(sd(mae.c1[,1]-mae.base[,1]),3)
dm.test(mae.c1[,1], mae.base[,1], alternative = c("two.sided"), h = 1, power = 1)
#c1 vs pca
round(mean(mae.c1[,1]-mae.pca[,1]),3)
round(sd(mae.c1[,1]-mae.pca[,1]),3)
dm.test(mae.c1[,1], mae.pca[,1], alternative = c("two.sided"), h = 1, power = 1)
#c1 vs fastICA
round(mean(mae.c1[,1]-mae.fastica[,1]),3)
round(sd(mae.c1[,1]-mae.fastica[,1]),3)
dm.test(mae.c1[,1], mae.fastica[,1], alternative = c("two.sided"), h = 1, power = 1)

#c2 vs baseline
round(mean(mae.c2[,1]-mae.base[,1]),3)
round(sd(mae.c2[,1]-mae.base[,1]),3)
dm.test(mae.c2[,1], mae.base[,1], alternative = c("two.sided"), h = 1, power = 1)
#c2 vs pca
round(mean(mae.c2[,1]-mae.pca[,1]),3)
round(sd(mae.c2[,1]-mae.pca[,1]),3)
dm.test(mae.c2[,1], mae.pca[,1], alternative = c("two.sided"), h = 1, power = 1)
#c2 vs fastICA
round(mean(mae.c2[,1]-mae.fastica[,1]),3)
round(sd(mae.c2[,1]-mae.fastica[,1]),3)
dm.test(mae.c2[,1], mae.fastica[,1], alternative = c("two.sided"), h = 1, power = 1)


#c3 vs baseline
round(mean(mae.c3[,1]-mae.base[,1]),3)
round(sd(mae.c3[,1]-mae.base[,1]),3)
dm.test(mae.c3[,1], mae.base[,1], alternative = c("two.sided"), h = 1, power = 1)
#c3 vs pca
round(mean(mae.c3[,1]-mae.pca[,1]),3)
round(sd(mae.c3[,1]-mae.pca[,1]),3)
dm.test(mae.c3[,1], mae.pca[,1], alternative = c("two.sided"), h = 1, power = 1)
#c3 vs fastICA
round(mean(mae.c3[,1]-mae.fastica[,1]),3)
round(sd(mae.c3[,1]-mae.fastica[,1]),3)
dm.test(mae.c3[,1], mae.fastica[,1], alternative = c("two.sided"), h = 1, power = 1)

#c1 vs c2
round(mean(mae.c1[,1]-mae.c2[,1]),3)
round(sd(mae.c1[,1]-mae.c2[,1]),3)
dm.test(mae.c1[,1], mae.c2[,1], alternative = c("two.sided"), h = 1, power = 1)
#c1 vs c3
round(mean(mae.c1[,1]-mae.c3[,1]),3)
round(sd(mae.c1[,1]-mae.c3[,1]),3)
dm.test(mae.c1[,1], mae.c3[,1], alternative = c("two.sided"), h = 1, power = 1)
#c3 vs fastICA
round(mean(mae.c2[,1]-mae.c3[,1]),3)
round(sd(mae.c2[,1]-mae.c3[,1]),3)
dm.test(mae.c2[,1], mae.c3[,1], alternative = c("two.sided"), h = 1, power = 1)


###PLOTS###
ref.1=data.frame(Sim.Feeders[,1], ref[,1], "r1")
colnames(ref.1)=c("ts", "Signal", "ref")
ref.2=data.frame(Sim.Feeders[,1], ref[,2], "r2")
colnames(ref.2)=c("ts", "Signal", "ref")
ref.3=data.frame(Sim.Feeders[,1], ref[,3], "r3")
colnames(ref.3)=c("ts", "Signal", "ref")
ref.df=rbind.data.frame(ref.1, ref.2, ref.3)
colnames(ref.df)=c("ts","Signal", "Reference")

ref.df %>%
  ggplot(aes(x = ts, y = Signal, group = Reference)) +
  geom_line() +
  facet_grid(~Reference)



ggplot(ref.df, aes(x = ts, y = Signal, group = 1)) + 
  geom_point(aes(color = Reference), size = 1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_minimal()
  


#####################################TEST 6 HIERARCHICAL SIMULATED TRAIN/TEST SPLIT#####################################T
#Test matrixes Simulation 3: Random Frequencies and White Noise#
library("fastICA", lib.loc="~/R/win-library/3.5")

#Load Simulated Data Sets
Sim.Feeders=data.frame(read.csv("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulated_Feeder.csv", header = TRUE, sep = ","))
Sim.SSEE=data.frame(read.csv("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulated_SSEE.csv", header = TRUE, sep = ","))
Sim.Feeders=Sim.Feeders[,2:dim(Sim.Feeders)[2]]
Sim.SSEE=Sim.SSEE[,2:dim(Sim.SSEE)[2]]


#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1

#Train/Test Split ON ROWS! (Time dimension)
#2. BUILD BLIND SOURCE SETTING AND RESULTS DATA FRAMES
X.feeders=as.matrix(Sim.Feeders[,2:dim(Sim.Feeders)[2]])
week.num=dim(X.feeders)[1]/(24*7)
week.end.train=(week.num)*0.6*(24*7)

X.feeders.train=X.feeders[1:(week.end.train),]
X.feeders.test=X.feeders[((week.end.train)+1):dim(X.feeders)[1],]




###1. REFERENCE SIGNAL EXTRACTION FROM FEEDERS
library("dplyr", lib.loc="~/R/win-library/3.5")
library("tsoutliers", lib.loc="~/R/win-library/3.5")
library("TSclust", lib.loc="~/R/win-library/3.5")
library("factoextra", lib.loc="~/R/win-library/3.5")
library("dendextend", lib.loc="~/R/win-library/3.5")
library("cluster", lib.loc="C:/Program Files/R/R-3.5.1/library")
library("NbClust", lib.loc="~/R/win-library/3.5")
library("dtwclust", lib.loc="~/R/win-library/3.5")
library("hydroGOF", lib.loc="~/R/win-library/3.5")
#1.1 SELECT NUMBER OF CLUSTERS
#DoParallel routine
require("doParallel")
# Create parallel workers
workers <- makeCluster(6L)
# Preload dtwclust in each worker; not necessary but useful
invisible(clusterEvalQ(workers, library("factoextra")))
# Register the backend; this step MUST be done
registerDoParallel(workers)
set.seed(12345)
AvgSilhouette=fviz_nbclust(X.feeders.train, FUN = hcut, hc_func=c("hclust"), hc_method=c("complete"), method = "silhouette")
#gap_stat=clusGap(Final.Data.Set, FUN = hcut, hc_func=c("hclust"), hc_method=c("complete"), hc_metric = "pearson",
#                 K.max = round(0.05*dim(Output.Period.Outliers)[2]), B = 50, method="globalmax")
# Stop parallel workers
stopCluster(workers)
# Go back to sequential computation
registerDoSEQ()
#3. Extract Maximum Silhouette
Cluster.Sil=AvgSilhouette$data
CLUS.NUM=as.numeric(Cluster.Sil[which.max(Cluster.Sil$y),1])

#1.2 RUN AGGLOMERATIVE HCLUST COMPLETE LINKAGE
dist_ts=TSclust::diss(SERIES = t(X.feeders.train), METHOD = "INT.PER") # note the dataframe must be transposed
hcPER.complete=hclust(dist_ts, method="complete")
plot(hcPER.complete)
rect.hclust(hcPER.complete, k=CLUS.NUM, border=2:12) #Change Clus.Num
sub_grp.complete=cutree(hcPER.complete, k = CLUS.NUM) #Change Clus.Num
t=data.frame(sub_grp.complete)
Output.HClust=data.frame(row.names(t), t, row.names=NULL)
names(Output.HClust)=c("Alim_ID", "Cluster")

#1.3 Compute Cluster Centroids
#Extract Memberships & Centroids
Clust.Center=list()
for (i in (1:CLUS.NUM)){
  aux1=as.character(t(Output.HClust%>%dplyr::filter(Cluster==i)%>%dplyr::select(Alim_ID)))
  idx=match(aux1, names(Sim.Feeders))
  aux2=data.frame(Sim.Feeders[1:(week.end.train),2],Sim.Feeders[1:(week.end.train), idx])
  Clust.Center[[i]]=base::rowMeans(aux2[,2:dim(aux2)[2]])
}

#1.4 Reference Signals
ref=data.frame(Clust.Center[[1]])
for (i in 2:CLUS.NUM){
  ref=cbind.data.frame(ref,Clust.Center[[i]])
}

#2. BUILD BLIND SOURCE SETTING AND RESULTS DATA FRAMES
X.ssee=as.matrix(Sim.SSEE[,2:dim(Sim.SSEE)[2]])
X.ssee.train=X.ssee[1:(week.end.train),]
X.ssee.test=X.ssee[((week.end.train)+1):dim(X.ssee)[1],]



#2.1 Empty Data Frames
#Y Empty Data Frames
c1.y.df=data.frame(Sim.SSEE[1:(week.end.train),1])
c2.y.df=data.frame(Sim.SSEE[1:(week.end.train),1])
c3.y.df=data.frame(Sim.SSEE[1:(week.end.train),1])
fastICA.df=data.frame(Sim.SSEE[1:(week.end.train),1])

#W Empty Data Frames
c1.w.df=data.frame(colnames(X.ssee.train))
c2.w.df=data.frame(colnames(X.ssee.train))
c3.w.df=data.frame(colnames(X.ssee.train))

#3. RUN ICA ALGORITMS IN TRAINING SET
#Benchmarks
#Run fastICA
ic=fastICA(X.ssee.train, CLUS.NUM, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
Xhat.ica=t(t(ic$A)%*%t(ic$S))

#PCA
Xpca=prcomp(X.ssee.train)
nComp=CLUS.NUM
Xhat.pca.aux=Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
Xhat.pca=scale(Xhat.pca.aux, center = -colMeans(X.ssee.train), scale = FALSE)

#Iterate over all extracted references
for (i in 1:CLUS.NUM){
  #Run cICA Algorithms
  c1=ICA_R(X.ssee.train, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
  c2=Zhang_cICA(X.ssee.train, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
  c3=Lin_fastICAR(X.ssee.train, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
  #Extract ICs
  c1.y.df=cbind.data.frame(c1.y.df, as.numeric(c1[[1]]))
  c2.y.df=cbind.data.frame(c2.y.df, as.numeric(c2[[1]]))
  c3.y.df=cbind.data.frame(c3.y.df, as.numeric(c3[[1]]))
  #Extract W (only for cICA implementations)
  c1.w.df=cbind.data.frame(c1.w.df, c1[[2]])
  c2.w.df=cbind.data.frame(c2.w.df, c2[[2]])
  c3.w.df=cbind.data.frame(c3.w.df, c3[[2]])
  
}
#Compute Reconstruction Xhat
WhitenMatrix=whiten_PCA(X.ssee.train)[[2]]
#c1: Lu-Rajapakse
c1.y=as.matrix(c1.y.df[,-1])
c1.w=as.matrix(c1.w.df[,-1])
Xhat.c1.aux=t(ginv(t(c1.w))%*%t(c1.y))
Xhat.c1=t(inv(WhitenMatrix)%*%t(Xhat.c1.aux))

#c2: Zilling-Zhang
c2.y=as.matrix(c2.y.df[,-1])
c2.w=as.matrix(c2.w.df[,-1])
Xhat.c2.aux=t(ginv(t(c2.w))%*%t(c2.y))
Xhat.c2=t(inv(WhitenMatrix)%*%t(Xhat.c2.aux))

#c3: Lin's
c3.y=as.matrix(c3.y.df[,-1])
c3.w=as.matrix(c3.w.df[,-1])
Xhat.c3.aux=t(ginv(t(c3.w))%*%t(c3.y))
Xhat.c3=t(inv(WhitenMatrix)%*%t(Xhat.c3.aux))

####GENERATE MSE dataFrames####
MSE.df=data.frame(round(mse(Xhat.pca,X.ssee.train),5), round(mse(Xhat.ica,X.ssee.train),5), round(mse(Xhat.c1,X.ssee.train),5), 
                  round(mse(Xhat.c2,X.ssee.train),5), round(mse(Xhat.c3,X.ssee.train),5))
colnames(MSE.df)=c("PCA", "ICA", "Lu-Rajapakse", "Zhang", "Lin's")
Total.MSE=colMeans(MSE.df)



#TO DO
#1. Add Noice to the reference signals
#2. In case i dont know anything use the three classical load profiles curves. (generate functions with the approximate shape of the load curves)
#3. Pure signals of a particular type (data-driven approach). T°
#4. Data-Driven vs Domain Knowledge for the load profiles extraction (review domain driven load profile construction).
#5. review kernel-PCA (intractable with data growth)
#6. time aspect (how frequently are you gonna run these models)
#7. Common for walltime, number of flops IDEA: See how the algorithms scale 
#8. See what if's scenarios.
#9. Simulate the hierarchical data sets (how to specify the aggregation assume indepence and simulate)
#10. Dream-up arbitrarily scaleable situation (M references, N hidden signals)
#11. Construct the R package
#12. generate quantiles curves to see the CI for the simulations (Use commom random numbers)


