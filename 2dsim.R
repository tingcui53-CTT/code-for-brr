#simulation
#2d
computer_simulator<-function(x){
  m=nrow(x)
  x1=x[,1]
  x2=x[,2]
  t1=x1^2
  t2=x2^2
  j1<-rep(1,m)
  j2<-rep(2,m)
  e<-rnorm(m,0,0.01)
  y=j1*t1+j2*t2+e
  return(y)
}
computer_simulator<-function(x){
  m=nrow(x)
  x1=x[,1]
  x2=x[,2]
  e<-rnorm(m,0,0.01)
  y=sin(8*x1+x2)+e
  return(y)
}
#simulation function
computer_simulator<-function(x){
  m=nrow(x)
  x1=x[,1]
  x2=x[,2]
  t1=x1^2
  t2=x2^2
  e<-rnorm(m,0,0.01)
  y=-sum(x1+x2)*exp(-sum(t1+t2))+e
  return(y)
}
computer_simulator<-function(x){
  m=nrow(x)
  x1=x[,1]
  x2=x[,2]
  t1=x1^2
  t2=x2^2
  e<-rnorm(m,0,0.01)
  y=-20*exp(-sqrt((t1+t2)/2)/5)-exp(2*pi*(x1+x2)/2)+20+e
  return(y)
}
#informative prior
#sampling
d=2
m=10
RMSE1<-array(0,2)
RMSE1p<-array(0,2)
for(t in 1:2){
  set.seed(t)
  ll=5000
  NT=20000
  u<-runif(NT+ll,0,1)
  v<-runif(NT+ll,0,1)
  xx<-cbind(u,v)
  yy <-computer_simulator(xx)
  plot(yy)
  #training dataset
  utrain<-u[1:ll]
  vtrain<-v[1:ll]
  xtrain<-xx[1:ll,]
  ytrain<-yy[1:ll]
  #test dataset
  utest<-u[(ll+1):(NT+ll)]
  vtest<-v[(ll+1):(NT+ll)]
  xtest<-xx[(ll+1):(NT+ll),]
  ytest<-yy[(ll+1):(NT+ll)]
  R<-function(theta,h){
    RH<-sum(-h^2%*%theta)
    R<-exp(RH)
    return(R)
  }
  #choose m knots from n
  CA<-function(X){
    for(k in 1:ll){
      C0=ll
      d<-2
      n<-ll
      m=max(5*d,log(n))
      t<-sample(1:n,m)
      U<-X[t]
      l=1
      c<-NULL
      for(i in 1:(m-1)){
        for(j in (i+1):m){
          c[l] <- 1/abs(U[i]-U[j])
          l <- l+1
        }
        
      }
      CA<-max(c)
      if(!is.na(CA)) {
        if( CA <= C0){
          X<-U
          C0<-CA
          t0<-t
          
        }
      }
      
    }
    return(t0)
  }
  m0<-CA(xtrain)
  # n=200
  # m=80
  # m0<-sample(1:n,m)
  #interpolation knots
  x1<-xtrain[m0,]
  Ua<-x1[,1]
  Va<-x1[,2]
  y1<-ytrain[m0]
  #reconstruction paramaterization
  RA<- function(x,theta){
    m<-nrow(x)
    Ua<-x[,1]
    Va<-x[,2]
    RA <- matrix(0,m,m)
    for (i in 1:m) {
      for(j in 1:m){
        RA[i,j]= R(theta,c(Ua[i]-Ua[j],Va[i]-Va[j]))
      }
    }
    return(RA)
  }
  ra <- function(x,X,theta){
    m<-nrow(x)
    Ua<-x[,1]
    Va<-x[,2]
    n<-nrow(X)
    U<-X[,1]
    V<-X[,2]
    ra<-matrix(0,m,n)
    for(i in 1:m){
      h <- cbind(U-Ua[i],V-Va[i])
      for(j in 1:n){
        ra[i,j]=R(theta,h[j,])
      }
    }
    return(ra)
  }
  b<-function(x,X,theta){
    RAA<-RA(x,theta)
    rA<-ra(x,X,theta)
    m<-nrow(x)
    Ua<-x[,1]
    Va<-x[,2]
    n<-nrow(X)
    U<-X[,1]
    V<-X[,2]
    GA<-cbind(rep(1,m),Ua,Va)
    I<-rep(1,m)
    RAT<-solve(RAA+1e-8*diag(m))
    GRA<-solve(t(GA)%*%RAT%*%GA+1e-8*diag(3))
    P<-RAT%*%GA%*%GRA
    Q<-(diag(I)-RAT%*%GA%*%GRA%*%t(GA))%*%RAT
    b <- matrix(0,m,n)
    for (i in 1:n){
      b[,i]<-P%*%c(1,U[i],V[i])+Q%*%rA[,i]
    }
    return(b)
  }
  b<-function(x,X,theta){
    RAA<-RA(x,theta)
    rA<-ra(x,X,theta)
    m<-nrow(x)
    Ua<-x[,1]
    Va<-x[,2]
    n<-nrow(X)
    U<-X[,1]
    V<-X[,2]
    GA<-cbind(rep(1,m),Ua,Va)
    I<-rep(1,m)
    RAT<-solve(RAA+1e-8*diag(m))
    GRA<-solve(t(GA)%*%RAT%*%GA+1e-8*diag(3))
    P<-RAT%*%GA%*%GRA
    Q<-(diag(I)-RAT%*%GA%*%GRA%*%t(GA))%*%RAT
    b <- matrix(0,m,n)
    for (i in 1:n){
      b[,i]<-P%*%c(1,U[i],V[i])+Q%*%rA[,i]
    }
    return(b)
  }
  Logf<-function(par){
    m=nrow(x1)
    d=ncol(x1)
    n=length(ytrain)
    sigma2<-par[1]
    theta<-par[c(2:(1+d))]
    B<-b(x1,xtrain,theta)
    mum<-solve(B%*%t(B)+1e-8*diag(m))%*%B%*%ytrain
    value<-(-n/2*log(sigma2)-((ytrain-t(mum)%*%B)%*%t(ytrain-t(mum)%*%B))/(2*sigma2))
    return(-value)
  }
  optLogf<-optim(c(0.1,rep(1,d)),Logf,method ="L-BFGS-B",lower=c(1e-8,rep(1e-8,d)))
  Vgamma=RA(x1,optLogf$par[c(2:(1+d))])
  SVgamma=solve(Vgamma+1e-8*diag(m))
  B<-b(x1,xtrain,optLogf$par[c(2:(1+d))])
  mum<-solve(B%*%t(B)+1e-8*diag(m))%*%B%*%ytrain
  Sigma1=solve(B%*%t(B)+SVgamma)
  mugamma1=Sigma1%*%(B%*%ytrain+SVgamma%*%mum)
  SSigma1=B%*%t(B)+SVgamma
  Q1=t(ytrain)%*%ytrain+t(mum)%*%SVgamma%*%mum-t(mugamma1)%*%SSigma1%*%mugamma1
  ytestpred<-function(x){
    x=matrix(xtest[i,],ncol=2)
    bstar=b(x1,x,optLogf$par[c(2:(1+d))])
    ystar=computer_simulator(x)
    Sigma2=solve(SVgamma+bstar%*%t(bstar)+B%*%t(B))
    mugamma2=Sigma2%*%(SVgamma%*%mum+bstar%*%ystar++B%*%ytrain)
    SSigma2=SVgamma+bstar%*%t(bstar)+B%*%t(B)
    omega1=(1-t(bstar)%*%Sigma2%*%bstar)^(-1)
    mugamma3=omega1*t(mugamma1)%*%SSigma1%*%Sigma2%*%bstar
    Q3=t(ytrain)%*%ytrain+t(mum)%*%SVgamma%*%mum-t(mugamma1)%*%SSigma1%*%Sigma2%*%SSigma1%*%mugamma1-omega1^(-1)*mugamma3^2
    v1=omega1*Q3
    return(c(mugamma3,v1))
  }
  N=20000
  yt<-array(0,N)
  varp<-array(0,N)
  for(i in 1 :N){
    yt[i]=ytestpred(matrix(xtest[i,],ncol=2))[1]
    varp[i]=ytestpred(matrix(xtest[i,],ncol=2))[2]
  }
  yt<-as.numeric(yt)
  re1<-(ytest-yt)^2
  RMSE1[t]<-sqrt(mean(re1))
  ytp<-t(mugamma1)%*%b(x1,xtest,optLogf$par[c(2:(1+d))])
  RMSE1p[t]<-sqrt(mean((ytest-ytp)^2))
}
MRMSE1<-mean(RMSE1)
MRMSE1p<-mean(RMSE1p)
v1<-var(yt)
#noninformative prior
RMSE2<-array(0,2)
RMSE2p<-array(0,2)
for(t in 1:2){
  set.seed(t)
Logff<-function(par){
  d=ncol(x1)
  n=length(ytrain)
  sigma2<-par[1]
  theta<-par[c(2:(1+d))]
  B<-b(x1,xtrain,theta)
  gammah=solve(B%*%t(B)+1e-8*diag(m))%*%B%*%ytrain
  value<-(-n/2*log(sigma2)-((ytrain-t(gammah)%*%B)%*%t(ytrain-t(gammah)%*%B))/(2*sigma2))
  return(-value)
}
optLogff<-optim(c(0.1,rep(1,d)),Logff,method ="L-BFGS-B",lower=c(1e-8,rep(1e-8,d)))
B<-b(x1,xtrain,optLogff$par[c(2:(1+d))])
Sigmagam=solve(B%*%t(B)+1e-8*diag(m))
SSigmagam=B%*%t(B)
gammahat=Sigmagam%*%B%*%ytrain
SS1=t(ytrain)%*%ytrain-t(gammahat)%*%SSigmagam%*%gammahat
ytestpred1<-function(x){
  x=matrix(xtest[i,],ncol=2)
  bstar=b(x1,x,optLogff$par[c(2:(1+d))])
  ystar=computer_simulator(x)
  Sigma2=solve(bstar%*%t(bstar)+B%*%t(B)+1e-8*diag(m))
  mugamma2=Sigma2%*%(bstar%*%ystar++B%*%ytrain)
  SSigma2=bstar%*%t(bstar)+B%*%t(B)
  omega1=(1-t(bstar)%*%Sigma2%*%bstar)^(-1)
  mugamma3=omega1*t(gammahat)%*%SSigmagam%*%Sigma2%*%bstar
  SS3=t(ytrain)%*%ytrain-t(gammahat)%*%SSigmagam%*%Sigma2%*%SSigmagam%*%gammahat-omega1^(-1)*mugamma3^2
  v2=omega1*SS3
  return(c(mugamma3,v2))
}
N=20000
yt1<-array(0,N)
varp2<-array(0,N)
for(i in 1 :N){
  yt1[i]=ytestpred1(matrix(xtest[i,],ncol=2))[1]
  varp2[i]=ytestpred1(matrix(xtest[i,],ncol=2))[2]
}
yt1<-as.numeric(yt1)
re2<-(ytest-yt1)^2
RMSE2[t]<-sqrt(mean(re2))
yt2<-t(gammahat)%*%b(x1,xtest,optLogff$par[c(2:(1+d))])
RMSE2p[t]<-sqrt(mean((ytest-yt2)^2))
}
MRMSE2<-mean(RMSE2)
MRMSE2p<-mean(RMSE2p)
v2=var(yt1)
#nystrom
#sampling
d=2
RMSEny<-array(0,2)
RMSEnyp<-array(0,2)
for(t in 1:10){
  set.seed(t)
  ll=5000
  NT=20000
  u<-runif(NT+ll,0,1)
  v<-runif(NT+ll,0,1)
  xx<-cbind(u,v)
  yy <-computer_simulator(xx)
  plot(yy)
  #training dataset
  utrain<-u[1:ll]
  vtrain<-v[1:ll]
  xtrain<-xx[1:ll,]
  ytrain<-yy[1:ll]
  #test dataset
  utest<-u[(ll+1):(NT+ll)]
  vtest<-v[(ll+1):(NT+ll)]
  xtest<-xx[(ll+1):(NT+ll),]
  ytest<-yy[(ll+1):(NT+ll)]
  R<-function(theta,h){
    RH<-sum(-h^2%*%theta)
    R<-exp(RH)
    return(R)
  }

  RA<- function(x,theta){
    m<-nrow(x)
    Ua<-x[,1]
    Va<-x[,2]
    RA <- matrix(0,m,m)
    for (i in 1:m) {
      for(j in 1:m){
        RA[i,j]= R(theta,c(Ua[i]-Ua[j],Va[i]-Va[j]))
      }
    }
    return(RA)
  }
  ra <- function(x,X,theta){
    m<-nrow(x)
    Ua<-x[,1]
    Va<-x[,2]
    n<-nrow(X)
    U<-X[,1]
    V<-X[,2]
    ra<-matrix(0,m,n)
    for(i in 1:m){
      h <- cbind(U-Ua[i],V-Va[i])
      for(j in 1:n){
        ra[i,j]=R(theta,h[j,])
      }
    }
    return(ra)
  }
  
  nystrom_kernel_approximation <- function(X,m,theta) {
    # X: 输入数据矩阵，每一行代表一个数据点
    # kernel_func: 核函数
    # m: 选取的子样本数
    n <- nrow(X) # 数据点的总数
    # 随机选择m个数据点作为子样本
    idx <- sample(1:n, m)
    C <- X[idx,]
    # 计算子矩阵
    K_cc <-RA(C, theta)
    # 计算子矩阵与原始数据点之间的核函数值
    K_ic <-ra(X, C,theta)
    #对子矩阵 svd 分解
    P=eigen(K_cc,symmetric = FALSE)$vectors
    K_ccsv<-P%*%diag(eigen(K_cc,symmetric = FALSE)$values)%*%t(P)
    # 计算近似的核矩阵
    K_approx <- K_ic %*% solve(K_ccsv+1e-8*diag(nrow(C))) %*% t(K_ic)
    return(K_approx)
  }
  #computation reduction nystrom 
  d=2
  k=5*d
  n<-nrow(xtrain)
  U<-xtrain[,1]
  V<-xtrain[,2]
  H<-cbind(rep(1,n),U,V)
  me<-sample(1:n,round(k))
  floggy<-function(par){
    n<-nrow(xtrain)
    d=ncol(xtrain)
    delta<-par[1]
    theta<-par[c(2:(1+d))]
    rmn<-ra(xtrain[me,],xtrain,theta)
    rmm<-ra(xtrain[me,],xtrain[me,],theta)
    rnm<-ra(xtrain,xtrain[me,],theta)
    nySigma<-rnm%*%solve(rmm+1e-8*diag(k))%*%rmn+delta*diag(n)
    SN<-(1/delta)*(diag(n)-rnm%*%solve(delta*rmm+rmn%*%rnm+1e-8*diag(k))%*%rmn)
    betap<-solve(t(H)%*%SN%*%H)%*%t(H)%*%SN%*%ytrain
    sigma2<-(t(ytrain-H%*%betap)%*%SN%*%(ytrain-H%*%betap))/n
    floggy=-(-n*log(sigma2)/2-0.5*log(det(nySigma))
             -((t(ytrain-H%*%betap)%*%SN%*%(ytrain-H%*%betap))/(2*sigma2)))
    return(floggy)
  }
  opthyperpama1<-optim(c(1,rep(1,d)),floggy,method ="L-BFGS-B",lower=c(1,rep(1e-8,d)) ,control = list(maxit=10)) 
  k=10*d
  me<-sample(1:n,round(k))
  rmn<-ra(xtrain[me,],xtrain,opthyperpama1$par[c(2:(1+d))])
  rmm<-ra(xtrain[me,],xtrain[me,],opthyperpama1$par[c(2:(1+d))])
  rnm<-ra(xtrain,xtrain[me,],opthyperpama1$par[c(2:(1+d))])
  nySigma<-rnm%*%solve(rmm+1e-8*diag(k))%*%rmn+opthyperpama1$par[2]*diag(n)
  SN<-(1/opthyperpama1$par[1])*(diag(n)-rnm%*%solve(opthyperpama1$par[1]*rmm+rmn%*%rnm+1e-8*diag(k))%*%rmn)
  betahat<-solve(t(H)%*%SN%*%H)%*%t(H)%*%SN%*%ytrain
  Hte<-cbind(rep(1,nrow(xtest)),xtest)
  Rtetr<-ra(xtrain,xtest,opthyperpama1$par[c(2:(1+d))])
  #Rtr<-RA(xtrain,opthyperpama1$par[c(3:8)])
  ytny<-Hte%*%betahat+t(Rtetr)%*%SN%*%(ytrain-H%*%betahat)
  reny<-(ytest-ytny)^2
  RMSEny[t]<-sqrt(mean(reny))
  RMSEnyp[t]<-sqrt(mean(ytest-Hte%*%betahat+t(Rtetr)%*%SN%*%(ytrain-H%*%betahat)))
}
MRMSEny<-mean(RMSEny)
MRMSEnyp<-mean(RMSEnyp)
vny=var(ytny)
#ny
rbf_kernel <- function(X1, X2, theta){
  D2 <- outer(rowSums(X1^2), rowSums(X2^2), "+") - 2*X1 %*% t(X2)
  exp(-theta * D2)
}
nystrom_features <- function(X, Xm, theta){
  K_mm <- rbf_kernel(Xm, Xm, theta)
  K_nm <- rbf_kernel(X, Xm, theta)
  
  eig <- eigen(K_mm, symmetric = TRUE)
  idx <- eig$values > 1e-8
  U <- eig$vectors[, idx, drop=FALSE]
  Dinv <- diag(1 / sqrt(eig$values[idx]))
  
  Phi <- K_nm %*% U %*% Dinv
  list(Phi = Phi, Xm = Xm)
}

loglik_nystrom <- function(par, X, y, H, Xm){
  delta <- exp(par[1])        # 保证正
  theta <- exp(par[2])
  
  ns <- nystrom_features(X, Xm, theta)
  Phi <- ns$Phi
  n <- nrow(X)
  
  # Woodbury inverse
  A <- diag(delta, n)
  B <- Phi
  Cinv <- diag(ncol(B))
  
  Sigma_inv <- (1/delta) * (
    diag(n) -
      B %*% solve(t(B)%*%B + delta*Cinv) %*% t(B)
  )
  
  # beta(theta, delta)
  beta_hat <- solve(
    t(H)%*%Sigma_inv%*%H,
    t(H)%*%Sigma_inv%*%y
  )
  
  r <- y - H %*% beta_hat
  
  # log determinant
  logdet <- n*log(delta) +
    determinant(diag(ncol(B)) + (1/delta)*t(B)%*%B,
                logarithm = TRUE)$modulus
  
  val <- -0.5 * (
    logdet +
      t(r) %*% Sigma_inv %*% r
  )
  
  return(-as.numeric(val))  # optim 最小化
}
fit_nystrom_MLE <- function(X, y, m){
  n <- nrow(X)
  
  Xm <- X[sample(1:n, m), ]
  H <- cbind(1, X)
  
  opt <- optim(
    par = c(log(0.1), log(1)),
    fn = loglik_nystrom,
    X = X, y = y, H = H, Xm = Xm,
    method = "L-BFGS-B"
  )
  
  list(
    delta = exp(opt$par[1]),
    theta = exp(opt$par[2]),
    Xm = Xm,
    H = H
  )
}

predict_nystrom <- function(fit, Xtr, ytr, Xte){
  delta <- fit$delta
  theta <- fit$theta
  Xm <- fit$Xm
  Htr <- fit$H
  Hte <- cbind(1, Xte)
  
  ns <- nystrom_features(Xtr, Xm, theta)
  Phi <- ns$Phi
  
  Sigma_inv <- (1/delta) * (
    diag(nrow(Xtr)) -
      Phi %*% solve(t(Phi)%*%Phi + delta*diag(ncol(Phi))) %*% t(Phi)
  )
  
  beta_hat <- solve(
    t(Htr)%*%Sigma_inv%*%Htr,
    t(Htr)%*%Sigma_inv%*%ytr
  )
  
  r <- ytr - Htr %*% beta_hat
  K_te_m <- rbf_kernel(Xte, Xm, theta)
  K_tr_m <- rbf_kernel(Xtr, Xm, theta)
  
  yhat <- Hte %*% beta_hat +
    K_te_m %*% solve(t(K_tr_m)%*%K_tr_m + delta*diag(ncol(K_tr_m))) %*%
    t(K_tr_m) %*% r
  
  as.vector(yhat)
}
RMSE_nystrom <- numeric(2)

for(t in 1:2){
  
  set.seed(t)
  ll <- 5000
  NT <- 20000
  
  X <- matrix(runif((ll+NT)*2), ncol=2)
  y <- computer_simulator(X)
  
  Xtr <- X[1:ll,]; ytr <- y[1:ll]
  Xte <- X[(ll+1):(ll+NT),]
  yte <- y[(ll+1):(ll+NT)]
  
  fit <- fit_nystrom_MLE(Xtr, ytr, m = 20)
  
  yhat <- predict_nystrom(fit, Xtr, ytr, Xte)
  
  RMSE_nystrom[t] <- sqrt(mean((yte - yhat)^2))
}

mean(RMSE_nystrom)
sd(RMSE_nystrom)
#RFF
#rff
repmat <- function(A, N, M) {
  kronecker(matrix(1, N, M), A)
}
generateRFF <- function(eta, X, OmegaD) {
  # Implement the Algorithm 2
  # X : input data matrix
  # OmegaD : matrix with entries being i.i.d. N(0,1)
  # D : number of features
  #X<-xtrain
  #OmegaD<-s
  #eta<-c(opt$par[c(1:6)])
  n <- nrow(X)
  d <- ncol(X)
  D <- nrow(OmegaD)
  theta <- sqrt(2 * eta)
  Omega <- OmegaD %*% diag(theta)
  Z <- cbind(cos(X %*% t(Omega)), sin(X %*% t(Omega)))
  Z<-Z/sqrt(D)
  return(Z)
}
#computation reductiuon rff
d=2
k<-100*d
n<-nrow(xtrain)
U<-xtrain[,1]
V<-xtrain[,2]
H<-cbind(rep(1,n),U,V)
D1<-round(k)
s<-rmnorm(D1,rep(0,2),diag(2))
flogyy<-function(par){
  n<-nrow(xtrain)
  d=ncol(xtrain)
  delta<-par[1]
  theta<-par[c(2:(1+d))]
  Zgama<-generateRFF(theta,xtrain,s)
  SigmaRFF=Zgama%*%t(Zgama)+delta*diag(n)
  SS<-(1/delta)*(diag(n)-Zgama%*%solve(t(Zgama)%*%Zgama+delta*diag(2*k))%*%t(Zgama))
  betap<-solve(t(H)%*%SS%*%H+1e-8*diag(3))%*%t(H)%*%SS%*%ytrain
  sigma2<-(t(ytrain-H%*%betap)%*%SS%*%(ytrain-H%*%betap))/n
  flogyy=-(-n*log(sigma2)/2-0.5*log(det(SigmaRFF))
           -((t(ytrain-H%*%betap)%*%SS%*%(ytrain-H%*%betap))/(2*sigma2)))
  return(flogyy)
}
opthyperpama<-optim(c(1,rep(1,2)),flogyy,method ="L-BFGS-B",lower=c(1,rep(1e-8,2)),control=list(maxit=10)) 
k<-100*d
D1<-round(k)
s<-rmnorm (D1,rep(0,2),diag(2))
Zgama<-generateRFF(opthyperpama$par[c(2:(1+d))],xtrain,s)
SigmaRFF=Zgama%*%t(Zgama)+opthyperpama$par[1]*diag(n)
SS<-(1/opthyperpama$par[1])*(diag(n)-Zgama%*%solve(t(Zgama)%*%Zgama+opthyperpama$par[1]*diag(2*k))%*%t(Zgama))
Hte<-cbind(rep(1,nrow(xtest)),xtest)
betahat<-solve(t(H)%*%SS%*%H+1e-8*diag(3))%*%t(H)%*%SS%*%ytrain
Rtetr<-ra(xtrain,xtest,opthyperpama$par[c(2:(1+d))])
ytrff<-Hte%*%betahat+t(Rtetr)%*%SS%*%(ytrain-H%*%betahat)
# ytestpred1<-function(x){
#   x=matrix(xtest[i,],ncol=2)
#   rstar=ra(x,xtrain,opthyperpama$par[c(2:(1+d))])
#   ystar=computer_simulator(x)
#   y1=t(rstar)%*%SS%*%ytrain
#   Hte<-cbind(rep(1,nrow(x)),x)
#   h1=t(Hte)-t(rstar)%*%SS%*%H
#   omega=(1+opthyperpama$par[1]-t(rstar)%*%SS%*%rstar)^(-1)
#   mugamma=y1+h1%*%betahat
#   return(c(mugamma,omega))
# }
# ytrff<-array(0,NT)
# #varp2<-array(0,NT)
# for(i in 1 :NT){
#   ytrff[i]=ytestpred1(matrix(xtest[i,],ncol=2)[1]
#   #varp2[i]=ytestpred1(matrix(xtest[i,],ncol=6))[2]
# }
# ytrff<-as.numeric(ytrff)
rerff<-(ytest-ytrff)^2
RMSErff[t]<-sqrt(mean(rerff))
RMSErffp[t]<-sqrt(mean(ytest-Hte%*%betahat+t(Rtetr)%*%SS%*%(ytrain-H%*%betahat)))
}
MRMSErff<-mean(RMSErff)
MRMSErffp<-mean(RMSErffp)
vrff=var(ytrff)