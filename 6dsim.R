computer_simulator<-function(x){
  m=nrow(x)
  x1=x[,1]
  x2=x[,2]
  x3=x[,3]
  x4=x[,4]
  x5=x[,5]
  x6=x[,6]
  t1=x1^2
  t2=x2^2
  t3=x3^2
  t4=x4^2
  t5=x5^2
  t6=x6^2
  j1<-rep(1,m)
  j2<-rep(2,m)
  j3<-rep(3,m)
  j4<-rep(4,m)
  j5<-rep(5,m)
  j6<-rep(6,m)
  e<-rnorm(m,0,0.01)
  y=j1*t1+j2*t2+j3*t3+j4*t4+j5*t5+j6*t6+e
  return(y)
}
computer_simulator<-function(x){
  m=nrow(x)
  x1=x[,1]
  x2=x[,2]
  x3=x[,3]
  x4=x[,4]
  x5=x[,5]
  x6=x[,6]
  e<-rnorm(m,0,0.01)
  y=sin(8*x1+x2+x3+x4+x5+x6)+e
  return(y)
}
#simulation function
computer_simulator<-function(x){
  m=nrow(x)
  x1=x[,1]
  x2=x[,2]
  x3=x[,3]
  x4=x[,4]
  x5=x[,5]
  x6=x[,6]
  t1=x1^2
  t2=x2^2
  t3=x3^2
  t4=x4^2
  t5=x5^2
  t6=x6^2
  e<-rnorm(m,0,0.01)
  y=-sum(x1+x2+x3+x4+x5+x6)*exp(-sum(t1+t2+t3+t4+t5+t6))+e
  return(y)
}
R<-function(theta,h){
  RH<-sum(-h^2%*%theta)
  R<-exp(RH)
  return(R)
}
#reconstruction paramaterization
RA<- function(x,theta){
  m<-nrow(x)
  U1a<-x[,1]
  U2a<-x[,2]
  U3a<-x[,3]
  U4a<-x[,4]
  U5a<-x[,5]
  Va<-x[,6]
  RA <- matrix(0,m,m)
  for (i in 1:m) {
    for(j in 1:m){
      RA[i,j]= R(theta,c(U1a[i]-U1a[j],U2a[i]-U2a[j],U3a[i]-U3a[j],U4a[i]-U4a[j],U5a[i]-U5a[j],Va[i]-Va[j]))
    }
  }
  return(RA)
}
ra <- function(x,X,theta){
  m<-nrow(x)
  U1a<-x[,1]
  U2a<-x[,2]
  U3a<-x[,3]
  U4a<-x[,4]
  U5a<-x[,5]
  Va<-x[,6]
  n<-nrow(X)
  U1<-X[,1]
  U2<-X[,2]
  U3<-X[,3]
  U4<-X[,4]
  U5<-X[,5]
  V<-X[,6]
  ra<-matrix(0,m,n)
  for(i in 1:m){
    h <- cbind(U1-U1a[i],U2-U2a[i],U3-U3a[i],U4-U4a[i],U5-U5a[i],V-Va[i])
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
  U1a<-x[,1]
  U2a<-x[,2]
  U3a<-x[,3]
  U4a<-x[,4]
  U5a<-x[,5]
  Va<-x[,6]
  n<-nrow(X)
  U1<-X[,1]
  U2<-X[,2]
  U3<-X[,3]
  U4<-X[,4]
  U5<-X[,5]
  V<-X[,6]
  GA<-cbind(rep(1,m),U1a,U2a,U3a,U4a,U5a,Va)
  I<-rep(1,m)
  RAT<-solve(RAA+1e-8*diag(m))
  GRA<-solve(t(GA)%*%RAT%*%GA+1e-8*diag(7))
  P<-RAT%*%GA%*%GRA
  Q<-(diag(I)-RAT%*%GA%*%GRA%*%t(GA))%*%RAT
  b <- matrix(0,m,n)
  for (i in 1:n){
    b[,i]<-P%*%c(1,U1[i],U2[i],U3[i],U4[i],U5[i],V[i])+Q%*%rA[,i]
  }
  return(b)
}
CA <- function(X, m){
  n <- nrow(X)
  for(k in 1:n){
    C0 <- n
    t <- sample(1:n, m)
    U <- X[t,]
    l <- 1
    c_vec <- NULL
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        c_vec[l] <- 1/sqrt(sum((U[i,]-U[j,])^2) + 1e-8)
        l <- l+1
      }
    }
    CA_value <- max(c_vec)
    if(!is.na(CA_value) && CA_value <= C0){
      X_knots <- U
      C0 <- CA_value
      t0 <- t
    }
  }
  return(t0)
}

# Nystrom方法函数
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
  delta <- exp(par[1])
  theta <- exp(par[2])
  
  ns <- nystrom_features(X, Xm, theta)
  Phi <- ns$Phi
  n <- nrow(X)
  
  A <- diag(delta, n)
  B <- Phi
  Cinv <- diag(ncol(B))
  
  Sigma_inv <- (1/delta) * (
    diag(n) -
      B %*% solve(t(B)%*%B + delta*Cinv) %*% t(B)
  )
  
  beta_hat <- solve(
    t(H)%*%Sigma_inv%*%H,
    t(H)%*%Sigma_inv%*%y
  )
  
  r <- y - H %*% beta_hat
  
  logdet <- n*log(delta) +
    determinant(diag(ncol(B)) + (1/delta)*t(B)%*%B,
                logarithm = TRUE)$modulus
  
  val <- -0.5 * (
    logdet +
      t(r) %*% Sigma_inv %*% r
  )
  
  return(-as.numeric(val))
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

# RFF方法函数
rff_features <- function(X, W, b){
  n <- nrow(X)
  D <- nrow(W)
  
  X_norm <- X / sqrt(rowSums(X^2) + 1e-8)
  W_norm <- W / sqrt(rowSums(W^2) + 1e-8)
  
  inner_prod <- X_norm %*% t(W_norm)
  inner_prod <- pmin(pmax(inner_prod, -100), 100)
  
  Z <- sqrt(2/D) * cos(inner_prod + matrix(b, n, D, byrow=TRUE))
  Z
}

loglik_rff <- function(par, X, y, H, W0, b){
  delta <- exp(par[1])
  theta <- exp(par[2])
  
  W <- sqrt(2*theta) * W0
  Phi <- rff_features(X, W, b)
  n <- nrow(X)
  D <- ncol(Phi)
  
  Sigma_inv <- (1/delta) * (
    diag(n) -
      Phi %*% solve(t(Phi)%*%Phi + delta*diag(D)) %*% t(Phi)
  )
  
  beta_hat <- solve(
    t(H)%*%Sigma_inv%*%H,
    t(H)%*%Sigma_inv%*%y
  )
  
  r <- y - H %*% beta_hat
  
  logdet <- n*log(delta) +
    determinant(diag(D) + (1/delta)*t(Phi)%*%Phi,
                logarithm = TRUE)$modulus
  
  val <- -0.5 * (logdet + t(r)%*%Sigma_inv%*%r)
  return(-as.numeric(val))
}

fit_rff_MLE <- function(X, y, D){
  n <- nrow(X)
  d <- ncol(X)
  
  H <- cbind(1, X)
  W0 <- matrix(rnorm(D*d), D, d)
  b <- runif(D, 0, 2*pi)
  
  init_par <- c(log(0.01), log(1.0))
  lower_bounds <- c(log(1e-10), log(1e-3))
  upper_bounds <- c(log(100), log(1000))
  
  opt <- tryCatch({
    optim(
      par = init_par,
      fn = loglik_rff,
      X = X, y = y, H = H,
      W0 = W0, b = b,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(
        maxit = 1000,
        factr = 1e7,
        pgtol = 1e-4
      )
    )
  }, error = function(e) {
    message("优化失败: ", e$message)
    return(NULL)
  })
  
  if(is.null(opt)) {
    return(list(
      delta = 0.1,
      theta = 1.0,
      W0 = W0,
      b = b,
      H = H,
      convergence = -1
    ))
  }
  
  list(
    delta = exp(opt$par[1]),
    theta = exp(opt$par[2]),
    W0 = W0,
    b = b,
    H = H,
    convergence = opt$convergence,
    value = opt$value
  )
}

predict_rff <- function(fit, Xtr, ytr, Xte){
  delta <- fit$delta
  theta <- fit$theta
  W0 <- fit$W0
  b <- fit$b
  Htr <- fit$H
  Hte <- cbind(1, Xte)
  
  W <- sqrt(2*theta) * W0
  Phi_tr <- rff_features(Xtr, W, b)
  Phi_te <- rff_features(Xte, W, b)
  
  Sigma_inv <- (1/delta) * (
    diag(nrow(Xtr)) -
      Phi_tr %*% solve(t(Phi_tr)%*%Phi_tr + delta*diag(ncol(Phi_tr))) %*%
      t(Phi_tr)
  )
  
  beta_hat <- solve(
    t(Htr)%*%Sigma_inv%*%Htr,
    t(Htr)%*%Sigma_inv%*%ytr
  )
  
  r <- ytr - Htr %*% beta_hat
  
  yhat <- Hte %*% beta_hat +
    Phi_te %*% solve(t(Phi_tr)%*%Phi_tr + delta*diag(ncol(Phi_tr))) %*%
    t(Phi_tr) %*% r
  
  as.vector(yhat)
}

# 主循环：训练集大小和m值情况
d <- 6
ll_values <- c(5000, 10000, 15000, 20000)  # 训练集大小
m_values <- c(10*d, 15*d, 20*d)             # m值

# 存储结果的列表
results <- list()

for(ll_idx in 1:length(ll_values)){
  ll_current <- ll_values[ll_idx]
  cat(paste("\n=== 训练集大小 ll =", ll_current, "===\n"))
  
  results_ll <- list()  # 存储当前训练集大小的结果
  
  for(m_idx in 1:length(m_values)){
    m_current <- m_values[m_idx]
    cat(paste("\n--- m =", m_current, "---\n"))
    
    # 存储当前m值的所有结果
    current_results <- list()
    
    # 1. 信息先验方法
    RMSE1 <- array(0, 2)
    
    for(t in 1:2){
      set.seed(t)
      NT <- 20000
      
      # 生成数据
      u1 <- runif(NT + ll_current, 0, 1)
      u2 <- runif(NT + ll_current, 0, 1)
      u3 <- runif(NT + ll_current, 0, 1)
      u4 <- runif(NT + ll_current, 0, 1)
      u5 <- runif(NT + ll_current, 0, 1)
      v <- runif(NT + ll_current, 0, 1)
      xx <- cbind(u1, u2, u3,u4,u5, v)
      yy <- computer_simulator(xx)
      
      # 训练数据集
      xtrain <- xx[1:ll_current,]
      ytrain <- yy[1:ll_current]
      
      # 测试数据集
      xtest <- xx[(ll_current+1):(NT+ll_current),]
      ytest <- yy[(ll_current+1):(NT+ll_current)]
      
      # 选择m个节点
      m0 <- CA(xtrain, m_current)
      x1 <- xtrain[m0,]
      
      # 对数似然函数
      Logf <- function(par){
        m <- nrow(x1)
        d <- ncol(x1)
        n <- length(ytrain)
        sigma2 <- par[1]
        theta <- par[c(2:(1+d))]
        B <- b(x1, xtrain, theta)
        mum <- solve(B %*% t(B) + 1e-8*diag(m)) %*% B %*% ytrain
        value <- (-n/2*log(sigma2) - ((ytrain - t(mum) %*% B) %*% t(ytrain - t(mum) %*% B))/(2*sigma2))
        return(-value)
      }
      
      # 优化
      optLogf <- optim(c(0.1, rep(1, d)), Logf, method = "L-BFGS-B", 
                       lower = c(1e-8, rep(1e-8, d)))
      
      # 预测
      B <- b(x1, xtrain, optLogf$par[c(2:(1+d))])
      mum <- solve(B %*% t(B) + 1e-8*diag(m_current)) %*% B %*% ytrain
      ytp <- t(mum) %*% b(x1, xtest, optLogf$par[c(2:(1+d))])
      
      RMSE1[t] <- sqrt(mean((ytest - ytp)^2))
    }
    
    current_results$informative_RMSE <- mean(RMSE1)
    
    # 2. 无信息先验方法
    RMSE2 <- array(0, 2)
    
    for(t in 1:2){
      set.seed(t)
      NT <- 20000
      
      # 生成数据
      u1 <- runif(NT + ll_current, 0, 1)
      u2 <- runif(NT + ll_current, 0, 1)
      u3 <- runif(NT + ll_current, 0, 1)
      u4 <- runif(NT + ll_current, 0, 1)
      u5 <- runif(NT + ll_current, 0, 1)
      v <- runif(NT + ll_current, 0, 1)
      xx <- cbind(u1, u2, u3,u4,u5, v)
      yy <- computer_simulator(xx)
      
      # 训练数据集
      xtrain <- xx[1:ll_current,]
      ytrain <- yy[1:ll_current]
      
      # 测试数据集
      xtest <- xx[(ll_current+1):(NT+ll_current),]
      ytest <- yy[(ll_current+1):(NT+ll_current)]
      
      m0 <- CA(xtrain, m_current)
      x1 <- xtrain[m0,]
      
      Logff <- function(par){
        d <- ncol(x1)
        n <- length(ytrain)
        sigma2 <- par[1]
        theta <- par[c(2:(1+d))]
        B <- b(x1, xtrain, theta)
        gammah <- solve(B %*% t(B) + 1e-8*diag(m_current)) %*% B %*% ytrain
        value <- (-n/2*log(sigma2) - ((ytrain - t(gammah) %*% B) %*% t(ytrain - t(gammah) %*% B))/(2*sigma2))
        return(-value)
      }
      
      optLogff <- optim(c(0.1, rep(1, d)), Logff, method = "L-BFGS-B", 
                        lower = c(1e-8, rep(1e-8, d)))
      
      B <- b(x1, xtrain, optLogff$par[c(2:(1+d))])
      gammahat <- solve(B %*% t(B) + 1e-8*diag(m_current)) %*% B %*% ytrain
      yt2 <- t(gammahat) %*% b(x1, xtest, optLogff$par[c(2:(1+d))])
      
      RMSE2[t] <- sqrt(mean((ytest - yt2)^2))
    }
    
    current_results$noninformative_RMSE <- mean(RMSE2)
    
    # 3. Nystrom方法
    RMSE_nystrom <- numeric(2)
    
    for(t in 1:2){
      set.seed(t)
      NT <- 20000
      
      X <- matrix(runif((ll_current + NT) * d), ncol = d)
      y <- computer_simulator(X)
      
      Xtr <- X[1:ll_current,]
      ytr <- y[1:ll_current]
      Xte <- X[(ll_current+1):(ll_current+NT),]
      yte <- y[(ll_current+1):(ll_current+NT)]
      
      fit <- fit_nystrom_MLE(Xtr, ytr, m = m_current)
      yhat <- predict_nystrom(fit, Xtr, ytr, Xte)
      
      RMSE_nystrom[t] <- sqrt(mean((yte - yhat)^2))
    }
    
    current_results$nystrom_RMSE <- mean(RMSE_nystrom)
    
    # 4. RFF方法
    RMSE_RFF <- numeric(2)
    
    for(t in 1:2){
      set.seed(t)
      NT <- 20000
      
      X <- matrix(runif((ll_current + NT) * d), ncol = d)
      y <- computer_simulator(X)
      
      Xtr <- X[1:ll_current,]
      ytr <- y[1:ll_current]
      Xte <- X[(ll_current+1):(ll_current+NT),]
      yte <- y[(ll_current+1):(ll_current+NT)]
      
      fit <- fit_rff_MLE(Xtr, ytr, D = m_current)
      yhat <- predict_rff(fit, Xtr, ytr, Xte)
      
      RMSE_RFF[t] <- sqrt(mean((yte - yhat)^2))
    }
    
    current_results$rff_RMSE <- mean(RMSE_RFF)
    
    # 存储当前m值的结果
    results_ll[[paste0("m_", m_current)]] <- current_results
    
    # 打印当前m值的结果
    cat(paste("  m =", m_current, "的结果:\n"))
    cat(paste("    信息先验RMSE:", round(current_results$informative_RMSE, 6), "\n"))
    cat(paste("    无信息先验RMSE:", round(current_results$noninformative_RMSE, 6), "\n"))
    cat(paste("    Nystrom RMSE:", round(current_results$nystrom_RMSE, 6), "\n"))
    cat(paste("    RFF RMSE:", round(current_results$rff_RMSE, 6), "\n"))
  }
  
  # 存储当前训练集大小的结果
  results[[paste0("ll_", ll_current)]] <- results_ll
}

# 打印所有结果汇总
cat("\n=== 所有结果汇总 ===\n")
for(ll_val in ll_values){
  cat(paste("\n训练集大小 ll =", ll_val, ":\n"))
  res_ll <- results[[paste0("ll_", ll_val)]]
  
  for(m_val in m_values){
    cat(paste("\n  m =", m_val, ":\n"))
    res <- res_ll[[paste0("m_", m_val)]]
    cat(paste("    信息先验方法 RMSE:", round(res$informative_RMSE, 6), "\n"))
    cat(paste("    无信息先验方法 RMSE:", round(res$noninformative_RMSE, 6), "\n"))
    cat(paste("    Nystrom方法 RMSE:", round(res$nystrom_RMSE, 6), "\n"))
    cat(paste("    RFF方法 RMSE:", round(res$rff_RMSE, 6), "\n"))
  }
}

# 创建结果表格以便分析
cat("\n=== 结果表格 ===\n")
cat("训练集大小\tm值\t信息先验\t无信息先验\tNystrom\t\tRFF\n")
cat("----------------------------------------------------------------------------\n")

for(ll_val in ll_values){
  res_ll <- results[[paste0("ll_", ll_val)]]
  
  for(m_val in m_values){
    res <- res_ll[[paste0("m_", m_val)]]
    cat(sprintf("%d\t\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n",
                ll_val, m_val,
                res$informative_RMSE,
                res$noninformative_RMSE,
                res$nystrom_RMSE,
                res$rff_RMSE))
  }
  if(ll_val != ll_values[length(ll_values)]){
    cat("----------------------------------------------------------------------------\n")
  }
}