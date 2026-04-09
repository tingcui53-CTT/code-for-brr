# ============================================================================
# 完整的MCMC诊断代码 - 四种方法比较
# Nyström, RFF, Informative Prior, Non-Informative Prior
# ============================================================================

# 加载必要的包
library(MASS)
library(invgamma)
library(LearnBayes)
library(actuar)
library(coda)
library(ggplot2)

# ============================================================================
# 数据准备部分
# ============================================================================

# 假设数据已加载
# Folds5x2_pp 应该是您的数据框
d <- 4
m <- 5*d

# 提取变量
u1 <- Folds5x2_pp$AT
u2 <- Folds5x2_pp$V
u3 <- Folds5x2_pp$AP
v <- Folds5x2_pp$RH
xx <- cbind(u1, u2, u3, v)
yy <- Folds5x2_pp$PE

# 划分训练集和测试集
ll <- sample(9568, 8000)
xtrain <- xx[ll, ]
xtest <- xx[-ll, ]
ytrain <- yy[ll]
ytest <- yy[-ll]

# ============================================================================
# 基础函数定义
# ============================================================================

# 协方差函数
R <- function(theta, h) {
  RH <- sum(-h^2 %*% theta)
  R <- exp(RH)
  return(R)
}

# 节点选择函数
CA <- function(X) {
  for(k in 1:8000) {
    C0 <- 8000
    d <- 4
    n <- 8000
    m <- max(5*d, log(n))
    t <- sample(1:n, m)
    U <- X[t]
    l <- 1
    c <- NULL
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        c[l] <- 1/abs(U[i] - U[j])
        l <- l + 1
      }
    }
    CA <- max(c)
    if(!is.na(CA)) {
      if(CA <= C0) {
        X <- U
        C0 <- CA
        t0 <- t
      }
    }
  }
  return(t0)
}

# 选择节点
m0 <- CA(xtrain)
x1 <- xtrain[m0, ]
U1a <- x1[, 1]
U2a <- x1[, 2]
U3a <- x1[, 3]
Va <- x1[, 4]
y1 <- ytrain[m0]

# 相关矩阵函数
RA <- function(x, theta) {
  m <- nrow(x)
  U1a <- x[, 1]
  U2a <- x[, 2]
  U3a <- x[, 3]
  Va <- x[, 4]
  RA <- matrix(0, m, m)
  for(i in 1:m) {
    for(j in 1:m) {
      RA[i, j] <- R(theta, c(U1a[i] - U1a[j], U2a[i] - U2a[j], 
                             U3a[i] - U3a[j], Va[i] - Va[j]))
    }
  }
  return(RA)
}

ra <- function(x, X, theta) {
  m <- nrow(x)
  n <- nrow(X)
  ra <- matrix(0, m, n)
  for(i in 1:m) {
    h <- cbind(X[, 1] - x[i, 1], X[, 2] - x[i, 2], 
               X[, 3] - x[i, 3], X[, 4] - x[i, 4])
    for(j in 1:n) {
      ra[i, j] <- R(theta, h[j, ])
    }
  }
  return(ra)
}

b <- function(x, X, theta) {
  RAA <- RA(x, theta)
  rA <- ra(x, X, theta)
  m <- nrow(x)
  n <- nrow(X)
  U1a <- x[, 1]
  U2a <- x[, 2]
  U3a <- x[, 3]
  Va <- x[, 4]
  GA <- cbind(rep(1, m), U1a, U2a, U3a, Va)
  I <- rep(1, m)
  RAT <- solve(RAA + 1e-8 * diag(m))
  GRA <- solve(t(GA) %*% RAT %*% GA + 1e-8 * diag(5))
  P <- RAT %*% GA %*% GRA
  Q <- (diag(I) - RAT %*% GA %*% GRA %*% t(GA)) %*% RAT
  b <- matrix(0, m, n)
  for(i in 1:n) {
    b[, i] <- P %*% c(1, X[i, 1], X[i, 2], X[i, 3], X[i, 4]) + Q %*% rA[, i]
  }
  return(b)
}

# ============================================================================
# 信息性先验优化
# ============================================================================

Logf <- function(par) {
  m <- nrow(x1)
  d <- ncol(x1)
  n <- length(ytrain)
  theta <- par[c(1:d)]
  B <- b(x1, xtrain, theta)
  mum <- solve(B %*% t(B) + 1e-8 * diag(m)) %*% B %*% ytrain
  sigma2 <- (ytrain - t(mum) %*% B) %*% t(ytrain - t(mum) %*% B)
  value <- (-n/2 * log(sigma2) - n/2)
  return(-value)
}

optLogf <- optim(rep(1, d), Logf, method = "L-BFGS-B", lower = rep(1e-8, d))
Vgamma <- RA(x1, optLogf$par[c(1:d)])
SVgamma <- solve(Vgamma + 1e-8 * diag(m))
B <- b(x1, xtrain, optLogf$par[c(1:d)])
mum <- solve(B %*% t(B) + 1e-8 * diag(m)) %*% B %*% ytrain

# ============================================================================
# 非信息性先验优化
# ============================================================================

Logff <- function(par) {
  d <- ncol(x1)
  n <- length(ytrain)
  theta <- par[c(1:d)]
  B <- b(x1, xtrain, theta)
  gammah <- solve(B %*% t(B) + 1e-8 * diag(m)) %*% B %*% ytrain
  sigma2 <- (ytrain - t(gammah) %*% B) %*% t(ytrain - t(gammah) %*% B)
  value <- (-n/2 * log(sigma2) - n/2)
  return(-value)
}

optLogff <- optim(rep(1, d), Logff, method = "L-BFGS-B", lower = rep(1e-8, d))

# ============================================================================
# Nyström方法
# ============================================================================

cat("\n========== Nyström方法 ==========\n")

n <- nrow(xtrain)
k <- 5 * d
U1 <- xtrain[, 1]
U2 <- xtrain[, 2]
U3 <- xtrain[, 3]
V <- xtrain[, 4]
H <- cbind(rep(1, n), U1, U2, U3, V)
me <- sample(1:n, round(k))

# 优化超参数
floggy <- function(par) {
  n <- nrow(xtrain)
  d <- ncol(xtrain)
  delta <- par[1]
  theta <- par[c(2:(1+d))]
  rmn <- ra(xtrain[me, ], xtrain, theta)
  rmm <- ra(xtrain[me, ], xtrain[me, ], theta)
  rnm <- ra(xtrain, xtrain[me, ], theta)
  nySigma <- rnm %*% solve(rmm + 1e-8 * diag(k)) %*% rmn + delta * diag(n)
  SN <- (1/delta) * (diag(n) - rnm %*% solve(delta * rmm + rmn %*% rnm + 1e-8 * diag(k)) %*% rmn)
  betap <- solve(t(H) %*% SN %*% H + 1e-8 * diag(d+1)) %*% t(H) %*% SN %*% ytrain
  sigma2 <- (t(ytrain - H %*% betap) %*% SN %*% (ytrain - H %*% betap)) / n
  floggy <- -(-n * log(sigma2)/2 - 0.5 * log(det(nySigma)) -
                ((t(ytrain - H %*% betap) %*% SN %*% (ytrain - H %*% betap)) / (2 * sigma2)))
  return(floggy)
}

opthyperpama1 <- optim(c(1, rep(1, d)), floggy, method = "L-BFGS-B", 
                       lower = c(1e-2, rep(1e-8, d)), upper = c(1, rep(10, d)), 
                       control = list(maxit = 10))

# 计算预测
k <- 5 * d
me <- sample(1:n, round(k))
rmn <- ra(xtrain[me, ], xtrain, opthyperpama1$par[c(2:(1+d))])
rmm <- ra(xtrain[me, ], xtrain[me, ], opthyperpama1$par[c(2:(1+d))])
rnm <- ra(xtrain, xtrain[me, ], opthyperpama1$par[c(2:(1+d))])
nySigma <- rnm %*% solve(rmm) %*% rmn + opthyperpama1$par[2] * diag(n)
SN <- (1/opthyperpama1$par[1]) * (diag(n) - rnm %*% solve(opthyperpama1$par[1] * rmm + rmn %*% rnm + 1e-8 * diag(k)) %*% rmn)
betahat <- solve(t(H) %*% SN %*% H) %*% t(H) %*% SN %*% ytrain
Hte <- cbind(rep(1, nrow(xtest)), xtest)
Rtetr <- ra(xtrain, xtest, opthyperpama1$par[c(2:(1+d))])
ytny <- Hte %*% betahat + t(Rtetr) %*% SN %*% (ytrain - H %*% betahat)
RMSEny <- sqrt(mean((ytest - ytny)^2))

# Nyström MCMC
Lthetany <- function(thetap, beta, sigma2) {
  n <- length(ytrain)
  value <- (-(n/2) * log(sigma2) - 0.5 * log(det(nySigma)) -
              ((t(ytrain - H %*% beta) %*% SN %*% (ytrain - H %*% beta)) / (2 * sigma2)) + 
              2 * log(2) - 0.5 * sum(thetap))
  return(value)
}

cat("\n运行Nyström MCMC...\n")
N <- 10000
p <- d + 1
theta_nystrom <- matrix(0, N, d)
sigma2_nystrom <- array(0, N)
beta_nystrom <- matrix(0, N, p)
sigma2_nystrom[1] <- 0.1
theta_nystrom[1, ] <- opthyperpama1$par[c(2:(1+d))]

accept_count_nystrom <- 0
cum_mean_sigma2_nystrom <- numeric(N)
start_time_nystrom <- Sys.time()

for(k in 2:N) {
  betahat <- solve(t(H) %*% SN %*% H) %*% t(H) %*% SN %*% ytrain
  beta_nystrom[k, ] <- rmnorm(1, betahat, sigma2_nystrom[k-1] * solve(t(H) %*% SN %*% H))
  
  sigma2_nystrom[k] <- rinvgamma(1, n/2, 
                                 scale = (t(ytrain - H %*% beta_nystrom[k, ]) %*% 
                                            SN %*% (ytrain - H %*% beta_nystrom[k, ])) / 2)
  
  lu <- rmnorm(1, log(abs(theta_nystrom[k-1, ])), 0.2 * diag(d))
  u <- exp(lu)
  pitheta <- Lthetany(theta_nystrom[k-1, ], beta_nystrom[k, ], sigma2_nystrom[k])
  piu <- Lthetany(u, beta_nystrom[k, ], sigma2_nystrom[k])
  minvalue <- min(1, piu / pitheta)
  l <- runif(1, 0, 1)
  
  if(!is.na(l < minvalue)) {
    theta_nystrom[k, ] <- u
    accept_count_nystrom <- accept_count_nystrom + 1
  } else {
    theta_nystrom[k, ] <- theta_nystrom[k-1, ]
  }
  
  cum_mean_sigma2_nystrom[k] <- mean(sigma2_nystrom[1:k])
  
  if(k %% 1000 == 0) cat("  Nyström迭代:", k, "/", N, "\n")
}

end_time_nystrom <- Sys.time()
comp_time_nystrom <- as.numeric(difftime(end_time_nystrom, start_time_nystrom, units = "secs"))
accept_rate_nystrom <- accept_count_nystrom / (N - 1)

# 计算ESS
burn_in <- 2001
sigma2_mcmc_nystrom <- mcmc(sigma2_nystrom[burn_in:N])
theta_mcmc_nystrom <- mcmc(theta_nystrom[burn_in:N, ])
ess_sigma2_nystrom <- effectiveSize(sigma2_mcmc_nystrom)
ess_theta_nystrom <- mean(effectiveSize(theta_mcmc_nystrom))

cat("\nNyström诊断结果:\n")
cat("  接受率:", accept_rate_nystrom, "\n")
cat("  计算时间:", comp_time_nystrom, "秒\n")
cat("  ESS (σ²):", ess_sigma2_nystrom, "\n")
cat("  ESS (θ):", ess_theta_nystrom, "\n")
cat("  RMSE:", RMSEny, "\n")

# ============================================================================
# RFF方法
# ============================================================================

cat("\n========== RFF方法 ==========\n")

# RFF生成函数
generateRFF <- function(eta, X, OmegaD) {
  n <- nrow(X)
  d <- ncol(X)
  D <- nrow(OmegaD)
  thetaa <- sqrt(2 * eta)
  Omega <- OmegaD %*% diag(thetaa)
  Z <- cbind(cos(X %*% t(Omega)), sin(X %*% t(Omega)))
  Z <- Z / sqrt(D)
  return(Z)
}

d <- 4
k <- 5 * d
n <- nrow(xtrain)
H <- cbind(rep(1, n), xtrain[, 1], xtrain[, 2], xtrain[, 3], xtrain[, 4])
D1 <- round(k)
s <- rmnorm(D1, rep(0, d), diag(d))

flogyy <- function(par) {
  n <- nrow(xtrain)
  d <- ncol(xtrain)
  delta <- par[1]
  theta <- par[c(2:(1+d))]
  Zgama <- generateRFF(theta, xtrain, s)
  SigmaRFF <- Zgama %*% t(Zgama) + delta * diag(n)
  SS <- (1/delta) * (diag(n) - Zgama %*% solve(t(Zgama) %*% Zgama + delta * diag(2*k)) %*% t(Zgama))
  betap <- solve(t(H) %*% SS %*% H + 1e-8 * diag(d+1)) %*% t(H) %*% SS %*% ytrain
  sigma2 <- (t(ytrain - H %*% betap) %*% SS %*% (ytrain - H %*% betap)) / n
  flogyy <- -(-n * log(sigma2)/2 - 0.5 * log(det(SigmaRFF)) -
                ((t(ytrain - H %*% betap) %*% SS %*% (ytrain - H %*% betap)) / (2 * sigma2)))
  return(flogyy)
}

opthyperpama <- optim(c(1, rep(1, d)), flogyy, method = "L-BFGS-B",
                      lower = c(1e-2, rep(1e-8, d)), upper = c(1, rep(12.5, d)),
                      control = list(maxit = 10))

k <- 5 * d
D1 <- round(k)
s <- rmnorm(D1, rep(0, d), diag(d))
Zgama <- generateRFF(opthyperpama$par[c(2:(1+d))], xtrain, s)
SigmaRFF <- Zgama %*% t(Zgama) + opthyperpama$par[1] * diag(n)
SS <- (1/opthyperpama$par[1]) * (diag(n) - Zgama %*% solve(t(Zgama) %*% Zgama + opthyperpama$par[1] * diag(2*k)) %*% t(Zgama))
betahat <- solve(t(H) %*% SS %*% H + 1e-8 * diag(d+1)) %*% t(H) %*% SS %*% ytrain
Rtetr <- ra(xtrain, xtest, opthyperpama$par[c(2:(1+d))])
ytrff <- Hte %*% betahat + t(Rtetr) %*% SS %*% (ytrain - H %*% betahat)
RMSErff <- sqrt(mean((ytest - ytrff)^2))

# RFF MCMC
Lthetarff <- function(thetap, beta, sigma2) {
  n <- length(ytrain)
  value <- (-(n/2) * log(sigma2) - 0.5 * log(det(SigmaRFF)) -
              ((t(ytrain - H %*% beta) %*% SS %*% (ytrain - H %*% beta)) / (2 * sigma2)) +
              2 * log(2) - 0.5 * sum(thetap))
  return(value)
}

cat("\n运行RFF MCMC...\n")
p=5
theta_rff <- matrix(0, N, d)
sigma2_rff <- array(0, N)
beta_rff <- matrix(0, N, p)
sigma2_rff[1] <- 0.1
theta_rff[1, ] <- opthyperpama$par[c(2:(1+d))]

accept_count_rff <- 0
cum_mean_sigma2_rff <- numeric(N)
start_time_rff <- Sys.time()

for(k in 2:N) {
  betahat <- solve(t(H) %*% SS %*% H) %*% t(H) %*% SS %*% ytrain
  beta_rff[k, ] <- rmnorm(1, betahat, sigma2_rff[k-1] * solve(t(H) %*% SS %*% H))
  
  sigma2_rff[k] <- rinvgamma(1, n/2, 
                             scale = (t(ytrain - H %*% beta_rff[k, ]) %*% 
                                        SS %*% (ytrain - H %*% beta_rff[k, ])) / 2)
  
  lu <- rmnorm(1, log(abs(theta_rff[k-1, ])), 0.2 * diag(d))
  u <- exp(lu)
  pitheta <- Lthetarff(theta_rff[k-1, ], beta_rff[k, ], sigma2_rff[k])
  piu <- Lthetarff(u, beta_rff[k, ], sigma2_rff[k])
  minvalue <- min(1, piu / pitheta)
  l <- runif(1, 0, 1)
  
  if(!is.na(l < minvalue)) {
    theta_rff[k, ] <- u
    accept_count_rff <- accept_count_rff + 1
  } else {
    theta_rff[k, ] <- theta_rff[k-1, ]
  }
  
  cum_mean_sigma2_rff[k] <- mean(sigma2_rff[1:k])
  
  if(k %% 1000 == 0) cat("  RFF迭代:", k, "/", N, "\n")
}

end_time_rff <- Sys.time()
comp_time_rff <- as.numeric(difftime(end_time_rff, start_time_rff, units = "secs"))
accept_rate_rff <- accept_count_rff / (N - 1)

sigma2_mcmc_rff <- mcmc(sigma2_rff[burn_in:N])
theta_mcmc_rff <- mcmc(theta_rff[burn_in:N, ])
ess_sigma2_rff <- effectiveSize(sigma2_mcmc_rff)
ess_theta_rff <- mean(effectiveSize(theta_mcmc_rff))

cat("\nRFF诊断结果:\n")
cat("  接受率:", accept_rate_rff, "\n")
cat("  计算时间:", comp_time_rff, "秒\n")
cat("  ESS (σ²):", ess_sigma2_rff, "\n")
cat("  ESS (θ):", ess_theta_rff, "\n")
cat("  RMSE:", RMSErff, "\n")

# ============================================================================
# 信息性先验 MCMC
# ============================================================================

cat("\n========== 信息性先验方法 ==========\n")

Ltheta <- function(thetap, gamma, sigma2) {
  n <- length(ytrain)
  B <- b(x1, xtrain, thetap)
  value <- (-(n/2) * log(sigma2) - 
              (((ytrain - t(gamma) %*% B) %*% t(ytrain - t(gamma) %*% B)) / (2 * sigma2)) +
              2 * log(2) - 0.5 * sum(thetap))
  return(value)
}

cat("\n运行信息性先验MCMC...\n")
theta_info <- matrix(0, N, d)
sigma2_info <- array(0, N)
gamma_info <- matrix(0, N, m)

theta_info[1, ] <- optLogf$par[c(2:(1+d))]
gamma_info[1, ] <- mum
sigma2_info[1] <- optLogf$par[1]

accept_count_info <- 0
cum_mean_sigma2_info <- numeric(N)
start_time_info <- Sys.time()

for(k in 2:N) {
  B <- b(x1, xtrain, theta_info[k-1, ])
  Sigmagamma <- solve(B %*% t(B) + solve(Vgamma + 1e-8 * diag(m)))
  mugamma <- Sigmagamma %*% (B %*% ytrain + SVgamma %*% mum)
  gamma_info[k, ] <- rmnorm(1, mugamma, sigma2_info[k-1] * Sigmagamma)
  
  sigma2_info[k] <- rinvgamma(1, (n + m)/2, 
                              scale = (t(ytrain) %*% ytrain + t(mum) %*% SVgamma %*% mum -
                                         2 * t(mugamma) %*% solve(Sigmagamma) %*% gamma_info[k, ] +
                                         t(gamma_info[k, ]) %*% solve(Sigmagamma) %*% gamma_info[k, ]) / 2)
  
  lu <- rmnorm(1, log(abs(theta_info[k-1, ])), 0.1 * diag(d))
  u <- exp(lu)
  pitheta <- Ltheta(theta_info[k-1, ], gamma_info[k, ], sigma2_info[k])
  piu <- Ltheta(u, gamma_info[k, ], sigma2_info[k])
  minvalue <- min(1, piu / pitheta)
  l <- runif(1, 0, 1)
  
  if(!is.na(l < minvalue)) {
    theta_info[k, ] <- u
    accept_count_info <- accept_count_info + 1
  } else {
    theta_info[k, ] <- theta_info[k-1, ]
  }
  
  cum_mean_sigma2_info[k] <- mean(sigma2_info[1:k])
  
  if(k %% 1000 == 0) cat("  信息性先验迭代:", k, "/", N, "\n")
}

end_time_info <- Sys.time()
comp_time_info <- as.numeric(difftime(end_time_info, start_time_info, units = "secs"))
accept_rate_info <- accept_count_info / (N - 1)

sigma2_mcmc_info <- mcmc(sigma2_info[burn_in:N])
theta_mcmc_info <- mcmc(theta_info[burn_in:N, ])
ess_sigma2_info <- effectiveSize(sigma2_mcmc_info)
ess_theta_info <- mean(effectiveSize(theta_mcmc_info))

# 计算信息性先验RMSE
B <- b(x1, xtest, colMeans(theta_info[burn_in:N, ]))
ystarpred <- t(colMeans(gamma_info[burn_in:N, ])) %*% B
RMSEF1 <- sqrt(mean((ytest - as.numeric(ystarpred))^2))

cat("\n信息性先验诊断结果:\n")
cat("  接受率:", accept_rate_info, "\n")
cat("  计算时间:", comp_time_info, "秒\n")
cat("  ESS (σ²):", ess_sigma2_info, "\n")
cat("  ESS (θ):", ess_theta_info, "\n")
cat("  RMSE:", RMSEF1, "\n")

# ============================================================================
# 非信息性先验 MCMC
# ============================================================================

cat("\n========== 非信息性先验方法 ==========\n")

cat("\n运行非信息性先验MCMC...\n")
theta_noninfo <- matrix(0, N, d)
sigma2_noninfo <- array(0, N)
gamma_noninfo <- matrix(0, N, m)

theta_noninfo[1, ] <- optLogff$par[c(2:(1+d))]
gamma_noninfo[1, ] <- y1
sigma2_noninfo[1] <- optLogff$par[1]

accept_count_noninfo <- 0
cum_mean_sigma2_noninfo <- numeric(N)
start_time_noninfo <- Sys.time()

for(k in 2:N) {
  B <- b(x1, xtrain, theta_noninfo[k-1, ])
  gammahat <- solve(B %*% t(B) + 1e-6 * diag(m)) %*% B %*% ytrain
  Sigmagamma <- solve(B %*% t(B) + 1e-6 * diag(m))
  gamma_noninfo[k, ] <- rmnorm(1, gammahat, sigma2_noninfo[k-1] * Sigmagamma)
  
  sigma2_noninfo[k] <- rinvgamma(1, n/2, 
                                 scale = t(ytrain - t(gamma_noninfo[k, ]) %*% B) %*% 
                                   (ytrain - t(gamma_noninfo[k, ]) %*% B) / 2)
  
  lu <- rmnorm(1, log(abs(theta_noninfo[k-1, ])), 0.2 * diag(d))
  u <- exp(lu)
  pitheta <- Ltheta(theta_noninfo[k-1, ], gamma_noninfo[k, ], sigma2_noninfo[k])
  piu <- Ltheta(u, gamma_noninfo[k, ], sigma2_noninfo[k])
  minvalue <- min(1, piu / pitheta)
  l <- runif(1, 0, 1)
  
  if(!is.na(l < minvalue)) {
    theta_noninfo[k, ] <- u
    accept_count_noninfo <- accept_count_noninfo + 1
  } else {
    theta_noninfo[k, ] <- theta_noninfo[k-1, ]
  }
  
  cum_mean_sigma2_noninfo[k] <- mean(sigma2_noninfo[1:k])
  
  if(k %% 1000 == 0) cat("  非信息性先验迭代:", k, "/", N, "\n")
}

end_time_noninfo <- Sys.time()
comp_time_noninfo <- as.numeric(difftime(end_time_noninfo, start_time_noninfo, units = "secs"))
accept_rate_noninfo <- accept_count_noninfo / (N - 1)

sigma2_mcmc_noninfo <- mcmc(sigma2_noninfo[burn_in:N])
theta_mcmc_noninfo <- mcmc(theta_noninfo[burn_in:N, ])
ess_sigma2_noninfo <- effectiveSize(sigma2_mcmc_noninfo)
ess_theta_noninfo <- mean(effectiveSize(theta_mcmc_noninfo))

# 计算非信息性先验RMSE
B <- b(x1, xtest, colMeans(theta_noninfo[burn_in:N, ]))
ystarpred1 <- t(colMeans(gamma_noninfo[burn_in:N, ])) %*% B
RMSEF2 <- sqrt(mean((ytest - as.numeric(ystarpred1))^2))

cat("\n非信息性先验诊断结果:\n")
cat("  接受率:", accept_rate_noninfo, "\n")
cat("  计算时间:", comp_time_noninfo, "秒\n")
cat("  ESS (σ²):", ess_sigma2_noninfo, "\n")
cat("  ESS (θ):", ess_theta_noninfo, "\n")
cat("  RMSE:", RMSEF2, "\n")

# ============================================================================
# 生成诊断图
# ============================================================================

cat("\n========== 生成诊断图 ==========\n")

# Trace Plots
par(mfrow = c(2, 2))
plot(sigma2_nystrom[burn_in:N], type = "l", main = "Nyström: Trace Plot of σ²",
     xlab = "Iteration", ylab = "σ²", col = "blue")
abline(h = mean(sigma2_nystrom[burn_in:N]), col = "red", lty = 2)

plot(sigma2_rff[burn_in:N], type = "l", main = "RFF: Trace Plot of σ²",
     xlab = "Iteration", ylab = "σ²", col = "red")
abline(h = mean(sigma2_rff[burn_in:N]), col = "red", lty = 2)

plot(sigma2_info[burn_in:N], type = "l", main = "Informative Prior: Trace Plot of σ²",
     xlab = "Iteration", ylab = "σ²", col = "green")
abline(h = mean(sigma2_info[burn_in:N]), col = "red", lty = 2)

plot(sigma2_noninfo[burn_in:N], type = "l", main = "Non-Informative Prior: Trace Plot of σ²",
     xlab = "Iteration", ylab = "σ²", col = "purple")
abline(h = mean(sigma2_noninfo[burn_in:N]), col = "red", lty = 2)

# Cumulative Mean Plots
par(mfrow = c(2, 2))
plot(cum_mean_sigma2_nystrom[burn_in:N], type = "l", main = "Nyström: Cumulative Mean of σ²",
     xlab = "Iteration", ylab = "Cumulative Mean", col = "blue")

plot(cum_mean_sigma2_rff[burn_in:N], type = "l", main = "RFF: Cumulative Mean of σ²",
     xlab = "Iteration", ylab = "Cumulative Mean", col = "red")

plot(cum_mean_sigma2_info[burn_in:N], type = "l", main = "Informative Prior: Cumulative Mean of σ²",
     xlab = "Iteration", ylab = "Cumulative Mean", col = "green")

plot(cum_mean_sigma2_noninfo[burn_in:N], type = "l", main = "Non-Informative Prior: Cumulative Mean of σ²",
     xlab = "Iteration", ylab = "Cumulative Mean", col = "purple")

par(mfrow = c(1, 1))

# ============================================================================
# 性能汇总表
# ============================================================================

cat("\n========== 性能汇总表 ==========\n")

performance_summary <- data.frame(
  Method = c("Nyström", "RFF", "Informative Prior", "Non-Informative Prior"),
  Acceptance_Rate = round(c(accept_rate_nystrom, accept_rate_rff, 
                            accept_rate_info, accept_rate_noninfo), 4),
  Time_Seconds = round(c(comp_time_nystrom, comp_time_rff, 
                         comp_time_info, comp_time_noninfo), 2),
  ESS_sigma2 = round(c(ess_sigma2_nystrom, ess_sigma2_rff, 
                       ess_sigma2_info, ess_sigma2_noninfo), 0),
  ESS_theta = round(c(ess_theta_nystrom, ess_theta_rff, 
                      ess_theta_info, ess_theta_noninfo), 0),
  RMSE = round(c(RMSEny, RMSErff, RMSEF1, RMSEF2), 4)
)

print(performance_summary)

# ============================================================================
# 生成LaTeX表格
# ============================================================================

cat("\n========== LaTeX表格 ==========\n")

cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Performance summary of MCMC algorithms for four methods.}\n")
cat("\\label{tab:mcmc_performance}\n")
cat("\\begin{tabular}{lccccc}\n")
cat("\\toprule\n")
cat("Method & Acceptance Rate & Time (s) & ESS ($\\sigma^2$) & ESS ($\\theta$) & RMSE \\\\\n")
cat("\\midrule\n")
for(i in 1:4) {
  cat(sprintf("%s & %.4f & %.2f & %.0f & %.0f & %.4f \\\\\n",
              performance_summary$Method[i],
              performance_summary$Acceptance_Rate[i],
              performance_summary$Time_Seconds[i],
              performance_summary$ESS_sigma2[i],
              performance_summary$ESS_theta[i],
              performance_summary$RMSE[i]))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")

# ============================================================================
# 保存结果
# ============================================================================

saveRDS(list(
  nystrom = list(sigma2 = sigma2_nystrom, theta = theta_nystrom, 
                 acceptance_rate = accept_rate_nystrom, time = comp_time_nystrom,
                 ess_sigma2 = ess_sigma2_nystrom, ess_theta = ess_theta_nystrom,
                 cum_mean = cum_mean_sigma2_nystrom),
  rff = list(sigma2 = sigma2_rff, theta = theta_rff,
             acceptance_rate = accept_rate_rff, time = comp_time_rff,
             ess_sigma2 = ess_sigma2_rff, ess_theta = ess_theta_rff,
             cum_mean = cum_mean_sigma2_rff),
  informative = list(sigma2 = sigma2_info, theta = theta_info,
                     acceptance_rate = accept_rate_info, time = comp_time_info,
                     ess_sigma2 = ess_sigma2_info, ess_theta = ess_theta_info,
                     cum_mean = cum_mean_sigma2_info),
  noninformative = list(sigma2 = sigma2_noninfo, theta = theta_noninfo,
                        acceptance_rate = accept_rate_noninfo, time = comp_time_noninfo,
                        ess_sigma2 = ess_sigma2_noninfo, ess_theta = ess_theta_noninfo,
                        cum_mean = cum_mean_sigma2_noninfo),
  performance_table = performance_summary
), "mcmc_diagnostics_results.rds")

cat("\n诊断结果已保存到 mcmc_diagnostics_results.rds\n")
cat("\n========== 诊断完成 ==========\n")