makeSim <- function(n, m, x, beta, meanDispPairs, sf = rep(1, m)) {
  idx <- sample(nrow(meanDispPairs), n, replace = TRUE)
  mu0 <- meanDispPairs[idx, 1]
  disp <- meanDispPairs[idx, 2]
  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  muMat <- matrix(rep(mu, times = m) * rep(sf, each = n), ncol = m)
  list(
    mat = matrix(rnbinom(n * m, mu = muMat, size = 1 / disp), ncol = m),
    disp = disp,
    mu0 = mu0
  )
}

makeSim_GP <- function(n, m, x, beta, meanDispPairs_gp, sf = rep(1, m)) {
  idx <- sample(nrow(meanDispPairs), n, replace = TRUE)
  mu0 <- meanDispPairs[idx, 1]

  # disp <- meanDispPairs[idx,2]
  gp_lam2 <- meanDispPairs[idx, 3]

  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  # muMat <- matrix(rep(mu, times=m) * rep(sf, each=n), ncol=m)

  res_mat <- matrix(0, nrow = n, ncol = m)
  for (gi in 1:n) {
    lam1_g <- mu[gi, ] * (1 - gp_lam2[gi])
    res_mat[gi, ] <- sapply(lam1_g, RMKdiscrete::rLGP, n = 1, lambda = gp_lam2[gi])
  }

  list(
    mat = res_mat,
    gp_lam = gp_lam2,
    mu0 = mu0
  )
}


makeSim_GP_v1 <- function(n, m, x, beta, MDP_gp, sf = rep(1, m)) {
  idx <- sample(nrow(meanDispPairs), n, replace = TRUE)
  mu0 <- meanDispPairs[idx, 1]

  # disp <- meanDispPairs[idx,2]
  var0 <- meanDispPairs[idx, 4]

  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  # muMat <- matrix(rep(mu, times=m) * rep(sf, each=n), ncol=m)

  res_mat <- matrix(0, nrow = n, ncol = m)
  for (gi in 1:n) {
    mvpair <- RMKdiscrete::LGPMVP(mu = mu[gi, ], sigma2 = var0[gi])
    RMKdiscrete::rLGP(n = 1, theta = mvpair[, 1], lambda = mvpair[, 2])
    res_mat[gi, ] <- RMKdiscrete::rLGP(n = 1, theta = mvpair[, 1], lambda = mvpair[, 2])
  }

  list(
    mat = res_mat,
    mu0 = mu0,
    var0 = var0
  )
}

makeSim_linear_var <- function(n, m, x, beta, MDP_gp, sf = rep(1, m)) {
  idx <- sample(nrow(meanDispPairs), n, replace = TRUE)
  mu0 <- meanDispPairs[idx, 1]
  # disp <- meanDispPairs[idx,2]
  mul <- rlnorm(n, meanlog = 0, sdlog = 1.7)
  vf <- mu0 * 3
  var0 <- vf * mul

  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  # muMat <- matrix(rep(mu, times=m) * rep(sf, each=n), ncol=m)

  res_mat <- round(matrix(rgamma(n * m, as.numeric(mu)^2 / rep(mul, m) / vf, as.numeric(mu) / rep(mul, m) / vf), n, m))
  list(
    mat = res_mat,
    mu0 = mu0,
    var0 = var0
  )
}


makeSim_under_linear_var <- function(n, m, x, beta, MDP_gp, sf = rep(1, m)) {
  idx <- sample(nrow(meanDispPairs), n, replace = TRUE)
  mu0 <- meanDispPairs[idx, 1]
  # disp <- meanDispPairs[idx,2]
  mul <- rlnorm(n, meanlog = -1, sdlog = 1.2)
  vf <- mu0 * 1 / 3
  var0 <- vf * mul

  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  # muMat <- matrix(rep(mu, times=m) * rep(sf, each=n), ncol=m)

  res_mat <- round(matrix(rgamma(n * m, as.numeric(mu)^2 / rep(mul, m) / vf, as.numeric(mu) / rep(mul, m) / vf), n, m))
  list(
    mat = res_mat,
    mu0 = mu0,
    var0 = var0
  )
}


makeSim_log_var <- function(n, m, x, beta, MDP_gp, sf = rep(1, m)) {
  idx <- sample(nrow(meanDispPairs), n, replace = TRUE)
  mu0 <- meanDispPairs[idx, 1]
  # disp <- meanDispPairs[idx,2]
  mul <- rlnorm(n, meanlog = 0, sdlog = 1.2)
  vf <- mu0 + log(mu0 + 1)
  var0 <- vf * mul

  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  # muMat <- matrix(rep(mu, times=m) * rep(sf, each=n), ncol=m)

  res_mat <- round(matrix(rgamma(n * m, as.numeric(mu)^2 / rep(mul, m) / vf, as.numeric(mu) / rep(mul, m) / vf), n, m))
  list(
    mat = res_mat,
    mu0 = mu0,
    var0 = var0
  )
}
