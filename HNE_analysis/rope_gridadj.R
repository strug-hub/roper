# Robust adjustment factors -----------------------------------------------
A <-
  function(x, pred) {
    sum(pred * x ^ 2) - ((sum(pred * x)) ^ 2) / (sum(pred))
  }
B <-
  function(x, pred, res) {
    sum((res ^ 2) * (x - (sum(pred * x)) / (sum(pred))) ^ 2)
  }


#' Rope alternative full likelihood
#'
#' computing adj.profile likelihood over auto-search grid
rope_gridadj <- function(datmat, ngi, X_model, accuracy = 1000) {

  require(foreach)

  # Define parallel actions

  Ni.adj <- edgeR::calcNormFactors(datmat) * colSums(datmat)
  fit.mglm <-
    edgeR::mglmLevenberg(datmat,
      X_model,
      dispersion = 0,
      offset = log(Ni.adj)
    )

  ygi <- as.numeric(datmat[ngi, ])
  if (sum(ygi) == 0) {
    list(
      MLE = NA,
      grid = NA,
      log.norm.lik = NA,
      adjfact = NA,
      mlevszero = NA,
      extreme_control = NA,
      allzero = TRUE
    )
  } else {
    data_w <- data.frame(ygi, X_model, Ni.adj)

    X_reduced <- X_model[, -dim(X_model)[2]]

    # Define function parameters
    profile.theta <- "x_PI"
    theta.off <- data_w[, names(data_w) == profile.theta]
    offset.glm <- "Ni.adj"
    glm.off <- data_w[, names(data_w) == offset.glm]

    muyhat <- fit.mglm$fitted.values[ngi, ]
    resy <- datmat[ngi, ] - fit.mglm$fitted.values[ngi, ]

    A.est <- A(X_model[, "x_PI"], muyhat)
    B.est <- B(X_model[, "x_PI"], muyhat, resy)

    pre_mle <- fit.mglm$coefficients[, "x_PI"][ngi]

    # new grid based on SE of LRT
    se_l <- sqrt(A.est / B.est)
    mul.base <- 3

    # Add muliplying factors
    mul.fac <- 0
    MIN_0 <- log(1)
    # Control extreme adjustment
    extreme_control <- FALSE

    # pre-mm
    fit.p.i <-
      edgeR::mglmLevenberg(
        t(as.matrix(datmat[ngi, ])),
        X_reduced,
        dispersion = 0,
        offset = (pre_mle * theta.off + log(glm.off))
      )
    pre.mm <-
      sum(dpois(
        x = datmat[ngi, ],
        lambda = fit.p.i$fitted.values[1, ],
        log = T
      ))

    # Control extreme adjustment
    extreme_control <- FALSE
    while (MIN_0 > log(1 / 1000)) {
      beta_grid <-
        sort(c(
          seq(
            pre_mle - mul.base * se_l * (2^mul.fac),
            pre_mle + mul.base * se_l * (2^mul.fac),
            length.out = accuracy
          ),
          pre_mle
        ))
      # beta_grid <- sort(c(seq(pre_mle - mul.adj * se_l,pre_mle + mul.adj * se_l,length.out =accuracy),0,pre_mle))

      # left boundary
      p.l <- beta_grid[1]
      fit.p.i <-
        edgeR::mglmLevenberg(
          t(as.matrix(datmat[ngi, ])),
          X_reduced,
          dispersion = 0,
          offset = (p.l * theta.off + log(glm.off))
        )
      b.l <-
        sum(dpois(
          x = datmat[ngi, ],
          lambda = fit.p.i$fitted.values[1, ],
          log = T
        ))

      # right boundary
      p.r <- beta_grid[length(beta_grid)]
      fit.p.i <-
        edgeR::mglmLevenberg(
          t(as.matrix(datmat[ngi, ])),
          X_reduced,
          dispersion = 0,
          offset = (p.r * theta.off + log(glm.off))
        )
      b.r <-
        sum(dpois(
          x = datmat[ngi, ],
          lambda = fit.p.i$fitted.values[1, ],
          log = T
        ))

      MIN_0 <- min(c(b.l - pre.mm, b.r - pre.mm) * (A.est / B.est))
      mul.fac <- mul.fac + 1
    }

    beta_grid <-
      sort(c(
        seq(
          pre_mle - mul.base * se_l * (2^mul.fac),
          pre_mle + mul.base * se_l * (2^mul.fac),
          length.out = accuracy
        ),
        pre_mle,
        0
      ))

    if ((beta_grid[1] < -50) | (beta_grid[length(beta_grid)] > 50)) {
      extreme_control <- TRUE
      beta_grid <-
        sort(c(seq(-50, 50, length.out = accuracy), 0, pre_mle))
    }

    ##
    log.lik.t1 <-
      foreach(
        i = iterators::icount(length(beta_grid)),
        .packages = "edgeR",
        .combine = c
      ) %dopar% {
        pi <- beta_grid[i]
        fit.p.i <-
          edgeR::mglmLevenberg(
            t(as.matrix(datmat[ngi, ])),
            X_reduced,
            dispersion = 0,
            offset = (pi * theta.off + log(glm.off))
          )
        sum(dpois(
          x = datmat[ngi, ],
          lambda = fit.p.i$fitted.values[1, ],
          log = T
        ))
      }

    ### Create output
    theta <- beta_grid[is.na(log.lik.t1) != 1]
    log.lik <- log.lik.t1[is.na(log.lik.t1) != 1]
    profile.lik <- exp(log.lik)
    mm <- max(log.lik, na.rm = TRUE)
    log.norm.lik <- log.lik - mm
    profile.lik.norm <- exp(log.norm.lik)
    pmle <- max(theta[profile.lik.norm == max(profile.lik.norm)])

    ### Robust Adjustment
    log.norm.lik.adj <- log.norm.lik * (A.est / B.est)
    # MLE vs 0 Ratio
    # mlevsz <- exp(0 - log.norm.lik.adj[match(0,beta_grid)])
    lmlevsz <- 0 - log.norm.lik.adj[match(0, beta_grid)]
    list(
      MLE = pmle,
      grid = theta,
      log.norm.lik = log.norm.lik,
      adjfact = A.est / B.est,
      log.norm.lik.adj = log.norm.lik.adj,
      log_mlevszero = lmlevsz,
      extreme_control = extreme_control,
      allzero = FALSE,
      mul.iter = mul.fac
    )
  }
}
