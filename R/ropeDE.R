# Main function: RoPE_fast ---------------------------------------------------------------
geneDE_fast <- function(datmat, X_model, x_PI_idx) {
  require(methods)
  require(edgeR)
  # Use external mglm fit to save time
  # induced parameters
  Ni.adj <- calcNormFactors(datmat) * colSums(datmat)
  fit.mglm <-
    mglmLevenberg(datmat,
      X_model,
      dispersion = 0,
      offset = log(Ni.adj)
    )
  X_reduced <- X_model[, -x_PI_idx]
  # Define function parameters
  muyhat <- fit.mglm$fitted.values
  resy <- datmat - fit.mglm$fitted.values
  A.est <-
    apply(
      X = muyhat,
      MARGIN = 1,
      FUN = A,
      x = X_model[, x_PI_idx]
    )
  B.est <-
    mapply(
      FUN = B,
      pred = as.data.frame(t(muyhat)),
      res = as.data.frame(t(resy)),
      MoreArgs = list(x = X_model[, x_PI_idx])
    )

  pre_mle <- fit.mglm$coefficients[, x_PI_idx]
  loglikemle <- rep(NA, nrow(datmat))
  loglike0 <- rep(NA, nrow(datmat))

  for (ngi in (1:nrow(datmat))[rowSums(datmat) > 0]) {
    fit.p.i <-
      mglmLevenberg(
        t(as.matrix(datmat[ngi, ])),
        X_reduced,
        dispersion = 0,
        offset = (pre_mle[ngi] * X_model[, x_PI_idx] + log(Ni.adj))
      )
    loglikemle[ngi] <-
      sum(dpois(
        x = datmat[ngi, ],
        lambda = fit.p.i$fitted.values[1, ],
        log = T
      ))
  }

  for (ngi in (1:nrow(datmat))[rowSums(datmat) > 0]) {
    fit.p.i <-
      mglmLevenberg(t(as.matrix(datmat[ngi, ])),
        X_reduced,
        dispersion = 0,
        offset = log(Ni.adj)
      )
    loglike0[ngi] <-
      sum(dpois(
        x = datmat[ngi, ],
        lambda = fit.p.i$fitted.values[1, ],
        log = T
      ))
  }

  log_mlevszero_gene_all <- loglikemle - loglike0
  adj <- (A.est / B.est)
  adj_logLR <- adj * log_mlevszero_gene_all
  pvals <-
    sapply(
      2 * adj_logLR,
      pchisq,
      df = 1,
      ncp = 0,
      lower.tail = F,
      log.p = FALSE
    )
  padj <- p.adjust(pvals, method = "BH")
  data.frame(
    logFC = pre_mle,
    logLR = log_mlevszero_gene_all,
    adj = A.est / B.est,
    adj_logLR = adj_logLR,
    pvals = pvals,
    q = padj
  )
}
