#' Robust profile likelihood-based differential expression analysis
#'
#' Fit a robust profile likelihood ratio test for a working Poisson log-linear model.
#'
#' @param datmat A matrix of filtered RNA-seq counts.
#' @param X_model A numerical design matrix contains DE experiment group (parameter of interest) and other covariates
#' @param x_PI_idx An integer indicates the column index of the parameter of interest, usually the DE experiment group.
#' @return A data frame with columns contain the DE result of log fold change, unadjusted and adjusted log likelihood ratio, adjustment factors, p-values and Benjamini & Hochberg adjusted p-values.
#' @export
#' @importFrom edgeR calcNormFactors mglmLevenberg
#' @importFrom stats dpois p.adjust pchisq
#' @examples
#' count <- matrix(rnbinom(n = 5e4, mu = 100, size = 1 / 0.5), ncol = 50)
#' mod <- model.matrix(~ gl(n = 2, k = 25))
#' res <- rope(datmat = count, X_model = mod, x_PI_idx = dim(mod)[2])
#' res
rope <- function(datmat, X_model, x_PI_idx) {
  Ni.adj <- edgeR::calcNormFactors(datmat) * colSums(datmat)
  fit.mglm <-
    edgeR::mglmLevenberg(datmat,
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
      edgeR::mglmLevenberg(
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
      edgeR::mglmLevenberg(t(as.matrix(datmat[ngi, ])),
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
    padj = padj
  )
}
