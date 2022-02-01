runDESeq2 <- function(e, retDDS=FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}

runDESeq2Outliers <- function(e, retDDS=FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  ddsDefault <- DESeq(dds, quiet=TRUE)
  ddsNoRepl <- ddsDefault
  if (ncol(e) >= 14) {
    # insert original maximum Cook's distances
    # so the rows with replacement will be filtered
    # this avoid re-running with minReplicateForReplace=Inf
    mcols(ddsNoRepl)$maxCooks <- apply(assays(ddsNoRepl)[["cooks"]], 1, max)
  }
  resDefault <- results(ddsDefault)
  resNoFilt <- results(ddsDefault, cooksCutoff=FALSE)
  resNoRepl <- results(ddsNoRepl)
  resList <- list("DESeq2"=resDefault, "DESeq2-noFilt"=resNoFilt, "DESeq2-noRepl"=resNoRepl)
  resOut <- lapply(resList, function(res) {
    pvals <- res$pvalue
    padj <- res$padj
    pvals[is.na(pvals)] <- 1
    pvals[rowSums(exprs(e)) == 0] <- NA
    padj <- p.adjust(pvals,method="BH")
    padj[is.na(padj)] <- 1
    list(pvals=pvals, padj=padj)
  })
  return(resOut)
}

runEdgeR <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  dgel <- estimateGLMCommonDisp(dgel, design)
  dgel <- estimateGLMTrendedDisp(dgel, design)
  dgel <- estimateGLMTagwiseDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predbeta=predbeta[,"pData(e)$conditionB"], predbeta10=predbeta10[,"pData(e)$conditionB"])
}

runEdgeRRobust <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # settings for robust from robinson_lab/edgeR_robust/robust_simulation.R
  dgel <- estimateGLMRobustDisp(dgel, design, maxit=6)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predbeta=predbeta[,"pData(e)$conditionB"])
}

runDESeq <- function(e) {
  cds <- newCountDataSet(exprs(e), pData(e)$condition)
  cds <- DESeq::estimateSizeFactors(cds)
  cds <- DESeq::estimateDispersions(cds)
  suppressWarnings({capture.output({fit1 <- DESeq::fitNbinomGLMs(cds, count ~ condition)})})
  suppressWarnings({capture.output({fit0 <- DESeq::fitNbinomGLMs(cds, count ~ 1)})})
  pvals <- DESeq::nbinomGLMTest(fit1, fit0)
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=fit1$conditionB))
}

runDSS <- function(e) {
  X <- as.matrix(exprs(e))
  colnames(X) <- NULL
  designs <- as.character(pData(e)$condition)
  seqData <- newSeqCountSet(X, designs)
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData)
  result <- waldTest(seqData, "B", "A")
  result <- result[match(rownames(seqData),rownames(result)),]
  pvals <- result$pval
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=( log2(exp(1)) * result$lfc ))
}

runVoom <- function(e) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=tt$logFC)
}

runrope <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  rope.fit <- rope(datmat = dgel$counts, X_model = design,x_PI_idx = 2)
  pvals <- rope.fit$pvals
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta= rope.fit$logFC)
}

runrope.perm <- function(e) {
  design <- model.matrix( ~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  rope.fit <-
    rope(datmat = dgel$counts,
                X_model = design,
                x_PI_idx = 2)
  flag.na <- rowSums(exprs(e)) == 0

  logFC = ifelse(test = flag.na,
                 yes = 0,
                 no = rope.fit$logFC)
  logLR = ifelse(test = flag.na,
                 yes = 0,
                 no = rope.fit$logLR)
  adj = ifelse(test = flag.na,
               yes = 1,
               no = rope.fit$adj)
  adj.logLR = ifelse(test = flag.na,
                     yes = 0,
                     no = rope.fit$adj_logLR)

  list(
    logFC = logFC,
    logLR = logLR,
    adj = adj,
    adj.logLR = adj.logLR
  )
}

