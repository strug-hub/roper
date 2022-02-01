run_edgeR_exact<-function(object) {
  #A function that runs edgeR exact test
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(edgeR)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  y <- DGEList(counts=count, group=cond)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)

  et <- exactTest(y)
  et <- et$table
  et$q <- p.adjust(et$PValue, method="BH")
  et$q[is.na(et$q)] <- 1

  et$Genes <- rownames(et)
  et <- merge(et, de.gene, by="Genes")
  colnames(et) <- c("Genes", "LFC", "logCPM", "p", "q", "de.genes")
  return(list(result=et,  inputs = object, tool.name="edgeR exact"))
}

# LZ
run_edgeR_glm<-function(object) {
  #A function that runs edgeR GLM
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(edgeR)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator


  y <- DGEList(counts=count, group=cond)
  y <- calcNormFactors(y)
  design <- model.matrix(~cond)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  dFit <- glmFit(y, design)
  dLrt <- glmLRT(dFit, coef=2)
  dLrt <- dLrt$table
  dLrt$q <- p.adjust(dLrt$PValue,method="BH")
  dLrt$q[is.na(dLrt$q)] <- 1

  dLrt$Genes <- rownames(dLrt)
  dLrt <- merge(dLrt, de.gene, by="Genes")
  dLrt <- dLrt[, -4]
  colnames(dLrt) <- c("Genes", "LFC", "logCPM", "p", "q", "de.genes")

  return(list(result=dLrt, inputs = object, tool.name="edgeR GLM"))
}

run_edgeR_glm_l<-function(object) {
  #A function that runs edgeR GLM
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(edgeR)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator


  y <- DGEList(counts=count, group=cond)
  y <- calcNormFactors(y)
  design <- model.matrix(~cond)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  dFit <- glmFit(y, design)
  dLrt <- glmLRT(dFit, coef=2)
  dLrt <- dLrt$table
  dLrt$q <- p.adjust(dLrt$PValue,method="BH")
  dLrt$q[is.na(dLrt$q)] <- 1

  dLrt$Genes <- rownames(dLrt)
  dLrt <- merge(dLrt, de.gene, by="Genes")
  dLrt <- dLrt[, -4]
  colnames(dLrt) <- c("Genes", "LFC", "logCPM", "p", "q", "de.genes")

  return(list(result=dLrt, tool.name="edgeR GLM"))
}




run_edgeR_robust_10priorDF<-function(object) {
  #A function that runs robust edgeR GLM
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  #default prior degrees of freedm (estimateGLMRobustDisp(.,prior.df = 10))
  require(edgeR)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  y <- DGEList(counts=count, group=cond)
  y <- calcNormFactors(y)
  design <- model.matrix(~cond)
  y <- estimateGLMRobustDisp(y, design, prior.df = 10)

  ## DE test
  dFit <- glmFit(y, design)
  dLrt <- glmLRT(dFit,coef=2)                       # Test

  dLrt <- dLrt$table
  dLrt$q <- p.adjust(dLrt$PValue, method="BH")
  dLrt$q[is.na(dLrt$q)] <- 1

  dLrt$Genes <- rownames(dLrt)
  dLrt <- merge(dLrt, de.gene, by="Genes")
  dLrt <- dLrt[, -4]
  colnames(dLrt) <- c("Genes", "LFC", "logCPM", "p", "q", "de.genes")

  return(list(result=dLrt, inputs = object, tool.name="edgeR robust (prior.df=10)"))
}

run_edgeR_robust_10priorDF_l<-function(object) {
  #A function that runs robust edgeR GLM
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  #default prior degrees of freedm (estimateGLMRobustDisp(.,prior.df = 10))
  require(edgeR)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  y <- DGEList(counts=count, group=cond)
  y <- calcNormFactors(y)
  design <- model.matrix(~cond)
  y <- estimateGLMRobustDisp(y, design, prior.df = 10)

  ## DE test
  dFit <- glmFit(y, design)
  dLrt <- glmLRT(dFit,coef=2)                       # Test

  dLrt <- dLrt$table
  dLrt$q <- p.adjust(dLrt$PValue, method="BH")
  dLrt$q[is.na(dLrt$q)] <- 1

  dLrt$Genes <- rownames(dLrt)
  dLrt <- merge(dLrt, de.gene, by="Genes")
  dLrt <- dLrt[, -4]
  colnames(dLrt) <- c("Genes", "LFC", "logCPM", "p", "q", "de.genes")

  return(list(result=dLrt,tool.name="edgeR robust (prior.df=10)"))
}






run_edgeR_robust_auto.priorDF<-function(object) {
  #A function that runs robust edgeR GLM
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  #15 prior degrees of freedm (estimateGLMRobustDisp(.,prior.df = 15))
  require(edgeR)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  y <- DGEList(counts=count, group=cond)
  pDF <- getPriorN(y)  #Calculating the appropriate prior degrees of freedom
  y <- calcNormFactors(y)
  design <- model.matrix(~cond)
  y <- estimateGLMRobustDisp(y, design, prior.df = pDF)

  ## DE test
  dFit <- glmFit(y, design)
  dLrt <- glmLRT(dFit,coef=2)                       # Test

  dLrt <- dLrt$table
  dLrt$q <- p.adjust(dLrt$PValue, method="BH")
  dLrt$q[is.na(dLrt$q)] <- 1

  dLrt$Genes <- rownames(dLrt)
  dLrt <- merge(dLrt, de.gene, by="Genes")
  dLrt <- dLrt[, -4]
  colnames(dLrt) <- c("Genes", "LFC", "logCPM", "p", "q", "de.genes")

  return(list(result=dLrt, inputs = object, tool.name="edgeR robust (prior.df=auto)", prior.DF=pDF))
}


run_edgeR_ql<-function(object) {
  #A function that runs edgeR quasi-likelihood test
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(edgeR)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  y <- DGEList(counts=count, group=cond)
  y <- calcNormFactors(y)
  design <- model.matrix(~cond)
  y <- estimateDisp(y, design)

  fit <- glmQLFit(y, design, robust=TRUE)
  dFit <- glmQLFTest(fit)
  dLrt <- dFit$table
  dLrt$q <- p.adjust(dLrt$PValue, method="BH")
  dLrt$q[is.na(dLrt$q)] <- 1

  dLrt$Genes <- rownames(dLrt)
  dLrt <- merge(dLrt, de.gene, by="Genes")
  dLrt <- dLrt[, -4]
  colnames(dLrt) <- c("Genes", "LFC", "logCPM", "p", "q", "de.genes")

  return(list(result=dLrt, inputs = object, tool.name="edgeR QL"))
}

run_DESeq <- function(object){
  #A function that runs DESeq
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(DESeq)
  count <- object$counts                    # count matrixs
  cond  <- as.factor(object$group)                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  cds <- newCountDataSet(count, cond) ## Initialize new DESeq object
  cds <- estimateSizeFactors(cds) ## estimate size factor
  cds_norm <- cds

  ## estimate dispersion parameters
  cds  <- estimateDispersions(cds, sharingMode="maximum",
                                  method= "pooled", fitType='local')#all are default
  res <- nbinomTest(cds, levels(cond)[1], levels(cond)[2]) ## differential expression
  res <- res[, c(1, 2, 6, 7, 8)]
  colnames(res) <- c("Genes", "Mean", "LFC", "p", "q")
  res$q[is.na(res$q)] <- 1

  res <- merge(res, de.gene, by="Genes")
  return(list(result=res, inputs = object, tool.name="DESeq"))
}

run_DESeq2_indFilterDeactive<- function(object) {
  #A function that runs DESeq
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  #The eindependent filtering is diabled results(., independentFiltering = F)
  require(DESeq2)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator


  colData <- data.frame(cond)
  rownames(colData) <- colnames(count)
  d_deseq2 <- DESeqDataSetFromMatrix(countData=count, colData=colData, design=~cond)

  # Differential expression analysis
  d_deseq2 <- DESeq(d_deseq2)
  res  <- results(d_deseq2, independentFiltering = FALSE)
  res <- as.data.frame(res)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,5,6,7)]
  colnames(res) <- c("Mean", "LFC", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="DESeq2 (ind.filter=FALSE)"))
}


run_DESeq2_indFilterDeactive_l <- function(object) {
  #A function that runs DESeq
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  #The eindependent filtering is diabled results(., independentFiltering = F)
  require(DESeq2)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator


  colData <- data.frame(cond)
  rownames(colData) <- colnames(count)
  d_deseq2 <- DESeqDataSetFromMatrix(countData=count, colData=colData, design=~cond)

  # Differential expression analysis
  d_deseq2 <- DESeq(d_deseq2)
  res  <- results(d_deseq2, independentFiltering = FALSE)
  res <- as.data.frame(res)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,5,6,7)]
  colnames(res) <- c("Mean", "LFC", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, tool.name="DESeq2 (ind.filter=FALSE)"))
}




run_DESeq2_indFilterActive<- function(object) {
  #A function that runs DESeq
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  #The independent filtering is enabled results(., independentFiltering = T)
  require(DESeq2)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator


  colData <- data.frame(cond)
  rownames(colData) <- colnames(count)
  d_deseq2 <- DESeqDataSetFromMatrix(countData=count, colData=colData, design=~cond)

  # Differential expression analysis
  d_deseq2 <- DESeq(d_deseq2)
  res  <- results(d_deseq2)
  res <- as.data.frame(res)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,5,6,7)]
  colnames(res) <- c("Mean", "LFC", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="DESeq2 (ind.filter=TRUE)"))
}


run_DESeq2_indFilterDeactive_CooksOff<- function(object) {
  #A function that runs DESeq
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  #The eindependent filtering is diabled results(., independentFiltering = F)
  require(DESeq2)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator


  colData <- data.frame(cond)
  rownames(colData) <- colnames(count)
  d_deseq2 <- DESeqDataSetFromMatrix(countData=count, colData=colData, design=~cond)

  # Differential expression analysis
  d_deseq2 <- DESeq(d_deseq2)
  res  <- results(d_deseq2, independentFiltering = FALSE, cooksCutoff=FALSE)
  res <- as.data.frame(res)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,5,6,7)]
  colnames(res) <- c("Mean", "LFC", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="DESeq2 (ind.filter=FALSE, cooksCutoff=Inf)"))
}

run_limmaVoom<- function(object){
  #A function that runs limmaVoom
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(limma)
  require(edgeR)

  #limma voom with least square estimatetion lmFit(., method="ls")
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  # Create model matrix and set up dge object
  design <- model.matrix(~cond)
  dgel <- DGEList(counts=count)

  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)

  # Voom transformation - estimate dispersion
  v <- voom(dgel, design, plot=FALSE)

  # Test for differential expression
  fit <- lmFit(v,design, method="ls")  #default
  fit <- eBayes(fit, trend = FALSE, robust = FALSE)  #default

  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")

  res <- as.data.frame(de_table)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,4,5,7)]
  colnames(res) <- c("LFC", "Mean", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="limmaVoom (LSE)"))
}

run_limmaVoom_l<- function(object){
  #A function that runs limmaVoom
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(limma)
  require(edgeR)

  #limma voom with least square estimatetion lmFit(., method="ls")
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  # Create model matrix and set up dge object
  design <- model.matrix(~cond)
  dgel <- DGEList(counts=count)

  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)

  # Voom transformation - estimate dispersion
  v <- voom(dgel, design, plot=FALSE)

  # Test for differential expression
  fit <- lmFit(v,design, method="ls")  #default
  fit <- eBayes(fit, trend = FALSE, robust = FALSE)  #default

  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")

  res <- as.data.frame(de_table)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,4,5,7)]
  colnames(res) <- c("LFC", "Mean", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, tool.name="limmaVoom (LSE)"))
}


run_limmaVoom_robust<- function(object){
  #A function that runs limmaVoom
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(limma)
  require(edgeR)

  #limma voom with least square estimatetion lmFit(., method="robust")
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  # Create model matrix and set up dge object
  design <- model.matrix(~cond)
  dgel <- DGEList(counts=count)

  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)

  # Voom transformation - estimate dispersion
  v <- voom(dgel, design, plot=FALSE)

  # Test for differential expression
  fit <- lmFit(v,design, method="ls")  #default
  fit <- eBayes(fit, robust = TRUE, trend = FALSE)

  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")

  res <- as.data.frame(de_table)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,4,5,7)]
  colnames(res) <- c("LFC", "Mean", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="limmaVoom (Robust)"))
}

run_limmaTrended <- function(object){
  #A function that runs limmaVoom
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(limma)
  require(edgeR)

  #limma trended with least square estimatetion lmFit(., method="ls")

  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  # Create model matrix and set up dge object
  design <- model.matrix(~cond)
  dgel   <- DGEList(counts=count)

  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)

  y <- new("EList")
  y$E <- edgeR::cpm(dgel, log = TRUE, prior.count = 3)

  fit <- lmFit(y, design, method="ls")  #default
  fit <- eBayes(fit, trend = TRUE, robust = FALSE)

  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none", adjust.method = "BH")

  res <- as.data.frame(de_table)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,4,5,7)]
  colnames(res) <- c("LFC", "Mean", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="limmaTrended (LSE)"))
}
run_limmaTrended_robust <- function(object){
  #A function that runs limmaVoom
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(limma)
  require(edgeR)

  #limma trended with least square estimatetion lmFit(., method="ls")

  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  # Create model matrix and set up dge object
  design <- model.matrix(~cond)
  dgel   <- DGEList(counts=count)

  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)

  y <- new("EList")
  y$E <- edgeR::cpm(dgel, log = TRUE, prior.count = 3)

  fit <- lmFit(y, design, method="ls")  #default
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none", adjust.method = "BH")

  res <- as.data.frame(de_table)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,4,5,7)]
  colnames(res) <- c("LFC", "Mean", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="limmaTrended (Robust)"))
}

run_limmaVoom_QW <- function(object){
  #A function that runs limmaVoom + quality weight
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(limma)
  require(edgeR)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator

  # Create model matrix and set up dge object
  design <- model.matrix(~cond)
  dgel <- DGEList(counts=count)

  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)

  # Voom transformation with sample quality weights - estimate dispersion
  v <-  voomWithQualityWeights(dgel,design,normalization="none",plot=FALSE)

  # Test for differential expression
  fit <- lmFit(v,design)
  fit <- eBayes(fit)

  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")

  res <- as.data.frame(de_table)
  res$Genes <- rownames(res)
  res <- res[,c(1,2,4,5,7)]
  colnames(res) <- c("LFC", "Mean", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="limmaVoom QW"))
}

run_limmaVst <- function(object){
  #A function that runs limma+variance stablizing transformation
  # @param object contains non normalized count matrices, grouping variable and and indicator variable of true DE genes

  require(limma)
  require(DESeq)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator


  ## Initialize design matrix and new DESeq object
  design <- model.matrix(~cond)
  cds <- newCountDataSet(count, cond)
  cds <- estimateSizeFactors(cds) ## Estimate size factors

  ## Estimate dispersion parameters
  cds <- estimateDispersions(cds, method='blind', fitType='local') # Options set as in Soneson research

  ## Apply variance stabilizing transformation
  # own function that adapts getVarianceStabilizedData function
  vst <- DESeq::getVarianceStabilizedData(cds)

  # Limma test for differential expression
  fit <- lmFit(vst, design)
  fit <- eBayes(fit)
  res <- topTable(fit,coef=ncol(design),n=nrow(vst),sort.by="none")

  res <- res[, c(1, 2,4, 5)]
  res$Genes <- rownames(res)
  colnames(res) <- c("LFC", "Mean", "p", "q", "Genes")
  res$q[is.na(res$q)] <- 1
  res <- merge(res, de.gene, by="Genes")

  return(list(result=res, inputs = object, tool.name="limmaVST"))
}

calc.perf.metrics <- function(q, null,  thrld){
  #This function calculates the four performance comparison metrics: TPR, TNR, FPR, FNR, and FDR
  #TPR
  if(sum(null==1) > 0){
    tpr <- sum(q<thrld & null == 1)/sum(null == 1)
  }
  else { tpr <- NA }

  #TNR
  if(sum(null==0) > 0){
    tnr <- sum(q >= thrld & null == 0)/sum(null == 0)
  }
  else { tnr <- NA }

  #FPR
  if(sum(null==0) > 0){
    fpr <- sum(q<thrld & null == 0)/sum(null == 0)
  }
  else { fpr <- NA }

  #FNR
  if(sum(null == 1) >0){
    fnr <- sum(q>=thrld & null == 1)/sum(null == 1)
  }
  else {  fnr <- NA }

  #FDR
  if(sum(q<thrld) > 0){
    fdr <- sum(q<thrld & null == 0)/sum(q<thrld)
  }
  else { fdr <- 0 } #by definition of FDR (BH95)

  return(list(TPR=tpr, TNR=tnr, FPR=fpr, FNR=fnr, FDR=fdr))
}

run_rope <- function(object) {
  require(edgeR)
  library(roper)
  count <- object$counts                    # count matrixs
  cond  <- object$group                     # Grouping variable
  de.gene <- as.data.frame(object$de.gene)  # a data frame for DE genes abd biotype indicator
  design <- model.matrix(~cond)
  rope.fit <- rope(datmat = count, X_model = design,x_PI_idx = 2)
  rope.res <- rope.fit
  rope.res$q[is.na(rope.res$q)] <- 1

  rope.res$Genes <- rownames(rope.res)
  rope.res <- merge(rope.res,de.gene,by = "Genes")
  rope.res <- rope.res[,-c(3,5)]
  colnames(rope.res) <- c("Genes", "LFC", "adj", "p", "q", "de.genes")
  return(list(result=rope.res, tool.name="rope"))
}

