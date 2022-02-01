generate.count <- function(par.list){
  require(SimSeq)

  if(is.null(rownames(par.list$counts.source))) {stop("Genes ID/name is not provided.")}
  if(length(rownames(par.list$counts.source)) != length(unique(rownames(par.list$counts.source))))
  {stop("Genes ID/name are not unique.")}

  #Preparing inputs
  counts.s   <- par.list$counts.source
  group.s    <- par.list$group.source
  k.ind      <- par.list$k.ind
  n.genes    <- par.list$n.genes
  n.initial  <- par.list$n.initial
  n.diff     <- round(par.list$p.diff*n.initial)
  nf         <- par.list$nf
  wghts      <- par.list$wghts

  #Generating counts
  counts.simseq.list <- SimData(counts = counts.s, treatment = group.s, replic = NULL,
                                sort.method = "unpaired", k.ind = k.ind, n.genes = n.initial,
                                n.diff = n.diff, norm.factors = nf, weights = wghts, switch.trt = F)

  counts.simseq <- counts.simseq.list$counts    # Simulated Count matrix from SimSeq
  genes.samp <- counts.simseq.list$genes.subset # Genes sampled from source matrix
  de.genes <- counts.simseq.list$DE.genes       # DE genes sampled from source matrix
  ee.genes <- genes.samp[ ! (genes.samp %in% de.genes) ] # EE genes sampled from source matrix
  samp.col <- counts.simseq.list$col # Columns sampled in SimSeq algorithm
  de.genes.sim <- counts.simseq.list$genes.subset %in% de.genes # logical vector giving which genes are DE in simulted matrix
  de.genes.sim <- data.frame(Genes = rownames(counts.simseq), de.genes = ifelse(de.genes.sim,1,0))
  trt <- counts.simseq.list$treatment

  #filtering genes with no expression (at least 1 expression each group)
  keep <- names(which(rowSums(counts.simseq[, c(1:k.ind)]) > 0 &
                        rowSums(counts.simseq[, c((k.ind+1):(2*k.ind))])>0))

  out.object <- list(counts=counts.simseq,
                     group=trt,
                     de.gene=de.genes.sim,
                     setting=list(k.ind = k.ind, n.genes=n.genes, p.diff=p.diff))

  return(out.object)
}
