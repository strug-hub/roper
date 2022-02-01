source("./generate_count_fun.R")
# source('./np_sim_gtex.R')
# processing_GTEx_data(g_size = 200)

require(SimSeq)
require(fdrtool)
library(edgeR)
library(doParallel)

n.itration = 50
read <- readRDS(file = "./dat_Sim/GTEx_Data_400and400.RData")


counts.source <-
  read$counts

group.source <- read$group
counts.source <-
  counts.source[which(rowSums(counts.source[, group.source == levels(group.source)[1]]) > 0 &
                        rowSums(counts.source[, group.source == levels(group.source)[2]]) > 0), ]
dim(counts.source)


lib.sizes <- apply(counts.source, 2, sum)
nf <-
  calcNormFactors(counts.source) * lib.sizes              #using edgeR package

## Compute weights to sample DE genes in SimData function using edgeR normalization method TMM
probs <-
  CalcPvalWilcox(
    counts.source,
    treatment = group.source,
    replic = NULL,
    sort.method = "unpaired",
    sorted = TRUE,
    nf,
    exact = FALSE
  )

wghts <-
  1 - fdrtool(probs,
              statistic = "pvalue",
              plot = FALSE,
              verbose = FALSE)$lfd

# saveRDS(list(probs=probs, weights=wghts), file="weights.RData", ascii=TRUE)


#Simulation Parameter setting

# PDE <- c(0, 0.1,0.2,0.3) #proportion of true DE genes
# PDE <- c(0,0.1,0.2,0.3)
PDE <- c(0.2)


for (pde in 1:length(PDE)) {
  n.genes <- 15000
  p.diff <- PDE[pde]  #0, 0.01,  0.05, 0.1, 0.2, 0.3
  n.initial <-
    18000 #number of genes to be simulated before filtration

  par.setting <-
    list(
      list(
        counts.sources = counts.source,
        group.source = group.source,
        nf = nf,
        k.ind = 30,
        n.initial = n.initial,
        n.genes = n.genes,
        p.diff = p.diff,
        wghts = wghts
      ),
      list(
        counts.sources = counts.source,
        group.source = group.source,
        nf = nf,
        k.ind = 50,
        n.initial = n.initial,
        n.genes = n.genes,
        p.diff = p.diff,
        wghts = wghts
      ),
      list(
        counts.sources = counts.source,
        group.source = group.source,
        nf = nf,
        k.ind = 100,
        n.initial = n.initial,
        n.genes = n.genes,
        p.diff = p.diff,
        wghts = wghts
      ),
      list(
        counts.sources = counts.source,
        group.source = group.source,
        nf = nf,
        k.ind = 150,
        n.initial = n.initial,
        n.genes = n.genes,
        p.diff = p.diff,
        wghts = wghts
      ),
      list(
        counts.sources = counts.source,
        group.source = group.source,
        nf = nf,
        k.ind = 200,
        n.initial = n.initial,
        n.genes = n.genes,
        p.diff = p.diff,
        wghts = wghts
      )
    )

  n.samples <- length(par.setting)
  par.list <- rep(par.setting, each = n.itration)

  set.seed(123)
  sim.counts <- lapply(par.list, FUN = generate.count)
  rm(par.list)
  # saveRDS(sim.counts, file=paste0("./sim_counts", 100*p.diff, "_L_study_PDE.RData"))


  source("./DE_tools_np.R")
  cl <- 3
  cl <- makeCluster(cl)

  res.edgeR.glm                         = parLapplyLB(cl, sim.counts, run_edgeR_glm_l)
  print(paste("... ... edgeR GLM...completed"))

  res.edgeR.robust.10priorDF            = parLapplyLB(cl, sim.counts, run_edgeR_robust_10priorDF_l)
  print(paste("... ... edgeR rob 10pDE...completed"))

  res.DESeq2.indFilterDeactive          = parLapplyLB(cl, sim.counts, run_DESeq2_indFilterDeactive_l)
  print(paste("... ... DESeq2 I ...completed"))

  res.limmaVoom                         = parLapplyLB(cl, sim.counts, run_limmaVoom_l)
  print(paste("... ... limma voom ...completed"))

  res.rope = parLapplyLB(cl, sim.counts, run_rope)
  print(paste("... ... rope ...completed"))
  stopCluster(cl)

  results <- list(
    res.edgeR.glm                         = res.edgeR.glm,
    res.edgeR.robust.10priorDF = res.edgeR.robust.10priorDF,
    res.DESeq2.indFilterDeactive = res.DESeq2.indFilterDeactive,
    res.limmaVoom = res.limmaVoom,
    res.rope = res.rope
  )
  saveRDS(results, file = paste0("./res_Sim/npsim_result.v15k.RData"))
}
