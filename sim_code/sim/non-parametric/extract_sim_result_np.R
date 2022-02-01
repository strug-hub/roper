require(SimSeq)
require(fdrtool)
library(edgeR)
library(doParallel)
library(DESeq2)
n.samples = 5
n.itration = 50

source("./DE_tools_np")

sim20_result <- readRDS("./res_Sim/npsim_result.v15k.RData")

results <- c(sim20_result)

# saveRDS(results, file=paste0("./rs3/final_DE_result", 100*p.diff, "_small_study_PDE.RData"))

#-------------------------------------------------------------------------
#Compute performance measures
print(paste(".... summarizing analyzed simulated counts"))

perf.metrics.list <- list()
fdr.thrld= seq(0, 1, 0.005)

DE.result <- results
perf.metrics.sub.list <- list()
for(j in 1:length(fdr.thrld)){
  print(j)
  thrld <- fdr.thrld[j]
  perf.metrics <- lapply(DE.result, function(x){
    t(sapply(x, function(y){
      res <- y$result
      sumr.mRNA <- calc.perf.metrics(res$q, res$de.genes, thrld = thrld)
      FNR.mRNA <- sumr.mRNA$FNR
      TPR.mRNA <- sumr.mRNA$TPR
      TNR.mRNA <- sumr.mRNA$TNR
      FPR.mRNA <- sumr.mRNA$FPR
      FDR.mRNA <- sumr.mRNA$FDR

      rep.size <- y$inputs$setting$k.ind
      prop.DE  <- y$inputs$setting$p.diff
      # biotype.composition <- unlist(y$inputs$gene.biotype.composition)
      c(FNR.mRNA=FNR.mRNA, TPR.mRNA=TPR.mRNA, TNR.mRNA=TNR.mRNA,FPR.mRNA=FPR.mRNA,
        FDR.mRNA=FDR.mRNA,
        rep.size=rep.size, prop.DE=prop.DE)
    }))
  })
  perf.metrics2 <- as.data.frame(do.call("rbind", perf.metrics))
  perf.metrics2$DE.tool  <- rep(names(perf.metrics), each=n.itration*n.samples)
  perf.metrics2$alpha    <- thrld
  perf.metrics.sub.list[j] <- list(perf.metrics2)
}
perf.metrics.sub.list <- do.call("rbind", perf.metrics.sub.list)


# saveRDS(perf.metrics.sub.list, paste0("./rs3/performance.results", 100*p.diff, "_small_study_PDE.RData"))
print(paste(".... done."))



# Create Figures

sim.res.all <- perf.metrics.sub.list
# sim.res.all[sim.res.all$DE.tool == "",]$DE.tool = "rope"
sim.res.all[sim.res.all$DE.tool == "res.DESeq2.indFilterDeactive",]$DE.tool = "DESeq2 (setting1)"
sim.res.all[sim.res.all$DE.tool == "res.edgeR.glm",]$DE.tool = "edgeR"
sim.res.all[sim.res.all$DE.tool == "res.edgeR.robust.10priorDF",]$DE.tool = "edgeR (Robust)"
sim.res.all[sim.res.all$DE.tool == "res.limmaVoom",]$DE.tool = "limma (Voom)"
sim.res.all[sim.res.all$DE.tool == "res.rope",]$DE.tool = "rope"


table(sim.res.all$DE.tool)

saveRDS(sim.res.all,file = "./res_Sim/extract_npsim_result.v15k.RData")
