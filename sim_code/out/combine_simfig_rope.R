library(tidyverse)
library(ggsci)

load("../sim/parametric/results_simulateOutliers_nb_v6.RData")

# check_sample_sizes_parametric ------------------------------------------------------
unique(res_nb_v6$m)
# parametric --------------------------------------------------------------

res <- res_nb_v6
res$algorithm <-
  factor(res$algorithm,
    levels = c("RoPE", "edgeR", "edgeR-robust", "DESeq2", "limma-voom")
  )

res$m <- factor(res$m)
levels(res$m) <- paste0("m=", levels(res$m))
res$percentOutlier <- 100 * res$percentOutlier
res$percentOutlier <- factor(res$percentOutlier)
levels(res$percentOutlier) <-
  paste0(levels(res$percentOutlier), "% outlier")


resSensPadj <- res[res$senspadj < 0.1, ]
resSensPadj <- resSensPadj[nrow(resSensPadj):1, ]
resSensPadj <-
  resSensPadj[!duplicated(with(resSensPadj, paste(
    algorithm, m,
    percentOutlier
  ))), ]

plot_resSensPadj <-
  resSensPadj[resSensPadj$percentOutlier == "0% outlier", ]

res <- droplevels(res)
plot_dat_para <- res[res$percentOutlier == "0% outlier", ]
plot_dat_para <- subset(plot_dat_para, select = -percentOutlier)
plot_dat_para$sim <- rep("parametric", nrow(plot_dat_para))
plot_dat_para

# non-paramatric ----------------------------------------------------------

sim.res.all <-
  readRDS("../sim/non-parametric/res_Sim/extract_npsim_result.v15k.RData")

it.1 <- c(rep(30, 50), rep(50, 50), rep(100, 50), rep(150, 50), rep(200, 50))
sim.res.all$rep.size <- rep(it.1, nrow(sim.res.all) / length(it.1))


tools0 <- unique(sim.res.all$DE.tool)

tools <- tools0
tools.name <-
  c(
    "edgeR",
    "edgeR (Robust)",
    "DESeq2 (setting1)",
    "limma (Voom)",
    "RoPE"
  )
cols <- c("black", "blue", "forestgreen", "darkorange3", "red")

head(sim.res.all)

sim.res.nonp <- as_tibble(sim.res.all)
sim.res.nonp
sim.res.nonp.summary <-
  sim.res.nonp %>%
  group_by(alpha, DE.tool, rep.size) %>%
  summarise(
    FNR = median(FNR.mRNA),
    TPR = median(TPR.mRNA),
    TNR = median(TNR.mRNA),
    FPR = median(FPR.mRNA),
    FDR = median(FDR.mRNA)
  )
sim.res.nonp.summary


# combine_par_np ----------------------------------------------------------
sim.p <- as_tibble(plot_dat_para)
sim.np <- ungroup(sim.res.nonp.summary)
# sample size
levels(sim.p$m)[levels(sim.p$m) == "m=60"] <- "n/2=30"
levels(sim.p$m)[levels(sim.p$m) == "m=100"] <- "n/2=50"
levels(sim.p$m)[levels(sim.p$m) == "m=200"] <- "n/2=100"
levels(sim.p$m)[levels(sim.p$m) == "m=300"] <- "n/2=150"
levels(sim.p$m)[levels(sim.p$m) == "m=400"] <- "n/2=200"


colnames(sim.p)[colnames(sim.p) == "m"] <- "n.size"
sim.p

sim.np$rep.size <- factor(sim.np$rep.size)

levels(sim.np$rep.size)[levels(sim.np$rep.size) == "30"] <- "n/2=30"
levels(sim.np$rep.size)[levels(sim.np$rep.size) == "50"] <- "n/2=50"
levels(sim.np$rep.size)[levels(sim.np$rep.size) == "100"] <-
  "n/2=100"
levels(sim.np$rep.size)[levels(sim.np$rep.size) == "150"] <-
  "n/2=150"
levels(sim.np$rep.size)[levels(sim.np$rep.size) == "200"] <-
  "n/2=200"

colnames(sim.np)[colnames(sim.np) == "rep.size"] <- "n.size"

sim.np

# algorithm
sim.p$algorithm
unique(sim.np$DE.tool)
sim.np[sim.np$DE.tool == "DESeq2 (setting1)", ]$DE.tool <- "DESeq2"
sim.np[sim.np$DE.tool == "edgeR", ]$DE.tool <- "edgeR"
sim.np[sim.np$DE.tool == "edgeR (Robust)", ]$DE.tool <-
  "edgeR-robust"
sim.np[sim.np$DE.tool == "limma (Voom)", ]$DE.tool <- "limma-voom"
sim.np[sim.np$DE.tool == "RoPE", ]$DE.tool <- "RoPE"
sim.np$DE.tool <-
  factor(sim.np$DE.tool,
    levels = c("RoPE", "edgeR", "edgeR-robust", "DESeq2", "limma-voom")
  )
sim.np$DE.tool

sim.np$sim <- "non-parametric"

names(sim.p)
names(sim.np)

sim.np <-
  sim.np %>% dplyr::rename(
    algorithm = DE.tool,
    sensitivity = TPR,
    oneminusspec = FPR,
    oneminusprec = FDR
  )

sim.np$senspadj <- sim.np$alpha
sim.np$precpadj <- sim.np$alpha
sim.np <- sim.np %>% select(-FNR, -TNR, -alpha)
sim.p <- sim.p %>% select(-pvals)

sim.all <- dplyr::union(sim.p, sim.np)
sim.all$sim <-
  factor(sim.all$sim, levels = c("parametric", "non-parametric"))
sim.all$sim
sim.all

resSensPadj <- sim.all[sim.all$senspadj < 0.1, ]
resSensPadj <- resSensPadj[nrow(resSensPadj):1, ]
resSensPadj <-
  resSensPadj[!duplicated(with(resSensPadj, paste(
    algorithm, n.size,
    sim
  ))), ]
resSensPadj


# plot --------------------------------------------------------------------

tp2 <- ggplot(
  sim.all,
  aes(x = oneminusprec, y = sensitivity, color = algorithm)
)

tp2.p <- tp2 + scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3)) + scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 1)) + geom_line(size = 0.8) + theme_bw() + facet_grid(sim ~
n.size) +
  xlab("1 - precision (FDR)") + ylab("sensitivity (TPR)") + coord_cartesian(xlim = c(-0.03, 0.35), ylim = c(
    0.69,
    1
  )) + geom_point(
    aes(x = oneminusprec, y = sensitivity, shape = algorithm),
    size = 1.8,
    data = sim.all[sim.all$precpadj == 0.1, ]
  ) + geom_vline(xintercept = 0.1, linetype = "dashed") + theme(
    legend.position =
      "bottom"
  )



tp2.p + scale_color_npg()

tp1 <-
  ggplot(sim.all, aes(x = precpadj, y = oneminusprec, color = algorithm))

tp1.p <- tp1 + scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3)) + scale_y_continuous(breaks = c(
  0,
  0.1, 0.2, 0.3
)) + geom_line(size = 0.8) + theme_bw() + facet_grid(sim ~
n.size) +
  geom_abline(intercept = 0, slope = 1) + xlab("adjusted p-value") + ylab("1 - precision (FDR)") +
  coord_cartesian(xlim = c(-0.03, 0.35), ylim = c(-0.03, 0.35)) + geom_point(
    aes(
      x = precpadj,
      y = oneminusprec, shape = algorithm
    ),
    size = 1.8,
    data = sim.all[sim.all$precpadj == 0.1, ]
  ) + theme(legend.position = "bottom")


tp1.p + scale_color_npg()
