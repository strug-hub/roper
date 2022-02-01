library(tidyverse)
library(edgeR)
library(DESeq2)

# Load pre-processed clinical and RNA-Seq counts
load("../vold_dat/dat_list_v1.RData")
Count.v0 <-
  read.table(
    "../active_infection_dge_cre20201119/active_infection_dge_count.csv",
    header = T,
    sep = ",",
    row.names = 1
  )


# HNE analysis: DESeq2 ----------------------------------------------------

obj.ana3 <- dat.list.v1$gf
# Design matrix
design(obj.ana3) <-
  ~ pc_1 + pc_2 + pc_3 + unc + female + rin + ptprc_tmm + imm_cell_prop +
    longi_pa_binary2

obj.ana3.de2 <- DESeq(obj.ana3)
res.de2.v0 <- results(obj.ana3.de2)
res.de2 <- res.de2.v0[order(res.de2.v0$pvalue), ]
# DESeq2 Results
res.de2$Genes <-
  as.character(Count.v0[match(rownames(res.de2), Count.v0$Name), ]$Description)
res.de2_m3 <- res.de2

# HNE analysis: edgeR ----------------------------------------------------
y <-
  DGEList(
    counts = counts(obj.ana3),
    group = obj.ana3$longi_pa_binary2
  )
y <- calcNormFactors(y)
design <-
  model.matrix(
    ~ obj.ana3$pc_1 + obj.ana3$pc_2 + obj.ana3$pc_3 + obj.ana3$unc +
      obj.ana3$female + obj.ana3$rin + obj.ana3$ptprc_tmm +
      obj.ana3$imm_cell_prop + obj.ana3$longi_pa_binary2
  )
y <- estimateDisp(y, design)

dFit <- glmFit(y, design)
dLrt_0 <- glmLRT(dFit)

dLrt <- dLrt_0$table
dLrt$q <- p.adjust(dLrt$PValue, method = "BH")
dLrt$q[is.na(dLrt$q)] <- 1

dLrt <- dLrt[order(dLrt$PValue), ]
dLrt$Genes <-
  as.character(Count.v0[match(rownames(dLrt), Count.v0$Name), ]$Description)
edgeR_m3 <- dLrt

# HNE analysis: voom ----------------------------------------------------
count <- counts(obj.ana3)
cond <-
  obj.ana3$longi_pa_binary3

design <-
  model.matrix(
    ~ obj.ana3$pc_1 + obj.ana3$pc_2 + obj.ana3$pc_3 + obj.ana3$unc +
      obj.ana3$female + obj.ana3$rin + obj.ana3$ptprc_tmm +
      obj.ana3$imm_cell_prop + obj.ana3$longi_pa_binary2
  )

dgel <- DGEList(counts = count)

# Calculate normalization factors - TMM
dgel <- calcNormFactors(dgel)

# Voom transformation - estimate dispersion
v <- voom(dgel, design, plot = FALSE)

# Test for differential expression
fit <- lmFit(v, design, method = "ls") # default
fit <- eBayes(fit, trend = FALSE, robust = FALSE) # default

# Return top table
de_table <-
  topTable(fit,
    coef = ncol(design),
    n = nrow(dgel),
    sort.by = "P"
  )
de_table$Genes <-
  as.character(Count.v0[match(rownames(de_table), Count.v0$Name), ]$Description)
head(de_table)
lv_m3 <- de_table

# save(
#   res.de2_m3,
#   edgeR_m3,
#   lv_m3,
#   file = "../Results/PsAresult_othV0.RData"
# )

# HNE analysis: RoPE ----------------------------------------------------
library(roper)
obj.ana.rplr <- dat.list.v1$gf
raw_count <- counts(obj.ana.rplr)

col_var <- as.data.frame(obj.ana.rplr@colData)
x_PI <- as.numeric(col_var$longi_pa_binary2) - 1
col_var$x_PI <- x_PI

X_model_m3 <-
  model.matrix(~ pc_1 + pc_2 + pc_3 + unc + female + rin + ptprc_tmm + imm_cell_prop +
    x_PI,
  data = col_var
  )

RoPE_raw_res <-
  rope(
    datmat = raw_count,
    X_model = X_model_m3,
    x_PI_idx = dim(X_model_m3)[2]
  )

table_rplr_m3 <- RoPE_raw_res[order(RoPE_raw_res$pvals), ]
table_rplr_m3$Genes <-
  as.character(Count.v0[match(rownames(table_rplr_m3), Count.v0$Name), ]$Description)
RoPE_m3 <- table_rplr_m3

# save(RoPE_m3, file = "../Results/PsAresult_RoPE_m3.RData")


# Results in manuscript ---------------------------------------------------
# Figure 3: Standardized Profile Likelihood for SLC9A3----------------------------------------------------------------
source("./rope_gridadj.R")
SLC9A3_RPLR <- rope_gridadj(
  datmat = raw_count, ngi = match("ENSG00000066230.11", rownames(raw_count)),
  X_model = X_model_m3
)


g_dat_9A3 <-
  data.frame(
    x = SLC9A3_RPLR$grid,
    y = exp(SLC9A3_RPLR$log.norm.lik.adj),
    yori = exp(SLC9A3_RPLR$log.norm.lik)
  )

col_new <- c("Original" = "black", "Adjusted" = "blue")

g_9A3 <-
  ggplot(data = g_dat_9A3, aes(x = x, y = y)) +
  geom_line(aes(col = "Adjusted")) +
  geom_line(aes(x = x, y = yori, col = "Original")) +
  theme_bw() +
  labs(x = "log fold change (PsA infection vs no infection)", y = "Std Profile Likelihood") +
  geom_text(x = -2.5, y = 1, label = "MLE = -1.2705") +
  geom_text(x = -2.5, y = 0.9, label = "adjfact fact = 0.0011") +
  geom_text(x = -2.5, y = 0.8, label = "LR = exp(11.1033)") +
  geom_vline(xintercept = 0, col = "red") +
  scale_colour_manual(name = "Std Likelihood: ", values = col_new) +
  theme(legend.position = "bottom", text = element_text(size = 11))
g_9A3


# Figure 4: Venn diagram of the top 20 ranked DE genes----------------------------------------------------------------
deg_de3 <- res.de2_m3$Genes[1:20]
deg_er3 <- edgeR_m3$Genes[1:20]
deg_lv3 <- lv_m3$Genes[1:20]
deg_rplr3 <- table_rplr_m3$Genes[1:20]

library(RAM)

group.venn(
  list(
    DESeq2 = deg_de3,
    edgeR = deg_er3,
    voom = deg_lv3,
    RoPE = deg_rplr3
  ),
  label = TRUE,
  cat.pos = c(0, 0, 0, 0),
  lab.cex = 1.1
)


# Figure S9 Boxplot of protein-coding DE gene expression--------
library(grex)
data("gtexv7")

id <- gtexv7[]
df <- grex(id)

rplr.sig.tab <- table_rplr_m3 %>% filter(padj < 0.05)
# Check RopE genes
rplr.sig.tab$Genes %in% df$hgnc_symbol

rplr.sig.tab$Genes[!(rplr.sig.tab$Genes %in% df$hgnc_symbol)]

rplr.sig.tab.pc <- rplr.sig.tab[rplr.sig.tab$Genes %in% df$hgnc_symbol, ]

df.RoPE <- df[match(rplr.sig.tab.pc$Genes, df$hgnc_symbol), ]

# Choose protein coding genes
RoPE.top.pc <- df.RoPE[df.RoPE$gene_biotype == "protein_coding", ]$hgnc_symbol

# sig DE boxplots ---------------------------------------------------------

rplr.sig.tab <- rplr.sig.tab[rplr.sig.tab$Genes %in% RoPE.top.pc, ]



count_sig <- log(raw_count[match(rownames(rplr.sig.tab), rownames(raw_count)), ] + 0.1) %>%
  t() %>%
  as_tibble() %>%
  rename_with(~ rplr.sig.tab$Genes) %>%
  mutate(PsA = obj.ana.rplr$longi_pa_binary2)


dat_group_box <-
  count_sig %>%
  pivot_longer(!PsA, names_to = "gene", values_to = "log.count") %>%
  mutate(gene = factor(gene, levels = rplr.sig.tab$Genes)) %>%
  mutate(Rank = factor(match(gene, table_rplr_m3$Genes)))



pp.n <- ggplot(data = dat_group_box, aes(x = gene, y = log.count, fill = PsA)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, size = 7.5, vjust = 0.5, hjust = 1))
pp.n.p <- pp.n + labs(fill = "Pseudomonas Infection") + scale_fill_manual(labels = c("0" = "No", "1" = "Yes"), values = c("blue", "red"))


pp.n.p + ggsci::scale_fill_npg()
