library(tidyverse)
library(DESeq2)

# Read HNE Data
var.tab.v0 <- read.table("../active_infection_dge_cre20201119/active_infection_dge_covar.csv", header = T, sep = ",", row.names = 1)
var.tab.v0$rna_seq_id <- make.names(as.character(var.tab.v0$rna_seq_id))
Count.v0 <- read.table("../active_infection_dge_cre20201119/active_infection_dge_count.csv", header = T, sep = ",", row.names = 1)

# decreased by 71 - 65 = 6 samples. New size = 65
table(var.tab.v0$longi_pa_binary2)

# Sample filtering by active infection status:
var.tab.v1 <- var.tab.v0[!is.na(var.tab.v0$longi_pa_binary2), ]

# Create data object
raw.count.v0 <- Count.v0[, colnames(Count.v0) %in% var.tab.v1$rna_seq_id]
rownames(raw.count.v0) <- Count.v0$Name

# Count.v0[,names(which(colSums(raw.count.v0) == 0))]
# Remove all genes of all zeros: 9952
raw.count.v1 <- raw.count.v0[rowSums(raw.count.v0) > 0, ]
Count.v1 <- Count.v0[rowSums(raw.count.v0) > 0, ]

nrow(raw.count.v0) - nrow(raw.count.v1)

# Remove outliers
counts.source0 <- raw.count.v1

# Arrange X by infection status
var.tab.v2 <- var.tab.v1 %>% arrange(longi_pa_binary2)
counts.source <- raw.count.v1 <- counts.source0[, var.tab.v2$rna_seq_id]


# Change variable types to factors
var.tab.v3 <- var.tab.v2
fac.col <- c("unc", "female", "longi_pa_binary2")
var.tab.v3[fac.col] <- lapply(var.tab.v3[fac.col], factor)
group.source <- var.tab.v2$longi_pa_binary2


#
rownames(var.tab.v3) <- make.names(rownames(var.tab.v3))


dat.obj <- DESeqDataSetFromMatrix(countData = raw.count.v1, colData = var.tab.v3, design = ~longi_pa_binary2)

# gene filtering
counts.source2 <- counts.source[apply(counts.source, 1, function(x) {
  all(tapply(x, group.source, function(y) sum(y > 0)) > 10)
}), ]

dim(counts.source)
dim(counts.source2)

dat.obj.gf <- DESeqDataSetFromMatrix(countData = counts.source2, colData = var.tab.v3, design = ~longi_pa_binary2)

dat.list.v1 <- list(original = dat.obj, gf = dat.obj.gf)

save(dat.list.v1, file = "../vold_dat/dat_list_v1.RData")

# library(countsimQC)
# tempDir <- "../QC_report/"
#
# countsimQCReport(ddsList = dat.list.v1, outputFile = "countsim_report_v1.html",
#                  outputDir = tempDir, outputFormat = "html_document",
#                  showCode = FALSE, forceOverwrite = TRUE,
#                  savePlots = FALSE, description = "This is my QC summary.",
#                  maxNForCorr = 27, maxNForDisp = Inf,
#                  calculateStatistics = F, subsampleSize = 27,
#                  kfrac = 0.01, kmin = 5,
#                  permutationPvalues = FALSE, nPermutations = NULL)
