library(data.table)
processing_GTEx_data <- function(g_size = 50, seed = 123) {
  ###### Read in sample attributes

  temp_samp_attr <- fread("./data/GTEx_v7_Annotations_SampleAttributesDS.txt",
    header = TRUE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE, quote = "\""
  )
  ###### Read in count data

  # Read file
  input <- fread("./data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct", header = TRUE, stringsAsFactors = FALSE, sep = "\t", nrows = -1, skip = 2)


  # Replace '.' in colnames by '-'
  colnames(input) <- gsub(".", "-", colnames(input), fixed = TRUE)

  temp_samp_attr <- temp_samp_attr[temp_samp_attr$SAMPID %in% colnames(input), ]

  # Check technology used for all samples in input matrix
  table(temp_samp_attr$SMGEBTCHT) # -> ok, all samples contain TrueSeq data


  samp_sel <- temp_samp_attr[
    temp_samp_attr$SMTSD %in% c("Esophagus - Mucosa", "Lung"),
    c("SAMPID", "SMTSD")
  ]


  samp_sel <- samp_sel[order(samp_sel$SMTSD), ]
  table(samp_sel$SMTSD)

  set.seed(seed)
  indices_c1 <- which(samp_sel$SMTSD == "Esophagus - Mucosa")
  indices_c2 <- which(samp_sel$SMTSD == "Lung")

  # Select sample size per group
  # g_size <- 50
  g_size_rm <- g_size * 2 + 3
  id_samp_sel <- c(sample(indices_c1, g_size), sample(indices_c2, g_size))

  samp_sel <- samp_sel[id_samp_sel, ]

  # Assign simpler names to samples
  names(samp_sel) <- c("sampleID", "tissue")
  samp_sel$samp_name <- paste(ifelse(grepl("Esophagus", samp_sel$tissue), "A", "B"),
    c(
      1:sum(grepl("Esophagus", samp_sel$tissue)),
      1:sum(grepl("Lung", samp_sel$tissue))
    ),
    sep = ""
  )
  samp_sel$condition <- ifelse(grepl("Esophagus", samp_sel$tissue), "condA", "condB")

  input <- input[, c("Name", "Description", samp_sel$sampleID), with = F]
  colnames(input)[-c(1, 2)] <- samp_sel$samp_name[match(colnames(input)[-c(1, 2)], samp_sel$sampleID)]

  # Check whether all identifiers are unique
  input$ID <- substr(input$Name, 1, 15)
  length(input$ID) == length(unique(input$ID))

  # -> all identifiers are unique
  # Check whether all GeneIDs have a unique description
  length(input$ID) == length(unique(input$Description))


  library(grex)
  data("gtexv7")
  id <- gtexv7[]
  df <- grex(id)
  tail(df)

  # one to one matching of input and df
  input$ID != df$ensembl_id

  df$Description <- input$Description
  temp_genes <- df[, c("ensembl_id", "gene_biotype", "Description")]
  colnames(temp_genes) <- c("V1", "V2", "Description")
  # Check whether all genes are in the annotation file :: YES
  table(input$ID %in% temp_genes$V1)


  length(unique(input$Description[input$ID %in% temp_genes$V1]))
  # Subset gene data and count to only include genes that are in the annotation file
  #   AND that are protein coding
  input2 <- input[input$ID %in% temp_genes[temp_genes$V2 == "protein_coding", ]$V1, ]
  table(temp_genes$V2)

  length(unique(input2$Description))

  # Aggregate rows of genes with multiple IDs for the same description
  input3 <- aggregate(. ~ Description, input2[, -c(1, get("g_size_rm"))], sum)

  # Adding the gene type column to the counts data frame
  temp_genes2 <- temp_genes[temp_genes$V2 == "protein_coding", ]
  table(temp_genes2$V2)
  temp_genes2 <- temp_genes2[, -1]
  colnames(temp_genes2) <- c("Gene_Type", "Description")
  head(temp_genes2)
  dim(temp_genes2)
  temp_genes2 <- aggregate(temp_genes2, by = list(temp_genes2$Description), FUN = max)
  temp_genes2 <- temp_genes2[, -1]
  dim(temp_genes2)
  head(temp_genes2)

  input4 <- merge(input3, temp_genes2, by = "Description")
  head(input4)
  table(input4$Gene_Type)
  length(input4$Description) == length(unique(input4$Description))
  # -> 18770 genes

  # Create count matrix
  counts_full <- input3[, -1]
  rownames(counts_full) <- input3$Description
  genedata_full <- input3$Description

  ###### Further subset to only include genes with total rowcount >=1

  # Select IDs with average count>1
  # sel_id <- rowMeans(counts_full)>=1
  sel_id <- rowSums(counts_full) >= 1
  table(sel_id)
  # 17663 genes with average count larger than 1

  # Subset count matrix and genedata vector
  counts <- counts_full[sel_id, ]
  genedata <- genedata_full[sel_id]

  # Transform count dataframes to count matrices
  counts_full <- as.matrix(counts_full)
  counts <- as.matrix(counts)

  ###### Create some new variables

  nSamp <- ncol(counts)
  group <- as.factor(samp_sel$condition[match(colnames(counts), samp_sel$samp_name)])

  ###### Remove redundant objects
  rm(input, input2, input3, samp_sel, temp_genes, temp_samp_attr, sel_id)
  GTEx_data_full <- list(counts = counts, group = group, counts_full = counts_full)
  outpath <- paste("./dat_Sim/GTEx_Data_", g_size, "and", g_size, ".RData", sep = "")
  saveRDS(GTEx_data_full,
    file = outpath
  )
  rm(list = ls())
}

processing_GTEx_data(g_size = 400)
