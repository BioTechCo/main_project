#!/usr/bin/env Rscript

# this line is to mute 🐛 from vscode
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 1) {
  # default output file
  args[2] <- "out.txt"
}

# result/gene_list.txt
# INPUT_FILE = args[1]
INPUT_FILE = "result/gene_list.txt"
# CLUSTER_NUM = args[2]
CLUSTER_NUM = 5
# OUTPUT_PDF = args[3]
OUTPUT_PDF = "dendrogram.pdf"
# OUTPUT_CSV = args[4]
OUTPUT_CSV = "similarity_matrix.csv"


library(GOSemSim)
library(org.Hs.eg.db)
library(dendextend)
data <- read.csv(INPUT_FILE)
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = data$gene,
  keytype = "SYMBOL",
  column = "ENTREZID"
)
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

d_all <- lapply(c("BP", "CC", "MF"), function(ont) {
  godata("org.Hs.eg.db", ont = ont, computeIC = FALSE)
})

# calculate semantic similarity matrices
sim_matrices <- lapply(d_all, function(d) {
  mgeneSim(entrez_ids, semData = d, drop = "", measure = "Wang", verbose = TRUE)
})
sim_matrix_bp <- sim_matrices[[1]] # nolint
sim_matrix_cc <- sim_matrices[[2]]
sim_matrix_mf <- sim_matrices[[3]]

# gather all unique genes
all_genes <- unique(
  c(
    row.names(sim_matrix_BP),
    row.names(sim_matrix_CC),
    row.names(sim_matrix_MF)
  )
)
# create a combined matrix
combined_matrix <- matrix(NA, length(all_genes), length(all_genes))
rownames(combined_matrix) <- all_genes
colnames(combined_matrix) <- all_genes
fill_matrix <- function(main_matrix, sub_matrix) {
  rows <- rownames(sub_matrix)
  cols <- colnames(sub_matrix)
  main_matrix[rows, cols] <- sub_matrix
  return(main_matrix)
}
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_BP)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_CC)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_MF)

# transform Entrez IDs to gene symbols using mapIds
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = row.names(combined_matrix),
  keytype = "ENTREZID",
  column = "SYMBOL"
)
# update row and column names of combined_matrix
rownames(combined_matrix) <- gene_symbols
colnames(combined_matrix) <- gene_symbols

# write the combined matrix to a CSV file
write.csv(combined_matrix, file = OUTPUT_CSV, row.names = TRUE)

distance_matrix <- as.dist(1 - combined_matrix)

# below is the code for hierarchical clustering
cat("Number of NA values: ", sum(is.na(distance_matrix)), "\n")
cat("Number of NaN values: ", sum(is.nan(distance_matrix)), "\n")
cat("Number of Infinite values: ", sum(is.infinite(distance_matrix)), "\n")

distance_matrix[is.na(distance_matrix) | is.nan(distance_matrix) |
    is.infinite(distance_matrix)
] <- 1
# hierarchical clustering using ward.D2 method
hc <- hclust(distance_matrix, method = "ward.D2")

# convert the hclust object to dendrogram
dend <- as.dendrogram(hc)

# plot the dendrogram
number_of_cluster <- CLUSTER_NUM
clusters <- cutree(hc, number_of_cluster)
colors_vector <- c("red", "blue", "green", "yellow", "purple")

pdf(OUTPUT_PDF, width = 5, height = 10)
colored_dend <- dend %>%
  set("labels_col", colors_vector[clusters]) %>%
  set("branches_k_color", k = number_of_cluster) %>%
  set("labels_cex", 0.5) %>%
  plot(horiz = TRUE, main = "Gene Distance Tree with Colored Clusters")
dev.off()
