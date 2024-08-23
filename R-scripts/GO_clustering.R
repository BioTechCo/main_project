library(GOSemSim)
library(org.Hs.eg.db)

input_file <- "breast/result/GSE89093_nc/train100/dbeta_hyper_TSS_0.02.csv"
cluster_num <- 5
output_pdf <- "breast/result/GSE89093_nc/train100/dendrogram.pdf"
output_csv <- "breast/result/GSE89093_nc/train100/distance_matrix.csv"

# read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = data$gene,
  keytype = "SYMBOL",
  column = "ENTREZID"
)

# remove NA values
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# load GO data for biological process (BP), cellular component (CC), and molecular function (MF)
d_all <- lapply(c("BP", "CC", "MF"), function(ont) {
  godata(annoDb = "org.Hs.eg.db", ont = ont, computeIC = FALSE)
})

# calculate semantic similarity matrices
sim_matrices <- lapply(d_all, function(d) {
  mgeneSim(entrez_ids, semData = d, drop = "", measure = "Wang", verbose = TRUE)
})
sim_matrix_bp <- sim_matrices[[1]]
sim_matrix_cc <- sim_matrices[[2]]
sim_matrix_mf <- sim_matrices[[3]]

# gather all unique genes
all_genes <- unique(
  c(
    row.names(sim_matrix_bp),
    row.names(sim_matrix_cc),
    row.names(sim_matrix_mf)
  )
)

# create a combined matrix
combined_matrix <- matrix(NA, length(all_genes), length(all_genes))
rownames(combined_matrix) <- all_genes
colnames(combined_matrix) <- all_genes

# fill the combined matrix with the similarity matrices
fill_matrix <- function(main_matrix, sub_matrix) {
  rows <- rownames(sub_matrix)
  cols <- colnames(sub_matrix)
  main_matrix[rows, cols] <- sub_matrix
  return(main_matrix)
}
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_bp)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_cc)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_mf)

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

# calculate the distance matrix
distance_matrix <- as.dist(1 - combined_matrix)

# check for NA, NaN, and Infinite values
cat("Number of NA values: ", sum(is.na(distance_matrix)), "\n")
cat("Number of NaN values: ", sum(is.nan(distance_matrix)), "\n")
cat("Number of Infinite values: ", sum(is.infinite(distance_matrix)), "\n")

# replace NA, NaN, and Infinite values with 1
distance_matrix[is.na(distance_matrix) | is.nan(distance_matrix) |
    is.infinite(distance_matrix)
] <- 1

# write the combined matrix to a CSV file
write.csv(1 - combined_matrix, file = output_csv, row.names = TRUE)
