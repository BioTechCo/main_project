library(GOSemSim)
library(org.Hs.eg.db)
library(rentrez)

input_file <- "breast/result/GSE89093_nc/train100/dbeta_hyper_TSS_0.02.csv"
output_csv <- "breast/result/GSE89093_nc/train100/distance_matrix.csv"

# Define utility functions
# -------------------------
# inspect the number of NA and NaN values in a vector
nan_count <- function(x, data_name) {
  cat(
    "Number of (NA | NaN | non-NA/NaN) values in ", data_name, ":",
    sum(is.na(x)), " | ", sum(is.nan(x)), " | ", sum(!is.na(x) & !is.nan(x)),
    "\n"
  )
}

# fetch official gene names from NCBI Gene
fetch_official_gene_names <- function(gene_symbols) {
  official_names <- sapply(gene_symbols, function(gene) {
    search_result <- entrez_search(
      db = "gene",
      term = paste0(gene, "[sym]"),
      retmax = 1
    )
    if (length(search_result$ids) > 0) {
      gene_details <- entrez_summary(
        db = "gene",
        id = search_result$ids[1]
      )
      return(gene_details$name)
    } else {
      return(NA)
    }
  })
  return(official_names)
}

# map gene symbols to Entrez IDs or vice versa
gene_mapper <- function(keys, keytype, column) {
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = keys,
    keytype = keytype,
    column = column
  )
  return(entrez_ids)
}

# calculate similarity between genes
calculate_similarity <- function(d, ids) {
  mgeneSim(ids, semData = d, drop = "", measure = "Wang", verbose = TRUE)
}

# fill the combined matrix with the similarity matrices
fill_matrix <- function(main_matrix, sub_matrix) {
  rows <- rownames(sub_matrix)
  cols <- colnames(sub_matrix)
  main_matrix[rows, cols] <- sub_matrix
  return(main_matrix)
}

# Main script
# -----------
# read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
gene_list <- data$gene
entrez_ids <- gene_mapper(gene_list, "SYMBOL", "ENTREZID")
nan_count(entrez_ids, "Entrez IDs")

# replace missing Entrez IDs with official gene names if any
official_names <- fetch_official_gene_names(data$gene[is.na(entrez_ids)])
cat("Official names for missing Entrez IDs(source: NCBI Gene): \n")
print(official_names)
gene_list_modified <- gene_list[!is.na(entrez_ids)]
gene_list_modified <- c(gene_list_modified, unname(official_names))
gene_list_modified <- gene_list_modified[!is.na(gene_list_modified)]

# map gene symbols to Entrez IDs again
entrez_ids_modified <- gene_mapper(gene_list_modified, "SYMBOL", "ENTREZID")
nan_count(entrez_ids_modified, "Modified Entrez IDs")

# load GO data for BP, CC, and MF
d_all <- lapply(c("BP", "CC", "MF"), function(ont) {
  godata(annoDb = "org.Hs.eg.db", ont = ont, computeIC = FALSE)
})

# calculate similarity matrices
sim_matrices <- lapply(d_all, calculate_similarity, ids = entrez_ids)
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

combined_matrix <- fill_matrix(combined_matrix, sim_matrix_bp)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_cc)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_mf)

# transform Entrez IDs to gene symbols using mapIds
gene_symbols <- gene_mapper(row.names(combined_matrix), "ENTREZID", "SYMBOL")

# replace official gene names with missing Entrez IDs if any
official_names_vector <- official_names
gene_symbols_restored <- ifelse(
  gene_symbols %in% official_names_vector,
  names(official_names_vector)[match(gene_symbols, official_names_vector)],
  gene_symbols
)
gene_symbols_with_name <- setNames(gene_symbols_restored, names(gene_symbols))

# update row and column names of combined_matrix
rownames(combined_matrix) <- gene_symbols_with_name
colnames(combined_matrix) <- gene_symbols_with_name

# calculate the distance matrix
distance_matrix <- 1 - combined_matrix

nan_count(entrez_ids_modified, "Modified Entrez IDs")
nan_count(gene_symbols, "Gene symbols")
swapped_entrez_ids <- setNames(names(entrez_ids_modified), entrez_ids_modified)
diff <- setdiff(swapped_entrez_ids, gene_symbols)
print(diff)

# print columns of entrez_ids
write.csv(distance_matrix, file = output_csv, row.names = TRUE)
