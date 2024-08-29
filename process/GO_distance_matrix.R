input_file <- "breast/result/GSE89093_nc/train100/dbeta_hyper_TSS_0.02.csv"
output_csv <- "breast/result/GSE89093_nc/train100/distance_matrix.csv"


source("R-scripts/utils.R")


# read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
gene_list <- data$gene
entrez_ids <- gene_mapper(gene_list, "SYMBOL", "ENTREZID")
nan_count(entrez_ids, "Entrez IDs")


# replace missing Entrez IDs with official gene names if any
official_names <- fetch_official_gene_names(data$gene[is.na(entrez_ids)])
cat("Official names for missing Entrez IDs(source: NCBI Gene): \n")
print(official_names)
gene_list_modified <- modify_gene_list(gene_list, entrez_ids, official_names)


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
dimnames(combined_matrix) <- list(all_genes, all_genes)


# fill the combined matrix with the similarity matrices
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_bp)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_cc)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_mf)


# transform Entrez IDs to gene symbols using mapIds
gene_symbols <- gene_mapper(row.names(combined_matrix), "ENTREZID", "SYMBOL")


# replace official gene names with missing Entrez IDs if any
gene_symbols_with_name <- restore_gene_symbols(gene_symbols, official_names)


# update row and column names of combined_matrix
dimnames(combined_matrix) <- list(
    gene_symbols_with_name,
    gene_symbols_with_name
)


# calculate the distance matrix
distance_matrix <- 1 - combined_matrix


# print the number of missing values before and after calculating similarity
nan_count(entrez_ids_modified, "Modified Entrez IDs")
nan_count(gene_symbols, "Gene symbols")
swapped_entrez_ids <- setNames(names(entrez_ids_modified), entrez_ids_modified)
diff <- setdiff(swapped_entrez_ids, gene_symbols)
print(diff)


# print columns of entrez_ids
write.csv(distance_matrix, file = output_csv, row.names = TRUE)
