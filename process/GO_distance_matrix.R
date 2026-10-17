# Define the base directory
library(glue)
library(blogdown)
cancer_type <- readline("Enter the cancer type: ")
config <- read_toml(glue("config/{cancer_type}.toml"))
base_out_dir <- config$base_out_dir

# Construct the file paths using the base directory
input_file <- file.path(config$input_file)
output_bp_csv <- file.path(base_out_dir, "distance_matrix_bp.csv")
output_cc_csv <- file.path(base_out_dir, "distance_matrix_cc.csv")
output_mf_csv <- file.path(base_out_dir, "distance_matrix_mf.csv")
output_terms_count_csv <- file.path(base_out_dir, "terms_count.csv")

dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)

source("R-scripts/utils.R")


# read the input file
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
gene_list <- data$gene
entrez_ids <- gene_mapper(gene_list, "SYMBOL", "ENTREZID")
nan_and_na_count <- nan_count(entrez_ids, "Entrez IDs")


if (nan_and_na_count > 0) {
    # replace missing Entrez IDs with official gene names if any
    official_names <- fetch_official_gene_names(data$gene[is.na(entrez_ids)])
    cat("Official names for missing Entrez IDs(source: NCBI Gene): \n")
    print(official_names)
    gene_list_modified <- modify_gene_list(
        gene_list, entrez_ids, official_names
    )

    # map gene symbols to Entrez IDs again
    entrez_ids_modified <- gene_mapper(gene_list_modified, "SYMBOL", "ENTREZID")
    nan_count(entrez_ids_modified, "Modified Entrez IDs")
}

# load GO data for BP, CC, and MF
d_all <- lapply(c("BP", "CC", "MF"), function(ont) {
    godata(OrgDb  = "org.Hs.eg.db", ont = ont, computeIC = TRUE)
})

# a list of count of GO terms for each ontology
count <- lapply(d_all, function(godata_obj) {
    go_terms <- names(godata_obj@IC) # Extract GO terms from the IC slot
    return(length(go_terms)) # Return the number of GO terms
})

# Write count of count of GO terms for each ontology to a CSV file
count_vector <- c("BP" = count[[1]], "CC" = count[[2]], "MF" = count[[3]])
count_df <- data.frame(term = names(count_vector), count = count_vector)
write.csv(count_df, file = output_terms_count_csv, row.names = FALSE)

# calculate similarity matrices
if (nan_and_na_count > 0) {
    sim_matrices <- lapply(
        d_all,
        calculate_similarity,
        ids = entrez_ids_modified
    )
} else {
    sim_matrices <- lapply(d_all, calculate_similarity, ids = entrez_ids)
}
sim_matrix_bp <- sim_matrices[[1]]
sim_matrix_cc <- sim_matrices[[2]]
sim_matrix_mf <- sim_matrices[[3]]

# transform Entrez IDs to gene symbols using mapIds
gene_symbols_bp <- gene_mapper(colnames(sim_matrix_bp), "ENTREZID", "SYMBOL")
gene_symbols_cc <- gene_mapper(colnames(sim_matrix_cc), "ENTREZID", "SYMBOL")
gene_symbols_mf <- gene_mapper(colnames(sim_matrix_mf), "ENTREZID", "SYMBOL")


# replace official gene names with missing Entrez IDs if any
if (nan_and_na_count > 0) {
    gene_symbols_bp_restored <- restore_gene_symbols(
        gene_symbols_bp, official_names
    )
    gene_symbols_cc_restored <- restore_gene_symbols(
        gene_symbols_cc, official_names
    )
    gene_symbols_mf_restored <- restore_gene_symbols(
        gene_symbols_mf, official_names
    )
    # update row and column names of similarity matrices
    dimnames(sim_matrix_bp) <- list(
        gene_symbols_bp_restored,
        gene_symbols_bp_restored
    )
    dimnames(sim_matrix_cc) <- list(
        gene_symbols_cc_restored,
        gene_symbols_cc_restored
    )
    dimnames(sim_matrix_mf) <- list(
        gene_symbols_mf_restored,
        gene_symbols_mf_restored
    )
} else {
    # update row and column names of similarity matrices
    dimnames(sim_matrix_bp) <- list(gene_symbols_bp, gene_symbols_bp)
    dimnames(sim_matrix_cc) <- list(gene_symbols_cc, gene_symbols_cc)
    dimnames(sim_matrix_mf) <- list(gene_symbols_mf, gene_symbols_mf)
}

# calculate the distance matrices
distance_matrix_bp <- 1 - sim_matrix_bp
distance_matrix_cc <- 1 - sim_matrix_cc
distance_matrix_mf <- 1 - sim_matrix_mf


# print the number of missing values before and after calculating similarity
if (nan_and_na_count > 0) {
    nan_count(entrez_ids_modified, "Modified Entrez IDs")
    nan_count(gene_symbols_bp, "Gene symbols (BP)")
    nan_count(gene_symbols_cc, "Gene symbols (CC)")
    nan_count(gene_symbols_mf, "Gene symbols (MF)")
    swapped_entrez_ids <- setNames(
        names(entrez_ids_modified),
        entrez_ids_modified
    )

    # difference between swapped Entrez IDs and restored gene symbols
    diff_bp <- setdiff(swapped_entrez_ids, gene_symbols_bp)
    diff_cc <- setdiff(swapped_entrez_ids, gene_symbols_cc)
    diff_mf <- setdiff(swapped_entrez_ids, gene_symbols_mf)
    print(diff_bp)
    print(diff_cc)
    print(diff_mf)
} else {
    nan_count(entrez_ids, "Entrez IDs")
    nan_count(gene_symbols_bp, "Gene symbols (BP)")
    nan_count(gene_symbols_cc, "Gene symbols (CC)")
    nan_count(gene_symbols_mf, "Gene symbols (MF)")
}

# print columns of entrez_ids
write.csv(distance_matrix_bp, file = output_bp_csv, row.names = TRUE)
write.csv(distance_matrix_cc, file = output_cc_csv, row.names = TRUE)
write.csv(distance_matrix_mf, file = output_mf_csv, row.names = TRUE)
