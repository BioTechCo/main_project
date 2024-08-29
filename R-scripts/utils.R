library(GOSemSim)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(rentrez)


nan_count <- function(x, data_name) {
    cat(
        "Number of (NA | NaN | non-NA/NaN) values in ", data_name, ":",
        sum(is.na(x)), " | ",
        sum(is.nan(x)), " | ", sum(!is.na(x) & !is.nan(x)),
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


# Remove genes without Entrez IDs and append official names
modify_gene_list <- function(gene_list, entrez_ids, official_names) {
    gene_list_modified <- gene_list[!is.na(entrez_ids)]
    gene_list_modified <- c(gene_list_modified, unname(official_names))
    gene_list_modified <- gene_list_modified[!is.na(gene_list_modified)]
    return(gene_list_modified)
}


# Restore gene symbols with official names if available
restore_gene_symbols <- function(gene_symbols, official_names) {
    gene_symbols_restored <- ifelse(
        gene_symbols %in% official_names,
        names(official_names)[match(gene_symbols, official_names)],
        gene_symbols
    )
    gene_symbols_with_name <- setNames(
        gene_symbols_restored,
        names(gene_symbols)
    )
    return(gene_symbols_with_name)
}
