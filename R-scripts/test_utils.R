# to run the tests, run `testthat::test_dir("R-scripts")`
source("utils.R")


library(testthat)


# Test fetch_official_gene_names
test_that("fetch_official_gene_names retrieves official names correctly", {
    gene_symbols <- c("TP53", "BRCA1", "INVALIDGENE")
    official_names <- fetch_official_gene_names(gene_symbols)
    expect_length(official_names, length(gene_symbols))
    expect_false(is.na(official_names[1])) # TP53 should have an official name
    expect_false(is.na(official_names[2])) # BRCA1 should have an official name
    expect_true(is.na(official_names[3])) # INVALIDGENE should return NA
})


# Test gene_mapper
test_that("gene_mapper maps gene symbols to Entrez IDs correctly", {
    gene_symbols <- c("TP53", "BRCA1")
    entrez_ids <- gene_mapper(gene_symbols, "SYMBOL", "ENTREZID")
    expect_length(entrez_ids, length(gene_symbols))
    expect_true(all(!is.na(entrez_ids)))
    expect_equal(names(entrez_ids), gene_symbols)
})


# Test fill_matrix
test_that("fill_matrix fills the main matrix correctly", {
    main_matrix <- matrix(0, nrow = 3, ncol = 3)
    rownames(main_matrix) <- c("gene1", "gene2", "gene3")
    colnames(main_matrix) <- c("gene1", "gene2", "gene3")
    sub_matrix <- matrix(1, nrow = 2, ncol = 2)
    rownames(sub_matrix) <- c("gene1", "gene2")
    colnames(sub_matrix) <- c("gene1", "gene2")
    filled_matrix <- fill_matrix(main_matrix, sub_matrix)
    for (i in 1:3) {
        for (j in 1:3) {
            if (i <= 2 & j <= 2) {
                expect_equal(filled_matrix[i, j], 1)
            } else {
                expect_equal(filled_matrix[i, j], 0)
            }
        }
    }
})

# Test calculate_similarity
test_that("fill_matrix fills the main matrix correctly", {
    main_matrix1 <- matrix(NA, nrow = 4, ncol = 4)
    main_matrix2 <- matrix(NA, nrow = 4, ncol = 4)
    main_matrix3 <- matrix(NA, nrow = 4, ncol = 4)
    genes <- c("gene1", "gene2", "gene3", "gene4")
    gene1 <- c("gene1", "gene2")
    gene2 <- c("gene1", "gene3")
    gene3 <- c("gene1", "gene2", "gene4")

    dimnames(main_matrix1) <- list(genes, genes)
    sub_matrix1 <- matrix(1, nrow = 2, ncol = 2)
    dimnames(sub_matrix1) <- list(gene1, gene1)
    main_matrix1 <- fill_matrix(main_matrix1, sub_matrix1)
    dimnames(main_matrix2) <- list(genes, genes)
    sub_matrix2 <- matrix(2, nrow = 2, ncol = 2)
    dimnames(sub_matrix2) <- list(gene2, gene2)
    main_matrix2 <- fill_matrix(main_matrix2, sub_matrix2)
    dimnames(main_matrix3) <- list(genes, genes)
    sub_matrix3 <- matrix(3, nrow = 3, ncol = 3)
    dimnames(sub_matrix3) <- list(gene3, gene3)
    main_matrix3 <- fill_matrix(main_matrix3, sub_matrix3)
    stacked_array <- abind(main_matrix1, main_matrix2, main_matrix3, along = 3)
    average_matrix <- average_non_zero(stacked_array)
    result_matrix <- matrix(c(
        2, 2, 2, 3, 2, 2, 0, 3, 2, 0, 2, 0, 3, 3, 0, 3
    ), nrow = 4, ncol = 4)
    for (i in 1:4) {
        for (j in 1:4) {
            expect_equal(average_matrix[i, j], result_matrix[i, j])
        }
    }
})

# Test modify_gene_list
test_that("modify_gene_list modifies gene list correctly", {
    gene_list <- c("ATF7IP", "C17orf97")
    entrez_ids <- c("ATF7IP" = "55729", "C17orf97" = NA, "ASXL2" = "0000")
    official_names <- c("C17orf97" = "LIAT1")
    modified_list <- modify_gene_list(gene_list, entrez_ids, official_names)
    expect_equal(modified_list, c(
        "ATF7IP", "LIAT1"
    ))
})


# Test restore_gene_symbols
test_that("restore_gene_symbols restores gene symbols correctly", {
    gene_symbols <- c("71" = "ACTG1", "155" = "ADRB3", "55821" = "ALLC")
    official_names <- c(
        "ACTG1_name_missing_in_db" = "ACTG1",
        "ADRB3_name_missing_in_db" = "ADRB3"
    )
    gene_symbols
    official_names
    restored_symbols <- restore_gene_symbols(gene_symbols, official_names)
    restored_symbols
    expect_equal(restored_symbols, c(
        "71" = "ACTG1_name_missing_in_db",
        "155" = "ADRB3_name_missing_in_db",
        "55821" = "ALLC"
    ))
})
