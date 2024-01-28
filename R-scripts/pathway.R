

library(BiocManager)
library(clusterProfiler)
library(org.Hs.eg.db)


data <- read.csv("../result/result_basic_dbeta_0.35.csv")

gene_symbols <- data[,3]
print(gene_symbols)
entrez_ids <- mapIds(org.Hs.eg.db, keys=gene_symbols, keytype="SYMBOL", column="ENTREZID")
sorted_gene_list <- sort(entrez_ids, decreasing = TRUE)

# 進行 KEGG Pathway 富集分析  
# test = bitr_kegg(entrez_ids, fromType = "kegg", toType = "Path", organism = "hsa")
kk <- enrichKEGG(gene         = entrez_ids,
                 organism     = 'hsa')

# save kk as a flie
 write.csv(kk, file = "kegg_enrichment.csv")
 new = read.csv("kegg_enrichment.csv")
    new$URL = paste0("https://www.kegg.jp/kegg-bin/show_pathway?", new$ID, "/", new$geneID)
    write.csv(new, file = "kegg_enrichment.csv")
print("done")