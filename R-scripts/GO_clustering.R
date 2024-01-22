# nolint start
# 載入需要的套件
library(GOSemSim)
library(org.Hs.eg.db)
library(stats)
library(dendextend)
# 讀取CSV檔案並將基因名稱(symbol)轉換成Entrez Gene IDs
data <- read.csv("result/gene_list.txt")
entrez_ids <- mapIds(org.Hs.eg.db, keys=data$gene, keytype="SYMBOL", column="ENTREZID")
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# 創建GO資料物件，考慮所有三個GO分支
d_all <- lapply(c("BP", "CC", "MF"), function(ont) {
    godata('org.Hs.eg.db', ont=ont, computeIC=FALSE)
})

# 計算基於所有三個GO分支的語義相似度
sim_matrices <- lapply(d_all, function(d) {
    mgeneSim(entrez_ids, semData=d,drop="", measure="Wang",verbose = TRUE)
})

# 現在有三個矩陣：基於BP、CC和MF的語義相似度
sim_matrix_BP <- sim_matrices[[1]]
sim_matrix_CC <- sim_matrices[[2]]
sim_matrix_MF <- sim_matrices[[3]]

# 獲取所有的基因名稱
all_genes <- unique(c(row.names(sim_matrix_BP), row.names(sim_matrix_CC), row.names(sim_matrix_MF)))

# 初始化一個全為NA的矩陣
combined_matrix <- matrix(NA, length(all_genes), length(all_genes))
rownames(combined_matrix) <- all_genes
colnames(combined_matrix) <- all_genes

# 填充每個相似度矩陣的值
fill_matrix <- function(main_matrix, sub_matrix) {
  rows <- rownames(sub_matrix)
  cols <- colnames(sub_matrix)
  main_matrix[rows, cols] <- sub_matrix
  return(main_matrix)
}

combined_matrix <- fill_matrix(combined_matrix, sim_matrix_BP)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_CC)
combined_matrix <- fill_matrix(combined_matrix, sim_matrix_MF)

# 如果你希望將NA值替換為0：
# combined_matrix[is.na(combined_matrix)] <- 0

# 使用mapIds轉換Entrez IDs為基因符號
gene_symbols <- mapIds(org.Hs.eg.db, keys=row.names(combined_matrix), keytype="ENTREZID", column="SYMBOL")

# 更新combined_matrix的行名和列名
rownames(combined_matrix) <- gene_symbols
colnames(combined_matrix) <- gene_symbols

# 如果你想保存這個矩陣：
write.csv(combined_matrix, file="combined_matrix_symbols.csv", row.names=TRUE)


# 計算距離矩陣
distance_matrix <- as.dist(1-combined_matrix)

# 檢查是否存在NA、NaN或Inf的值
sum(is.na(distance_matrix))
sum(is.nan(distance_matrix))
sum(is.infinite(distance_matrix))

# 將非正常值替換為1（表示最大距離）
distance_matrix[is.na(distance_matrix) | is.nan(distance_matrix) | is.infinite(distance_matrix)] <- 1

# 進行層次聚類
hc <- hclust(distance_matrix, method="ward.D2") # 這裡我使用的是average linkage方法，但你可以根據需要選擇其他方法

# 將hclust物件轉換為dendrogram物件
dend <- as.dendrogram(hc)

# 使用cutree將層次聚類結果切成k個集群
k <- 5 # 例如，我們切成3個集群

clusters <- cutree(hc, k)# clusters矩陣現在包含每個基因的集群分配

# 定義三種顏色
colors_vector <- c("red", "blue", "green", "yellow", "purple")

# 使用color_branches著色dendrogram
colored_dend <- dend %>% 
    set("labels_col", colors_vector[clusters]) %>%
    set("branches_k_color", k = k) %>%
    set("labels_cex", 0.5)

# 打開一個更高的圖形裝置
pdf("Dendrogram.pdf", width=10, height=80)

# 繪製樹狀圖
colored_dend %>%
  set("labels_cex", 0.5) %>%  # 減小標籤大小
  plot(horiz = TRUE, main = "Gene Distance Tree with Colored Clusters") 

#如果pdf印不出東西表示這個寫錯了，ylim = c(0, 0.5)去掉後就可以跑出來了

#colored_dend %>%
#set("labels_cex", 0.5) %>%  # 減小標籤大小
#plot(horiz = TRUE, main = "Gene Distance Tree with Colored Clusters")
#目前使用這個跑得出來

dev.off()

# nolint end

