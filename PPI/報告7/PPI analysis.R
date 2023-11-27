BiocManager::install(c("STRINGdb","igraph"),ask = F,update = F)
library(STRINGdb)
######################### 選擇STRINGdb類型 #########################
string_db <- STRINGdb$new( version="11.5",
                           species=9606,   #人9606 
                           score_threshold=700, #蛋白互作的得分 默認400, 低150，高700，極高900
                           input_directory="") #可自己導入數據


data <- read.csv("1926Bio_pvalue_ttest.csv")
dat_map <- string_db$map(my_data_frame=data, 
                         my_data_frame_id_col_names="gene", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = TRUE )
hits <- dat_map$STRING_id

## PPI
png("string_PPI.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits)
dev.off()
# 存檔
img <- png::readPNG("string_PPI.png")
png::writePNG(img, "string_PPI.png")


## PPI_halo  #给PPI添加上下調信息
# 過濾p值並添加颜色列（下調基因為綠色，上調基因為红色）
dat_map_color <- string_db$add_diff_exp_color(subset(dat_map, P.Value<0.01),
                                              logFcColStr="logFC" )
payload_id <- string_db$post_payload(dat_map_color$STRING_id,
                                     colors=dat_map_color$color)
png("string_PPI_halo.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits, payload_id=payload_id )
dev.off()
# 存檔
img <- png::readPNG("string_PPI_halo.png")
png::writePNG(img, "string_PPI_halo.png")

## iGraph clustering 互作網絡分簇
#algorithm: fastgreedy(默認), walktrap, edge.betweenness
clustersList <- string_db$get_clusters(string_ids = hits ,
                                       algorithm  = "fastgreedy" ) 
# 分6群
png("string_PPI_iGraph_cluster.png",units="in",width = 15,height = 10,res=400)
par(mfrow=c(2,5))
for(i in 1:10){
  string_db$plot_network(clustersList[[i]])
}
dev.off()
# 存檔
img <- png::readPNG("string_PPI_iGraph_cluster.png")
png::writePNG(img, "string_PPI_iGraph_cluster.png")

############################## 獲取蛋白互作信息用於後續可視化 ###############
library(magrittr)
library(dplyr)
dat_link <- string_db$get_interactions(hits)
# 轉換 stringID 為 gene symbol
dat_link$from <- dat_map[match(dat_link$from,dat_map$STRING_id),'gene']
dat_link$to <- dat_map[match(dat_link$to,dat_map$STRING_id),'gene']  
colnames(dat_link) <- c('node1','node2','combined_score')
# 去除重複
dat_link <- dat_link %>% distinct(node1, node2, .keep_all = T)
write.csv(dat_link,'string_link.csv',row.names = F,quote = F)
