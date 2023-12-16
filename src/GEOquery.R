# nolint start
library(GEOquery)
file450 <- "GPL13534_series_matrix.txt.gz"
file850 <- "GSE207998-GPL21145_series_matrix.txt.gz"

# grab the data
gse450 <- getGEO(GEO="GSE207998", filename=file450,GSEMatrix = TRUE,getGPL = FALSE)
x450 <- exprs(object = gse)
write.csv(x = x, file = "450series.csv", quote = F, row.names = F)

# can not grab the data? download it from website and read the gz file
gse850=getGEO(filename=file.path("GSE207998-GPL21145_series_matrix.txt.gz"),GSEMatrix = TRUE,getGPL = FALSE)
experimentData(gse850)
x850 <- exprs(object = gse850)
write.csv(x = x850, file = "850series.csv", quote = F, row.names = F)
# nolint end
