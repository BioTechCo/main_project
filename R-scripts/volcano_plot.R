# nolint start
library(ggplot2)
library(ggrepel)
dTN<-read.csv("result/DMP_max.csv")
threshold <- 0.35
# dTN$Change = as.factor(
#     ifelse(dTN$adj.P.Val < 0.05 & abs(dTN$logFC) > threshold,
#         ifelse(dTN$logFC > threshold , ifelse(dTN$feature == 'TSS200','TSS200',ifelse(dTN$feature == 'TSS1500','TSS1500','Hyper')), 'Hypo'),
#     'NS')
# )
dTN$Change = as.factor(
    ifelse(dTN$adj.P.Val < 0.05 & abs(dTN$logFC) > threshold,
        ifelse(dTN$logFC > threshold , 'Hyper', 'Hypo'),
    'NS')
)
write.csv(dTN, "result/DMP_logFC.csv", row.names = FALSE)
table(dTN$Change) 
dTN$label <- ifelse(dTN$adj.P.Val< 0.0005 & abs(dTN$logFC) >= 1,dTN$Gene.symbol,"")

ggplot(dTN, aes(x=logFC, y=-log10(adj.P.Val),color=Change)) + 
    geom_point(alpha=0.4, size=2) +  
    theme_bw(base_size = 12) +  
    xlab("logFC") + 
    ylab("-Log10(adj.P.adj)") + 
    scale_colour_manual(values = c('steelblue','brown','gray','green','#9b50a0')) + 
    geom_hline(yintercept = -log10(0.05), lty = 4) + 
    geom_vline(xintercept = c(-threshold, threshold), lty = 4) + 
    geom_label_repel(data = dTN, aes(label = label),
                   size = 3,box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10000) +
    theme(plot.margin = margin(t = 1, r = 1, b = 10, l = 1, unit = "in"),
    plot.background = element_rect(fill = "#ffffff"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 150,hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right")
# nolint end
