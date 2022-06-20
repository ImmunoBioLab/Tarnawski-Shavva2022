#Graph-based clustering
#Original authors were using Seurat, so they probably used this as well.
ssData$infomap <- clusterCells(ssData, use.dimred = "PCA", BLUSPARAM = NNGraphParam(cluster.fun = "infomap")) %>% factor()
ssData$walktrap <- clusterCells(ssData, use.dimred = "PCA", BLUSPARAM = SNNGraphParam(k = 10, type = "rank", cluster.fun = "walktrap")) %>% factor()
ssData$louvain <- clusterCells(ssData, use.dimred = "PCA", BLUSPARAM = NNGraphParam(k = 7, type = "jaccard", cluster.fun = "louvain")) %>% factor()
ssData$eigen <- clusterCells(ssData, use.dimred = "PCA", BLUSPARAM = NNGraphParam(cluster.fun = "leading_eigen")) %>% factor()

set.seed(100)
ssData$kmeans <- clusterCells(ssData, use.dimred="PCA", BLUSPARAM = KmeansParam(centers = 8))

#Display all clustering ways on UMAP
pdf(file.path(techDir,  "UMAP_clusters.pdf"), width = 12, height = 12)
gridExtra::grid.arrange(
  plotReducedDim(ssData, "UMAP", colour_by = "subtype"),
  plotReducedDim(ssData, "UMAP", colour_by = "infomap"),
  plotReducedDim(ssData, "UMAP", colour_by = "walktrap"),
  plotReducedDim(ssData, "UMAP", colour_by = "louvain"),
  plotReducedDim(ssData, "UMAP", colour_by = "eigen"),
  plotReducedDim(ssData, "UMAP", colour_by = "kmeans"),
  ncol = 3, nrow = 2
)
dev.off()

#Show ChATness distribution in Louvain clusters
pdf(file.path(figDir,  "ClusterChatnessDimRed.pdf"), width = 10, height = 5)
gridExtra::grid.arrange(
  plotReducedDim(ssData, dimred = "UMAP", colour_by = "louvain"),
  plotReducedDim(ssData, dimred = "UMAP", colour_by = "ChATness"),
  nrow = 1, ncol = 2
)
dev.off()

#Show ChAT expression per cells in Louvain clusters
ssData$ChATExpr <- logcounts(ssData) %>% .["CHAT",]

pdf(file.path(figDir,  "ClusterExprChatnessDimRed.pdf"), width = 10, height = 5)
gridExtra::grid.arrange(
  plotReducedDim(ssData, dimred = "UMAP", colour_by = "louvain"),
  plotReducedDim(ssData, dimred = "UMAP", colour_by = "ChATExpr"),
  nrow = 1, ncol = 2
)
dev.off()

#Cluster spread by celltype
celltypeByCluster <- table(ssData$subtype, ssData$louvain)

pdf(file.path(techDir, "CelltypeDistributionCluster.pdf"), width = 12, height = 9)
pheatmap(log2(celltypeByCluster+1), color=viridis::viridis(101))[[4]]
dev.off()

pdf(file.path(techDir, "CelltypeDistributionPatient.pdf"), width = 12, height = 9)
pheatmap(log2(celltypeByPatient+1), color=viridis::viridis(101))[[4]]
dev.off()

#Calculate percentage of ChAT+ cells per cluster
chatCelltype <- colData(ssData) %>% .[, c("louvain", "ChATness")] %>%
  split(., .$louvain) %>%
  lapply(., function(cluster) {
    data.frame(Cluster = unique(cluster$louvain),
               ChATness = c("Positive", "Negative"),
               Count = c(sum(cluster$ChATness), nrow(cluster) - sum(cluster$ChATness)),
               Total = nrow(cluster)
    )
  }) %>%
  do.call("rbind", .)
rownames(chatCelltype) <- NULL

totalCells <- apply(celltypeByCluster, 2, sum) %>%
  data.frame(Cluster = names(.),
             Count = .)
rownames(totalCells) <- NULL

totalCells %<>% .[match(unique(chatCelltype$Cluster), totalCells$Cluster),]
totalCells$ChATPercentage <- round(chatCelltype[seq(1, nrow(chatCelltype), 2), "Count"]/chatCelltype[seq(1, nrow(chatCelltype), 2), "Total"]*100, 0)
totalCells %<>% dplyr::arrange(., desc(ChATPercentage))
totalCells$Cluster %<>% factor(., levels = unique(.))

#Plot it as a bargraph
pdf(file.path(figDir, "barChATnessCluster.pdf"), width = 8, height = 6)
ggplot() +
  geom_bar(data = totalCells, mapping = aes(x = Cluster, y = ChATPercentage), stat = "identity", fill = viridis::viridis(9)) +
  geom_text(data = totalCells, mapping = aes(x = Cluster, y = ChATPercentage, label = ChATPercentage), vjust = -0.5, size = 5) +
  ylab("% of ChAT+ cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, size = 12, hjust = 1),
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14),
        axis.title.x = element_blank())
dev.off()

#Violin plot of ChAT expression per cluster
pdf(file.path(figDir, "violinChATnessCluster.pdf"), width = 6, height = 3)
plotExpression(ssData, exprs_values = "tpms", features = "CHAT", x = "louvain", colour_by = "louvain") +
  theme(axis.title.x = element_blank())
dev.off()

#Remove unneceesary objects
rm(totalCells)
rm(chatCelltype)
rm(celltypeByCluster)
rm(celltypeByPatient)

#Marker detection
markers <- scoreMarkers(ssData, ssData$louvain)

#Select top 5 markers per cluster (ordered by mean AUC and expressed in every cell of the cluster)
allMarkers <- lapply(markers, function(cluster) {
  genes <- cluster %>% .[order(.$mean.AUC, decreasing = TRUE),] %>% .[.$self.detected == 1,] %>% rownames(.)
  genes <- genes[1:5]
}) 

annotation_row <- data.frame(Cluster = as.character(rep(c(1:9), each = 5)))
rownames(annotation_row) <- unlist(allMarkers)

#Plot a heatmap of top 5 markers of each sluter
pdf(file.path(figDir, "clusterMarkersHeatmap.pdf"), width = 5, height = 8)
plotGroupedHeatmap(ssData, features = unlist(allMarkers), group = "louvain", center = TRUE, zlim = c(-3, 3), 
                   cluster_rows = FALSE, annotation_row = annotation_row)[[4]]
dev.off()

#Save sheets with cluster markers per cluster
wb <- createWorkbook()
for(i in seq_along(markers)) {
  addWorksheet(wb, str_c("Cluster", names(markers)[i], sep = "_"))
  writeData(wb, i, as.data.frame(markers[[i]][order(markers[[i]]$mean.AUC, decreasing = TRUE),]), startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
}
saveWorkbook(wb, file.path(getwd(), "VSS", "Clusters_&_Markers.xlsx"), TRUE)

wb <- createWorkbook()
for(i in seq_along(markers)) {
  geneNames <- allMarkers[[i]]
  
  geneDf <- markers[[i]] %>% .[rownames(.) %in% geneNames,] %>% .[order(.$mean.AUC, decreasing = TRUE),]
  
  addWorksheet(wb, str_c("Cluster", names(markers)[i], sep = "_"))
  
  writeData(wb, i, as.data.frame(geneDf), startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
}
saveWorkbook(wb, file.path(getwd(), "VSS", "Clusters_&_Markers_Heatmap.xlsx"), TRUE)
rm(wb)
rm(annotation_row)
rm(allMarkers)
rm(markers)