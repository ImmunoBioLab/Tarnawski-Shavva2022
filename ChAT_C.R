#Correlate ChAT with everything and select top genes
corDir <- file.path(getwd(), "VSS", "CorFigures")

chatCor <- lapply(rownames(ssData), function(gene) {
  corTmp <- cor.test(tpms(ssData)["CHAT",], tpms(ssData)[gene,], method = "pearson")
  
  data.frame(Gene = gene,
             rVal = corTmp$estimate,
             pVal = corTmp$p.value
  )
}) %>% do.call("rbind", .)
chatCor %<>% dplyr::arrange(., desc(rVal))
chatCor$fdr <- p.adjust(chatCor$pVal, method = "fdr", n = nrow(chatCor)-1)
chatCor$logFDR <- log(chatCor$fdr, 10) * -1

#How many cells are top genes expressed in?
chatCor %<>% dplyr::mutate(., PosCells = tpms(ssData)[.$Gene,] %>% apply(., 1, function(gene) sum(gene > 0)))
chatCor %<>% dplyr::mutate(., PosPercent = round(.$PosCells/ncol(tpms(ssData))*100, 0))

#How many ChAT+  and ChAT- cells are top genes expressed in?
chatCor %<>% dplyr::mutate(., ChPosCells = tpms(ssData)[.$Gene, ssData$ChATness == TRUE] %>% apply(., 1, function(gene) sum(gene > 0)))
chatCor %<>% dplyr::mutate(., ChPosPercent = round(.$ChPosCells/ncol(tpms(ssData)[, ssData$ChATness == TRUE])*100, 0))
chatCor %<>% dplyr::mutate(., ChNegCells = tpms(ssData)[.$Gene, ssData$ChATness == FALSE] %>% apply(., 1, function(gene) sum(gene > 0)))
chatCor %<>% dplyr::mutate(., ChNegPercent = round(.$ChNegCells/ncol(tpms(ssData)[, ssData$ChATness == FALSE])*100, 0))
chatCor %<>% dplyr::mutate(., ChDiff = .$ChPosPercent - .$ChNegPercent)

#Save
wb <- openxlsx::createWorkbook()
addWorksheet(wb, "Correlation")
writeData(wb, 1, chatCor, startRow = 1, startCol = 1, colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb, file.path(getwd(), "VSS", "ChAT_correlation.xlsx"), TRUE)
rm(wb)

#Best cor genes: r < 0.1, posPercent < 5. negPercent < 95
#Worst cor genes: 
topCorGenes <- chatCor %>% .[.$rVal > 0.5 & .$ChPosPercent > 60 & .$ChNegPercent  < 30, "Gene"]
botCorGenes <- chatCor %>% .[.$rVal < -0.3, "Gene"] %>% .[!is.na(.)] %>% rev() %>% .[c(1:length(topCorGenes)-1)]  

intGenes <- c(topCorGenes, botCorGenes) %>% .[!. %in% "CHAT"]

#Plot correlation with ChAT for all top genes in all cells
corPlots <- lapply(seq_along(intGenes), function(index) {
  gene <- data.frame(CHATexp = tpms(ssData)["CHAT",],
                     GeneExp = tpms(ssData)[intGenes[index],]
  )
  
  pVal <- chatCor %>% .[.$Gene == intGenes[index], "fdr"] %>% formatC(., format = "e", digits = 2)
  rVal <- chatCor %>% .[.$Gene == intGenes[index], "rVal"] %>% round(., 2)
  
  ggplot(data = gene, aes(y = CHATexp, x = GeneExp)) +
    geom_point(color = "black", alpha = 0.5) +
    geom_text(data = gene, mapping = aes(x = Inf, y = 2.5, vjust = 1, hjust = 1,
                                         label = paste0("r == ", rVal)), inherit.aes = FALSE, parse = TRUE) +
    geom_text(data = gene, mapping = aes(x = Inf, y = 2.5, vjust = 2.65, hjust = 1,
                                         label = paste0("p == ", pVal)), inherit.aes = FALSE, parse = TRUE) +
    xlab(intGenes[index]) +
    ylab("CHAT") +
    theme_classic() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
          panel.grid.major = element_line(color = "gray80", size = 0.2), panel.grid.minor = element_line(color = "gray90", size = 0.1))
})

pdf(file.path(corDir, "CorGenes.pdf"), width = 15, height = 15)
ggpubr::ggarrange(plotlist = corPlots)
dev.off()

#Violin plot
pdf(file.path(corDir, "corGenesViolin.pdf"), width = 14, height = 4)
plotExpression(ssData, exprs_values = "tpms" , features = intGenes, x = "ChATness", colour_by = "louvain", ncol = 8, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), legend.position = "none")
dev.off()

#Dot plot per cell
pdf(file.path(corDir, "corGenesDotPlot.pdf"), width = 8, height = 4)
geneDotPlot(ssData, countAssay = "tpms", geneNames = c("CHAT", intGenes), orderRow = "CHAT", logTransform = TRUE, adjust = TRUE,
            intCol = "louvain", low = viridis::viridis(101)[1], high = viridis::viridis(101)[101])
dev.off()

#Psudobulk DE
summed <- ssData %>% aggregateAcrossCells(., id = colData(.)[, c("Patient", "ChATness")], use.assay.type = "logcounts")

#Creating up a DGEList object for use in edgeR:
y <- DGEList(logcounts(summed), samples = colData(summed))

#Discarding due to low cell count in a sample
discarded <- summed$ncells < 5

y <- y[,!discarded]

y <- calcNormFactors(y)

par(mfrow=c(2,3))
pdf(file.path(techDir, "DistributionExpression.pdf"), width = 4, height = 4)
for (i in seq_len(ncol(y))) {
  plotMD(y, column=i)
}
dev.off()
par(mfrow=c(1,1))

pdf(file.path(techDir, "SampleDistribution.pdf"), width = 8, height = 8)
plotMDS(cpm(y, log=TRUE), 
        col=ifelse(y$samples$ChATness, "red", "blue"))
dev.off()

design <- model.matrix(~factor(ChATness), y$samples)

y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust = TRUE)

pdf(file.path(techDir, "VariationExpression.pdf"), width = 6, height = 6)
plotBCV(y)
plotQLDisp(fit)
dev.off()

res <- glmQLFTest(fit, coef = ncol(design))

topTags(res)

res$table$log10pValue <- log10(res$table$PValue)*-1
res$table$Gene <- rownames(res$table)
res$table %<>% dplyr::arrange(., desc(log10pValue))

#Remove unnecessary objects
rm(y)
rm(fit)
rm(i)
rm(design)

topTags(res)

res$table$log10pValue <- log10(res$table$PValue)*-1
res$table$Gene <- rownames(res$table)
res$table %<>% dplyr::arrange(., desc(log10pValue))

#Volcano plot
pdf(file.path(deDir, "VolcanoPlot.pdf"), width = 8, height = 8)
EnhancedVolcano(res$table,
                lab = rownames(res$table),
                x = 'logFC',
                y = 'PValue',
                cutoffLineWidth = 0,
                FCcutoff = log2(res$table[1, "logFC"]),
                pCutoff = res$table[1, "PValue"],
                col = c("grey30", "grey30", "grey30", "red2"),
                selectLab = c("CHAT", intGenes), max.overlaps = 30, pointSize = 1,
                labSize = 2.0, drawConnectors = TRUE) +
  geom_point(data = res$table[res$table$Gene %in% c("CHAT", intGenes),], aes(x = logFC, y = log10pValue), color = "tomato")
dev.off()

#Remove unneessary objects
rm(res)
rm(summed)

#Get genes that reasonably predict ChAT
geneNames <- intGenes[1:8]

#Mark cells that have at least one of the genes
colData(ssData)$ChMarkerIn <- tpms(ssData)[geneNames, ] %>% apply(., 2, function(cell) sum(cell) > 1)
colData(ssData)$ChMarkerSum <- tpms(ssData)[geneNames, ] %>% apply(., 2, function(cell) sum(cell))

#Get a violin plot of ChAT in cells marked by top correlated genes
pdf(file.path(deDir, "MarkerChAT.pdf"), width = 3, height = 5)
plotExpression(ssData, exprs_values = "tpms", features = "CHAT", x = "ChMarkerIn", colour_by = "ChMarkerIn") +
  scale_x_discrete(labels = c("Negative", "Positive")) +
  scale_color_manual(values = c("black", "tomato")) +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_blank(), legend.position = "none")
dev.off()