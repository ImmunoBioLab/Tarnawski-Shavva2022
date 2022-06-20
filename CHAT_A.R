library(openxlsx)
library(magrittr)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(pheatmap)
library(SingleR)
library(edgeR)
library(EnhancedVolcano)
library(biomaRt)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(bluster)

#DotPlot for Patients/Genes
#Show differential expression of genes per patient (while indicating number of cells expressing the gene)
setGeneric("geneDotPlot", function(object, countAssay = "logcounts", geneNames = NULL, intCol = "subtype", orderRow = NULL, logTransform = FALSE, adjust = TRUE,
                                   low = "#132B43", high = "#56B1F7") standardGeneric("geneDotPlot"))
setMethod("geneDotPlot", signature(object = "SingleCellExperiment"), function(object, countAssay = "logcounts", geneNames = NULL, intCol = "subtype", 
                                                                              orderRow = NULL, logTransform = FALSE, adjust = TRUE, low = "#132B43", high = "#56B1F7") {
  countMat <- assay(object, countAssay)
  condDf <- colData(object)
  
  if(is.factor(condDf[[intCol]]) == FALSE) {
    condDf[[intCol]] %<>% factor(., levels = unique(condDf[[intCol]] ))
  }
  
  if(is.factor(geneNames) == TRUE) {
    geneNames <- as.character(geneNames)
  } 
  
  idList <- condDf %>% split(., .[[intCol]]) %>%
    lapply(., function(group) rownames(group))
  
  if("All" %in% geneNames || is.null(geneNames) == TRUE) {
    geneNames <- rownames(countMat)
  } else if(!all(geneNames %in% rownames(countMat))) {
    stop("Desired genes are not present in the oject!")
  }
  
  countMat %<>% .[geneNames,]
  geneNames %<>% factor(., levels = .)
  
  if(logTransform == TRUE) {
    cat("Performing log transformation...", "\n")
    countMat %<>% log()
    countMat[is.infinite(countMat)] <- 0
  }
  
  #Count percentage of cells that express the gene per group
  posCells <- lapply(seq_along(idList), function(index) {
    #Split by group
    ids <- idList[[index]]
    expMat <- countMat[,ids]
    
    #Count cells positive for the gene
    cellCounts <- apply(expMat, 1, function(gene) sum(gene > 0, na.rm = TRUE)/length(gene)*100)
    
    data.frame(Group = names(idList)[index],
               Gene = names(cellCounts),
               Percentage = cellCounts)
  }) %>% do.call("rbind", .)
  
  #Analyze differential expression per group
  countMat[is.na(countMat) == TRUE] <- 0
  
  #Scale per gene in all groups combined (to see differences between groups)
  countMat %<>% apply(., 1, function(gene) gene/mean(gene)) %>% t()
  
  #Calculate mean per cluster per gene and assemble into a dataframe
  meanExp <- lapply(seq_along(idList), function(index) {
    ids <- idList[[index]]
    
    expMat <- countMat[,ids]
    
    geneExp <- apply(expMat, 1, function(gene) mean(gene))
    
    data.frame(Group = names(idList)[index],
               Gene = names(geneExp),
               Expression = geneExp
    )
  }) %>% do.call("rbind", .)
  
  posCells %<>% dplyr::left_join(., meanExp, c("Group", "Gene"))
  posCells$Gene %<>% factor(., levels = geneNames)
  posCells$Group %<>% factor(., levels = levels(condDf[[intCol]]))
  
  #Adjust to 100 (100 in highest)
  if(adjust == TRUE) {
    if(length(adjust) > 1) {
      stop("Limits needs to be a logical vector with length 1!")
    } else if(is.logical(adjust) == FALSE) {
      stop("Limits needs to be a logical vector with length 1!")
    } else {
      cat("Performing adjustment...", "\n")
      
      posCells %<>% split(., .$Gene) %>%
        lapply(., function(gene) {
          gene$Expression <- gene$Expression/max(gene$Expression)*100
          
          return(gene)
        }) %>%
        do.call("rbind", .)
    }
  }
  
  if(is.null(orderRow) == FALSE) {
    if(!orderRow %in% geneNames) {
      stop("Cannot order by requested gene - its not in the geneNames!")
    } else {
      cat("Ordering by row...", "\n")
      intColLevels <- posCells %>% .[.$Gene == orderRow,] %>% dplyr::arrange(., Expression) %>% .[["Group"]]
      
      posCells$Group %<>% factor(., levels = intColLevels)
      posCells %<>% dplyr::arrange(., Group)
    }
  }
  
  Graph <- ggplot() +
    geom_point(data = posCells, mapping = aes(x = Gene, y = Group, size = Percentage, color = Expression)) +
    scale_color_gradient(low = low, high = high) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())
  
  if(adjust == TRUE) {
    Graph <- Graph +
      guides(size = guide_legend(title = "Pos. Cells %"),
             color = guide_colorbar(title = "Exp. %")
      )
  } else {
    Graph <- Graph +
      guides(size = guide_legend(title = "Pos. Cells %"),
             color = guide_colorbar(title = "Exp.")
      )
  }
  
  return(Graph)
})

#Accessor for tpms spot
setGeneric("tpms", function(object) standardGeneric("tpms"))
setMethod("tpms", signature(object = "SingleCellExperiment"), function(object) assay(object, "tpms"))
setGeneric("tpms<-", function(x, value) standardGeneric("tpms<-"))
setMethod("tpms<-", signature(x = "SingleCellExperiment"), function(x, value) {
  assay(x, "tpms") <- value
  
  return(x)
})


#Setting up the directories
setwd("Z:/Groups/Peder Olofsson/VSS/Bioconductor/hChAT/")

dir.create(file.path(getwd(), "VSS"))
dir.create(file.path(getwd(), "VSS", "Figures"))
dir.create(file.path(getwd(), "VSS", "QC"))
dir.create(file.path(getwd(), "VSS", "Tech"))
dir.create(file.path(getwd(), "VSS", "CorFigures"))
dir.create(file.path(getwd(), "VSS", "DeFigures"))
dir.create(file.path(getwd(), "VSS", "ClustFigures"))
figDir <- file.path(getwd(), "VSS", "Figures")
techDir <- file.path(getwd(), "VSS", "Tech")
fileDir <- file.path(getwd(), "VSS")
deDir <-  file.path(getwd(), "VSS", "DeFigures")
clDir <-  file.path(getwd(), "VSS", "ClustFigures")

#Load processed and normalized logcounts from Eric Kort
ssData <- readRDS(file.path(getwd(), "Data", "2020", "PBMC_merged_filtered_alra.rds"))

#Make data into a SingleCellExperiment for convenience
metadata <- rownames(ssData)

metadata <- data.frame(Id = str_c("Cell", "_", c(1:length(metadata))),
                       Patient = unlist(str_extract_all(metadata, "Patient\\_\\d+")),
                       CellLabel = unlist(str_extract_all(metadata, "bc[A-Z]+")),
                       RawId = metadata)
rownames(metadata) <- metadata$Id
metadata$Patient %<>% factor(., levels = unique(.))
metadata$CellLabel %<>% factor(., levels = unique(.))
metadata$RawId %<>% factor(., levels = unique(.))

ssData %<>% t()
colnames(ssData) <- metadata$Id

geneNames <- rownames(ssData)

ssData <- SingleCellExperiment(list(logcounts = ssData),
                               colData = metadata)
rm(metadata)

#Label ChATness
colData(ssData)$ChATness <- logcounts(ssData)["CHAT",] > 0

#Get CD4 signature (selection of CD4 cells is described in https://github.com/VanAndelInstitute/chat_paper/blob/master/analysis.Rmd)
selectCD4 <- openxlsx::read.xlsx(file.path(getwd(), "VSS", "Chat_all_cells true and false cd4.xlsx"))
selectCD4 %<>% .[.$Cell_ID %in% colData(ssData)$RawId,]
colData(ssData)$selectCD4 <- selectCD4 %>% .[match(colData(ssData)$RawId, .$Cell_ID), "CD4_cell"]
rm(selectCD4)

#Get raw values from transformed ones
expMat <- logcounts(ssData) %>% exp()
expMat <- expMat - 1
expMat[expMat == 0] <- NA
expMat <- expMat/10000
expMat <- expMat * sum(expMat, na.rm = TRUE)

#Select only CD4 cells
ssData <- SingleCellExperiment(list(counts = expMat),
                               colData = colData(ssData))
ssData %<>% .[, colData(.)$selectCD4 == TRUE]

#Exclude patients 10, 21, 22, 23, 25 - they have less then 10 CD4 cells per patient
ssData %<>% .[, !colData(.)$Patient %in% c("Patient_10", "Patient_21", "Patient_22", "Patient_23", "Patient_25")]
ssData$Patient %<>% factor()

#Get genes with mean reads above 5 (in cells where genes are expressed)
expMat <- counts(ssData)
geneMeans <- rowMeans(expMat, na.rm = TRUE)
geneMeans %<>% .[!is.na(.)]
expMat[is.na(expMat) == TRUE] <- 0

sum(geneMeans < 5, na.rm = TRUE)

nsGenes <- geneMeans[geneMeans > 5]
nsGenes <- names(geneMeans)[geneMeans > 5]

#Get ensembl ids for genes through HGNC symbol
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Obtain transcript length for all transcripts
ensembl_id <- getBM(attributes = c('external_gene_name', "ensembl_transcript_id", "transcript_biotype",
                                   "ensembl_exon_id", "exon_chrom_start",
                                   "exon_chrom_end", "cds_length"),
                    filters = 'external_gene_name',
                    values = nsGenes,
                    mart = ensembl)

all(nsGenes %in% ensembl_id$external_gene_name)

#Replace outdated gene names with new ones (split from https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks)
#Select genes that were not found in ensembl and split into chunks (serve connection fails at times and, as it take a silly amount of time to retrieve data, this allows to save data lost to crashes)
nsGenes %<>% .[!. %in% ensembl_id$external_gene_name] %>% split(., ceiling(seq_along(.)/10))

#Look up new gene names
newNames <- list()
for(i in seq_along(nsGenes)) {
  print(i*10)
  newNames[[i]] <- Seurat::UpdateSymbolList(nsGenes[[i]], timeout = 10, several.ok = FALSE, verbose = TRUE)
}

#First run: 1:133
newNames %<>% .[c(1:133)]

for(i in c(134:length(nsGenes))) {
  print(i*10)
  newNames[[i]] <- Seurat::UpdateSymbolList(nsGenes[[i]], timeout = 10, several.ok = FALSE, verbose = TRUE)
}

newGeneNames <- data.frame(OldNames = unlist(nsGenes),
                           NewNames = unlist(newNames)
)

#Save gene data
wb <- openxlsx::createWorkbook()
addWorksheet(wb, "NewGeneNames")
writeData(wb, 1, newGeneNames, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file.path(getwd(), "VSS", "NewGeneNames.xlsx"), TRUE)
rm(wb)

#Join all genes
ensembl_id_new <- getBM(attributes = c('external_gene_name', "ensembl_transcript_id", "transcript_biotype",
                                       "ensembl_exon_id", "exon_chrom_start",
                                       "exon_chrom_end", "cds_length"),
                        filters = 'external_gene_name',
                        values = newGeneNames$NewNames[apply(newGeneNames, 1, function(gene) length(unique(gene)) == 2)],
                        mart = ensembl)
ensembl_id %<>% rbind(., ensembl_id_new)
rm(ensembl_id_new)
rm(nsGenes)

#Save gene data
wb <- openxlsx::createWorkbook()
addWorksheet(wb, "Genes")
writeData(wb, 1, ensembl_id, startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file.path(getwd(), "VSS", "GeneData.xlsx"), TRUE)
rm(wb)

#Calculate lengths of transcripts
transLengths <- ensembl_id %>% split(., .$external_gene_name) %>%
  lapply(., function(gene) {
    #If protein-coding isoforms exist - remove non-coding isoforms
    if(sum(is.na(gene$cds_length)) < nrow(gene)) {
      gene %<>% .[!is.na(gene$cds_length),]
    }
    
    #Calculate length of each transcript
    gene %<>% split(., .$ensembl_transcript_id) %>%
      lapply(., function(transcript) {
        GRanges(
          seqnames = Rle(rep("chr1", times = nrow(transcript))),
          ranges = IRanges(start = transcript$exon_chrom_start,
                           end = transcript$exon_chrom_end),
          strand = Rle(strand(rep("*", times = nrow(transcript)))),
          transcript_id = transcript$ensembl_transcript_id,
          exon_id = transcript$ensembl_exon_id
        ) %>% width() %>% sum()
      }) %>% unlist()
    
    #Select the longest transcript
    max(gene, na.rm = TRUE)
  }) %>% unlist()
length(transLengths) == nrow(rawMatrix)

transLengths <- data.frame(SYMBOL = names(transLengths), transLength = transLengths)
rownames(transLengths) <- transLengths$SYMBOL

#Save gene data
wb <- openxlsx::createWorkbook()
addWorksheet(wb, "TransLengths")
writeData(wb, 1, transLengths, startRow = 1, startCol = 1, colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb, file.path(getwd(), "VSS", "TransLengths.xlsx"), TRUE)
rm(wb)

#Replace old names with renewed ones
metadata <- colData(ssData)

rownames(expMat)[which(rownames(expMat) %in% newGeneNames$OldNames)] <- newGeneNames$NewNames
rm(newGeneNames)

#Select genes with high enough mean reads and the ones whose transcripts have data in Ensembl
expMat %<>% .[rownames(transLengths),]


#Make a new scExperiment with remade raw values, metadata and gene data
ssData <- SingleCellExperiment(list(counts = expMat),
                               colData = metadata, 
                               rowData = transLengths)

#Remove unneeded objects
rm(metadata)
rm(transLengths)
rm(expMat)
rm(geneMeans)
rm(nsGenes)

#Calculate TPM - adopted from https://www.biostars.org/p/171766/
#Doing this before removing genes in next step - to take them into account in terms of sequencing depth
tpms <- do.call("cbind", lapply(1:ncol(counts(ssData)), function(i) {
  rpk <- counts(ssData)[, i]/rowData(ssData)$transLength
  scale_factor <- sum(rpk, na.rm = TRUE)
  tpm <- round(ceiling(rpk/scale_factor*1E6), 0)
})
)
colnames(tpms) <- colnames(counts(ssData))
tpms(ssData) <- tpms
rm(tpms)

#How many cells are genes expressed in? What percantage?
rowData(ssData)$PositiveCells <- counts(ssData) %>% apply(., 1, function(gene) sum(gene > 0))
rowData(ssData)$PositivePercentage <- counts(ssData) %>% apply(., 1, function(gene) sum(gene > 0)/length(gene)*100)

#A lot of genes which are only expressed by a couple of cells show up as DE in this analysis. It is unlikely that they
#determine biological difference between ChAT+ and ChAT- cells, thus we are removing them.
rareGenes <- rownames(rowData(ssData))[rowData(ssData)$PositivePercentage < 10]

ssData %<>% .[!rownames(counts(.)) %in% rareGenes,]

#Remove cells that have less then 50 genes expressed
ssData %<>% .[, !apply(counts(.), 2, function(cell) sum(cell) < 50)]

#QC for abundance of mitochonrial transcripts
location <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(counts(ssData)),
                   column = "SEQNAME", keytype = "SYMBOL")
is.mito <- which(location == "MT")

ssData %<>% addPerCellQC(., subsets = list(Mito = is.mito))


sum(colData(ssData)$subsets_Mito_percent > 5)
#28 cells have more than 5% of its transcripts as mitochonrdial genes

#Removing one cell and recording mitoQC
reasons <- data.frame(
  UMIcount = colData(ssData)$sum < 1e3,
  GeneCount = colData(ssData)$detected < 2e2, 
  Mito = colData(ssData)$subsets_Mito_percent > 5
)
reasons$discard <- rowSums(reasons) > 0 
colSums(as.matrix(reasons))
ssData$discard <- reasons$discard

pdf(file.path(getwd(), "VSS", "QC", "QC.pdf"), width = 30, height = 30)
gridExtra::grid.arrange(
  plotColData(ssData, x = "Patient", y = "sum", colour_by = "discard") + 
    scale_y_log10() + ggtitle("Total UMI count"),
  plotColData(ssData, x = "Patient", y = "detected", colour_by = "discard") + 
    scale_y_log10() + ggtitle("Detected genes"),
  plotColData(ssData, x = "Patient",  y = "subsets_Mito_percent", colour_by = "discard") + 
    ggtitle("Mito percent"),
  nrow = 3, ncol = 1
)
dev.off()

lost <- calculateAverage(counts(ssData)[,!reasons$discard])
kept <- calculateAverage(counts(ssData)[,reasons$discard])

logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

pdf(file.path(getwd(), "VSS", "QC", "Mito_distrubution.pdf"), width = 5, height = 5)
plotColData(ssData, x = "sum", y = "subsets_Mito_percent", 
            colour_by = "discard") +
  theme(panel.border = element_rect(color = "grey"))
dev.off()

pdf(file.path(getwd(), "VSS", "QC", "Discarded.pdf"), width = 5, height = 5)
plot(abundance, logFC, xlab = "Average count", ylab = "Log-FC (lost/kept)", pch = 16)
points(abundance[is.mito], logFC[is.mito], col = "dodgerblue", pch = 16)
dev.off()

#Remove unneeded objects
rm(reasons)
rm(lost)
rm(kept)
rm(logged)
rm(abundance)
rm(logFC)
rm(is.mito)
rm(location)
rm(rareGenes)

#Removing cells
ssData <- ssData[,!colData(ssData)$discard]

#Normalization by deconvolution
set.seed(100)
clust <- quickCluster(ssData)
ssData <- computeSumFactors(ssData, cluster = clust, min.mean = 0.1)
ssData <- logNormCounts(ssData)

#Looking for most variable genes for PCA
dec <- modelGeneVar(ssData, block = ssData$Patient)

#Visualizing the fit:
fit <- metadata(dec)
par(mfrow=c(1,12))
blocked.stats <- dec$per.block

pdf(file.path(techDir, "VariableGenesFit.pdf"), width = 5, height = 5)
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main = i, pch = 16, cex = 0.5,
       xlab = "Mean of log-expression", ylab = "Variance of log-expression")
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col = "red", pch=16)
  curve(curfit$trend(x), col = 'dodgerblue', add = TRUE, lwd = 2) 
}
dev.off()

dec[order(dec$bio, decreasing = TRUE), 1:6]

#Most variable genes
hvg <- getTopHVGs(dec, n = 2000)

#PCA
set.seed(100)
ssData <- runPCA(ssData, subset_row = hvg)
reducedDim(ssData, "PCA") <- reducedDim(ssData, "PCA")[, 1:50]

set.seed(00101001101)
ssData <- runTSNE(ssData, dimred = "PCA")

set.seed(1100101001)
ssData <- runUMAP(ssData, dimred = "PCA")

pdf(file.path(techDir,  "DimRed.pdf"), width = 9, height = 9)
gridExtra::grid.arrange(
  plotReducedDim(ssData, dimred = "PCA", colour_by = "Patient"),
  plotReducedDim(ssData, dimred = "PCA", colour_by = "ChATness"),
  plotReducedDim(ssData, dimred="TSNE", colour_by = "Patient"),
  plotReducedDim(ssData, dimred="TSNE", colour_by = "ChATness"),
  plotReducedDim(ssData, dimred = "UMAP", colour_by = "Patient"),
  plotReducedDim(ssData, dimred = "UMAP", colour_by = "ChATness"),
  nrow = 3, ncol = 2
)
dev.off()

#Remove unnecessary objects
rm(blocked.stats)
rm(dec)
rm(fit)
rm(current)
rm(curfit)
rm(hvg)
rm(i)
rm(clust)

#Annotation
ref <- DatabaseImmuneCellExpressionData(cell.ont = "all")
pred <- SingleR(test = ssData, ref = ref, labels = ref$label.main)
to.remove <- pruneScores(pred)

#Transfer labels to colData
ssData$celltype <- pred$labels

#Fine labels
pred2 <- SingleR(test = ssData, ref = ref, labels = ref$label.fine)

#Transfer labels to colData
ssData$subtype <- pred2$labels 

#Filter out non-T cells
ssData %<>% .[, !to.remove]
ssData %<>% .[, colData(.)$celltype %in% c("T cells, CD4+")]
ssData %<>% .[, !colData(.)$subtype %in% c("B cells, naive", "Monocytes, CD14+", "T cells, CD8+, naive, stimulated", "Monocytes, CD16+", "NK cells")]

#Remove unnecessary objects
rm(ref)
rm(pred)
rm(pred2)
rm(to.remove)