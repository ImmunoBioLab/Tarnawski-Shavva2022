# Tarnawski-Shavva2022

Contains code for SingleCellExperiment pipleline used to analyse CD4+ T cells. Selection of CD4+ cells is described in https://github.com/VanAndelInstitute/chat_paper/blob/master/analysis.Rmd.

Content:
CHAT_A:
  1) Loading data;
  2) Selection of genes and cells for downstream analysis;
  3) TPM calculation;
  4) Normalization;
  5) Dimensional reduction;
  6) Filtering CD4+ cells based on automated celltype assignment.

CHAT_B:
  1) Clustering;
  2) Marker gene detection.

CHAT_C:
  1) Correlation of ChAT expression vs other genes and selection of top correlated genes.
  2) Pseudobulk DE analysis of ChAT+ vs ChAT- CD4+ T cells;
  3) Filtering cells based on expression of top correlated genes.

Chat_all_cells true and false cd4.xlsx:
  Excel sheet containing cell IDs and their CD4 (TRUE/FALSE) status - necessary for selecting cells for downstream analysis.
