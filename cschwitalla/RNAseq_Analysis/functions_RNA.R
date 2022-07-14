################################################################################
###                                Purpose                                   ###
################################################################################
## Title: Functions for RNAseq Analysis
## Author: Carolin Schwitalla
##
## Description: This script contains all functions that are used for the RNAseq
##              analysis for my master thesis
################################################################################
###                                                                          ###
################################################################################

# create tx2gene df that maps transcripts to genenames for tximport
# according to DESeq2 documentation

## DESCRIPTION: function that creates a tabel, mapping transcript IDs to
##              gene names
## Requirements: load EnsDb.Hsapiens.v86
##
##
create_tx2gene <- function() {
  edb <- EnsDb.Hsapiens.v86
  k <- keys(edb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(edb, k, "GENENAME", "TXNAME")
  return(tx2gene)
}


## DESCRIPTION: function that performs DE analysis
## PARAMETERS:
##             - deseq2_object: DESeq(dds)
##             - contrast: vector of strings, containing condition column
##                         and 2 particular levels that should be compared
##                         [vector]
##             - lfc: threshold for log2FoldChange [numeric]
##             - padj: threshold for p adjusted value [numeric]
##
DE_analysis <- function(deseq2_object, contrast, lfc, padj) {
  # get shrunken results of pairwise comparison
  # shrinkage reduces noisy lfcs for low count genes
  # contrast is used because of more than 2 conditions -->
  # need to specify which conditions will be compared
  # ashr type of shrinkage because by using contrast it's the only option
  DE_df <- lfcShrink(deseq2_object, contrast = contrast, type = "ashr")
  # include columns for UP , DOWN , SIG , expressed genes
  DE_df$significant <- "NO"
  DE_df$significant[abs(DE_df$log2FoldChange) > lfc & DE_df$padj < padj] <- "SIG"
  DE_df$diffexp[DE_df$log2FoldChange > 0 & DE_df$significant == "SIG"] <- "UP"
  DE_df$diffexp[DE_df$log2FoldChange < 0 & DE_df$significant == "SIG"] <- "DOWN"
  DE_df$Names <- DE_df@rownames
  return(DE_df)
}

## DESCRIPTION: function that makes Volcano plot
## PARAMETERS:
##             - res_df: DE dataframe from pairwise comparison [data frame]
##             - y_val: set either padju or pval for y axis [string]
##             - title: title of volcano plot [string]
##             - lfc: log2FoldChange threshold [numeric]
##             - padj: threshold for p adjusted value [numeric]
## OUTPUT: plots a volcano plot of genes compared between two conditions
##
make_volcano <- function(res_df, y_val, title, lfc, padju) {
  EnhancedVolcano(res_df,
    lab = rownames(res_df),
    x = "log2FoldChange",
    y = y_val,
    FCcutoff = lfc,
    pCutoff = padju,
    ylim = c(0, max(log10(res_df$padj))),
    xlim = c(
      min(res_df$log2FoldChange),
      max(res_df$log2FoldChange)
    ),
    ylab = bquote(paste(-Log[10], .(y_val))),
    selectLab = NA, # if TRUE: gene names will be visible as labels
    labSize = 3, # gene label size
    drawConnectors = TRUE, # gene and corresponding labels are
    # connected with arrows
    title = title,
    legendPosition = "bottom",
    legendLabels = c(
      "NS", "Log2 FC",
      y_val,
      paste(y_val, "& Log2 FC")
    )
  )
}



## DESCRIPTION: mapping meta data to colors for annotation color in oncoplots
## PARAMETERS:
##            - region_column: matadata df column were tumor region /conditions
##                             are listed
##            - patient_column: metadata df column were patient ids are listed
##            - sex_column: metadata df column were the patients sex is listed
##            - meth_column: metadata df column were the patients
##                           MGMT methylation status is listed
## OUTPUT: returns a list of lists with colors for each factor level
##
create_annotation_color <- function(patient_column,
                                    region_column,
                                    sex_column,
                                    meth_column) {
  patient_scale <- colorRampPalette(c("#543005", "#f5f5f5", "#003c30"))
  region_color <- c("#4B6C22", "#74add1", "#9e0142", "#fdae61")
  sex_color <- c("#D5BD9E", "#565E71")
  meth_color <- c("#4B6C22", "#74add1")

  # get the number of uniqe patient ids to extract colors from color scale
  sample_num <- length(unique(patient_column))
  # get a specific color palette with num of patients
  patient_color <- patient_scale(sample_num)
  # set names to asign for each level the right color
  names(patient_color) <- unique(patient_column)
  names(region_color) <- sort(unique(region_column))
  names(sex_color) <- sort(unique(sex_column))
  names(meth_color) <- sort(unique(meth_column))
  annotation_color <- list(
    Tumor_region = region_color,
    Patient_ID = patient_color,
    Sex = sex_color,
    MGMT_methylation = meth_color
  )
  return(annotation_color)
}






## DESCRIPTION: function that plots heatmap for a gene selection
## PARAMETERS:
##             - gene_selection: gene names that will be plotted [vector]
##             - vsd: transformed count data object from DESeq2 [DESeq2_obj]
##             - batch: column in vsd df for which batch correction
##                      should be performed [vector]
## OUTPUT: heatmap of selected genes with metadata column annotation
make_heatmap <- function(gene_selection, vsd, batch, annotation_color, k) {
  # get vst normalized counts df for gene selection
  vsd_mat <- as.data.frame(assay(vsd)) %>%
    dplyr::filter(row.names(assay(vsd)) %in% gene_selection)
# remove batch effects with limma because it was suggested by deseq2 authors
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  if (!is.null(batch)) {
    design_matrix <- model.matrix(~Tumor_region, colData(vsd))
    vsd_rm_batch <- limma::removeBatchEffect(vsd_mat,
                                             batch,
                                             design = design_matrix)
    vsd_mat <- vsd_rm_batch
  }
  # annotation for heatmap
  annotation_df <- data.frame(
    "Tumor_region" = vsd@colData@listData$Tumor_region,
    "Patient_ID" = vsd@colData@listData$Patient_ID,
    "Sex" = vsd@colData@listData$Sex,
    "MGMT_methylation" = vsd@colData@listData$MGMT
  )
  if (!is.null(k)) {
    pheatmap(
      vsd_mat,
      scale = "row", # genes = rows are z-score transformed
      show_rownames = FALSE,
      annotation_col = annotation_df,
      annotation_colors = annotation_color,
      kmeans_k = k
    )
  } else {
    pheatmap(
      vsd_mat,
      scale = "row", # genes = rows are z-score transformed
      show_rownames = FALSE,
      annotation_col = annotation_df,
      annotation_colors = annotation_color,
      cutree_cols = 4
    )
  }
# 
#   # draw heatmap
#   pheatmap(
#     vsd_mat,
#     scale = "row", # genes = rows are z-score transformed
#     show_rownames = FALSE,
#     annotation_col = annotation_df,
#     annotation_colors = annotation_color
#   )
}



## DESCRIPTION: function that plots PCA with or without batch correction
## PARAMETERS:
##             - dds_default: dds object from DESeq2 [dds object]
##             - batch: vsd column of the batch [vsd column]
##
## OUTPUT: PCA plot
plot_pca <- function(dds_default, batch) {
  vsd <- vst(dds_default, blind = FALSE)
  vsd_mat <- assay(vsd)
  if (!is.null(batch)) {
    design_matrix <- model.matrix(~Tumor_region, colData(vsd))
    vsd_rm_batch <- limma::removeBatchEffect(vsd_mat,
                                             batch,
                                             design = design_matrix)
    vsd_mat <- vsd_rm_batch
  }
  assay(vsd) <- vsd_mat
  plotPCA(vsd, intgroup = "Tumor_region")
}



## DESCRIPTION: function that plots historams of raw counts for each sample
## PARAMETERS:
##             - txi_data: txi object with raw counts data [txi object]
##             - col_row_num: set panels for plotting window [vector]
##             - sample_num: number of samples that should be plotted [numeric]
##             - xlim: x axis limits (lower & upper) for each histogram [vector]
##             - ylim: y axis limits (lower & upper) for each histogram [vector]
##             - breaks: breaks for each histogram [numeric]
## OUTPUT: for each sample a histogram of raw counts distribution is plotted
plot_rawcounts_hist <- function(txi_data,
                                col_row_num,
                                sample_num,
                                xlim,
                                ylim,
                                breaks) {
  # creat df of rounded raw count from the txi object
  data <- txi_data$counts %>%
    round() %>%
    data.frame()
  # plot raw counts histograms for each sample side by side
  par(mfrow = col_row_num) # define panels
  for (sample in 1:sample_num) {
    hist(data[, sample],
      breaks = breaks,
      xlim = xlim,
      ylim = ylim,
      main = sample
    )
  }
}
