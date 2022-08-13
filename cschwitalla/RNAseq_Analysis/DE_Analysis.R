################################################################################
###                                Purpose                                   ###
################################################################################
## DESCRIPTIOM: This script performs differential expression analysis on Salmon
##              quantification files, based on the DESeq2 package.
##
##
## INPUT: Salmon quantification files
## OUTPUT: for each pairwise comparison: - count table with DE expression
##                                         Annotation --> SIG, DOWN, UP
##                                       - Volcano plot
##                                       - MA plot
################################################################################
###                             References                                   ###
################################################################################
## Differential expression analysis and testing is based on the DEseq2 package
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html



################################################################################
###                             Manual Steps                                 ###
################################################################################
# clear R environment
rm(list = ls())
# set directory of the input files for the analysis
input_dir <- "/Users/cschwitalla/Documents/transcriptomics_results/Pipeline_running/transcriptomics_results/Quant_files/"
# set output directory for results table and plots
output_dir <- "/Users/cschwitalla/Documents/RNAseq_analysis/"
# set directory of the metadata file
metadata_file <- "/Users/cschwitalla/Documents/RNAseq_analysis/Metadata_GB.tsv"
# set directory of the script that contains RNAseq analysis functions
source("/Users/cschwitalla/git/students/cschwitalla/RNAseq_Analysis/functions_RNA.R")
# set lfc and padju values for significance thresholds
lfc <- 3
padju <- 0.01
################################################################################
###                           Load libraries                                 ###
################################################################################
required_libs <- c(
  "DESeq2", "tximport", "ashr", "ggVennDiagram", "ggplot2",
  "dplyr", "readxl", "tidyr", "pheatmap", "RColorBrewer",
  "EnsDb.Hsapiens.v86", "ggrepel", "utils", "EnhancedVolcano", "limma"
)


suppressMessages(invisible(lapply(required_libs, library, character.only = TRUE)))


################################################################################
###                            Load Data                                     ###
################################################################################
# Load meta data --> Metadata_GB.tsv in workdir
metadata <- read.table(file = metadata_file, sep = "\t", header = TRUE)
metadata2 <- metadata[-grep(("QATLV129AQ|QATLV139AX|QATLV162AW|QATLV171AV|QATLV188AQ"),metadata$QBiC.Code),]
# get filenames of inputDir
file_names <- list.files(path = input_dir)
# files without ben + outlier sample
filnames_excl <- grep(("NEC|INF|T1"), file_names, value = TRUE)
filnames_excl <- filnames_excl[c(1:7,9:45)]
# TO DO: create/ import tx2gene dataframe
tx2gene <- create_tx2gene()
# import Salmon quant files with tximport and tx2gene
txi <- tximport(paste0(input_dir, file_names), type = "salmon", tx2gene = tx2gene)
txi2 <- tximport(paste0(input_dir, filnames_excl), type = "salmon", tx2gene = tx2gene)

# make Deseq2 object of raw read counts and specify metadata and project design
dds <- DESeqDataSetFromTximport(txi,
  colData = metadata,
  design = ~Tumor_region
)
dds2 <- DESeqDataSetFromTximport(txi2,
                                colData = metadata2,
                                design = ~Tumor_region
)
# pre filtering, keep only rows with at least 10 reads total
# WHY:reduces computation time and memory usage
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds2 <- dds2[keep, ]
# model read counts, estimate LFCs and test for DE expressed genes
# with default parameters
dds_default <- DESeq(dds)
dds_minor <- DESeq(dds2)
################################################################################
###                            DE_Analysis                                   ###
################################################################################
# list of genes that are diff expressed
DE_genes <- c()
# get each pairwise comparison for condition that should be tested
for (com in apply(combn(levels(dds_default$Tumor_region), 2), 2, paste, collapse = " vs ")) {
  # get condition vector for DE analysis function
  comparison <- unlist(strsplit(com, " vs "))
  # make results data frame
  res_df <- DE_analysis(dds_default, c("Tumor_region", comparison), lfc, padju)
  # get DE genes and append them into DE_genes list
  DE_genes <- append(DE_genes, res_df$Names[res_df$significant == "SIG"])
  # write results as tsv into the workdir to stor results
  write.table(res_df, file = paste(output_dir, com, ".tsv"), sep = "\t", quote = FALSE)
  # make MA plot and save it as png in output dir
  pdf(paste0(output_dir, "MAplot_", com, ".pdf"))
  DESeq2::plotMA(res_df, xlim = c(1, 1e5), ylim = c(-4, 4), main = com)
  dev.off()
  # make Volcano plot and save it as pdf in output dir
  pdf(paste0(output_dir, "Volcano_", com, ".pdf"), width = 10, height = 8)
  make_volcano(res_df, "padj", com, lfc, padju)
  dev.off()
}


################################################################################
###                              Heatmap                                     ###
################################################################################
set.seed(8)
# varianze stabilizing transformation of the data to plot as heatmap
vsd <- vst(dds_minor, blind = FALSE)
# create list of annotation colors for heatmap plotting
annotation_color <- create_annotation_color(
  metadata$Patient_ID,
  metadata$Tumor_region,
  metadata$Sex,
  metadata$MGMT
)
# plot heatmap with bacth corection for Patientens
make_heatmap(DE_genes, vsd, vsd$Patient_ID, annotation_color, k = NULL)
# plot without batch correction
make_heatmap(DE_genes, vsd, batch = NULL, annotation_color)

#make_heatmap(DE_genes, vsd, vsd$Patient_ID, annotation_color, k = 4)
################################################################################
###                               PCA                                        ###
################################################################################
# plot pca with batch correction
plot_pca(dds_default = dds_default, batch = vsd$Patient_ID)
# plot oca without batch correction
plot_pca(dds_default = dds_default, batch = NULL)



# Inspection--------------------------------------
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_default)[["cooks"]]), range=0, las=2)

res <- results(dds_default)
plotDispEsts(dds_default)

# TODO: make a function out of it  ---------------------------------------------
sampleDists <- dist(t(assay(vsd)))


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Tumor_region, vsd$Patient_ID, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)

# test for enrichment analyiss and stuff----------------------------------------
# ensembla and org db sind nicht vergleichbar und deshlab lieber nicht mischen
# muss hier ne gute lösung finden
