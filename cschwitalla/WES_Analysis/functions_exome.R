################################################################################
###                                Purpose                                   ###
################################################################################
## Title: Functions for WES Analysis
## Author: Carolin Schwitalla
##
## Description: This script contains all functions that are used for the WES
##              analysis for my master thesis
################################################################################
###                             Hard coded Variables                         ###
################################################################################
# colors for oncoplot itself
vc_cols <- setNames(
  c("#A00202", "#438EA9", "#436211", "#323949", "#EE6D17", "#FFC61A"),
  c(
    "Frame_Shift_Del", "Missense_Mutation", "Nonsense_Mutation",
    "Multi_Hit", "Frame_Shift_Ins", "In_Frame_Del"
  )
)


################################################################################
###                                Functions                                 ###
################################################################################
#oncogenic pathways
onco_pathways <- function(maf, output_dir, dir) {
  pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)],"_onco_pathway.pdf"),
      width = 7,
      height = 8)
  OncogenicPathways(maf = maf)
  dev.off()
}


# tcga compare
tcga_comapre <- function(maf, dir, kit_size, output_dir) {
  pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)],"_mutload.pdf"),
      width = 10,
      height = 5)
  tcgaCompare(maf = maf, cohortName = names(region_dir)[match(dir, region_dir)], logscale = TRUE, capture_size = kit_size/1000000)
  dev.off()
}



#make oncoplot pdf
make_oncoplot <- function(output_dir, region_dir, maf, anno_color, title, genes, sample_order) {

    pdf(paste0(output_dir, names(region_dir)[match(dir, region_dir)], title,"_oncoplot.pdf"),
        width = 10,
        height = 8)
    # make oncoplot for current maf that will exported as pdf in the output dir
    oncoplot(
      maf = maf,
      genes = genes,
      draw_titv = FALSE,
      clinicalFeatures = c("Tumor_region", "Patient_ID"),
      leftBarData = compute_vaf(maf), # see function script
      annotationColor = anno_color, # see function script
      colors = vc_cols,
      sampleOrder = sample_order,
      sortByAnnotation = TRUE,
      gene_mar = 6
    )# see function script
    dev.off()

    
  
}
# Variant allel Frequency
# ----- from Rike
# Calculate Varaint allele frequency
compute_vaf <- function(maf_obj) {
  genes_vaf <- subsetMaf(
    maf = maf_obj, query = "t_depth > 0",
    fields = c("t_depth", "t_alt_count"), mafObj = FALSE
  )
  genes_vaf$VAF <- as.numeric(genes_vaf$t_alt_count) / as.numeric(genes_vaf$t_depth)
  genes_vaf <- genes_vaf[, mean(VAF), Hugo_Symbol]
  colnames(genes_vaf)[2] <- "VAF"
  genes_vaf <- genes_vaf[which(genes_vaf$VAF > 0.02)]
  return(genes_vaf)
}



## DESCRIPTION: mapping meta data to colors for annotation color in oncoplots
## PARAMETERS:
##             - dregion_column: matadata df column were tumor region / conditions
##                               are listed
##             - patient_column: metadata df column were patien ids are listed
## OUTPUT: returns a list of lists
##
create_annotation_color <- function(patient_column, region_column) {
  patient_scale <- colorRampPalette(c("#543005", "#f5f5f5", "#003c30"))
  region_color <- c("#4B6C22", "#74add1", "#9e0142", "#fdae61")
  #get the number of unique patient ids to extract colors from color scale
  sample_num <- length(unique(patient_column))
  # get a specific color palette with num of patients
  patient_color <- patient_scale(sample_num)
  # set names to asign for each level the right color
  names(patient_color) <- unique(patient_column)
  names(region_color) <- sort(unique(region_column))
  annotation_color <- list(Tumor_region = region_color, Patient_ID = patient_color)
  return(annotation_color)
}
