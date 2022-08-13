################################################################################
###                                Purpose                                   ###
################################################################################
## DESCRIPTIOM: This script performs WES analysis by using mainly the maftools
##              package
##
## INPUT: maf files
## OUTPUT: Oncoplots
################################################################################
###                               TODOs                                      ###
################################################################################
# TODO: affected pathways from maftools 
# TODO: appropriate filtering for oncoplot and in general 
# TODO: known mutations in GB 

################################################################################
###                               Manual steps                               ###
################################################################################
## Clear the R environment
rm(list = ls())
# set working directory (top folder were the maf fiels are)
setwd("/Users/cschwitalla/Documents/WES_data/sarek3_annotation/strelka/")

# list of directories were maffiles for each region are +
# set regionName as name in same order like dirs to know which maf file
# belogs to which region

region_dir <- setNames(
  c(
    "./NEC//maf_converted/",
    "./T1/maf_converted/",
    "./INF/maf_converted/",
    "./BEN/maf_converted/"
  ),
  c("NEC", "T1", "INF", "BEN")
)
# define output directory
output_dir <- "/Users/cschwitalla/Documents/WES_analysis/"
# set path to the metadata file
metadata_file <- "/Users/cschwitalla/Documents/WES_analysis/WES_metadata.tsv"
# bed file for capture kit size
bed_file_dir <- "../../hglft_genome_298d9_1b60f0.sort.merged.bed"

# set path to functions file required for this script
source("/Users/cschwitalla/git/students/cschwitalla/WES_Analysis/functions_exome.R")
source("/Users/cschwitalla/git/students/cschwitalla/Ligandomics_Analysis/functions_ligandomics.R")
################################################################################
###                             Load libraries                               ###
################################################################################
# Load the libraries
required_libs <- c("maftools", "ggplot2", "circlize", "UpSetR", "dplyr", "tuple",
                   "ggVennDiagram")

suppressMessages(invisible(lapply(required_libs, library, character.only = T)))

################################################################################
###                            Load Data                                     ###
################################################################################
# load metadata file
metadata <- read.table(file = metadata_file, sep = "\t", header = TRUE)
metadata <- metadata %>% arrange(factor(metadata$Tumor_region, levels = c("BEN", "INF", "T1", "NEC")))
#calc capture kit size 
bed_file <- read.csv(file = bed_file_dir, 
                     sep = "\t",
                     header = FALSE)
kit_size = sum(bed_file$V3 - bed_file$V2)

# list to store the merged maf files for each region
regionmaf_list <- list()
all_files <- c()
# looooooop --------------------------------------------------------------------
for (dir in region_dir) {
  # get file paths
  files <- list.files(dir, full.names = TRUE, pattern = ".maf")
  all_files <- append(all_files, files)
  # merge maf files + metadata
  maf <- merge_mafs(files, clinicalData = metadata)
  # get clinical significance 
  clin_sig_region <- maf@data[,c("Hugo_Symbol", "CLIN_SIG")] %>%
    filter(!grepl("likely_benign|benign", CLIN_SIG) & CLIN_SIG != "")

  # append maf file into list
  regionmaf_list <- append(regionmaf_list, maf)
  clinsig_list <- append(clinsig_list,clin_sig_region)
  #make annotation list with metadata
  anno_color <- create_annotation_color(patient_column = metadata$Patient_ID,
                                        region_column = metadata$Tumor_region)
  # make pdf
  make_oncoplot(output_dir, region_dir, maf, anno_color, "_all", genes = NULL, sample_order = NULL)
  # make ocoplot with only clinical sig genes
  make_oncoplot(output_dir, region_dir, maf, anno_color, "_ClinSig_genes", c(unique(clin_sig_region$Hugo_Symbol)),sample_order = NULL)
  # mutload compare tcga
  tcga_comapre(maf, dir, kit_size, output_dir)
  #onco pathways
  onco_pathways(maf, output_dir, dir)
}
# set names for merged maf files in the list, to see which region is it
names(regionmaf_list) <- c("NEC", "T1", "INF", "BEN")

# make venn diagram including all genes from all regions
all_genes_set <- setNames(vector("list", length = length(names(regionmaf_list))),
                          c(names(regionmaf_list))
)
for (i in names(regionmaf_list)) {
  all_genes_set[i] <- list(regionmaf_list[[i]]@data$Hugo_Symbol)
}
all_genes_venn_data <- plot_custom_venn(all_genes_set, "all genes")

# merge all tumor region mafs + metadata
all_maf <- merge_mafs(all_files, clinicalData = metadata)
# make oncoplot with all mafs
sample_order <- unique(metadata$Tumor_Sample_Barcode)
make_oncoplot(output_dir, "all", all_maf, anno_color, "_all", genes = NULL, sample_order = sample_order)


#nix bei raus gekommen :(
test <- clinicalEnrichment(all_maf, clinicalFeature = "Tumor_region")
plotEnrichmentResults(enrich_res = test, pVal = 0.3)



# -- gene selection befor plotting----------------------------------------------
NEC_maf <- regionmaf_list$NEC
NEC_maf@data[,c("Hugo_Symbol", "CLIN_SIG")]

oncoplot(
  maf = nec_maf,
  genes = unique(genes$Hugo_Symbol),
  draw_titv = FALSE,
  clinicalFeatures = c("Patient_ID", "Tumor_region"),
  leftBarData = compute_vaf(nec_maf), # see function script
  annotationColor = anno_color, # see function script
  colors = vc_cols,
  sortByAnnotation = TRUE,
  gene_mar = 6
)# see function script

genes <- NEC_maf@data[,c("Hugo_Symbol", "CLIN_SIG", "PUBMED", "dbSNP_RS", "Tumor_Sample_Barcode", "PolyPhen")] %>%
  filter(CLIN_SIG != "")



getGeneSummary(NEC_maf)

nec_files <- list.files(path = "NEC/maf_converted/", full.names = TRUE, pattern = ".maf")
nec_maf <- maf <- merge_mafs(nec_files, clinicalData = metadata)
nec_anno <- nec_maf@data 
nec_gene <- nec_maf@gene.summary
#subset maf 
library("vcfR")
vcf <- read.vcfR("/Users/cschwitalla/Documents/WES_data/sarek3_annotation/manta/NEC/filtered/Manta_QATLV067AW_vs_QATLV069AE.diploidSV_VEP.ann_VEP.ann_filtered.vcf", verbose = FALSE)
test <- vcf@fix


#MY VARIANT TEST----------------------------------------------------------------



library("readxl")
library("myvariant")
stuf <- read_xlsx("/Users/cschwitalla/Downloads/Binders Medulloblastoma.xlsx")
test_vcf <- readVcf("/Users/cschwitalla/Documents/exome_data_indexed/INF/raw_data/Strelka/somatic/filtered/concat_snvs_filtered.vcf")
hgvs <- formatHgvs(test_vcf, variant_type = "snp")

head(hgvs)


test_clinvar <- getVariants(hgvs, fields = "clinvar.type")
test_cadd <- getVariants(hgvs, fields = "cadd.consequence")

# include lollipopplot maybe
# circouise plots
#circos.initializeWithIdeogram(species = "hg19")

# circos plits ---------------
bed_file <- read.table("/Users/cschwitalla/Documents/exome_data_indexed/INF/raw_data/Strelka/somatic/filtered/all_unsorted.bed", sep = "\t", header = FALSE)
bed_file <- bed_file[c(1:3, 6:7)]
names(bed_file) <- c("chr", "start", "end", "ref", "alt")
circos.initializeWithIdeogram(bed_file)
circos.genomicTrack(test,)
test_split <- split(bed_file,bed_file$chr)


chr1 <- test_split$chr1
test_chr1 <- table(chr1)

test_freq <- apply(chr1[c("start", "end", "ref", "alt")], 2, table)



circos.initializeWithIdeogram(plotType = NULL)
circos.track(ylim = c(0,1), panel.fun = function(x,y){
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)  




circos.genomicTrackPlotRegion(bed_file, panel.fun = function(region, value , ...) {
  if(CELL_META$sector.index == "chr1") {
    print(head(region, n = 2))
    print(head(value, n = 2))
  }
})
  




