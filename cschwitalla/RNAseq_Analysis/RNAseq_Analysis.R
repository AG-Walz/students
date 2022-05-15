## I would suggest to move input content always to the top, this way you can
## find it faster
setwd("/Users/cschwitalla/Documents/transcriptomics_results/Pipeline_running/transcriptomics_results/Quant_files/")
## And I was wondering why you define the filedir with the same repository?
## After setwd, you can also mention:
##  filedir <- getwd()
## Or define the filedir first and then use that for
##  setwd(filedir)
filedir = "/Users/cschwitalla/Documents/transcriptomics_results/Pipeline_running/transcriptomics_results/Quant_files/"

# library("DESeq2")
# library("tximport")
# library("tximportData")
# library("readr")
# library("EnsDb.Hsapiens.v86") 
# library("utils")
# library("devtools")
# library("ashr")
# library("ggrepel")
# library("dplyr")
# library("ggplot2")
# library("ggVennDiagram")
# library("VennDiagram")
# library("RColorBrewer")
# library("tibble")
# library("pheatmap")
# library("viridis")
# library("clusterProfiler")
# library("tidyr")
# library("pcaExplorer")
# library("readxl")
# library("stringr")

## Suggestion Marissa
## ---------------------
necessaryLibs <- c("DESeq2", "tximport", "tximportData", "readr",
                   "EnsDb.Hsapiens.v86", "utils", "devtools", "ashr", "ggrepel",
                   "dplyr", "ggplot2", "ggVennDiagram", "VennDiagram",
                   "RColorBrewer", "tibble", "pheatmap", "viridis",
                   "clusterProfiler", "tidyr", "pcaExplorer", "readxl",
                   "stringr")
# invisible(lapply(necessaryLibs, BiocManager::install, update = F, ask = F))
suppressMessages(invisible(lapply(necessaryLibs, library, character.only = T)))
#list of all file names -----------------------------------
## Here, try to reduce the use of "=" in R when you assign a vector, this can 
## potentially raise errors
## Additionally, try to write your though process down

## Listing all files in the file directory
file_names <- list.files(path = filedir)
## Obtaining the absolute path of files that are matching the pattern
files <- dir(filedir, recursive= TRUE, pattern = "quant.sf", full.names = TRUE)
## Assign the relative path as the name
names(files) <- dir(filedir)

# CALC TPM FOR ALL MEIN REGIONS -------------------------------------

BEN_files = subset(files, )





# Metadata ------------------------------------
## Suggestion: all manual parameters on top
actuel_etadata = read_xlsx("/Users/cschwitalla/Downloads/Metadata patients P-141L.xlsx", col_names = TRUE)
metadata = read_xlsx("/Users/cschwitalla/Documents/QATLV_sample_preparations.xlsx", col_names = TRUE)

## If you have a stack overflow solution write down what this function does and
## include the link where you found it, you never know when this comes in handy

# stack over flow function
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

## I was wondering though, why include this step?
metadata = header.true(metadata)

## Maybe we need to sit down to discuss what you are trying to do here
t1_num = metadata[which(metadata$`Condition: tissues` == "T1 zone"),]
nec_num = metadata[which(metadata$`Condition: tissues` == "Necrotic zone"),]
inf_num = metadata[which(metadata$`Condition: tissues` == "Peritumoral infiltrating zone"),]
ben_num = metadata[which(metadata$`Condition: tissues` == "Benign"),]

methylation_status = actuel_etadata[1:15,c(2,13)]
patients_zh = methylation_status$`ZH number`
new_patients = c()
for(i in patients_zh){
  str1 = strsplit(i, "/")
  #print(str1)
  str2 = str1[[1]][1]
  str_replace_all(str2, pattern = " ", replacement = "")
  new_patients = append(new_patients, str2)
  
}
new_patients = new_patients[1:15]
new_patients[15] = "ZH1041"
new_patients[1] = "ZH0984"
new_patients[2] = "ZH0991"
new_patients[3] = "ZH0997"
  
methylation_status$`ZH number`= new_patients
names(metadata)[6] = "ZH number"

pat_meta=metadata$`ZH number`
mgmt_sta = c()
for(i in pat_meta){
  
}

metadata_2 = merge( methylation_status, metadata, by = "ZH number")

#construct samples table with condition information for each file ---------------------------

sample_names = c()

tumor_region = c()


for(i in file_names){
  string1 = strsplit(i, "-" )
  sample_names = append(sample_names, string1[[1]][1])
  n = length(string1[[1]])
  string2 = strsplit(string1[[1]][n], "_")
  
  tumor_region = append(tumor_region, string2[[1]][1])
}

samples = data.frame("Tumor_region"= tumor_region , "Sample_name" = sample_names, "File_name" = names(files))
samples$File_names <- names(files)
samples$run <- sample_names
samples$TumorRegion <- tumor_region



#tx2gen generation --> linking transcripts to genes

edb <- EnsDb.Hsapiens.v86
k <- keys(edb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(edb, k, "GENENAME", "TXNAME")
name2ID <- AnnotationDbi::select(edb,k, "GENEID", "TXNAME")
mapping_df = merge(tx2gene, name2ID, by = c("TXNAME","TXID"))

# import salmon files 
txi <- tximport(file_names, type = "salmon", tx2gene = tx2gene)

names(txi)
#deseq dataset from  txi object and tx2gene
dds <- DESeqDataSetFromTximport(txi, 
                                colData = samples, 
                                design = ~ Tumor_region)

#pre filtering --> keep only rows with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# RESULTS---------------------------------------

dds <- DESeq(dds)

## All of the code here (where you are comparing the conditions) can be reduced
## significantly. See example: 
## https://github.com/AG-Walz/Projects/blob/main/QBiC_coops/QLVMB/DEA/transcriptomics%20analysis.R#L137-L155

# https://support.bioconductor.org/p/98346/
# RESULTS----------------------------------------
res_NEC_vs_BEN = results(dds, contrast = c("Tumor_region", "NEC", "BEN"))
res_T1_vs_BEN = results(dds, contrast = c("Tumor_region", "T1", "BEN"))
res_INF_vs_BEN = results(dds, contrast = c("Tumor_region", "INF", "BEN"))

res_NEC_vs_T1 = results(dds, contrast = c("Tumor_region", "NEC", "T1"))
res_NEC_vs_INF = results(dds, contrast = c("Tumor_region", "NEC", "INF"))
res_T1_vs_INF = results(dds, contrast = c("Tumor_region", "T1", "INF"))


# SHRINKAGE ------------------------------------
# ashr for shrinkage because apeglm only for coef not for contast, and better than normal --> said DESeq2
resShrink_NEC_vs_BEN = lfcShrink(dds, contrast = c("Tumor_region", "NEC", "BEN"), type = "ashr")
resShrink_T1_vs_BEN = lfcShrink(dds, contrast = c("Tumor_region", "T1", "BEN"), type = "ashr")
resShrink_INF_vs_BEN = lfcShrink(dds, contrast = c("Tumor_region", "INF", "BEN"), type = "ashr")

resShrink_NEC_vs_T1 = lfcShrink(dds, contrast = c("Tumor_region", "NEC", "T1"), type = "ashr")
resShrink_NEC_vs_INF = lfcShrink(dds, contrast = c("Tumor_region", "NEC", "INF"), type = "ashr")
resShrink_T1_vs_INF = lfcShrink(dds, contrast = c("Tumor_region", "T1", "INF"), type = "ashr")





# MA-PLOT---------------------------------
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-4,4)
plotMA(resShrink_BEN_vs_NEC, xlim=xlim, ylim=ylim, main = "BEN vs NEC")
plotMA(resShrink_BEN_vs_T1, xlim=xlim, ylim=ylim, main = "BEN vs T1")
plotMA(resShrink_BEN_vs_INF, xlim=xlim, ylim=ylim, main = "BEN vs INF")

par(mfrow=c(1,3), mar=c(4,4,2,1))
ylim <- c(-4,4)
plotMA(resShrink_NEC_vs_T1,  ylim=ylim, main = "NEC vs T1")
plotMA(resShrink_NEC_vs_INF, ylim=ylim, main = "NEC vs INF")
plotMA(resShrink_T1_vs_INF, ylim=ylim, main = "T1 vs INF")




# GET SIGNIFICANT DIFF EXP ---------------------------------------
#include info for significance of diff expression 

## Try to put functions on the top as well or in another file with other 
## functions. This way you make sure that you collect them in one overview
## I do this myself as well, see:
## https://github.com/AG-Walz/Projects/blob/main/QBiC_coops/QLVMB/DEA/functions.R
getDiffexpressed <- function(Deseq_results_dataset){
  Deseq_results_dataset$significant <- "NO"
  Deseq_results_dataset$diffexp <- "NO"
  Deseq_results_dataset$significant[abs(Deseq_results_dataset$log2FoldChange) > 3 & Deseq_results_dataset$padj < 0.01] <- "SIG"
  Deseq_results_dataset$diffexp[Deseq_results_dataset$log2FoldChange > 0 & Deseq_results_dataset$significant == "SIG"] <- "UP"
  Deseq_results_dataset$diffexp[Deseq_results_dataset$log2FoldChange < 0 & Deseq_results_dataset$significant == "SIG"] <- "DOWN"
  Deseq_results_dataset$Names <- Deseq_results_dataset@rownames
  return(Deseq_results_dataset)
}
#results===========================================
res_NEC_vs_BEN = getDiffexpressed(res_NEC_vs_BEN)
res_T1_vs_BEN = getDiffexpressed(res_T1_vs_BEN)
res_INF_vs_BEN = getDiffexpressed(res_INF_vs_BEN)

res_NEC_vs_T1 = getDiffexpressed(res_NEC_vs_T1)
res_NEC_vs_INF = getDiffexpressed(res_NEC_vs_INF)
res_T1_vs_INF = getDiffexpressed(res_T1_vs_INF)




#shrink============================================
resShrink_NEC_vs_T1 = getDiffexpressed(resShrink_NEC_vs_T1)
resShrink_NEC_vs_INF = getDiffexpressed(resShrink_NEC_vs_INF)
resShrink_T1_vs_INF = getDiffexpressed(resShrink_T1_vs_INF)
resShrink_NEC_vs_BEN = getDiffexpressed(resShrink_NEC_vs_BEN)
resShrink_INF_vs_BEN = getDiffexpressed(resShrink_INF_vs_BEN)
resShrink_T1_vs_BEN= getDiffexpressed(resShrink_T1_vs_BEN)

# testing if the genes that found are in the same reihenfolge 



# DATA SETS OF DIFF EXPRESSED GENES ---------------------------------

# all vs BEN
necben_up = rownames(resShrink_NEC_vs_BEN[which(resShrink_NEC_vs_BEN$diffexp == "UP"),])
necben_down = rownames(resShrink_NEC_vs_BEN[which(resShrink_NEC_vs_BEN$diffexp == "DOWN"),])
infben_up = rownames(resShrink_INF_vs_BEN[which(resShrink_INF_vs_BEN$diffexp == "UP"),])
infben_down = rownames(resShrink_INF_vs_BEN[which(resShrink_INF_vs_BEN$diffexp == "DOWN"),])
t1ben_up = rownames(resShrink_T1_vs_BEN[which(resShrink_T1_vs_BEN$diffexp == "UP"),])  
t1ben_down = rownames(resShrink_T1_vs_BEN[which(resShrink_T1_vs_BEN$diffexp == "DOWN"),])  

# all vs INF
infnec_up = rownames(resShrink_INF_vs_NEC[which(resShrink_INF_vs_NEC$diffexp == "UP"),])
infnec_down = rownames(resShrink_INF_vs_NEC[which(resShrink_INF_vs_NEC$diffexp == "DOWN"),])
inft1_up = rownames(resShrink_INF_vs_T1[which(resShrink_INF_vs_T1$diffexp == "UP"),])
inft1_down = rownames(resShrink_INF_vs_T1[which(resShrink_INF_vs_T1$diffexp == "DOWN"),])

#all vs T1 
nect1_up = rownames(resShrink_NEC_vs_T1[which(resShrink_NEC_vs_T1$diffexp == "UP"),])
nect1_down = rownames(resShrink_NEC_vs_T1[which(resShrink_NEC_vs_T1$diffexp == "DOWN"),])

## combinde datasets + write to csv ==========================

#up regulated genes for each region 
NEC_UP = unique(c(necben_up,nect1_up, infnec_down))
T1_UP = unique(c(t1ben_up, inft1_down, nect1_down))
INF_UP = unique(c(infben_up,inft1_up, infnec_down))
BEN_UP = unique(c(necben_down, t1ben_down, infben_down))

write.csv(as.data.frame(NEC_UP), file = "/Users/cschwitalla/Documents/Intergration/NEC_UP.csv", row.names = FALSE, quote = FALSE, )
write.csv(as.data.frame(T1_UP), file = "/Users/cschwitalla/Documents/Intergration/T1_UP.csv", row.names = FALSE, quote = FALSE, )
write.csv(as.data.frame(INF_UP), file = "/Users/cschwitalla/Documents/Intergration/INF_UP.csv", row.names = FALSE, quote = FALSE, )
write.csv(as.data.frame(BEN_UP), file = "/Users/cschwitalla/Documents/Intergration/BEN_UP.csv", row.names = FALSE, quote = FALSE, )


# higly expressed genes for each region 

counts(dds, normalized=TRUE)
vsd <- vst(dds, blind = FALSE)
vsd_df = as.data.frame(assay(vsd))
coll_ann_list = colData(dds)[,c("Tumor_region")]
c(coll_ann_list)
names(vsd_df) = coll_ann_list
vsd_df = vsd_df[,order(colnames(vsd_df))]


BEN_df = vsd_df[which(startsWith(colnames(vsd_df), "BEN"))]
BEN_df$means <- rowMeans(BEN_df)
BEN_df = BEN_df[order(BEN_df$means, decreasing = TRUE),]
NEC_df = vsd_df[which(startsWith(colnames(vsd_df), "NEC"))]
NEC_df$means <- rowMeans(NEC_df)
NEC_df = NEC_df[order(NEC_df$means, decreasing = TRUE),]
INF_df = vsd_df[which(startsWith(colnames(vsd_df), "INF"))]
INF_df$means <- rowMeans(INF_df)
INF_df = INF_df[order(INF_df$means, decreasing = TRUE),]
T1_df =vsd_df[which(startsWith(colnames(vsd_df), "T1"))]
T1_df$means <- rowMeans(T1_df)
T1_df = T1_df[order(T1_df$means, decreasing = TRUE),]


write.csv(as.data.frame(rownames(BEN_df)), file = "/Users/cschwitalla/Documents/Intergration/BEN_highexpr.csv", row.names = FALSE, quote = FALSE, )
write.csv(as.data.frame(rownames(NEC_df)), file = "/Users/cschwitalla/Documents/Intergration/NEC_highexpr.csv", row.names = FALSE, quote = FALSE, )
write.csv(as.data.frame(rownames(T1_df)), file = "/Users/cschwitalla/Documents/Intergration/T1_highexpr.csv", row.names = FALSE, quote = FALSE, )
write.csv(as.data.frame(rownames(INF_df)), file = "/Users/cschwitalla/Documents/Intergration/INF_highexpr.csv", row.names = FALSE, quote = FALSE, )



# VENN DIAGRAMME ----------------------------

myCol <- brewer.pal(3, "Pastel2")

#diagram function
drawVenn <- function(set_list, titel){
  grid.newpage()
  venn <- venn.diagram(
    x = set_list,
    filename = NULL,
    lwd = 2,
    lty = "blank",
    fill = myCol,
    main = titel,
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
  )
  grid.draw(venn)
}



######
up_vs_BEN <- list(NEC = nec_up, INF = inf_up, T1 = t1_up)
down_vs_BEN <- list(NEC = rownames(resShrink_NEC_vs_BEN[which(resShrink_NEC_vs_BEN$diffexp == "DOWN"),]),
                    INF = rownames(resShrink_INF_vs_BEN[which(resShrink_INF_vs_BEN$diffexp == "DOWN"),]),
                    T1 = rownames(resShrink_T1_vs_BEN[which(resShrink_T1_vs_BEN$diffexp == "DOWN"),]))
up_vs_ben = c(up_vs_BEN)


time_up <- list(NEC = rownames(resShrink_NEC_vs_T1[which(resShrink_NEC_vs_T1$diffexp == "UP"),]),
                T1 = rownames(resShrink_T1_vs_INF[which(resShrink_T1_vs_INF$diffexp == "UP"),]),
                INF = rownames(resShrink_INF_vs_BEN[which(resShrink_INF_vs_BEN$diffexp == "UP"),]))
time_down <- list(NEC = rownames(resShrink_NEC_vs_T1[which(resShrink_NEC_vs_T1$diffexp == "DOWN"),]),
                  T1 = rownames(resShrink_T1_vs_INF[which(resShrink_T1_vs_INF$diffexp == "DOWN"),]),
                  INF = rownames(resShrink_INF_vs_BEN[which(resShrink_INF_vs_BEN$diffexp == "DOWN"),]))

# up regulation in the INfiltration zone vs nec vs T1 vs Bening
INF_up = list(INFvsNEC = rownames(resShrink_INF_vs_NEC[which(resShrink_INF_vs_NEC$diffexp == "UP"),]),
              INFvsT1 = rownames(resShrink_INF_vs_T1[which(resShrink_INF_vs_T1$diffexp == "UP"),]),
              INFvsBEN = rownames(resShrink_INF_vs_BEN[which(resShrink_INF_vs_BEN$diffexp == "UP"),]))
INF_down = list(INFvsNEC = rownames(resShrink_INF_vs_NEC[which(resShrink_INF_vs_NEC$diffexp == "DOWN"),]),
                INFvsT1 = rownames(resShrink_INF_vs_T1[which(resShrink_INF_vs_T1$diffexp == "DOWN"),]),
                INFvsBEN = rownames(resShrink_INF_vs_BEN[which(resShrink_INF_vs_BEN$diffexp == "DOWN"),]))


drawVenn(up_vs_BEN, "UP vs BEN")
drawVenn(down_vs_BEN, "DOWN vs BEN")
drawVenn(time_up, "Time UP")
drawVenn(time_down, "Time DOWN")

drawVenn(INF_up)

UP_INF_res = calculate.overlap(INF_up)

UP_INF_res
help("calculate.overlap")
displayVenn(INF_down)





venn <- Venn(up_vs_BEN)
dat = process_region_data(venn)
dat$item[7]

venn <- Venn(INF_up)
dat = process_region_data(venn)
#write.csv(dat$item[3], "INF_UP_vs_BEN.txt", row.names = FALSE)
dat$item[[3]]
small_INF_vs_BEN = resShrink_INF_vs_BEN[dat$item[[3]],]



NEC_vs_T1_up = rownames(resShrink_NEC_vs_T1)[which(resShrink_NEC_vs_T1$diffexp == "UP")]
NEC_vs_INF_up = rownames(resShrink_NEC_vs_INF)[which(resShrink_NEC_vs_INF$diffexp == "UP")]
T1_vs_INF_up = rownames(resShrink_T1_vs_INF)[which(resShrink_T1_vs_INF$diffexp == "UP")]


# GET LISTS of DIFF EXP GENES ----------------------------------------

#shrink
NEC_diff = c(rownames(resShrink_NEC_vs_BEN[which(resShrink_NEC_vs_BEN$significant == "SIG"),]), rownames(resShrink_NEC_vs_T1[which(resShrink_NEC_vs_T1$significant == "SIG"),]), rownames(resShrink_NEC_vs_INF[which(resShrink_NEC_vs_INF$significant == "SIG"),]) ) 
NEC_diff =unique(NEC_diff)
NEC_mapp = mapping_df %>% filter(mapping_df$GENENAME %in% NEC_diff)
NEC_up = c(rownames(resShrink_NEC_vs_BEN[which(resShrink_NEC_vs_BEN$diffexp == "UP"),]), rownames(resShrink_NEC_vs_T1[which(resShrink_NEC_vs_T1$diffexp == "UP"),]), rownames(resShrink_NEC_vs_INF[which(resShrink_NEC_vs_INF$diffexp == "UP"),]) ) 
NEC_up = unique(NEC_up)

NEC_UP_INF = c(rownames(resShrink_NEC_vs_INF[which(resShrink_NEC_vs_INF$diffexp == "UP"),]))
NEC_UP_T1 = c(rownames(resShrink_NEC_vs_T1[which(resShrink_NEC_vs_T1$diffexp == "UP"),]))
write.csv(as.data.frame(NEC_UP_INF), file = "NEC_up_inf.csv", row.names = FALSE, quote = FALSE, )



# VOLCANO PLOT----------------------------------------
# to do : color adjust that NO UP DOWN have color not as scale 
color_volcano = c("#e7a982", "#5b5b5b", "#82c0e7")
names(color_volcano) = c("UP", "NO", "DOWN")


volcanoPlot3 <- function(Deseq_results_dataset, title){
  ggplot(data = as.data.frame(Deseq_results_dataset), aes(x = log2FoldChange, y = -log10(padj), col = diffexp)) +
    geom_point(aes(alpha=ifelse(significant == "N", 0.1,0.7), color=as.data.frame(Deseq_results_dataset)$diffexp), show.legend = F) +
    geom_vline(xintercept = c(-3,3), col = "red") +
    geom_hline(yintercept = -log10(0.01), col = "red") +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +
    scale_color_manual(values = c(color_volcano)) +
    #geom_text_repel(aes(label=ifelse(padj <0.001 | abs(log2FoldChange) > 10, Names, "" )),
                   # color = "black", 
                    #box.padding = unit(0.3, "lines"),
                  #  point.padding = unit(0.5, "lines"),
                  #  segment.color = '#5b5b5b', alpha=0.5,
                   # max.overlaps = 10) +
    labs(x = "log fold change", y = "-log10(p adjusted value)", color="") +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.border = element_rect(color = "grey", fill=NA, size = 0.5),
          legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(fill = "transparent")) 
}


volcanoPlot3(resShrink_NEC_vs_INF, "NEC vs INF")
volcanoPlot3(resShrink_NEC_vs_T1, "NEC vs T1")
volcanoPlot3(resShrink_T1_vs_INF, "T1 vs INF")

table(resShrink_NEC_vs_T1$diffexp)


volcanoPlot3(resShrink_INF_vs_BEN, "INF vs BEN")
volcanoPlot3(resShrink_T1_vs_BEN, "T1 vs BEN")
volcanoPlot3(resShrink_NEC_vs_BEN, "NEC vs BEN")
library(EnhancedVolcano)

EnhancedVolcano(resShrink_NEC_vs_INF,lab = rownames(resShrink_NEC_vs_INF),
                x = "log2FoldChange",
                y = "padj",
                FCcutoff = 3,
                pCutoff = 0.01,
                selectLab = NA,
                legendPosition = "right",
                title = "NEC vs INF",
                legendLabels = c("NS", "Log2 FC", "Adjusted p-value","Adjusted p-value & Log2 FC"))

NT_vol = EnhancedVolcano(resShrink_NEC_vs_T1,lab = rownames(resShrink_NEC_vs_T1),
                x = "log2FoldChange",
                y = "pvalue",
                FCcutoff = 3,
                pCutoff = 0.01,
                selectLab = NA,
                title = "NEC vs T1",
                legendLabels = c("NS", "Log2 FC", "Adjusted p-value","Adjusted p-value & Log2 FC"))

TI_vol = EnhancedVolcano(resShrink_T1_vs_INF,lab = rownames(resShrink_T1_vs_INF),
                x = "log2FoldChange",
                y = "pvalue",
                FCcutoff = 3,
                pCutoff = 0.01,
                selectLab = NA,
                title = "T1 vs INF",
                legendLabels = c("NS", "Log2 FC", "Adjusted p-value","Adjusted p-value & Log2 FC"))

NB_vol = EnhancedVolcano(resShrink_NEC_vs_BEN,lab = rownames(resShrink_NEC_vs_BEN),
                x = "log2FoldChange",
                y = "pvalue",
                FCcutoff = 3,
                pCutoff = 0.01,
                selectLab = NA,
                title = "NEC vs BEN",
                legendLabels = c("NS", "Log2 FC", "Adjusted p-value","Adjusted p-value & Log2 FC"))
TB_vol = EnhancedVolcano(resShrink_T1_vs_BEN,lab = rownames(resShrink_T1_vs_BEN),
                x = "log2FoldChange",
                y = "pvalue",
                FCcutoff = 3,
                pCutoff = 0.01,
                selectLab = NA,
                title = "T1 vs BEN",
                legendLabels = c("NS", "Log2 FC", "Adjusted p-value","Adjusted p-value & Log2 FC"))
IB_vol = EnhancedVolcano(resShrink_INF_vs_BEN, 
                lab = rownames(resShrink_INF_vs_BEN),
                x = "log2FoldChange",
                y = "pvalue",
                FCcutoff = 3,pCutoff = 0.01,
                selectLab = NA,
                title = "INF vs BEN",
                legendLabels = c("NS", "Log2 FC", "Adjusted p-value","Adjusted p-value & Log2 FC"))



library(gridExtra)
library(grid)

grid.arrange(NT_vol,NI_vol,TI_vol ,ncol= 3)
grid.rect(gp = gpar(fill=NA))

grid.arrange(NB_vol,TB_vol,IB_vol ,ncol= 3)
grid.rect(gp = gpar(fill=NA))


#CLUSTERING -----------------------

#first get distances 
# Pearson correlation 

names(all_diff_test) =  metadata$MGMT#colData(dds)[,c("Tumor_region")]
names(all_diff_rlog) = colData(dds)[,c("Tumor_region")]

#dist euclidean
#column clustering 
euc_dist_vst = dist(t(all_diff_df), method = "euclidean")
euc_dist_vst_exben = dist(t(without_ben_df), method = "euclidean")
#row clustering 
euc_dist_vst_row = dist(all_diff_df, method = "euclidean")
euc_dist_vst_row_exben = dist(without_ben_df, method = "euclidean")
#euc_dist_rlog = dist(t(all_diff_rlog), method = "euclidean")
# correlation between columns
#pearson_all_diff_col = cor(all_diff_test[1:49], method = "pearson") # all diff test--> df from variancen stabilizing transformation function ( not vst() function)

#correlation between rows
#pearson_all_diff_row = cor(t(all_diff_test[1:49]), method = "pearson") # all diff test--> df from variancen stabilizing transformation function ( not vst() function)


#unsupervised hierarchical clustering
cluster_all_diff_col = hclust(as.dist(1-pearson_all_diff_col), method = "complete")
cluster_all_diff_row = hclust(as.dist(1-pearson_all_diff_row), method = "complete")

#columns
cluster_all_vst = hclust(euc_dist_vst, method = "complete")
cluster_exben_vst = hclust(euc_dist_vst_exben, method = "complete")
#cluster_all_rlog = hclust(euc_dist_rlog, method = "complete")
#rows
cluster_all_vst_rows = hclust(euc_dist_vst_row, method = "complete")
cluster_exben_vst_rows = hclust(euc_dist_vst_row_exben, method = "complete")
#plot
plot(cluster_all_vst)
plot(cluster_all_rlog)

cluster_all_diff


#calc distance of candidate genes

#to do means  + scale muss rein und condition 
distmatrix = data.frame(assay(vsd)[select,] )
names(distmatrix) <- df$`colData(dds)[, c("Tumor_region")]`
gene_dist <- dist(distmatrix)
hclust_matrix <- hclust(gene_dist, method = "complete")

plot(hclust_matrix, labels= FALSE)

gene_cluster <- cutree(hclust_matrix, k= 5) %>% enframe() %>% rename(gene = name, cluster = value)

distmatrix$gene <- rownames(distmatrix)
distmatrix$cluster <- gene_cluster$cluster

#profile plot 
distmatrix %>% ggplot(aes())





var





#TEST------------------------------------------
plotSparsity(dds)
#michael takes this as a diagnostic to see for when the majority of genes with high count have most of the row sum
#coming from individual samples
#not fit a negative binomial assumption: genes with many zeros and 
#individual very large counts are difficult to model with the negative binomial distribution.



# HEATMAP ------------------------



#normalizing the counts data with vst , varaince stabilizing transformation 
#rld <- rlog(dds, blind = FALSE)#funktioniert nicht 
vsd <- vst(dds, blind = FALSE)

#rlog_data <- assay(rld)
write.csv(as.data.frame(rlog_data), file = "/Users/cschwitalla/Documents/RNAseq_analysis/rlog_data.csv",quote = FALSE, sep = "  ", row.names = TRUE)



#test_vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

raw_data = counts(dds)
logtransf_raw = normTransform(dds)
# Prepare Column annotations ==================
#get only patient IDs without qbic code
Sample_names = colData(dds)[,c("Sample_name")]
Patient_ids = c()
for (name in Sample_names ){
  str = strsplit(name , "_")
  Patient_ids = append(Patient_ids, str[[1]][2])
}

# generate column annotation df
collumn_annotation = data.frame(Region = colData(dds)[,c("Tumor_region")], Patient = Patient_ids, Sex = metadata$Sex, MGMT_methylation = metadata$MGMT)

# get column annotation without BEN + get row number right 
collumn_annotation_without_ben = collumn_annotation[collumn_annotation$Region != "BEN",]
rownames(collumn_annotation_without_ben) <- 1:45

# color coding for annotations: ToDO!
column_colors = list(Region = c(NEC = "indianred2", T1 = "lightsalmon2", INF = "skyblue3", BEN = "darkolivegreen4"),
                     Sex = c(F = "wheat2", M = "lightsteelblue2"), MGMT_methylation = c(Methylated = "darkolivegreen3", Unmethylated = "darkorange1"))
                     #Patient = c(ZH0984 = "lightblue1" , ZH0991 = "lightcyan1", ZH0997 = "slategray1",ZH1002 = "lightsteelblue1" ,ZH1006 = "lightskyblue1",	
                      #           ZH1007 = "paleturquoise1",ZH1008 = "azure2",ZH1015 = "seashell1",ZH1019 = "honeydew1",ZH1022 = "ivory1", ZH1028 = "lemonchiffon1", 	
                      #           ZH1032 =  "bisque1", ZH1033 =  "mistyrose1", ZH1039 = "lavenderblush1" , ZH1041 =  "thistle1" )

# get the column numbers vor BEN to exclude them in the df without the bening samples
number_ben = which(collumn_annotation$Region == "BEN", arr.ind = TRUE)


#Prepare df for heatmap =========================
# get list of diff expressed genes unshrink
diff_NEC_BEN_us = res_NEC_vs_BEN$Names[res_NEC_vs_BEN$significant == "SIG"]
diff_T1_BEN_us = res_T1_vs_BEN$Names[res_T1_vs_BEN$significant == "SIG"]
diff_INF_BEN_us = res_INF_vs_BEN$Names[res_INF_vs_BEN$significant == "SIG"]
diff_NEC_INF_us = res_NEC_vs_INF$Names[res_NEC_vs_INF$significant == "SIG"]
diff_T1_INF_us = res_T1_vs_INF$Names[res_T1_vs_INF$significant == "SIG"]
diff_NEC_T1_us = res_NEC_vs_T1$Names[res_NEC_vs_T1$significant == "SIG"]


# get list of diff expressed genes shrink
diff_NEC_BEN = resShrink_NEC_vs_BEN$Names[resShrink_NEC_vs_BEN$significant == "SIG"]
diff_T1_BEN = resShrink_T1_vs_BEN$Names[resShrink_T1_vs_BEN$significant == "SIG"]
diff_INF_BEN = resShrink_INF_vs_BEN$Names[resShrink_INF_vs_BEN$significant == "SIG"]

diff_NEC_INF = resShrink_NEC_vs_INF$Names[resShrink_NEC_vs_INF$significant == "SIG"]
diff_T1_INF = resShrink_T1_vs_INF$Names[resShrink_T1_vs_INF$significant == "SIG"]
diff_NEC_T1 = resShrink_NEC_vs_T1$Names[resShrink_NEC_vs_T1$significant == "SIG"]

# all sig diff expressed genes from NEC, T1, INF, BEN
all_diff_genes = c(diff_NEC_BEN, diff_T1_BEN, diff_INF_BEN, diff_NEC_T1, diff_NEC_INF, diff_T1_INF)
#all_diff_genes_us = c(diff_NEC_BEN_us, diff_T1_BEN_us, diff_INF_BEN_us, diff_NEC_T1_us, diff_NEC_INF_us, diff_T1_INF_us)

#all sig diff expressed genes from NEC T1 INF
NTI_diff_genes = c(diff_NEC_T1, diff_T1_INF, diff_NEC_INF)
#NTI_diff_genes_us = c(diff_NEC_T1_us, diff_T1_INF_us, diff_NEC_INF_us)
# get the vsd normalized counts for diff expressed genes
all_diff_df <- as.data.frame(assay(vsd)) %>% filter(row.names(assay(vsd)) %in% all_diff_genes)
all_diff_test <- as.data.frame(assay(test_vsd)) %>% filter(row.names(assay(test_vsd)) %in% all_diff_genes)
raw_all_diff <- as.data.frame(raw_data[which(row.names(raw_data) %in% all_diff_genes),])
norm_all_diff <- as.data.frame(assay(logtransf_raw)) %>% filter(row.names(assay(logtransf_raw)) %in% all_diff_genes)
#all_diff_us <- as.data.frame(assay(test_vsd)) %>% filter(row.names(assay(test_vsd)) %in% all_diff_genes_us)
all_diff_rlog <- as.data.frame(assay(rld)) %>% filter(row.names(assay(rld)) %in% all_diff_genes)

#subset dataframe to exclude all BEN columns
NTI_diff_df <- as.data.frame(assay(vsd)) %>% filter(row.names(assay(vsd)) %in% NTI_diff_genes) 
#NTI_diff_us <- as.data.frame(assay(vsd)) %>% filter(row.names(assay(vsd)) %in% NTI_diff_genes_us) 

without_ben_df = subset(NTI_diff_df, select= -number_ben)
#without_ben_df_us = subset(NTI_diff_us, select= -number_ben)



# for each gene get the mean expressions
all_diff_df$gene_means <- rowMeans(all_diff_df)
all_diff_test$gene_means <- rowMeans(all_diff_test)
all_diff_rlog$gene_means <- rowMeans(all_diff_rlog)
without_ben_df$gene_means <- rowMeans(without_ben_df)

#substract gene means from all columns
all_diff_df[1:(ncol(all_diff_df)-1)] <- all_diff_df[1:(ncol(all_diff_df)-1)] - all_diff_df[,ncol(all_diff_df)]
all_diff_test[1:(ncol(all_diff_test)-1)] <- all_diff_test[1:(ncol(all_diff_test)-1)] - all_diff_test[,ncol(all_diff_test)]
without_ben_df[1:(ncol(without_ben_df)-1)] <- without_ben_df[1:(ncol(without_ben_df)-1)] - without_ben_df[,ncol(without_ben_df)]
all_diff_rlog[1:(ncol(all_diff_rlog)-1)] <- all_diff_rlog[1:(ncol(all_diff_rlog)-1)]-all_diff_rlog[,ncol(all_diff_rlog)]
# exclude last column and get column number right again
without_ben_df <- without_ben_df[1:45]
names(without_ben_df) <- 1:45

without_ben_df_us <- without_ben_df_us[1:45]
names(without_ben_df_us) <- 1:45




#draw heatmap function ==============================

drawHeatmap <- function(df, titel, color_list, annotation_df){
  pheatmap(df, 
           color = colorRampPalette(brewer.pal(9,"YlOrRd"))(999),
           border_color = NA,
           show_rownames = FALSE,
           show_colnames = FALSE,
           annotation_colors = color_list,
           main = titel,
           annotation_col = annotation_df,
           cluster_cols = TRUE,
           cluster_rows = TRUE
          # clustering_distance_cols = "correlation",
           #clustering_method = 
  )
}


#draw heatmap =================================0
drawHeatmap(without_ben_df, "All diff expressed genes between NEC, T1, INF", column_colors, collumn_annotation_without_ben)
drawHeatmap(without_ben_df_us, "All diff expressed genes between NEC, T1, INF", column_colors, collumn_annotation_without_ben)


drawHeatmap(all_diff_test[1:49], "All diff expressed genes between NEC, T1, INF, BEN", column_colors, collumn_annotation)
drawHeatmap(all_diff_us[1:49], "All diff expressed genes between NEC, T1, INF, BEN", column_colors, collumn_annotation)




#complex heatmap ==================
library("ComplexHeatmap")

# pheatmp for clustering --> 
library(circlize)

col_fun = colorRamp2(c(0, 11, 14), c("green", "white", "red"))

Ha = HeatmapAnnotation(Region = colData(dds)[,c("Tumor_region")], Patient = Patient_ids, Sex = metadata$Sex, MGMT_methylation = metadata$MGMT, col = column_colors)
Heatmap(all_diff_df[1:49], row_title = "genes", top_annotation = Ha,
        show_row_names = FALSE, show_column_names = FALSE, name = "counts",
        cluster_columns = cluster_all_vst, cluster_rows = cluster_all_vst_rows)#,clustering_distance_columns = "pearson")

#without ben 

metadata_exben = metadata %>% filter(!as.numeric(rownames(metadata)) %in% number_ben) 


Ha_exben = HeatmapAnnotation(Region = collumn_annotation_without_ben$Region, 
                       Patient = metadata_exben$`Source Name(s)`, 
                       Sex = metadata_exben$Sex,
                       MGMT_methylation = metadata_exben$MGMT, 
                       col = column_colors)


Heatmap(without_ben_df[1:45], row_title = "genes", top_annotation = Ha_exben,
        show_row_names = FALSE, show_column_names = FALSE, name = "counts",
        cluster_columns = cluster_exben_vst, cluster_rows = cluster_exben_vst_rows)


# res without shrinkage
Heatmap(all_diff_us[1:49], row_title = "genes", top_annotation = Ha, show_row_names = FALSE, show_column_names = FALSE, name = "counts", clustering_distance_rows = "pearson")#,clustering_distance_columns = "pearson")

#raw 
Heatmap(raw_all_diff[1:49], row_title = "genes", top_annotation = Ha, show_row_names = FALSE, show_column_names = FALSE, name = "counts", row_dend_reorder = FALSE, column_dend_reorder = FALSE, column_order = 1:49, column_title = "raw")


#unclusterd + vst
Heatmap(all_diff_test[1:49], row_title = "genes", top_annotation = Ha, show_row_names = FALSE, show_column_names = FALSE, name = "counts", row_dend_reorder = FALSE, column_dend_reorder = FALSE, column_order = 1:49, column_title =  "vst")

# only normTrabsformed -> log2(n+1)
Heatmap(norm_all_diff, row_title = "genes", top_annotation = Ha, show_row_names = FALSE, show_column_names = FALSE, name = "counts", row_dend_reorder = FALSE, column_dend_reorder = FALSE, column_order = 1:49, column_title =  "normTransformed")
#rlog transform 
#tod do --> besseres color coding 
Heatmap(all_diff_rlog, row_title = "genes", top_annotation = Ha, show_row_names = FALSE, show_column_names = FALSE, name = "counts", row_dend_reorder =FALSE, column_dend_reorder = FALSE, column_order = 1:49, column_title =  "rlog")


library("fpc")
# Marissas code reproduced
#calc number of clusters
nb <- pamk(all_diff_test[1:49], krange = 3:8, criterion = "ch", usepam = T)
nc <- nb$nc

first_heatmap <- pheatmap(raw_all_diff[1:49], scale = "row", silent = T)
#first_heat_clusters <- as.data.frame(cutree(first_heatmap@, k = nc), row.names = names(cutree(first_heatmap$tree_row, k=nc)))



#marissas code 
# Generate a heatmap the creates a cluster based on the Calinski-Harabasz index.
generateHeatmaps <- function(geneList, fileprefix, annotation, organismOrg) {
  ## target <- target[match(colnames(geneList), target$Sample.ID.Marissa),]
  # Make sure that an annotation df is created, an example is shown below.
  ## annotation <- data.frame(Group = groupID, sampleType = paste0(target$Region, "_", target$Population))
  ## rownames(annotation) <- colnames(geneList)
  # calculate the number of clusters using the Calinski-Harabasz index
  nb <- pamk(geneList, krange = 3:8, criterion = "ch", usepam = T)
  nc <- nb$nc
  
  # Generate the first heatmap with pheatmap to obtan the information to generate the second heatmap with.
  pres <- pheatmap(geneList, scale = "row", col = hmcol, fontsize_row = 0.05, fontsize_col = 6, silent = T)
  # Obtain the different clusters that are found in this heatmap by cutting the tree created by the
  # previous heatmap generation and the number of clusters that identified with pamk.
  pres.clust <- as.data.frame( cutree(pres$tree_row, k=nc), row.names = names(cutree(pres$tree_row, k=nc)))
  # Add the name of the cluster dataframe
  colnames(pres.clust) <- "Cluster"
  pdf(paste0(fileprefix, ".pdf"))
  # Create the second heatmap with the different clusters that are determined.
  pheatmap(geneList, main = "Heatmap Unsupervised Clustering ", scale = "row", show_rownames = F, show_colnames = F, cutree_rows = nc, annotation = annotation,
           annotation_row = pres.clust, annotation_names_row = F, col = hmcol)
  dev.off()
  
  # Make a matrix with the cluster number where it was defined and the beloning genesymbol.
  geneList <- cbind(geneList, cluster=pres.clust[match(rownames(geneList), rownames(pres.clust)),])
  geneList <- cbind(geneList, geneSym=mapIds(organismOrg, keys=rownames(geneList), column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first"))
  # Write the information of the matrix to a .csv file for further use. 
  write.table(geneList, file = paste0(fileprefix, ".csv"), quote = F, sep = ";", dec = ".", row.names = F)
}



library("edgeR")



# PCA -------------------------------
names(all_diff_df) = colData(dds)[,c("Tumor_region")]
plotPCA(test_vsd, c("Tumor_region"))





# GO ENRICHMENT ANALYSIS ------------------------------



# ich bracuhe df id name 


#BEN vs INF
geneList_ib = resShrink_INF_vs_BEN$log2FoldChange
names(geneList_ib) = rownames(resShrink_INF_vs_BEN)
geneList_ib = sort(geneList_ib, decreasing = TRUE)


geneList_ib_us = res_INF_vs_BEN$log2FoldChange
names(geneList_ib_us) = rownames(res_INF_vs_BEN)
geneList_ib_us = sort(geneList_ib_us, decreasing = TRUE)


us_genes = as.data.frame(geneList_ib_us)
us_genes$shrink

#NEC vs T1
geneList_nt1 = resShrink_NEC_vs_T1$log2FoldChange
names(geneList_nt1) = rownames(resShrink_NEC_vs_T1)
geneList_nt1 = sort(geneList_nt1, decreasing = TRUE)

#NEC vs INF
geneList_ni = resShrink_NEC_vs_INF$log2FoldChange
names(geneList_ni) = rownames(resShrink_NEC_vs_INF)
geneList_ni = sort(geneList_ni, decreasing = TRUE)

#T1 vs INF
geneList_t1i = resShrink_T1_vs_INF$log2FoldChange
names(geneList_t1i) = rownames(resShrink_T1_vs_INF)
geneList_t1i = sort(geneList_t1i, decreasing = TRUE)

#GO Profile for set of genes
## http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#annotating-and-exporting-results
## damit vielleicht mal annotieren 
genenames = names(geneList_nt1)

library("org.Hs.eg.db")
keytypes(org.Hs.eg.db)
#brauch ich garnicht lol
ggo <- groupGO(gene = genenames, "org.Hs.eg.db", keyType ="GENENAME", ont = "MF", readable = TRUE )
head(ggo)

##GENE SET ENRICHMENT ------------------

gse_nt1 = gseGO(geneList = geneList_nt1, "org.Hs.eg.db", ont = "MF", keyType ="GENENAME", pvalueCutoff = 0.05, pAdjustMethod = "BH")

gse_ni = gseGO(geneList = geneList_ni, "org.Hs.eg.db", ont = "BP", keyType ="GENENAME", pvalueCutoff = 0.05, pAdjustMethod = "BH")

gse_t1i = gseGO(geneList = geneList_t1i, "org.Hs.eg.db", ont = "BP", keyType ="GENENAME", pvalueCutoff = 0.05, pAdjustMethod = "BH")

gse_ib = gseGO(geneList = geneList_ib, "org.Hs.eg.db", ont = "BP", keyType ="GENENAME", pvalueCutoff = 0.05, pAdjustMethod = "BH")


# no gen set enrichment :( not for CC, MF or BP 

##GO over representation 

# gene --> list of overrepresentd genes or under --< list of diffexpgenes
#extract this list 
#T1_vs_INF_up = rownames(resShrink_T1_vs_INF)[which(resShrink_T1_vs_INF$diffexp == "UP")]
ib_up = rownames(resShrink_INF_vs_BEN)[which(resShrink_INF_vs_BEN$diffexp == "UP")]
nt1_up = rownames(resShrink_NEC_vs_T1)[which(resShrink_NEC_vs_T1$diffexp == "UP")]
ni_up = rownames(resShrink_NEC_vs_INF)[which(resShrink_NEC_vs_INF$diffexp == "UP")]
t1i_up = rownames(resShrink_T1_vs_INF)[which(resShrink_T1_vs_INF$diffexp == "UP")]



ib_df = bitr(ib_up, fromType = "SYMBOL", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)
ni_df = bitr(ni_up, fromType = "SYMBOL", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)
nt1_df = bitr(ni_up, fromType = "SYMBOL", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)
t1i_df = bitr(t1i_up, fromType = "SYMBOL", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)

ego_ib = enrichGO(gene = ib_df$ENSEMBL, universe = rownames(resShrink_INF_vs_BEN), "org.Hs.eg.db", ont = "MF", keyType ="ENSEMBL" ,pAdjustMethod = "BH", pvalueCutoff = 0.01, readable = TRUE )

ego_ni = enrichGO(gene = ni_df$ENSEMBL, universe = rownames(resShrink_NEC_vs_INF), "org.Hs.eg.db", ont = "MF", keyType ="ENSEMBL" ,pAdjustMethod = "BH", pvalueCutoff = 0.01, readable = TRUE )

ego_nt1 = enrichGO(gene = nt1_df$ENSEMBL, universe = rownames(resShrink_NEC_vs_T1), "org.Hs.eg.db", ont = "MF", keyType ="ENSEMBL" ,pAdjustMethod = "BH", pvalueCutoff = 0.01, readable = TRUE )

ego_t1i = enrichGO(gene = t1i_df$ENSEMBL, universe = rownames(resShrink_T1_vs_INF), "org.Hs.eg.db", ont = "MF", keyType ="ENSEMBL" ,pAdjustMethod = "BH", pvalueCutoff = 0.01, readable = TRUE )


##gprofiler2 ==================
library("gprofiler2")

gostres <- gost(query = up_vs_ben, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "bonferroni", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gostres, capped = FALSE, interactive = FALSE)
pp <- publish_gostplot(p,h)
publish_gosttable(gostres, highlight_terms = gostres$result[c(1:2,10,120),],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)
# PROFILE PLOTS ------------------------------

# create normalised counts tables for each tumor reagion 
vsd <- vst(dds, blind = FALSE)

ind = c(vsd@colData@listData$Tumor_region)
ben_num = c(which( ind %in% "BEN"))
t1_num = c(which( ind %in% "T1"))
nec_num = c(which( ind %in% "NEC"))
inf_num = c(which( ind %in% "INF"))

ben_norm_counts = assay(vsd)[,ben_num]
t1_norm_counts = assay(vsd)[,t1_num]
nec_norm_counts = assay(vsd)[,nec_num]
inf_norm_counts = assay(vsd)[,inf_num]

mean_counts_df <- data.frame("BEN" = rowMeans(ben_norm_counts), "NEC" = rowMeans(nec_norm_counts), "T1" = rowMeans(t1_norm_counts), "INF" = rowMeans(inf_norm_counts))


##count data transformation 
### variance stabilizing transformation vst
###regularized logarithm rlog
### remove dependece of the var on the mean 



#normalising values with vst 

vsd <- vst(dds, blind = FALSE)

set.seed(42)
rownames(as.data.frame(vsd@assays@data@listData))

geneprofiler(vsd, genelist = rownames(as.data.frame(vsd@assays@data@listData))  , intgroup = "Tumor_region", plotZ = FALSE )

pcaExplorer(dds = dds, dst = vsd)

# gene list --> all most mutated from exome sequencing
most_mut <- c("OR8U1","LNP1","CHIT1","TP53","TTN","HLA-DQB1","NF1","PTEN","EGFR","RYR3","ZNF880","HMCN1","LAMC3","PTOV1","KRT4","MROH2B","RYR2","SYNE1","ACAN","MUC16")

geneprofiler(mean_counts_df, most_mut, intgroup = "Tumor_region", plotZ = FALSE)
#genelist --> all upregulated between INF vs BEN
up_inf = c("ABCC3","ADM","BCL3","BLOC1S5-TXNDC5","CA12","COL4A1","COL5A1","ESM1","F8A2","GPRC5A","HAS2-AS1","KIF24","LIF","MMEL1","MMP19","MYO1G","PLA2G2A","PLAU","SPOCD1","STC2","TGFBI","TREM1","VEGFA")
geneprofiler(vsd, up_inf, intgroup = "Tumor_region", plotZ = FALSE)



mean_counts_df = t(mean_counts_df)
mean_counts_df[,most_mut]
rownames(mean_counts_df)
# with expression mean + ggplot



mean_counts_df = tibble::rownames_to_column(mean_counts_df,"Genes") #rownase (genes) make to column called genes
mean_c_df <- pivot_longer(mean_counts_df, cols= 2:5, names_to = "Tumor_region", values_to = "Mean_counts")


sub_df = filter(mean_c_df, mean_c_df$Genes %in% most_mut)
level_order = c("BEN", "INF", "T1", "NEC")

ggplot(sub_df, aes(x = factor(Tumor_region, level = level_order ), y = Mean_counts))+ geom_line(aes(group=Genes), size=0.5, alpha=0.3, color= brewer.pal("Set2")) + theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(color = "grey", fill=NA, size = 0.5),
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent")) 







## not important
countsdata = as.data.frame(counts(dds, normalized = FALSE, replaced = FALSE))
countsdata_norm = as.data.frame(counts(dds, normalized = TRUE, replaced = FALSE))

names(countsdata) = dds$Tumor_region

by_NEC <-countsdata[which(names(countsdata) == "NEC")]
by_T1 <- countsdata[which(names(countsdata) == "T1")]
by_INF <- countsdata[which(names(countsdata) == "INF")]
by_BEN<- countsdata[which(names(countsdata) == "BEN")]

rowMeans(by_NEC)
by_NEC$Means <- rowMeans(by_NEC)
by_T1$Means <- rowMeans(by_T1)
by_INF$Means <- rowMeans(by_INF)
by_BEN$Means <- rowMeans(by_BEN)

meanCounts <- data.frame(NEC = by_NEC$Means, T1 = by_T1$Means, INF = by_INF$Means, BEN = by_BEN$Means)

ggplot(meanCounts[1:20,])

