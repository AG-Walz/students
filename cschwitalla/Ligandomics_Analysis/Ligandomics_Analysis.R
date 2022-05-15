

setwd("/Users/cschwitalla/Documents/Immunopeptidomics/")

# library("stringr")
# library("RColorBrewer")
# library("tidyr")
# library("dplyr")
# library("tibble")
# library("ggplot2")
# library("viridis")
# library("VennDiagram")
# library("tuple")
# library("org.Hs.eg.db") 
# library("ggVennDiagram")
# library("readxl")

## Suggestion Marissa
## ---------------------
## Define the libraries that you want to use
necessaryLibs <- c("stringr", "RColorBrewer", "tidyr", "dplyr", "tibble",
                   "ggplot2", "viridis", "VennDiagram", "tuple", "org.Hs.eg.db",
                   "ggVennDiagram", "readxl")
## Install the libraries if necessary
# invisible(lapply(necessaryLibs, BiocManager::install, update = F, ask = F))
## Load the libraries
suppressMessages(invisible(lapply(necessaryLibs, library, character.only = T)))


# KNOWN BENIGN PEPTIDE DATA--------------------------------
#HLA-typing form Marcel Immunology================================
GB_HLA_types = read_xlsx("/Users/cschwitalla/Documents/Ligandomics_analysis/HLA-Typisierung_GBM.xlsx", col_names = TRUE)
#get list of unique HLA types
HLA_Types = GB_HLA_types$`HLA-Typing`

## These steps require more comments
HLA_Types_list = c()
for(string in HLA_Types){
  str1 = strsplit(string, ";", fixed = TRUE)
  for(splits in str1){
    HLA_Types_list = append(HLA_Types_list, splits)
  }
  
}
uniqe_HLA_types = unique(HLA_Types_list)
uniqe_HLA_types = as.data.frame(uniqe_HLA_types[2:41])

rewritten_hlatypes = c()
for(i in uniqe_HLA_types$`uniqe_HLA_types[2:41]`){
  str2 = strsplit(i, "*", fixed = TRUE)
  
  rewritten_hlatypes = append(rewritten_hlatypes, paste("HLA-", str2[[1]][1],str2[[1]][2], sep=""))
}

#write.table(rewritten_hlatypes, file= "/Users/cschwitalla/Documents/Ligandomics_analysis/Uniqe_HLA_Types.csv", sep = ",",quote = FALSE, row.names = FALSE, col.names = FALSE)

hla_str = ""
for(i in rewritten_hlatypes){
  hla_str = paste(hla_str, i, ",", sep="")
}
#write(hla_str, file= "/Users/cschwitalla/Documents/Ligandomics_analysis/Uniqe_HLA_Types.csv")

#Benign data Immunology -> from Marissa ==================================
#more specific
#less hits 
Benigniome_I = read.csv("/Users/cschwitalla/Documents/Ligandomics_analysis/newBenignmorespecific/Benign_class1.csv", header = FALSE, sep = ",", col.names = c("Peptide", "Acc", "V3", "Tissue", "V5"))
Benigniome_II = read.csv("/Users/cschwitalla/Documents/Ligandomics_analysis/newBenignmorespecific/Benign_class2.csv", header = FALSE, sep = ",", col.names = c("Peptide", "Acc", "V3", "Tissue", "V5","V6"))
#filter out all double mappers 
Benigniome_I_unique_acc = Benigniome_I[- grep(";",Benigniome_I$Acc),]
Benigniome_II_unique_acc = Benigniome_II[- grep(";", Benigniome_II$Acc),]


#HLA Ligand Atlas =================================================
## Try to reduce the length of your vector names, as this might lead to
## confusion later on if there is something going wrong
HLA_Ligand_atlas_df = read.csv("./hla_2020.12/HLA_aggregated.tsv", header = TRUE, sep= "\t")
HLA_Ligand_acc_df = read.csv("./hla_2020.12/HLA_protein_map.tsv", header = TRUE, sep = "\t")
#get mapping frequencie for every peptide id 
HLA_Ligand_mapperfreq = as.data.frame(table(HLA_Ligand_acc_df$peptide_sequence_id))
# get peptide ids with freq 1
Uniqe_mappers = HLA_Ligand_mapperfreq[which(HLA_Ligand_mapperfreq$Freq == 1),]
# get all acc that are uniqe mappers
HLA_Ligand_atlas_unique_acc = HLA_Ligand_acc_df %>% dplyr::filter(HLA_Ligand_acc_df$peptide_sequence_id  %in% Uniqe_mappers$Var1)
#split in class I nd class II 
HLA_Atlas_classI = HLA_Ligand_atlas_df[which((HLA_Ligand_atlas_df$hla_class == "HLA-I")|(HLA_Ligand_atlas_df$hla_class == "HLA-I+II")),]
HLA_Atlas_classII = HLA_Ligand_atlas_df[which((HLA_Ligand_atlas_df$hla_class == "HLA-II")|(HLA_Ligand_atlas_df$hla_class == "HLA-I+II")),]
# split protein acc in classI and classII 
HLA_Atlas_uniqueacc_classI = HLA_Atlas_classI %>% dplyr::filter(HLA_Atlas_classI$peptide_sequence_id %in% HLA_Ligand_atlas_unique_acc$peptide_sequence_id)
HLA_Atlas_uniqueacc_classII  = HLA_Atlas_classII %>% dplyr::filter(HLA_Atlas_classII$peptide_sequence_id %in%HLA_Ligand_atlas_unique_acc$peptide_sequence_id)

HLA_Atlas_I_uniq_acc = HLA_Ligand_atlas_unique_acc %>% dplyr::filter(HLA_Ligand_atlas_unique_acc$peptide_sequence_id %in% HLA_Atlas_uniqueacc_classI$peptide_sequence_id)
HLA_Atlas_II_uniq_acc = HLA_Ligand_atlas_unique_acc %>% dplyr::filter(HLA_Ligand_atlas_unique_acc$peptide_sequence_id %in% HLA_Atlas_uniqueacc_classII$peptide_sequence_id)

# add source protein acc to hla ligand atlas df 
HLA_Atlas_I_uniqaccs = merge(HLA_Atlas_uniqueacc_classI, HLA_Atlas_I_uniq_acc, by= "peptide_sequence_id")
HLA_Atlas_II_uniqaccs = merge(HLA_Atlas_uniqueacc_classII, HLA_Atlas_II_uniq_acc, by= "peptide_sequence_id")

#TUMOR REGION LIGANDOMICS DATA--------------------------------------

#file paths
classI_filedir = "./ClassI"
classII_filedir = "./ClassII"

#filenames with path 
files_cI <- dir(classI_filedir, recursive= TRUE, pattern = ".csv", full.names = TRUE)
files_cII <- dir(classII_filedir, recursive = TRUE, pattern = ".csv", full.names = TRUE)

#Creat dataframes for Class I and II Ligandomics data===========================

#function to read in ligandomics data from csv files 
createDataFrame <- function(filelist){
  #create empty classI dataframe 
  file_count = 0
  df <- data.frame(Patient_ID = character(),
                   Tumor_region = character(),
                   Sequence = character(),
                   Accessions = character())
  #loop over all files
  for(file in filelist){
    file_count = file_count +1
    tempfile = read.csv(file, header = TRUE, sep = "\t")
    patientid = strsplit(strsplit(file, split = "/")[[1]][3], split = "_")[[1]][1]
    tregion = strsplit(strsplit(strsplit(file, split = "/")[[1]][3], split = "_")[[1]][2], split= ".", fixed =  TRUE)[[1]][1]
    seq = tempfile$sequence
    acc = tempfile$accessions
    # temporaray df 
    temp_df = data.frame(Patient_ID = rep(patientid, times = length(seq)),
                         Sample_num = rep(paste("Sample", as.character(file_count), sep= "_"), times = length(seq) ),
                         Tumor_region = rep(tregion, times = length(seq)),
                         Sequence = seq,
                         Accessions = acc)
    # append to output df 
    df<- rbind(df, temp_df)
  }
  df$Sequence = str_replace_all(df$Sequence,"\\(Oxidation\\)", "") #get rid of the (Oxidation) string in sequences 
  df$Sequence = toupper(df$Sequence) # sequence all upper case
  return(df)
}

#read in data 
classI_df = createDataFrame(files_cI)
classII_df = createDataFrame(files_cII)



#DATA PREPERATION ----------------------------------------------

# exclude all intersections with bening db ======================
classI_filterd_df = subset(classI_df, !(Sequence %in% Benigniome_I$Peptide) & !(Sequence%in% HLA_Atlas_classI$peptide_sequence))
classII_filterd_df = subset(classII_df, !(Sequence %in% Benigniome_II$Peptide) & !(Sequence %in% HLA_Atlas_classII$peptide_sequence))


#Exclude peptides that mapp to multiple source proteins ====================================
#Benign filtered data
classI_filterd_uniqe_mappers = classI_filterd_df[- grep(";", classI_filterd_df$Accessions),]
classII_filterd_uniqe_mappers = classII_filterd_df[- grep(";", classII_filterd_df$Accessions),]
#Unfiltered data
classI_uniqe_mappers = classI_df[- grep(";", classI_df$Accessions),]
classII_uniqe_mappers = classII_df[- grep(";", classII_df$Accessions),]

#Get Tumorregion specific dataframes =======================================
#Unfilterd#############################
#multi mapper 
NEC_classI = classI_df[classI_df$Tumor_region == "NEC",]
NEC_classII= classII_df[classII_df$Tumor_region == "NEC",]

T1_classI = classI_df[classI_df$Tumor_region == "T1", ]
T1_classII = classII_df[classII_df$Tumor_region == "T1", ]

INF_classI = classI_df[classI_df$Tumor_region == "INF", ]
INF_classII = classII_df[classII_df$Tumor_region == "INF", ]

BEN_classI = classI_df[classI_df$Tumor_region == "BEN", ]
BEN_classII = classII_df[classII_df$Tumor_region == "BEN", ]

#unique mappers
## For loop with NEC, T1, INF, and BEN and use assign function
NEC_I_uniqe_acc= classI_uniqe_mappers[classI_uniqe_mappers$Tumor_region == "NEC",]
NEC_II_uniqe_acc = classII_uniqe_mappers[classII_uniqe_mappers$Tumor_region == "NEC",]

T1_I_uniqe_acc= classI_uniqe_mappers[classI_uniqe_mappers$Tumor_region == "T1", ]
T1_II_uniqe_acc= classII_uniqe_mappers[classII_uniqe_mappers$Tumor_region == "T1", ]

INF_I_uniqe_acc= classI_uniqe_mappers[classI_uniqe_mappers$Tumor_region == "INF", ]
INF_II_uniqe_acc= classII_uniqe_mappers[classII_uniqe_mappers$Tumor_region == "INF", ]

BEN_I_uniqe_acc= classI_uniqe_mappers[classI_uniqe_mappers$Tumor_region == "BEN", ]
BEN_II_uniqe_acc= classII_uniqe_mappers[classII_uniqe_mappers$Tumor_region == "BEN", ]



#Benigniome filtered #################################
#with multi mappers
NEC_classI_filterd = classI_filterd_df[classI_filterd_df$Tumor_region == "NEC",]
NEC_classII_filterd = classII_filterd_df[classII_filterd_df$Tumor_region == "NEC",]

T1_classI_filterd = classI_filterd_df[classI_filterd_df$Tumor_region == "T1", ]
T1_classII_filterd = classII_filterd_df[classII_filterd_df$Tumor_region == "T1", ]

INF_classI_filterd = classI_filterd_df[classI_filterd_df$Tumor_region == "INF", ]
INF_classII_filterd = classII_filterd_df[classII_filterd_df$Tumor_region == "INF", ]

BEN_classI_filterd = classI_filterd_df[classI_filterd_df$Tumor_region == "BEN", ]
BEN_classII_filterd = classII_filterd_df[classII_filterd_df$Tumor_region == "BEN", ]

#wunique mappers
NEC_I_uniqe_acc_filterd = classI_filterd_uniqe_mappers[classI_filterd_uniqe_mappers$Tumor_region == "NEC",]
NEC_II_uniqe_acc_filterd = classII_filterd_uniqe_mappers[classII_filterd_uniqe_mappers$Tumor_region == "NEC",]

T1_I_uniqe_acc_filterd= classI_filterd_uniqe_mappers[classI_filterd_uniqe_mappers$Tumor_region == "T1", ]
T1_II_uniqe_acc_filterd= classII_filterd_uniqe_mappers[classII_filterd_uniqe_mappers$Tumor_region == "T1", ]

INF_I_uniqe_acc_filterd= classI_filterd_uniqe_mappers[classI_filterd_uniqe_mappers$Tumor_region == "INF", ]
INF_II_uniqe_acc_filterd= classII_filterd_uniqe_mappers[classII_filterd_uniqe_mappers$Tumor_region == "INF", ]

BEN_I_uniqe_acc_filterd= classI_filterd_uniqe_mappers[classI_filterd_uniqe_mappers$Tumor_region == "BEN", ]
BEN_II_uniqe_acc_filterd= classII_filterd_uniqe_mappers[classII_filterd_uniqe_mappers$Tumor_region == "BEN", ]

Tumor_I_uniqe_acc_filterd = classI_filterd_uniqe_mappers[classI_filterd_uniqe_mappers$Tumor_region != "BEN", ]
Tumor_II_uniqe_acc_filterd = classII_filterd_uniqe_mappers[classII_filterd_uniqe_mappers$Tumor_region != "BEN", ]

#Get Tumorregion exclusive dataframes ==================================

#Class I ###################################
# multi mappers 
NEC_exclusive_classI_filterd  = subset(NEC_classI_filterd , !(Sequence %in% T1_classI_filterd$Sequence) & !(Sequence %in% INF_classI_filterd$Sequence) & !(Sequence %in% BEN_classI_filterd$Sequence))
T1_exclusive_classI_filterd = subset(T1_classI_filterd, !(Sequence %in% NEC_classI_filterd$Sequence) & !(Sequence %in% INF_classI_filterd$Sequence) & !(Sequence %in% BEN_classI_filterd$Sequence))
INF_exclusive_classI_filterd = subset(INF_classI_filterd, !(Sequence %in% T1_classI_filterd$Sequence) & !(Sequence %in% NEC_classI_filterd$Sequence) & !(Sequence %in% BEN_classI_filterd$Sequence))
BEN_exclusive_classI_filterd = subset(BEN_classI_filterd, !(Sequence %in% T1_classI_filterd$Sequence) & !(Sequence %in% INF_classI_filterd$Sequence) & !(Sequence %in% NEC_classI_filterd$Sequence))

#wunique mappers 
NEC_exclusive_I_uniqeAcc_filterd = subset(NEC_I_uniqe_acc_filterd, !(Sequence %in% T1_I_uniqe_acc_filterd$Sequence) & !(Sequence %in% INF_I_uniqe_acc_filterd$Sequence) & !(Sequence %in% BEN_I_uniqe_acc$Sequence))
T1_exclusive_I_uniqeAcc_filterd = subset(T1_I_uniqe_acc_filterd, !(Sequence %in% NEC_I_uniqe_acc_filterd$Sequence) & !(Sequence %in% INF_I_uniqe_acc_filterd$Sequence) & !(Sequence %in% BEN_I_uniqe_acc_filterd$Sequence))
INF_exclusive_I_uniqeAcc_filterd = subset(INF_I_uniqe_acc_filterd, !(Sequence %in% T1_I_uniqe_acc_filterd$Sequence) & !(Sequence %in% NEC_I_uniqe_acc_filterd$Sequence) & !(Sequence %in% BEN_I_uniqe_acc_filterd$Sequence))
BEN_exclusive_I_uniqeAcc_filterd = subset(BEN_I_uniqe_acc_filterd, !(Sequence %in% T1_I_uniqe_acc_filterd$Sequence) & !(Sequence %in% INF_I_uniqe_acc_filterd$Sequence) & !(Sequence %in% NEC_I_uniqe_acc_filterd$Sequence))


#Class II #######################################
#multi mappers
NEC_exclusive_classII_filterd = subset(NEC_classII_filterd, !(Sequence %in% T1_classII_filterd$Sequence) & !(Sequence %in% INF_classII_filterd$Sequence) & !(Sequence %in% BEN_classII_filterd$Sequence))
T1_exclusive_classII_filterd = subset(T1_classII_filterd, !(Sequence %in% NEC_classII_filterd$Sequence) & !(Sequence %in% INF_classII_filterd$Sequence) & !(Sequence %in% BEN_classII_filterd$Sequence))
INF_exclusive_classII_filterd = subset(INF_classII_filterd, !(Sequence %in% T1_classII_filterd$Sequence) & !(Sequence %in% NEC_classII_filterd$Sequence) & !(Sequence %in% BEN_classII_filterd$Sequence))
BEN_exclusive_classII_filterd = subset(BEN_classII_filterd, !(Sequence %in% T1_classII_filterd$Sequence) & !(Sequence %in% INF_classII_filterd$Sequence) & !(Sequence %in% NEC_classII_filterd$Sequence))


#unique mappers
NEC_exclusive_II_uniqeAcc_filterd = subset(NEC_II_uniqe_acc_filterd, !(Sequence %in% T1_II_uniqe_acc_filterd$Sequence) & !(Sequence %in% INF_II_uniqe_acc_filterd$Sequence) & !(Sequence %in% BEN_II_uniqe_acc_filterd$Sequence))
T1_exclusive_II_uniqeAcc_filterd = subset(T1_II_uniqe_acc_filterd, !(Sequence %in% NEC_II_uniqe_acc_filterd$Sequence) & !(Sequence %in% INF_II_uniqe_acc_filterd$Sequence) & !(Sequence %in% BEN_II_uniqe_acc_filterd$Sequence))
INF_exclusive_II_uniqeAcc_filterd = subset(INF_II_uniqe_acc_filterd, !(Sequence %in% T1_II_uniqe_acc_filterd$Sequence) & !(Sequence %in% NEC_II_uniqe_acc_filterd$Sequence) & !(Sequence %in% BEN_II_uniqe_acc_filterd$Sequence))
BEN_exclusive_I_uniqeAcc_filterd = subset(BEN_II_uniqe_acc_filterd, !(Sequence %in% T1_II_uniqe_acc_filterd$Sequence) & !(Sequence %in% INF_II_uniqe_acc_filterd$Sequence) & !(Sequence %in% NEC_II_uniqe_acc_filterd$Sequence))





#Get Protein ACC Functions--------------------------------
getProteinAcc <- function(list){
  Pacc = c()
  acc_only = c()
  
  for( i in list){
    str = strsplit(i, ";")
    Pacc = append(Pacc, str[[1]])
  }
  for (i in Pacc){
    str = strsplit(i, "|", fixed = TRUE)
    acc_only = append(acc_only, str[[1]][2])
  }
  return(acc_only)
}

getProteinAcc_uniqemappers <- function(list){
  acc_only = c()
  for (i in list){
    str = strsplit(i, "|", fixed = TRUE)
    acc_only = append(acc_only, str[[1]][2])
  }
  return(acc_only)
  
}

#PLOTS---------------------------------------------------------

#length distribution CLASS I + II ================

#TO DO: 
# - color adjust
# - % statt counts

classI_df$length = apply(classI_df,2,nchar)[,4]
table(classI_df$length)
counts_classI = aggregate(classI_df$length, by=list(Tumor_region = classI_df$Tumor_region, len = classI_df$length ), FUN=table )
table(classII_df$length)
counts_classI = counts_classI[order(counts_classI$Tumor_region),]
sum(counts_classI$x)


classII_df$length = apply(classII_df,2,nchar)[,4]
counts_classII = aggregate(classII_df$length, by=list(Tumor_region = classII_df$Tumor_region, len = classII_df$length ), FUN=table )
counts_classII = counts_classII[order(counts_classII$Tumor_region),]
sum(counts_classII$x)



#plot as bar plot 
ggplot(counts_classI, aes(x = len, y = x, group = Tumor_region, fill = Tumor_region)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = brewer.pal(4,"Set2")) +
  labs(x="Peptide length", y= "counts of HLA class I peptides") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))

ggplot(counts_classII, aes(x = len, y = x, fill = Tumor_region)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = brewer.pal(4,"Set2")) +
  labs(x="Peptide length", y= "counts of HLA class II peptides") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))


#plot 

ggplot(counts_classI ,aes(x = len, y = x, group = Tumor_region, color = Tumor_region)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(4,"Set2")) +
  labs(x="Peptide length", y= "counts of HLA class I peptides") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))

ggplot(counts_classII ,aes(x = len, y = x, group = Tumor_region, color = Tumor_region)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(4, "Set2")) +
  labs(x="Peptide length", y= "counts of HLA class II peptides") +
  theme_classic() +
  theme( axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1.2))

display.brewer.pal(n = 4,name = "Paired")


#saturation analysis ===================

# y axis = num of unique peptides
# x axis = num of samples

sat_classI <- classI_df[match(unique(classI_df$Sequence), classI_df$Sequence), ] %>%
  group_by(Sample_num) %>%
  summarise(num_unique = length(unique(Sequence))) %>%
  arrange(desc(num_unique)) %>%
  mutate(num_unique = cumsum(num_unique)) %>%
  mutate(num_samples = seq(1, length(Sample_num)))

sat_classII <- classII_df[match(unique(classII_df$Sequence), classII_df$Sequence), ] %>%
  group_by(Sample_num) %>%
  summarise(num_unique = length(unique(Sequence))) %>%
  arrange(desc(num_unique)) %>%
  mutate(num_unique = cumsum(num_unique)) %>%
  mutate(num_samples = seq(1, length(Sample_num)))

# fit curve 

fit_param_classI = NLSstAsymptotic(sortedXyData(c(0, sat_classI$num_samples), c(0, sat_classI$num_unique)))
fit_sat_classI = SSasymp(sat_classI$num_samples, abs(fit_param_classI[1]) + fit_param_classI[2], fit_param_classI[1], fit_param_classI[3])

fit_param_classII = NLSstAsymptotic(sortedXyData(c(0, sat_classII$num_samples), c(0, sat_classII$num_unique)))
fit_sat_classII = SSasymp(sat_classII$num_samples, abs(fit_param_classII[1]) + fit_param_classII[2], fit_param_classII[1], fit_param_classII[3])

fit_param_classI

#plot 
ggplot(sat_classI,aes(x = sat_classI$num_samples, y = fit_sat_classI)) +
  geom_smooth(method = "loess", se = F, color = "black",formula = y ~ x) +
  geom_hline(yintercept = abs(fit_param_classI[1])+fit_param_classI[2], linetype = "dashed", alpha = 0.6) +
  labs(x = "Number of Samples", y = "Number of unique Peptides") +
  theme(axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"))

ggplot(sat_classII,aes(x = sat_classII$num_samples, y = fit_sat_classII)) +
  geom_smooth(method = "loess", se = F, color = "black",formula = y ~ x) +
  geom_hline(yintercept = abs(fit_param_classII[1])+fit_param_classII[2], linetype = "dashed", alpha = 0.6) +
  labs(x = "Number of Samples", y = "Number of unique Peptides") +
  theme(axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 10, vjust = 1, hjust = 1, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"))









#VENN DIAGRAMS---------------------------------------------

#data sets for Venn diagrams =============================
#Unfiltered Sequences##################################
# class 1
BenDB_classI = c(Benigniome_I$Peptide, HLA_Atlas_classI$peptide_sequence)
seq_I = list(NEC = NEC_classI$Sequence, T1 =  T1_classI$Sequence, INF = INF_classI$Sequence, BEN = BEN_classI$Sequence, BENDB = BenDB_classI)
#class 2
BenDB_classII = c(Benigniome_II$Peptide, HLA_Atlas_classII$peptide_sequence)
seq_II = list(NEC = NEC_classII$Sequence, T1 =  T1_classII$Sequence, INF = INF_classII$Sequence, BEN = BEN_classII$Sequence, BENDB = BenDB_classII)

#Unfiltered unique mapper source protein accessions#############################
#class 1
setI_acc = list(NEC = getProteinAcc_uniqemappers(NEC_I_uniqe_acc$Accessions), T1 = getProteinAcc_uniqemappers(T1_I_uniqe_acc$Accessions), INF = getProteinAcc_uniqemappers(INF_I_uniqe_acc$Accessions), BEN = getProteinAcc_uniqemappers(BEN_I_uniqe_acc$Accessions), BENDB = c(Benigniome_I_unique_acc$Acc,HLA_Atlas_I_uniqaccs$uniprot_id)) #protein accessions
setI_seq = list(NEC = NEC_I_uniqe_acc$Sequence, T1 = T1_I_uniqe_acc$Sequence, INF = INF_I_uniqe_acc$Sequence, BEN = BEN_I_uniqe_acc$Sequence, BENDB = c(Benigniome_I_unique_acc$Peptide, HLA_Atlas_I_uniqaccs$peptide_sequence)) #peptide sequences from uniqe mappers

#class 2
setII_acc = list(NEC = getProteinAcc(NEC_II_uniqe_acc$Accessions), T1 = getProteinAcc(T1_II_uniqe_acc$Accessions), INF = getProteinAcc(INF_II_uniqe_acc$Accessions), BEN = getProteinAcc(BEN_II_uniqe_acc$Accessions), BENDB = c(Benigniome_II_unique_acc$Acc,HLA_Atlas_II_uniqaccs$uniprot_id)) #protein accessions
setII_seq = list(NEC = NEC_II_uniqe_acc$Sequence, T1 = T1_II_uniqe_acc$Sequence, INF = INF_II_uniqe_acc$Sequence, BEN = BEN_II_uniqe_acc$Sequence, BENDB = c(Benigniome_II_unique_acc$Peptide, HLA_Atlas_II_uniqaccs$peptide_sequence))#peptide sequences from uniqe mappers


#Draw Venn diagram ========================================

#class I acc and peptide sequences with uniqe protein mappers
ggVennDiagram(setI_acc, label_alpha = 0, label = "count", label_size = 3) + scale_fill_gradient(low = "papayawhip", high = "paleturquoise4")+ scale_color_manual(values = c(NEC = "grey", T1 = "grey", INF = "grey", BEN = "grey", BENDB = "grey")) + ggtitle("Class I source protein accessions")
ggVennDiagram(setI_seq, label_alpha = 0, label = "count", label_size = 3) + scale_fill_gradient(low = "papayawhip", high = "paleturquoise4")+ scale_color_manual(values = c(NEC = "grey", T1 = "grey", INF = "grey", BEN = "grey", BENDB = "grey")) + ggtitle("Class I peptide sequences")

#class II acc and peptide sequences with uniqe protein mappers
ggVennDiagram(setII_acc, label_alpha = 0, label = "count", label_size = 3) + scale_fill_gradient(low = "papayawhip", high = "paleturquoise4")+ scale_color_manual(values = c(NEC = "grey", T1 = "grey", INF = "grey", BEN = "grey", BENDB = "grey")) + ggtitle("Class II source protein accessions")
ggVennDiagram(setII_seq, label_alpha = 0, label = "count", label_size = 3) + scale_fill_gradient(low = "papayawhip", high = "paleturquoise4")+ scale_color_manual(values = c(NEC = "grey", T1 = "grey", INF = "grey", BEN = "grey", BENDB = "grey")) + ggtitle("Class II peptide sequences")


# class I and II peptide sequences with multi protein mappers
ggVennDiagram(seq_I, label_alpha = 0, label = "count", label_size = 3) + scale_fill_gradient(low = "papayawhip", high = "paleturquoise4")+ scale_color_manual(values = c(NEC = "grey", T1 = "grey", INF = "grey", BEN = "grey", BENDB = "grey")) + ggtitle("Class I peptide sequences; including multimappers")
ggVennDiagram(seq_II, label_alpha = 0, label = "count", label_size = 3) + scale_fill_gradient(low = "papayawhip", high = "paleturquoise4")+ scale_color_manual(values = c(NEC = "grey", T1 = "grey", INF = "grey", BEN = "grey", BENDB = "grey")) + ggtitle("Class II peptide sequences; including multimappers")

set_filter = list(NEC = NEC_I_uniqe_acc_filterd$Accessions, T1 = T1_I_uniqe_acc_filterd$Accessions, INF = INF_I_uniqe_acc_filterd$Accessions, BEN = Benigniome_I_unique_acc$Acc, BENDB = c(Benigniome_I_unique_acc$Acc,HLA_Atlas_I_uniqaccs$uniprot_id))
setII_filter = list(NEC = NEC_II_uniqe_acc_filterd$Accessions, T1 = T1_II_uniqe_acc_filterd$Accessions, INF = INF_II_uniqe_acc_filterd$Accessions, BEN = Benigniome_II_unique_acc$Acc, BENDB = c(Benigniome_II_unique_acc$Acc,HLA_Atlas_II_uniqaccs$uniprot_id))

ggVennDiagram(set_filter, label_alpha = 0, label = "count", label_size = 3) + scale_fill_gradient(low = "papayawhip", high = "paleturquoise4")+ scale_color_manual(values = c(NEC = "grey", T1 = "grey", INF = "grey", BEN = "grey")) + ggtitle("Class I protein acc, beningiome filterd")
ggVennDiagram(setII_filter, label_alpha = 0, label = "count", label_size = 3) + scale_fill_gradient(low = "papayawhip", high = "paleturquoise4")+ scale_color_manual(values = c(NEC = "grey", T1 = "grey", INF = "grey", BEN = "grey")) + ggtitle("Class II protein acc, beningiome filterd")


#WATERFALL PLOTS -----------------------------------------------------------------

# FUNCTIONS ======================================
#plots with protein acc
makeWaterfallDF_acc <- function(df1, df2) {
  
  #unique acc per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample 
  df1 = df1 %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
  df2 = df2 %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
  
  # union of unique acc --> x achse 
  uniq_acc <- unique(c(df1$Accessions, df2$Accessions))
  
  # neues df 
  waterfall_df = data.frame(uniq_acc)
  # ratio of peptide occurences für df1 und df2
  waterfall_df$set1 = sapply(uniq_acc, function(acc){
    mean(ifelse(is.na(matchAll(acc, df1$Accessions)), 0,
                (length(matchAll(acc, df1$Accessions))/length(unique(df1$Patient_ID)))*100))
  })
  
  waterfall_df$set2 = sapply(uniq_acc, function(acc){
    mean(ifelse(is.na(matchAll(acc, df2$Accessions)), 0,
                (length(matchAll(acc, df2$Accessions))/length(unique(df2$Patient_ID)))*100))
  })
  
  waterfall_df$Total = sapply(uniq_acc, function(acc){
    (length(matchAll(acc, c(df1$Accessions, df2$Accessions)))/
       length(c(unique(df1$Patient_ID), unique(df2$Patient_ID))))*100
  })
  
  # filtern des df 
  exclusive_1 = waterfall_df[which(waterfall_df$set2 == 0),]
  exclusive_2 = waterfall_df[which(waterfall_df$set1 == 0),]
  #ordnen
  exclusive_1 = exclusive_1[order(exclusive_1$set1, decreasing = T),]
  exclusive_2 = exclusive_2[order(exclusive_2$set2),]
  
  other = waterfall_df[which(waterfall_df$set1 != 0 & waterfall_df$set2 != 0),]
  other = other[order(other$set1, other$set2),]
  #create new waterfall df 
  waterfall_df = rbind(exclusive_1, other, exclusive_2)
  return(waterfall_df)
}
plotWaterfall_acc <- function(waterfall_df, name1, name2, title) {
  #color_legend = c( name1 = "salmon1", name2= "steelblue3")
  w = waterfall_df %>% mutate(uniq_acc = factor(uniq_acc, levels = unique(uniq_acc)))
  ggplot(w ,aes(x = uniq_acc)) +
    geom_bar( aes(y = set1, fill = name1), stat = "identity", width = 11) +
    geom_bar( aes(y = -set2, fill = name2), stat = "identity", width = 1) +
    scale_y_continuous(limits = c(-100,100), breaks = seq(-100, 100, 10)) +
    scale_fill_manual(values =  c("salmon1", "steelblue3")) +
    ylab("Frequency of positive ligandomes (%)") +
    xlab("Source Proteins") +
    labs(title = title) +       
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
}

# plots with peptide sequences
makeWaterfallDF <- function(df1, df2) {
  # union of unique sequencs --> x achse 
  uniq_seq <- unique(c(df1$Sequence, df2$Sequence))
  #unique sequence per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample 
  df1 = df1 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
  df2 = df2 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
  # neues df 
  waterfall_df = data.frame(uniq_seq)
  # ratio of peptide occurences für NEC und INF 
  waterfall_df$set1 = sapply(uniq_seq, function(seq){
    mean(ifelse(is.na(matchAll(seq, df1$Sequence)), 0,
                (length(matchAll(seq, df1$Sequence))/length(unique(df1$Patient_ID)))*100))
  })
  
  waterfall_df$set2 = sapply(uniq_seq, function(seq){
    mean(ifelse(is.na(matchAll(seq, df2$Sequence)), 0,
                (length(matchAll(seq, df2$Sequence))/length(unique(df2$Patient_ID)))*100))
  })
  
  waterfall_df$Total = sapply(uniq_seq, function(seq){
    (length(matchAll(seq, c(df1$Sequence, df2$Sequence)))/
       length(c(unique(df1$Patient_ID), unique(df2$Patient_ID))))*100
  })
  
  # filtern des df 
  exclusive_1 = waterfall_df[which(waterfall_df$set2 == 0),]
  exclusive_2 = waterfall_df[which(waterfall_df$set1 == 0),]
  #ordnen
  exclusive_1 = exclusive_1[order(exclusive_1$set1, decreasing = T),]
  exclusive_2 = exclusive_2[order(exclusive_2$set2),]
  
  other = waterfall_df[which(waterfall_df$set1 != 0 & waterfall_df$set2 != 0),]
  other = other[order(other$set1, other$set2),]
  #create new waterfall df 
  waterfall_df = rbind(exclusive_1, other, exclusive_2)
  return(waterfall_df)
}
plotWaterfall <- function(waterfall_df, name1, name2, title) {
  #color_legend = c( name1 = "salmon1", name2= "steelblue3")
  w = waterfall_df %>% mutate(uniq_seq = factor(uniq_seq, levels = unique(uniq_seq)))
  ggplot(w ,aes(x = uniq_seq)) +
    geom_bar( aes(y = set1, fill = name1), stat = "identity", width = 1) +
    geom_bar( aes(y = -set2, fill = name2), stat = "identity", width = 1) +
    scale_y_continuous(limits = c(-100,100), breaks = seq(-100, 100, 10)) +
    scale_fill_manual(values =  c("salmon1", "steelblue3")) +
    ylab("Frequency of positive ligandomes (%)") +
    xlab("Sequences") +
    labs(title = title) +       
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
}

#venn diagram for waterfall 
waterfall_venn <- function(listarea1,listarea2, category){
  grid.newpage()
  VennDiagram::draw.pairwise.venn(area1 = length(unique(listarea1)),area2 = length(unique(listarea2)),
                                  cross.area = length(Reduce(dplyr::intersect, list(unique(listarea1), unique(listarea2)))),
                                  #category = category, 
                                  fill = c("salmon1", "steelblue3"),
                                  cat.pos = c(0, 0), cat.dist = c(.05, .05),
                                  cex = 1.5,
                                  cat.cex = 1.5)
  
}

# get most frequent region exclusive peptides / proteins and most shared between all regions 
make_multi_waterfall_df <- function(df1, df2, df3){
  uniq_seq <- unique(c(df1$Sequence, df2$Sequence, df3$Sequence))
  
  #unique sequence per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample 
  df1 = df1 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
  df2 = df2 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
  df3 = df3 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
  # neues df 
  waterfall_df = data.frame(uniq_seq)
  # ratio of peptide occurences für NEC und INF 
  waterfall_df$set1 = sapply(uniq_seq, function(seq){
    mean(ifelse(is.na(matchAll(seq, df1$Sequence)), 0,
                (length(matchAll(seq, df1$Sequence))/length(unique(df1$Patient_ID)))*100))
  })
  
  waterfall_df$set2 = sapply(uniq_seq, function(seq){
    mean(ifelse(is.na(matchAll(seq, df2$Sequence)), 0,
                (length(matchAll(seq, df2$Sequence))/length(unique(df2$Patient_ID)))*100))
  })
  
  waterfall_df$set3 = sapply(uniq_seq, function(seq){
    mean(ifelse(is.na(matchAll(seq, df3$Sequence)), 0,
                (length(matchAll(seq, df3$Sequence))/length(unique(df3$Patient_ID)))*100))
  })
  

  
  waterfall_df$Total = sapply(uniq_seq, function(seq){
    (length(matchAll(seq, c(df1$Sequence, df2$Sequence, df3$Sequence)))/
       length(c(unique(df1$Patient_ID), unique(df2$Patient_ID), unique(df3$Patient_ID))))*100
  })
  
  
  # filtern des df 
  exclusive_1 = waterfall_df[which(waterfall_df$set2 == 0 & waterfall_df$set3 == 0 ),]
  exclusive_2 = waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set3 == 0 ),]
  exclusive_3 = waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set2 == 0 ),]
  #exclusive_4 = waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set2 == 0 & waterfall_df$set3 == 0),]
  
 
  shared = waterfall_df[which(waterfall_df$set1 != 0 & waterfall_df$set2 != 0 & waterfall_df$set3 != 0),]
  
  all = waterfall_df
  

  #create new waterfall df waterfall_df = rbind(exclusive_1, exclusive_2, exclusive_3, exclusive_4, shared)
  waterfall_df2 = list(Nec = exclusive_1, T1 = exclusive_2, INF = exclusive_3, Shared = shared, all = all)
  return(waterfall_df2)
  
}

make_multi_waterfall_acc <- function(df1, df2, df3){
  uniq_acc <- unique(c(df1$Accessions, df2$Accessions, df3$Accessions))
  
  #unique sequence per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample 
  df1 = df1 %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
  df2 = df2 %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
  df3 = df3 %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
  # neues df 
  waterfall_df = data.frame(uniq_acc)
  # ratio of peptide occurences für NEC und INF 
  waterfall_df$set1 = sapply(uniq_acc, function(acc){
    mean(ifelse(is.na(matchAll(acc, df1$Accessions)), 0,
                (length(matchAll(acc, df1$Accessions))/length(unique(df1$Patient_ID)))*100))
  })
  
  waterfall_df$set2 = sapply(uniq_acc, function(acc){
    mean(ifelse(is.na(matchAll(acc, df2$Accessions)), 0,
                (length(matchAll(acc, df2$Accessions))/length(unique(df2$Patient_ID)))*100))
  })
  
  waterfall_df$set3 = sapply(uniq_acc, function(acc){
    mean(ifelse(is.na(matchAll(acc, df3$Accessions)), 0,
                (length(matchAll(acc, df3$Accessions))/length(unique(df3$Patient_ID)))*100))
  })
  
  
  
  waterfall_df$Total = sapply(uniq_acc, function(acc){
    (length(matchAll(acc, c(df1$Accessions, df2$Accessions, df3$Accessions)))/
       length(c(unique(df1$Patient_ID), unique(df2$Patient_ID), unique(df3$Patient_ID))))*100
  })
  
  
  # filtern des df 
  exclusive_1 = waterfall_df[which(waterfall_df$set1 != 0 ),]
  exclusive_2 = waterfall_df[which(waterfall_df$set2 != 0),]
  exclusive_3 = waterfall_df[which(waterfall_df$set3 != 0),]
  #exclusive_4 = waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set2 == 0 & waterfall_df$set3 == 0),]
  
  
  #shared = waterfall_df[which(waterfall_df$set1 != 0 & waterfall_df$set2 != 0 & waterfall_df$set3 != 0),]
  
  #all = waterfall_df
  
  
  #create new waterfall df waterfall_df = rbind(exclusive_1, exclusive_2, exclusive_3, exclusive_4, shared)
  waterfall_df2 = list(Nec = exclusive_1, T1 =  exclusive_2, INF =  exclusive_3)
  return(waterfall_df2)
  
  
}




#MAKE WATERFALL DATAFRAMES==============================
#waterfallDF for sequences for uniqe protein mappers 
NEC_T1_w_df_uniqacc_seq = makeWaterfallDF(NEC_I_uniqe_acc_filterd,T1_I_uniqe_acc_filterd)
NEC_INF_w_df_uniqacc_seq = makeWaterfallDF(NEC_I_uniqe_acc_filterd,INF_I_uniqe_acc_filterd)
INF_T1_w_df_uniqacc_seq = makeWaterfallDF(INF_I_uniqe_acc_filterd, T1_I_uniqe_acc_filterd)

NEC_T1_w_df_uniqacc_seq_2 = makeWaterfallDF(NEC_II_uniqe_acc_filterd,T1_II_uniqe_acc_filterd)
NEC_INF_w_df_uniqacc_seq_2 = makeWaterfallDF(NEC_II_uniqe_acc_filterd,INF_II_uniqe_acc_filterd)
INF_T1_w_df_uniqacc_seq_2 = makeWaterfallDF(INF_II_uniqe_acc_filterd, T1_II_uniqe_acc_filterd)

#WaterfallDF for source proteins for uniqe protein mappers
NEC_T1_w_df_acc = makeWaterfallDF_acc(NEC_I_uniqe_acc_filterd,T1_I_uniqe_acc_filterd)
NEC_INF_w_df_acc = makeWaterfallDF_acc(NEC_I_uniqe_acc_filterd,INF_I_uniqe_acc_filterd)
INF_T1_w_df_acc = makeWaterfallDF_acc(INF_I_uniqe_acc_filterd, T1_I_uniqe_acc_filterd)

NEC_T1_w_df_acc_2 = makeWaterfallDF_acc(NEC_II_uniqe_acc_filterd,T1_II_uniqe_acc_filterd)
NEC_INF_w_df_acc_2 = makeWaterfallDF_acc(NEC_II_uniqe_acc_filterd,INF_II_uniqe_acc_filterd)
INF_T1_w_df_acc_2 = makeWaterfallDF_acc(INF_II_uniqe_acc_filterd, T1_II_uniqe_acc_filterd)


#PLOT WATERFALL================================
#PLOT waterfall for sequences for uniqe protein mappers 
plotWaterfall(NEC_T1_w_df_uniqacc_seq, name1 = "NEC", name2 =  "T1", "NEC vs T1 class I uniq acc, peptides")
plotWaterfall(NEC_INF_w_df_uniqacc_seq, name1 = "NEC", name2 =  "INF", "NEC vs INF class I uniq acc, peptides")
plotWaterfall(INF_T1_w_df_uniqacc_seq, name1 = "INF", name2 =  "T1", "INF vs T1 class I uniq acc, peptides")

plotWaterfall(NEC_T1_w_df_uniqacc_seq_2, name1 = "NEC", name2 =  "T1", "NEC vs T1 class II uniq acc, peptides")
plotWaterfall(NEC_INF_w_df_uniqacc_seq_2, name1 = "NEC", name2 =  "INF", "NEC vs INF class II uniq acc, peptides")
plotWaterfall(INF_T1_w_df_uniqacc_seq_2, name1 = "INF", name2 =  "T1", "INF vs T1 class II uniq acc, peptides")

#PLOT waterfall for source proteins for uniqe protein mappers
plotWaterfall_acc(NEC_T1_w_df_acc, name1 = "NEC", name2 =  "T1", "NEC vs T1 class I, uniq acc")
plotWaterfall_acc(NEC_INF_w_df_acc, name1 = "NEC", name2 =  "INF", "NEC vs INF class I, uniq acc")
plotWaterfall_acc(INF_T1_w_df_acc, name1 = "INF", name2 =  "T1", "INF vs T1 class I, uniq acc")


plotWaterfall_acc(NEC_T1_w_df_acc_2, name1 = "NEC", name2 =  "T1", "NEC vs T1 class II, uniq acc")
plotWaterfall_acc(NEC_INF_w_df_acc_2, name1 = "NEC", name2 =  "INF", "NEC vs INF class II, uniq acc")
plotWaterfall_acc(INF_T1_w_df_acc_2, name1 = "INF", name2 =  "T1", "INF vs T1 class II, uniq acc")




#VENN DIAGRAMS FOR WATERFALL==============================
#PLOT VENN for sequences for uniqe protein mappers 
waterfall_venn(NEC_I_uniqe_acc_filterd$Sequence, T1_I_uniqe_acc_filterd$Sequence, c("NEC", "T1"))
waterfall_venn(NEC_I_uniqe_acc_filterd$Sequence, INF_I_uniqe_acc_filterd$Sequence, c("NEC", "INF"))
waterfall_venn(INF_I_uniqe_acc_filterd$Sequence, T1_I_uniqe_acc_filterd$Sequence, c("INF", "T1"))

waterfall_venn(NEC_II_uniqe_acc_filterd$Sequence, T1_II_uniqe_acc_filterd$Sequence, c("NEC", "T1"))
waterfall_venn(NEC_II_uniqe_acc_filterd$Sequence, INF_II_uniqe_acc_filterd$Sequence, c("NEC", "INF"))
waterfall_venn(INF_II_uniqe_acc_filterd$Sequence, T1_II_uniqe_acc_filterd$Sequence, c("INF", "T1"))


#PLOT VENN for source proteins for uniqe protein mappers
waterfall_venn(NEC_I_uniqe_acc_filterd$Accessions, T1_I_uniqe_acc_filterd$Accessions, c("NEC", "T1"))
waterfall_venn(NEC_I_uniqe_acc_filterd$Accessions, INF_I_uniqe_acc_filterd$Accessions, c("NEC", "INF"))
waterfall_venn(INF_I_uniqe_acc_filterd$Accessions, T1_I_uniqe_acc_filterd$Accessions, c("INF", "T1"))

waterfall_venn(NEC_II_uniqe_acc_filterd$Accessions, T1_II_uniqe_acc_filterd$Accessions, c("NEC", "T1"))
waterfall_venn(NEC_II_uniqe_acc_filterd$Accessions, INF_II_uniqe_acc_filterd$Accessions, c("NEC", "INF"))
waterfall_venn(INF_II_uniqe_acc_filterd$Accessions, T1_II_uniqe_acc_filterd$Accessions, c("INF", "T1"))



#-MULTI WATERFALL ================================

#dataset for region exclusive peptides with highest frequency class I##############################
all_waterfall_withoutBEN = make_multi_waterfall_df(NEC_I_uniqe_acc_filterd, T1_I_uniqe_acc_filterd, INF_I_uniqe_acc_filterd)
NEC = all_waterfall_withoutBEN$Nec
T1 = all_waterfall_withoutBEN$T1
INF = all_waterfall_withoutBEN$INF 

all_waterfall_df_withoutBEn = make_multi_waterfall_acc(NEC_I_uniqe_acc_filterd, T1_I_uniqe_acc_filterd, INF_I_uniqe_acc_filterd)
NEC_acc = all_waterfall_df_withoutBEn$Nec
names(NEC_acc)[1]= "Accessions"
T1_acc = all_waterfall_df_withoutBEn$T1
names(T1_acc)[1] = "Accessions"
INF_acc = all_waterfall_df_withoutBEn$INF
names(INF_acc)[1] = "Accessions"


#classII region exclusive peptides with highest frequency 
all_waterfall_classII_without_BEN = make_multi_waterfall_df(NEC_II_uniqe_acc_filterd, T1_II_uniqe_acc_filterd, INF_II_uniqe_acc_filterd)
NEC_II = all_waterfall_classII_without_BEN$Nec
T1_II = all_waterfall_classII_without_BEN$T1
INF_II = all_waterfall_classII_without_BEN$INF 
shared_without_BEN_II_new = all_waterfall_classII_without_BEN$Shared
all_II = all_waterfall_classII_without_BEN$all

all_waterfall_classII_without_BEN_acc = make_multi_waterfall_acc(NEC_II_uniqe_acc_filterd, T1_II_uniqe_acc_filterd, INF_II_uniqe_acc_filterd)
NEC_acc_II = all_waterfall_classII_without_BEN_acc$Nec
names(NEC_acc_II)[1] = "Accessions"

T1_acc_II = all_waterfall_classII_without_BEN_acc$T1
names(T1_acc_II)[1] = "Accessions"

INF_acc_II = all_waterfall_classII_without_BEN_acc$INF
names(INF_acc_II)[1] = "Accessions"


#get Tumor peptide frequencies

uniqe_I = unique(Tumor_I_uniqe_acc_filterd$Sequence)
uniqe_II = unique(Tumor_II_uniqe_acc_filterd$Sequence)

Tumor_I = Tumor_I_uniqe_acc_filterd %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))
Tumor_II = Tumor_II_uniqe_acc_filterd %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))

Tumor_I_freq = data.frame(uniqe_I)
Tumor_II_freq = data.frame(uniqe_II)

Tumor_I_freq$Frequency = sapply(uniqe_I, function(seq){
  mean(ifelse(is.na(matchAll(seq, Tumor_I$Sequence)), 0,
              (length(matchAll(seq, Tumor_I$Sequence))/length(unique(Tumor_I$Patient_ID)))*100))
})

Tumor_II_freq$Frequency = sapply(uniqe_II, function(seq){
  mean(ifelse(is.na(matchAll(seq, Tumor_II$Sequence)), 0,
              (length(matchAll(seq, Tumor_II$Sequence))/ length(unique(Tumor_II$Patient_ID))) *100))
})

# get Tumor protein accesion frequencies
uniqe_I_acc = unique(Tumor_I_uniqe_acc_filterd$Accessions)
uniqe_II_acc = unique(Tumor_II_uniqe_acc_filterd$Accessions)

Tumor_I_acc = Tumor_I_uniqe_acc_filterd %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))
Tumor_II_acc = Tumor_II_uniqe_acc_filterd %>% group_by(Patient_ID) %>% summarise(Accessions = unique(Accessions))

Tumor_I_freq_acc = data.frame(uniqe_I_acc)
Tumor_II_freq_acc = data.frame(uniqe_II_acc)

Tumor_I_freq_acc$Frequency = sapply(uniqe_I_acc, function(acc){
  mean(ifelse(is.na(matchAll(acc, Tumor_I_acc$Accessions)), 0,
              (length(matchAll(acc, Tumor_I_acc$Accessions))/length(unique(Tumor_I_acc$Patient_ID)))*100))
})

Tumor_II_freq_acc$Frequency = sapply(uniqe_II_acc, function(acc){
  mean(ifelse(is.na(matchAll(acc, Tumor_II_acc$Accessions)), 0,
              (length(matchAll(acc, Tumor_II_acc$Accessions))/ length(unique(Tumor_II_acc$Patient_ID))) *100))
})

# GET PEPTIDES AND PROTEIN ACC and write to tsv###########################################
#class I 
NEC_peptide_list = NEC[which(NEC$set1 > 10),1:2]
names(NEC_peptide_list) = c("Sequence", "Peptide_frequency")
T1_peptide_list = T1[which(T1$set2 >10), c(1,3)]
names(T1_peptide_list) = c("Sequence", "Peptide_frequency")
INF_peptide_list = INF[which(INF$set3 > 10), c(1,4)]
names(INF_peptide_list) = c("Sequence", "Peptide_frequency")
Tumor_I_freq = Tumor_I_freq[which(Tumor_I_freq$Frequency > 10),]
names(Tumor_I_freq) = c("Sequence", "Peptide_frequency")

#class II
NEC_peptide_list_II = NEC_II[which(NEC_II$set1 > 10),1:2]
names(NEC_peptide_list_II) = c("Sequence", "Peptide_frequency")
T1_peptide_list_II = T1_II[which(T1_II$set2 >10), c(1,3)]
names(T1_peptide_list_II) = c("Sequence", "Peptide_frequency")
INF_peptide_list_II = INF_II[which(INF_II$set3 > 10), c(1,4)]
names(INF_peptide_list_II) = c("Sequence", "Peptide_frequency")
Tumor_II_freq = Tumor_II_freq[which(Tumor_II_freq$Frequency > 10),]
names(Tumor_II_freq) = c("Sequence", "Peptide_frequency")


names(Tumor_I_freq_acc) = c("Accessions", "Protein_frequency")
names(Tumor_II_freq_acc) = c("Accessions", "Protein_frequency")
#get protein accesions frequencies

NEC_peptide_list_acc = unique(subset(NEC_I_uniqe_acc_filterd[4:5], (NEC_I_uniqe_acc_filterd$Sequence %in% NEC_peptide_list$Sequence)))
T1_peptide_list_acc = unique(subset(T1_I_uniqe_acc_filterd[4:5], (T1_I_uniqe_acc_filterd$Sequence %in% T1_peptide_list$Sequence)))
INF_peptide_list_acc = unique(subset(INF_I_uniqe_acc_filterd[4:5], (INF_I_uniqe_acc_filterd$Sequence %in% INF_peptide_list$Sequence)))
Tumor_I_peptides_acc = unique(subset(Tumor_I_uniqe_acc_filterd[4:5], (Tumor_I_uniqe_acc_filterd$Sequence %in% Tumor_I_freq$Sequence)))

NEC_peptide_list_II_acc = unique(subset(NEC_II_uniqe_acc_filterd[4:5], (NEC_II_uniqe_acc_filterd$Sequence %in% NEC_peptide_list_II$Sequence)))
T1_peptide_list_II_acc = unique(subset(T1_II_uniqe_acc_filterd[4:5], (T1_II_uniqe_acc_filterd$Sequence %in% T1_peptide_list_II$Sequence)))
INF_peptide_list_II_acc = unique(subset(INF_II_uniqe_acc_filterd[4:5], (INF_II_uniqe_acc_filterd$Sequence %in% INF_peptide_list_II$Sequence)))
Tumor_II_peptides_acc = unique(subset(Tumor_II_uniqe_acc_filterd[4:5], (Tumor_II_uniqe_acc_filterd$Sequence %in% Tumor_II_freq$Sequence)))

#merge peptide frequencies and protein freaquencies 

NEC_I_frequend_peptides = merge(NEC_peptide_list, NEC_peptide_list_acc, by = "Sequence", all.x = TRUE)
NEC_I_frequend_peptides = merge(NEC_I_frequend_peptides, NEC_acc[1:2], by = "Accessions", all.x = TRUE)
NEC_I_frequend_peptides$Accessions = getProteinAcc_uniqemappers(NEC_I_frequend_peptides$Accessions)
names(NEC_I_frequend_peptides)[4] = "Protein_frequency"

T1_I_frequend_peptides = merge(T1_peptide_list, T1_peptide_list_acc, by = "Sequence", all.x = TRUE)
T1_I_frequend_peptides = merge(T1_I_frequend_peptides, T1_acc[c(1,3)], by = "Accessions", all.x = TRUE)
T1_I_frequend_peptides$Accessions = getProteinAcc_uniqemappers(T1_I_frequend_peptides$Accessions)
names(T1_I_frequend_peptides)[4] = "Protein_frequency"

INF_I_frequend_peptides = merge(INF_peptide_list, INF_peptide_list_acc, by = "Sequence", all.x = TRUE)
INF_I_frequend_peptides = merge(INF_I_frequend_peptides, INF_acc[c(1,4)], by = "Accessions", all.x = TRUE)
INF_I_frequend_peptides$Accessions = getProteinAcc_uniqemappers(INF_I_frequend_peptides$Accessions)
names(INF_I_frequend_peptides)[4] = "Protein_frequency"

NEC_II_frequend_peptides = merge(NEC_peptide_list_II, NEC_peptide_list_II_acc, by = "Sequence", all.x = TRUE)
NEC_II_frequend_peptides = merge(NEC_II_frequend_peptides, NEC_acc_II[1:2], by = "Accessions", all.x = TRUE)
NEC_II_frequend_peptides$Accessions = getProteinAcc_uniqemappers(NEC_II_frequend_peptides$Accessions)
names(NEC_II_frequend_peptides)[4] = "Protein_frequency"


T1_II_frequend_peptides = merge(T1_peptide_list_II, T1_peptide_list_II_acc, by = "Sequence", all.x = TRUE)
T1_II_frequend_peptides = merge(T1_II_frequend_peptides, T1_acc_II[c(1,3)], by = "Accessions", all.x = TRUE)
T1_II_frequend_peptides$Accessions = getProteinAcc_uniqemappers(T1_II_frequend_peptides$Accessions)
names(T1_II_frequend_peptides)[4] = "Protein_frequency"

INF_II_frequend_peptides = merge(INF_peptide_list_II, INF_peptide_list_II_acc, by = "Sequence", all.x = TRUE)
INF_II_frequend_peptides = merge(INF_II_frequend_peptides, INF_acc_II[c(1,4)], by = "Accessions", all.x = TRUE)
INF_II_frequend_peptides$Accessions = getProteinAcc_uniqemappers(INF_II_frequend_peptides$Accessions)
names(INF_II_frequend_peptides)[4] = "Protein_frequency"

Tumor_I_frequent_peptides = merge(Tumor_I_freq, unique(Tumor_I_uniqe_acc_filterd[4:5]), by = "Sequence", all.x = TRUE)
Tumor_I_frequent_peptides = merge(Tumor_I_frequent_peptides, Tumor_I_freq_acc, by = "Accessions", all.x = TRUE )
Tumor_I_frequent_peptides$Accessions = getProteinAcc_uniqemappers(Tumor_I_frequent_peptides$Accessions)

Tumor_II_frequent_peptides = merge(Tumor_II_freq, unique(Tumor_II_uniqe_acc_filterd[4:5]), by = "Sequence", all.x = TRUE)
Tumor_II_frequent_peptides = merge(Tumor_II_frequent_peptides, Tumor_II_freq_acc, by = "Accessions", all.x = TRUE )
Tumor_II_frequent_peptides$Accessions = getProteinAcc_uniqemappers(Tumor_II_frequent_peptides$Accessions)


write.table(NEC_I_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_I_frequend_peptides.tsv")
write.table(T1_I_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_I_frequend_peptides.tsv")
write.table(INF_I_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_I_frequend_peptides.tsv")
write.table(Tumor_I_frequent_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_I_frequend_peptides.tsv")
write.table(NEC_II_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_II_frequend_peptides.tsv")
write.table(T1_II_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_II_frequend_peptides.tsv")
write.table(INF_II_frequend_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_II_frequend_peptides.tsv")
write.table(Tumor_II_frequent_peptides,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_II_frequend_peptides.tsv")





# get protein acc for peptids and export it as csv or tsf 
write_region_pepacc_df <- function(region_peptide_list, region_excl_acc, filename){
  names(region_peptide_list)[1] = "Sequence"
  names(region_peptide_list)[2] = "Frequency"
  subdf = subset(region_excl_acc, (region_excl_acc$Sequence %in% region_peptide_list$Sequence))
  reduced_subdf = unique(subdf[,4:5])
  reduced_subdf$Accessions = getProteinAcc(reduced_subdf$Accessions)
  jointdf = merge(reduced_subdf, region_peptide_list, by = "Sequence")
  write.table(jointdf, sep = "\t", col.names = TRUE, file = paste("/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/", filename, sep = "" ), quote = FALSE, row.names = FALSE) 
  write.table(jointdf, sep = "\t", col.names = TRUE, file = paste("/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/", filename, sep = "" ), quote = FALSE, row.names = FALSE) 
  
  return(jointdf)
}

#classI
NEC_excl_acc_pep = write_region_pepacc_df(NEC_peptide_list, NEC_exclusive_I_uniqeAcc_filterd, "NEC_excl_pepacc.tsv")
T1_excl_acc_pep = write_region_pepacc_df(T1_peptide_list, T1_exclusive_I_uniqeAcc_filterd, "T1_excl_pepacc.tsv")
INF_excl_acc_pep = write_region_pepacc_df(INF_peptide_list, INF_exclusive_I_uniqeAcc_filterd, "INF_excl_pepacc.tsv")
shared_without_ben_acc_pep = write_region_pepacc_df(shared_withoutBEN_list, classI_uniqe_mappers, "shared_without_ben.tsv")
all_acc_pep = write_region_pepacc_df(all_I_list, classI_uniqe_mappers,"all_I_pepacc.tsv")
Tumor_I_acc_pep = write_region_pepacc_df(Tumor_I_freq, Tumor_I_uniqe_acc_filterd, "Tumor_I_pepacc.tsv")

#calssII
NEC_excl_acc_pep_II = write_region_pepacc_df(NEC_peptide_list_II, NEC_exclusive_II_uniqeAcc_filterd, "NEC_excl_pepacc_II.tsv")
T1_excl_acc_pep_II = write_region_pepacc_df(T1_peptide_list_II, T1_exclusive_II_uniqeAcc_filterd, "T1_excl_pepacc_II.tsv")
INF_excl_acc_pep_II = write_region_pepacc_df(INF_peptide_list_II, INF_exclusive_II_uniqeAcc_filterd, "INF_excl_pepacc_II.tsv")
shared_without_ben_acc_pep_II = write_region_pepacc_df(shared_withoutBEN_list_II_new, classII_uniqe_mappers, "shared_without_ben_II.tsv")
all_II_acc_pep = write_region_pepacc_df(all_II_list, classII_uniqe_mappers,"all_II_pepacc.tsv")

new_seq_II = c()
seq_II = Tumor_II_freq$uniqe_II
for(i in seq_II){
  #print(nchar(i))
  if(nchar(i) > 8){
    new_seq_II = append(new_seq_II, i )
  }
}

Tumor_II_freq = Tumor_II_freq %>% dplyr::filter(Tumor_II_freq$uniqe_II %in% new_seq_II)
table(Tumor_II_freq$uniqe_II)
Tumor_II_acc_pep = write_region_pepacc_df(Tumor_II_freq, Tumor_II_uniqe_acc_filterd, "Tumor_II_pepacc.tsv")


for( i in Tumor_II_acc_pep$Sequence){
  print(nchar(i))
}


check_len_II <- function(list){
  new_list = c()
  for(i in list){
    
    if(nchar(i) > 8){
      new_list = append(new_list, i )
    }
  }
  return(new_list)
}

#read in data for gioele 
NEC_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_excl_pepacc.tsv")
T1_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_excl_pepacc.tsv")
INF_excl_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_excl_pepacc.tsv")
shared_without_ben_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/shared_without_ben.tsv")
all_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/all_I_pepacc.tsv")


NEC_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ NEC_excl_pepacc_II.tsv")
T1_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ T1_excl_pepacc_II.tsv")
INF_excl_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ INF_excl_pepacc_II.tsv")
shared_without_ben_acc_pep_II = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ shared_without_ben_II.tsv")
all_II_acc_pep = read.csv(sep = "\t", "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/ all_II_pepacc.tsv")


#wirte only sequences to netMHCpan for class I 

write.table(NEC_excl_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/NEC_excl_peptides.tsv")
write.table(T1_excl_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/T1_excl_peptides.tsv")
write.table(INF_excl_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/INF_excl_peptides.tsv")
#write.table(shared_without_ben_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/shared_peptides.tsv")
write.table(Tumor_I_uniqe_acc_filterd$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCpan-4.1/tumor_regions/TumorI_peptides.tsv")
write.table(Tumor_II_acc_pep$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCIIpan-4.0/Tumor_regions//TumorII_peptides.tsv")

write.table(NEC_excl_acc_pep_II$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCIIpan-4.0/tumor_regions/NEC_excl_peptides_II.tsv")
write.table(check_len_II(T1_excl_acc_pep_II$Sequence),sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCIIpan-4.0/tumor_regions/T1_excl_peptides_II.tsv")
write.table(INF_excl_acc_pep_II$Sequence,sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/netMHCIIpan-4.0/tumor_regions/INF_excl_peptides_II.tsv")


#add tumor regions to Tumor_I and Tumor_II dfs 

nec_seq = NEC_I_uniqe_acc_filterd[3:4]
T1_seq = T1_I_uniqe_acc_filterd[3:4]
inf_seq = INF_I_uniqe_acc_filterd[3:4]
ben_seq = BEN_I_uniqe_acc_filterd[3:4]


nec_II_seq = NEC_II_uniqe_acc_filterd[3:4]
T1_II_seq = T1_II_uniqe_acc_filterd[3:4]
inf_II_seq = INF_II_uniqe_acc_filterd[3:4]
ben_II_seq = BEN_II_uniqe_acc_filterd[3:4]



write.table(nec_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_I_seq.tsv")
write.table(T1_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_I_seq.tsv")
write.table(inf_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_I_seq.tsv")
write.table(ben_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/BEN_I_seq.tsv")
write.table(nec_II_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_II_seq.tsv")
write.table(T1_II_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_II_seq.tsv")
write.table(inf_II_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_II_seq.tsv")
write.table(ben_II_seq,sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/BEN_II_seq.tsv")


test_join = unique(list(Tumor_I_acc_pep,nec_seq,T1_seq,inf_seq,ben_seq) %>% reduce(left_join, by = "Sequence"))
test_II_join = unique(list(Tumor_II_acc_pep,nec_II_seq,T1_II_seq,inf_II_seq,ben_II_seq) %>% reduce(left_join, by = "Sequence"))



library("purrr")

Tumor_I_acc_pep_o = merge(Tumor_I_acc_pep, nec_seq, by = "Sequence", all.x=TRUE)

Tumor_I_acc_pep_o  = unique(Tumor_I_acc_pep_o )

#PEPTIDE SELECTION FOR GIOELE -----------------------------------------------------------------------
#Tumor-exclusive peptides, presented on at least 2 patients --> xy-excl_acc_pep dfs enthalten das 



#Sorted based on ascending order of protein frequency in healthy tissue database --> in wie viele healthy tissues kommt es vor ?
#https://www.proteinatlas.org/about/download-->normal tissue data 
library("EnsDb.Hsapiens.v86")
library(ensembldb)

#read in normal tissue data
#Expression profiles for proteins in human tissues based on immunohistochemisty using tissue micro arrays

normal_tissue_proteom = read.csv("/Users/cschwitalla/Downloads/normal_tissue.tsv", header = TRUE, sep = "\t")
# filter for level--> all not detected will be removed

normal_tissue_proteom_filter = normal_tissue_proteom[-which(normal_tissue_proteom$Level == "Not detected"),]

normal_tissue_proteom_filter2 = normal_tissue_proteom_filter[-which(normal_tissue_proteom_filter$Reliability == "Uncertain"),]

normal_tissue_proteom_filter3 = normal_tissue_proteom_filter2 %>% dplyr::filter(normal_tissue_proteom_filter2$Gene.name %in% all_acc_pep$GeneName)
#genename tp protein acc
edb <- EnsDb.Hsapiens.v86
hasProteinData(edb)

#select(edb, keys = "P15502", keytype = "UNIPROTID", columns = "GENENAME")[1,][2]

all_acc_pep$GeneName = sapply(all_acc_pep$Accessions, function(gen){
  genname = select(edb, keys = gen, keytype = "UNIPROTID", columns = "GENENAME")[1,][2]
  
})



#normal_tissue_proteom_filter3 = normal_tissue_proteom_filter2 %>% filter(normal_tissue_proteom_filter2$)
#compute tissue abundance 
#3 = df3 %>% group_by(Patient_ID) %>% summarise(Sequence = unique(Sequence))

normal_tissue_df = normal_tissue_proteom_filter3 %>% group_by(Tissue) %>% summarise(Gene.name = unique(Gene.name))
uniq_genes = unique(normal_tissue_df$Gene.name)
normal_tissue_proteinf_df = data.frame(unique(normal_tissue_df$Gene.name))
# ratio of peptide occurences für NEC und INF 
normal_tissue_proteinf_df$Frequnecy = sapply(uniq_genes, function(gen){
  mean(ifelse(is.na(matchAll(gen, normal_tissue_df$Gene.name)), 0,
              (length(matchAll(gen, normal_tissue_df$Gene.name))/length(unique(normal_tissue_df$Tissue)))*100))
})

unique(normal_tissue_df$Tissue)
normal_tissue_df_protein_abundance = normal_tissue_proteinf_df

names(normal_tissue_df_protein_abundance)[1] = "GeneName"
names(normal_tissue_df_protein_abundance)[2] = "NormalTissueFrequency"

#combine normal tissue frequency 

all_acc_pep_2 = merge(all_acc_pep, normal_tissue_df_protein_abundance, by = "GeneName")


gNames = do.call(rbind.data.frame, all_acc_pep$GeneName)
names(gNames)[1] = "GeneName"
all_acc_pep$GeneName = gNames$GeneName

#Descending order of peptide frequency in glioblastoma tissue samples

#Descending order of protein source frequency in glioblastoma tissue samples

#protein frequency in non-glioblastoma malignant samples in descending order












#GET REGION EXCLUSIVE DATA FRAMES FOR INTEGRATION---------------------------------------------------------------
# files of tumor region exclusive peptides source protein accessions
# for integration of the different omics types

# veraltet NEU MACHEN UNBEDINGT !!======================
write.csv(Nec_excl_acc, file = "/Users/cschwitalla/Documents/Intergration//NEC_excl_acc.csv", row.names = FALSE, quote = FALSE)
write.csv(T1_excl_acc, file = "/Users/cschwitalla/Documents/Intergration//T1_excl_acc.csv", row.names = FALSE, quote = FALSE)
write.csv(INF_excl_acc, file = "/Users/cschwitalla/Documents/Intergration//INF_excl_acc.csv", row.names = FALSE, quote = FALSE)
write.csv(BEN_excl_acc, file = "/Users/cschwitalla/Documents/Intergration//BEN_excl_acc.csv", row.names = FALSE, quote = FALSE)
#=============================


test = t(as.data.frame(strsplit(files_cI, "_", fixed = TRUE)))
table(test)
