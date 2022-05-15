## So personally, I like it when even the loading of the libraries is as easy as possible.
## Therefore, you could reconsider doing the following for woud libraries:

## Suggestion Marissa
## ---------------------
## Define the libraries that you want to use
necessaryLibs <- c("maftools", "RColorBrewer", "VennDiagram", "ggplot2",
                   "devtools")
devtoolsLibs <- c("McGranahanLab/TcellExTRECT")
## Install the libraries if necessary
# invisible(lapply(necessaryLibs, BiocManager::install, update = F, ask = F))
## Load the libraries
suppressMessages(invisible(lapply(necessaryLibs, library, character.only = T)))
## Load libraries that need to be loaded in with devtools
suppressMessages(invisible(lapply(devtoolsLibs, install_github)))

## Code Carolin
## ---------------------
# library(maftools)
# library("RColorBrewer")
# library("VennDiagram")
# library("ggplot2")
#setwd("/Users/cschwitalla/Documents/WES_data") 
## I can suggest to move the setwd as one of the first lines in your script
setwd("/Users/dubbelaar/Downloads/data/") 

#library(TcellExTRECT)

## In order to go to the code here, I need to have the files 

# MAPPING OF PATIENT ID TO QBIC ID---------------------
NEC_meta = read.table(file = "./PBMC_vs_Necrotic.tsv", sep = "\t", header = FALSE)
T1_meta = read.table(file = "./PBMC_vs_T1.tsv", sep = "\t", header = FALSE)
INF_meta = read.table(file = "./PBMC_vs_Peritumoral.tsv", sep = "\t", header = FALSE)
BEN_meta = read.table(file = "./PBMC_Benign.tsv", sep = "\t", header = FALSE)

# das ist etwas anders wie ich dachte geht nicht einfach die bening samples reinzubringen 


#bening(germline) maf files
list_maf_ben = list.files(path = "./ben/germlineSarek/maf_converted/", full.names = TRUE)
list_maf_ben
ben_mafset = merge_mafs(list_maf_ben)
ben_names = list.files(path = "./ben/maf_converted/")

patientID_ben = vector()
for (i in ben_names){
  str_1 = strsplit(i, "_")
  patient_ben = str_1[[1]][2]
  print(patient_ben)
  patientID_ben = append(patientID_ben, patient_ben)
}
patientID_ben = patientID_ben[!duplicated(patientID_ben)]

clin_data_ben = data.frame(patientID_ben )
colnames(clin_data_ben) = "patient_ID"
region = rep("BEN",times=4)
clin_data_ben$Tumor_region = region
clin_data_ben$Tumor_Sample = patientID_ben

# BEN vs PBMC!!!
list_maf_benPBMC = list.files(path = "./ben/PBMC_vs_BEN/Annotation/VEP/maf_converted/", full.names = TRUE)
list_maf_benPBMC = list_maf_benPBMC[-1]
benPBMC_mafset = merge_mafs(list_maf_benPBMC)

benPBMC_names = list.files(path = "./ben/PBMC_vs_BEN/Annotation/VEP/maf_converted/")
benPBMC_names = benPBMC_names[-1]

tumorsamples_benPBMC = vector()
patientID_benPBMC = vector()
for (i in benPBMC_names){
  str_1 = strsplit(i, "_")
  tumor_benPBMC = str_1[[1]][2]
  patient_benPBMC = str_1[[1]][4]
  print(patient_benPBMC)
  tumorsamples_benPBMC = append(tumorsamples_benPBMC, tumor_benPBMC)
  patientID_benPBMC = append(patientID_benPBMC, patient_benPBMC)
}
tumorsamples_benPBMC = tumorsamples_benPBMC[!duplicated(tumorsamples_benPBMC)]
patientID_benPBMC = patientID_benPBMC[!duplicated(patientID_benPBMC)]

clin_data_benPBMC = data.frame(tumorsamples_benPBMC )
colnames(clin_data_benPBMC) = "Tumor_Sample"
clin_data_benPBMC$patient_ID = patientID_benPBMC
region = rep("BEN",times=4)
clin_data_benPBMC$Tumor_region = region


#read in necrotic maf files
list_maf_necrotic = list.files(path = "./necrotic/maf_converted/", full.names = TRUE)
list_maf_necrotic = list_maf_necrotic[-1]
necrotic_mafset = merge_mafs(list_maf_necrotic)


# clin data set Necrotic
necrotic_names = list.files(path = "./necrotic/maf_converted/")
necrotic_names = necrotic_names[-1]

tumorsamples_nec = vector()
patientID_nec = vector()
for (i in necrotic_names){
  str_1 = strsplit(i, "_")
  tumor_nec = str_1[[1]][2]
  patient_nec = str_1[[1]][4]
  print(patient_nec)
  tumorsamples_nec = append(tumorsamples_nec, tumor_nec)
  patientID_nec = append(patientID_nec, patient_nec)
}
tumorsamples_nec = tumorsamples_nec[!duplicated(tumorsamples_nec)]
patientID_nec = patientID_nec[!duplicated(patientID_nec)]

clin_data_nec = data.frame(tumorsamples_nec )
colnames(clin_data_nec) = "Tumor_Sample"
clin_data_nec$patient_ID = patientID_nec
region = rep("NEC",times=15)
clin_data_nec$Tumor_region = region



#read in t1 maf files
list_maf_t1 = list.files(path = "./t1/maf_converted/", full.names = TRUE)
list_maf_t1 = list_maf_t1[-1]
t1_mafset = merge_mafs(list_maf_t1)

#clin dataset t1

t1_names = list.files(path = "./t1/maf_converted/")
t1_names = t1_names[-1]

tumorsamples_t1 = vector()
patientID_t1 = vector()
for (i in t1_names){
  str_1 = strsplit(i, "_")
  tumor_t1 = str_1[[1]][2]
  patient_t1 = str_1[[1]][4]
  tumorsamples_t1 = append(tumorsamples_t1, tumor_t1)
  patientID_t1 = append(patientID_t1, patient_t1)
}

tumorsamples_t1 = tumorsamples_t1[!duplicated(tumorsamples_t1)]
patientID_t1 = patientID_t1[!duplicated(patientID_t1)]


clin_data_t1 = data.frame(tumorsamples_t1)
colnames(clin_data_t1) = "Tumor_Sample"
clin_data_t1$patient_ID = patientID_t1
region = rep("T1",times=15)
clin_data_t1$Tumor_region = region

#read in peritumoral maf files
list_maf_peritumoral = list.files(path = "./peritumoral/maf_converted/", full.names = TRUE)
list_maf_peritumoral = list_maf_peritumoral[-1]
peritumoral_mafset = merge_mafs(list_maf_peritumoral)

#clin data Infiltartion zone Peri( INF)


inf_names = list.files(path = "./peritumoral/maf_converted/")
inf_names = inf_names[-1]


tumorsamples_inf = vector()
patientID_inf = vector()
for (i in inf_names){
  str_1 = strsplit(i, "_")
  tumor_inf = str_1[[1]][2]
  patient_inf = str_1[[1]][4]
  tumorsamples_inf = append(tumorsamples_inf, tumor_inf)
  patientID_inf = append(patientID_inf, patient_inf)
}

tumorsamples_inf = tumorsamples_inf[!duplicated(tumorsamples_inf)]
patientID_inf = patientID_inf[!duplicated(patientID_inf)]


clin_data_inf = data.frame(tumorsamples_inf)
colnames(clin_data_inf) = "Tumor_Sample"
clin_data_inf$patient_ID = patientID_inf
region = rep("INF",times=14)
clin_data_inf$Tumor_region = region


#complete dataset 
clin_data_all = rbind(clin_data_nec,clin_data_t1, clin_data_inf,clin_data_benPBMC)
clin_data_excBEN = rbind(clin_data_nec,clin_data_t1, clin_data_inf)
names(clin_data_all)[names(clin_data_all) == "Tumor_Sample"] <- "Tumor_Sample_Barcode"
names(clin_data_excBEN)[names(clin_data_excBEN) == "Tumor_Sample"] <- "Tumor_Sample_Barcode"

list_maf_all = c(list_maf_necrotic,list_maf_peritumoral,list_maf_t1,list_maf_benPBMC)
list_all_maf_excBen = c(list_maf_necrotic,list_maf_peritumoral,list_maf_t1)
all_mafset = merge_mafs(list_maf_all, clinicalData = clin_data_all)
all_excBEN_mafset = merge_mafs(list_all_maf_excBen, clinicalData = clin_data_all)



#PATIENT MAF SETS-----------------
help("unique")
necp_IDs = unique(NEC_meta$V1)
nec_qb = unique(NEC_meta$V4)
nec_qb = nec_qb[seq(1,length(nec_qb), 2)]
NEC_meta = data.frame("Patient_ID" = necp_IDs, "QBIC_ID" = nec_qb)


t1p_IDs = unique(T1_meta$V1)
t1_qb = unique(T1_meta$V4)
t1_qb = t1_qb[seq(1,length(t1_qb), 2)]
T1_meta = data.frame("Patient_ID" = t1p_IDs, "QBIC_ID" = t1_qb)


infp_IDs = unique(INF_meta$V1)
inf_qb = unique(INF_meta$V4)
inf_qb = inf_qb[seq(1,length(inf_qb), 2)]
INF_meta = data.frame("Patient_ID" = infp_IDs, "QBIC_ID" = inf_qb)

patientID_list = NEC_meta$Patient_ID
#patientenlisten------------------------------

## This is a lot of code
## I would use an iterator 1:15 and then use that iterator to assign the values
## to a vector with the assign function 
p1 = c()
p2 = c()
p3 = c()
p4 = c()
p5 = c()
p6 = c()
p7 = c()
p8 = c()
p9 = c()
p10 = c()
p11 = c()
p12 = c()
p13 = c()
p14 = c()
p15 = c()


for (i in list_maf_necrotic){
  if(grepl(NEC_meta$QBIC_ID[1],i) == TRUE){
    p1 = append(p1,i)
  }
  if(grepl(NEC_meta$QBIC_ID[2],i) == TRUE){
    p2 = append(p2,i)
  }
  if(grepl(NEC_meta$QBIC_ID[3],i) == TRUE){
    p3 = append(p3,i)
  }
  if(grepl(NEC_meta$QBIC_ID[4],i) == TRUE){
    p4 = append(p4,i)
  }
  if(grepl(NEC_meta$QBIC_ID[5],i) == TRUE){
    p5 = append(p5,i)
  }
  if(grepl(NEC_meta$QBIC_ID[6],i) == TRUE){
    p6 = append(p6,i)
  }
  if(grepl(NEC_meta$QBIC_ID[7],i) == TRUE){
    p7 = append(p7,i)
  }
  if(grepl(NEC_meta$QBIC_ID[8],i) == TRUE){
    p8 = append(p8,i)
  }
  if(grepl(NEC_meta$QBIC_ID[9],i) == TRUE){
    p9 = append(p9,i)
  }
  if(grepl(NEC_meta$QBIC_ID[10],i) == TRUE){
    p10 = append(p10,i)
  }
  if(grepl(NEC_meta$QBIC_ID[11],i) == TRUE){
    p11 = append(p11,i)
  }
  if(grepl(NEC_meta$QBIC_ID[12],i) == TRUE){
    p12 = append(p12,i)
  }
  if(grepl(NEC_meta$QBIC_ID[13],i) == TRUE){
    p13 = append(p13,i)
  }
  if(grepl(NEC_meta$QBIC_ID[14],i) == TRUE){
    p14 = append(p14,i)
  }
  if(grepl(NEC_meta$QBIC_ID[15],i) == TRUE){
    p15 = append(p15,i)
  }
  
}

for (i in list_maf_t1){
  if(grepl(NEC_meta$QBIC_ID[1],i) == TRUE){
    p1 = append(p1,i)
  }
  if(grepl(NEC_meta$QBIC_ID[2],i) == TRUE){
    p2 = append(p2,i)
  }
  if(grepl(NEC_meta$QBIC_ID[3],i) == TRUE){
    p3 = append(p3,i)
  }
  if(grepl(NEC_meta$QBIC_ID[4],i) == TRUE){
    p4 = append(p4,i)
  }
  if(grepl(NEC_meta$QBIC_ID[5],i) == TRUE){
    p5 = append(p5,i)
  }
  if(grepl(NEC_meta$QBIC_ID[6],i) == TRUE){
    p6 = append(p6,i)
  }
  if(grepl(NEC_meta$QBIC_ID[7],i) == TRUE){
    p7 = append(p7,i)
  }
  if(grepl(NEC_meta$QBIC_ID[8],i) == TRUE){
    p8 = append(p8,i)
  }
  if(grepl(NEC_meta$QBIC_ID[9],i) == TRUE){
    p9 = append(p9,i)
  }
  if(grepl(NEC_meta$QBIC_ID[10],i) == TRUE){
    p10 = append(p10,i)
  }
  if(grepl(NEC_meta$QBIC_ID[11],i) == TRUE){
    p11 = append(p11,i)
  }
  if(grepl(NEC_meta$QBIC_ID[12],i) == TRUE){
    p12 = append(p12,i)
  }
  if(grepl(NEC_meta$QBIC_ID[13],i) == TRUE){
    p13 = append(p13,i)
  }
  if(grepl(NEC_meta$QBIC_ID[14],i) == TRUE){
    p14 = append(p14,i)
  }
  if(grepl(NEC_meta$QBIC_ID[15],i) == TRUE){
    p15 = append(p15,i)
  }
  
}
for (i in list_maf_peritumoral){
  if(grepl(NEC_meta$QBIC_ID[1],i) == TRUE){
    p1 = append(p1,i)
  }
  if(grepl(NEC_meta$QBIC_ID[2],i) == TRUE){
    p2 = append(p2,i)
  }
  if(grepl(NEC_meta$QBIC_ID[3],i) == TRUE){
    p3 = append(p3,i)
  }
  if(grepl(NEC_meta$QBIC_ID[4],i) == TRUE){
    p4 = append(p4,i)
  }
  if(grepl(NEC_meta$QBIC_ID[5],i) == TRUE){
    p5 = append(p5,i)
  }
  if(grepl(NEC_meta$QBIC_ID[6],i) == TRUE){
    p6 = append(p6,i)
  }
  if(grepl(NEC_meta$QBIC_ID[7],i) == TRUE){
    p7 = append(p7,i)
  }
  if(grepl(NEC_meta$QBIC_ID[8],i) == TRUE){
    p8 = append(p8,i)
  }
  if(grepl(NEC_meta$QBIC_ID[9],i) == TRUE){
    p9 = append(p9,i)
  }
  if(grepl(NEC_meta$QBIC_ID[10],i) == TRUE){
    p10 = append(p10,i)
  }
  if(grepl(NEC_meta$QBIC_ID[11],i) == TRUE){
    p11 = append(p11,i)
  }
  if(grepl(NEC_meta$QBIC_ID[12],i) == TRUE){
    p12 = append(p12,i)
  }
  if(grepl(NEC_meta$QBIC_ID[13],i) == TRUE){
    p13 = append(p13,i)
  }
  if(grepl(NEC_meta$QBIC_ID[14],i) == TRUE){
    p14 = append(p14,i)
  }
  if(grepl(NEC_meta$QBIC_ID[15],i) == TRUE){
    p15 = append(p15,i)
  }
  
}
## clin data for every patient 


creatClinData <- function(file_list,patientID){
  region_v = vector()
  tumor_v = vector()
  #normal_v = vector()
  for (i in file_list){
    str_1 = strsplit(i, "/")
    region = str_1[[1]][2]
    str_2 = strsplit(str_1[[1]][5],"_")
    region_v = append(region_v, region)
    tumor_v = append(tumor_v, str_2[[1]][2] )
    #normal_v = append(normal_v, str_2[[1]][4])
  }
  tumor_v = tumor_v[!duplicated(tumor_v)]
  #normal_v= normal_v[!duplicated(normal_v)]
  region_v = region_v[!duplicated(region_v)]
  
  clin_data_v = data.frame(tumor_v)
  colnames(clin_data_v) = "Tumor_Sample_Barcode"
  clin_data_v$patient_ID = rep(patientID,times=length(region_v))
  clin_data_v$Tumor_region = region_v
  
  return(clin_data_v)
}


clin_data_p1 = creatClinData(p1,patientID_list[1])
clin_data_p2 = creatClinData(p2,patientID_list[2])
clin_data_p3 = creatClinData(p3,patientID_list[3])
clin_data_p4 = creatClinData(p4,patientID_list[4])
clin_data_p5 = creatClinData(p5,patientID_list[5])
clin_data_p6 = creatClinData(p6,patientID_list[6])
clin_data_p7 = creatClinData(p7,patientID_list[7])
clin_data_p8 = creatClinData(p8,patientID_list[8])
clin_data_p9 = creatClinData(p9,patientID_list[9])
clin_data_p10 = creatClinData(p10,patientID_list[10])
clin_data_p11 = creatClinData(p11,patientID_list[11])
clin_data_p12 = creatClinData(p12,patientID_list[12])
clin_data_p13 = creatClinData(p13,patientID_list[13])
clin_data_p14 = creatClinData(p14,patientID_list[14])
clin_data_p15 = creatClinData(p15,patientID_list[15])

# merge maf files for all patients

p1_maf =  merge_mafs(p1, clinicalData = clin_data_p1)
p2_maf =  merge_mafs(p2, clinicalData = clin_data_p2)
p3_maf =  merge_mafs(p3, clinicalData = clin_data_p3)
p4_maf =  merge_mafs(p4, clinicalData = clin_data_p4)
p5_maf =  merge_mafs(p5, clinicalData = clin_data_p5)
p6_maf =  merge_mafs(p6, clinicalData = clin_data_p6)
p7_maf =  merge_mafs(p7, clinicalData = clin_data_p7)
p8_maf =  merge_mafs(p8, clinicalData = clin_data_p8)
p9_maf =  merge_mafs(p9, clinicalData = clin_data_p9)
p10_maf =  merge_mafs(p10, clinicalData = clin_data_p10)
p11_maf =  merge_mafs(p11, clinicalData = clin_data_p11)
p12_maf =  merge_mafs(p12, clinicalData = clin_data_p12)
p13_maf =  merge_mafs(p13, clinicalData = clin_data_p13)
p14_maf =  merge_mafs(p14, clinicalData = clin_data_p14)
p15_maf =  merge_mafs(p15, clinicalData = clin_data_p15)

# ONCOPLOTS----------------------------------------

#per patient ==========
par(mfrow = c(2,1))
oncoplot(maf = p1_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region", top = 40)
oncoplot(maf = p2_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p3_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p4_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p5_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p6_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p7_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p8_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p9_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p10_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p11_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p12_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p13_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p14_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")
oncoplot(maf = p15_maf, draw_titv = TRUE, clinicalFeatures = "Tumor_region")

#FILE EXPORT FOR INTEGRATION ---------------------------------
Nec_geneSumm = getGeneSummary(necrotic_mafset)
Nec_geneSumm_orderd = Nec_geneSumm[order(Nec_geneSumm$total, decreasing = TRUE),]

T1_geneSumm = getGeneSummary(t1_mafset)
T1_geneSumm_orderd = T1_geneSumm[order(T1_geneSumm$total, decreasing = TRUE),]

INF_geneSumm = getGeneSummary(peritumoral_mafset)
INF_geneSumm_orderd = INF_geneSumm[order(INF_geneSumm$total, decreasing = TRUE),]

BEN_geneSumm = getGeneSummary(benPBMC_mafset)
BEN_geneSumm_orderd = BEN_geneSumm[order(BEN_geneSumm$total, decreasing = TRUE),]

write.csv(Nec_geneSumm_orderd$Hugo_Symbol, file = "/Users/cschwitalla/Documents/Intergration//Nec_mutatedgenes.csv", row.names = FALSE, quote = FALSE)
write.csv(T1_geneSumm_orderd$Hugo_Symbol, file = "/Users/cschwitalla/Documents/Intergration//T1_mutatedgenes.csv", row.names = FALSE, quote = FALSE)
write.csv(INF_geneSumm_orderd$Hugo_Symbol, file = "/Users/cschwitalla/Documents/Intergration//INF_mutatedgenes.csv", row.names = FALSE, quote = FALSE)
write.csv(BEN_geneSumm_orderd$Hugo_Symbol, file = "/Users/cschwitalla/Documents/Intergration//BEN_mutatedgenes.csv", row.names = FALSE, quote = FALSE)

write.table(clin_data_nec, file = "/Users/cschwitalla/Documents/WES_data/necrotic/clin_data_nec.txt", row.names = FALSE, quote = FALSE)
write.table(clin_data_t1, file = "/Users/cschwitalla/Documents/WES_data/necrotic/clin_data_t1.txt", row.names = FALSE, quote = FALSE)
write.table(clin_data_inf, file = "/Users/cschwitalla/Documents/WES_data/necrotic/clin_data_INF.txt", row.names = FALSE, quote = FALSE)
write.table(clin_data_benPBMC, file = "/Users/cschwitalla/Documents/WES_data/necrotic/clin_data_benPBMC.txt", row.names = FALSE, quote = FALSE)

# PLOTS -------------------------------
#ONCOPLOTS ============================

#Variant allel Frequency
# ----- from Rike
# Calculate Varaint allele frequency
compute_vaf <- function(maf_obj){
  genes_vaf = subsetMaf(maf = maf_obj, query = "t_depth > 0", fields = c("t_depth","t_alt_count"), mafObj = FALSE)
  genes_vaf$VAF <- genes_vaf$t_alt_count/genes_vaf$t_depth
  genes_vaf <- genes_vaf[, mean(VAF), Hugo_Symbol]
  colnames(genes_vaf)[2] = "VAF"
  genes_vaf <- genes_vaf[which(genes_vaf$VAF > 0.02)]
  return(genes_vaf)
}






#necrotic
vaf_nec = compute_vaf(necrotic_mafset)
oncoplot(maf = necrotic_mafset, leftBarData = vaf_nec)


#t1
vaf_t1 = compute_vaf(t1_mafset)
oncoplot(maf = t1_mafset, leftBarData = vaf_t1)

#peritumoral
vaf_inf = compute_vaf(peritumoral_mafset)
oncoplot(maf = peritumoral_mafset, leftBarData = vaf_inf)

#ben 
oncoplot(maf = benPBMC_mafset, draw_titv = TRUE)



# combine all patient mafs+clindata
patien_clin_data = rbind(clin_data_p1,clin_data_p2, clin_data_p3,clin_data_p4,
                         clin_data_p5,clin_data_p6,clin_data_p7,clin_data_p8,
                         clin_data_p9,clin_data_p10,clin_data_p11, clin_data_p12,
                         clin_data_p13,clin_data_p14,clin_data_p15)
all_patient_list = c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15)
all_patient_mafset = merge_mafs(all_patient_list, clinicalData = patien_clin_data)

samples = c(patien_clin_data$Tumor_Sample_Barcode)
samples
oncoplot(maf = all_patient_mafset, draw_titv = TRUE,  clinicalFeatures = c("Tumor_region","patient_ID"), sortByAnnotation = TRUE, sampleOrder = samples,removeNonMutated = FALSE )

all_patient_mafset
#all n t1 p
library("ColourBrewer")
# hier kann man noch nach regionen clustern sozusagen 
brewer.pal(n = 9,"Set1")

vc_cols = RColorBrewer::brewer.pal(n = 6, name = "Set1")
#vc_cols = c("#3288bd", "#99d594", "#e6f598", "#fee08b", "#fc8d59", "#d53e4f")

#vc_cols= c("#FF7F00" ,"#1F78B4" , "#FB9A99", "#252525","#33A02C", "#6A3D9A" )
#vc_cols = c("#4575b4", "#91bfdb", "#e0f3f8", "#f33090", "#fc8d59", "#d73027")
vc_cols=c("#E41A1C", "#377EB8" ,"#4DAF4A", "#252525", "#FF7F00" ,"#FFFF33")
names(vc_cols) = c(  'Frame_Shift_Del','Missense_Mutation',
                     'Nonsense_Mutation',
                     'Multi_Hit',
                     'Frame_Shift_Ins',
                     'In_Frame_Del'
)
#ann_cols = 
library("svglite")


oncoplot(maf = all_mafset, draw_titv = TRUE, clinicalFeatures =  "Tumor_region", sortByAnnotation = TRUE, colors = vc_cols, removeNonMutated = FALSE)
oncoplot(maf = all_excBEN_mafset, clinicalFeatures = "Tumor_region", sortByAnnotation = TRUE, colors=vc_cols,removeNonMutated = FALSE )

oncoplot(maf = all_excBEN_mafset, clinicalFeatures = "Tumor_region", sortByAnnotation = TRUE, colors=vc_cols,removeNonMutated = FALSE )

ggsave("onco_marcel.pdf", width = 10, height = 5)


all_mafset@data$Tumor_Sample_Barcode

write.mafSummary(maf = all_mafset, basename = '/Users/cschwitalla/Documents/WES_analysis/all_mafset')

write(all_mafset, "/Users/cschwitalla/Documents/WES_data/all_mafset.maf")



library(dplyr)


# -----

oncoplot(maf = all_mafset, clinicalFeatures = "Tumor_region", sortByAnnotation = TRUE, colors=vc_cols,removeNonMutated = FALSE, leftBarData = compute_vaf(all_mafset))

vaf_all = compute_vaf(all_mafset)
# LNP1 kein vaf 

genes_vaf = subsetMaf(maf = all_mafset, query = "t_depth > 0", fields = c("t_depth","t_alt_count"), mafObj = FALSE)

LNPsubsetfor_vaf = subset(genes_vaf, Hugo_Symbol == "LNP1")
#exclude genes without VAF 
#get list of most muta genes --> filter out genes without VAF 
all_genes_summary = getGeneSummary(all_mafset)

all_genes_summary_filtered = all_genes_summary %>% filter(Hugo_Symbol %in% vaf_all$Hugo_Symbol)
top20genes = all_genes_summary_filtered$Hugo_Symbol[1:20]

sample_order = c(clin_data_all$Tumor_Sample_Barcode)
sample_order
pdf(file = "/Users/cschwitalla/Documents/WES_analysis/onopräsi_2.pdf", width = 10, height = 5 )
oncoplot(maf = all_mafset, 
         genes = top20genes, 
         #clinicalFeatures = c("patient_ID","Tumor_region"), 
         sortByAnnotation = TRUE, colors=vc_cols,
         removeNonMutated = FALSE, 
         leftBarData = compute_vaf(all_excBEN_mafset), 
         gene_mar = 6, 
         sampleOrder = sample_order)
dev.off()


oncoplot(maf = all_mafset, 
         genes = top20genes, 
         clinicalFeatures = c("patient_ID","Tumor_region"), 
         sortByAnnotation = TRUE, colors=vc_cols,
         removeNonMutated = FALSE, 
         leftBarData = compute_vaf(all_excBEN_mafset), 
         gene_mar = 6, 
         sampleOrder = sample_order)

#Capture kit size ----------------------------------

bed_file = read.csv(file = "hglft_genome_298d9_1b60f0.sort.merged.bed", sep = "\t", header = FALSE, labelPos = "all")

bed_file$Difference = bed_file$V3 -bed_file$V2
kit_size = sum(bed_file$Difference)
kit_size
#Lolipop plots -----------------------------

#make my own lollipop plot =============================
nec_patient_mapping = read.csv("/Users/cschwitalla/Documents/WES_data/NEC_PATIENT_MAPPING.tsv", header = FALSE, sep = "\t" )
t1_patient_mapping = read.csv("/Users/cschwitalla/Documents/WES_data/T1_PATIENT_MAPPING.tsv", header = FALSE, sep = "\t" )
inf_patient_mapping = read.csv("/Users/cschwitalla/Documents/WES_data/INF_PATIENT_MAPPING.tsv", header = FALSE, sep = "\t" )
ben_patient_mapping = read.csv("/Users/cschwitalla/Documents/WES_data/BEN_PATIENT_MAPPING.tsv", header = FALSE, sep = "\t" )

nec_patient_mapping = unique(nec_patient_mapping)
t1_patient_mapping = unique(t1_patient_mapping)
inf_patient_mapping = unique(inf_patient_mapping)
ben_patient_mapping = unique(ben_patient_mapping)

names(nec_patient_mapping) = c("Patient", "Tumor_Sample_Barcode")
names(t1_patient_mapping) = c("Patient", "Tumor_Sample_Barcode")
names(inf_patient_mapping) = c("Patient", "Tumor_Sample_Barcode")
names(ben_patient_mapping) = c("Patient", "Tumor_Sample_Barcode")


nec_maf_subset = necrotic_mafset@data
t1_maf_subset = t1_mafset@data
inf_maf_subset = peritumoral_mafset@data
ben_maf_subset = benPBMC_mafset@data

nec_maf_subset = select(nec_maf_subset,"Hugo_Symbol", "Tumor_Sample_Barcode","Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2","Allele", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "PolyPhen")
t1_maf_subset = select(t1_maf_subset,"Hugo_Symbol", "Tumor_Sample_Barcode","Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2","Allele", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "PolyPhen")
inf_maf_subset = select(inf_maf_subset,"Hugo_Symbol", "Tumor_Sample_Barcode","Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2","Allele", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "PolyPhen")
ben_maf_subset = select(ben_maf_subset,"Hugo_Symbol", "Tumor_Sample_Barcode","Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2","Allele", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "PolyPhen")

nec_maf_subset = merge(nec_maf_subset, nec_patient_mapping,by = "Tumor_Sample_Barcode", all.x = TRUE)
t1_maf_subset = merge(t1_maf_subset, nec_patient_mapping,by = "Tumor_Sample_Barcode", all.x = TRUE)
inf_maf_subset = merge(inf_maf_subset, nec_patient_mapping,by = "Tumor_Sample_Barcode", all.x = TRUE)
ben_maf_subset = merge(ben_maf_subset, nec_patient_mapping,by = "Tumor_Sample_Barcode", all.x = TRUE)

# make lolli plot 

# to do 
#adjust  axis labels 
#adjust text sizes and every size etc 
#add labels for domains etc 
Lollipop_plot <- function(genename, maffile_subset){
  
  gene_mafsubset = maffile_subset[which(maffile_subset$Hugo_Symbol == genename),]
  gene_protein_pos = gene_mafsubset$Protein_position
  total_len = c()
  pos = c()
  for(i in gene_protein_pos){
    str = strsplit(i, "/", fixed = TRUE)
    str2 = strsplit(str[[1]][1], "-", fixed = TRUE)
    pos = append(pos,str2[[1]][1])
    total_len = append(total_len, str[[1]][2])
  }
  gene_mafsubset$AA_mut = pos
  gene_mafsubset$Total_protein = total_len
  mut_occurance = ddply(gene_mafsubset,.(AA_mut,Amino_acids,Variant_Classification, Total_protein), nrow)
  
  geneID = genename
  #Protein domain source.
  gff = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
  gff = readRDS(file = gff)
  data.table::setDT(x = gff)
  prot = gff[HGNC %in% geneID]
  txs = unique(prot$refseq.ID)
  if(length(txs) > 1){
    prot = prot[which(prot$aa.length == max(prot$aa.length)),]
  }

  prot = prot[,domain_lenght := End - Start][order(domain_lenght, decreasing = TRUE)][,domain_lenght := NULL]
  domains = unique(prot[,Label])
  domain_cols = brewer.pal(n = 8, name = "Pastel1")
  
  if(length(domains) > length(domain_cols)){
    domain_cols = sample(colours(), size = length(domains), replace = FALSE)
  }
  domain_cols = domain_cols[1:length(domains)]
  names(domain_cols) = domains
  prot[, domainCol := domain_cols[prot[, Label]]]
  
  # Plot
  ggplot(mut_occurance, aes(x=as.numeric(AA_mut), y=V1, label = paste(Amino_acids, AA_mut, sep = ": "))) +
    geom_segment( aes(x=as.numeric(AA_mut), xend=as.numeric(AA_mut), y=0, yend=V1), color ="grey") +
    geom_point(aes(colour = factor(Variant_Classification)), size=4) + 
    geom_text(nudge_x = 0, nudge_y = 0.4)+
    scale_x_continuous(limits = c(1,as.numeric(total)))+
    annotate("rect", xmin = 1, xmax = as.numeric(total), ymin = -0.25, ymax = 0.25, alpha = 0.5, fill = "grey")+
    annotate("rect", xmin = prot[,Start], xmax = prot[,End], ymin = -0.4, ymax = 0.4, fill=prot[,domainCol], alpha = 0.75, color = "lightsteelblue4")+
    annotate("text", x = ((prot[,End]-prot[,Start])/2) +prot[,Start], y = 0, label = prot[,Label])+
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  
}

Lollipop_plot("OR8U1", nec_maf_subset[which(nec_maf_subset$Hugo_Symbol == "OR8U1"),])



#test stuff -------------------------


OR8U1_nec = nec_maf_subset[which(nec_maf_subset$Hugo_Symbol == "OR8U1"),]


TP53_nec = nec_maf_subset[which(nec_maf_subset$Hugo_Symbol == "TP53"),]

OR8U1_nec_pos = OR8U1_nec$Protein_position
total_len = c()
pos = c()
for(i in OR8U1_nec_pos){
  str = strsplit(i, "/", fixed = TRUE)
  str2 = strsplit(str[[1]][1], "-", fixed = TRUE)
  pos = append(pos,str2[[1]][1])
  total_len = append(total_len, str[[1]][2])
}

OR8U1_nec$AA_mut = pos
OR8U1_nec$Total_protein = total_len


#TP53_nec$num_patients = c(1,1,1,1,1,1)# das muss man besser machen
library(plyr)
test_df = ddply(OR8U1_nec,.(AA_mut,Amino_acids,Variant_Classification, Total_protein), nrow)
TP53_nec_2 = merge(TP53_nec, test_df, by = c( "AA_mut","Amino_acids"), all.x = TRUE)

tabel_summary = table(apply(select(TP53_nec, "AA_mut", "Patient","Amino_acids"),1,function(x) paste(sort(x), collapse = "_")))
tabel_summary

data <- data.frame(x=seq(1,30), y=abs(rnorm(30)))
vc_cols=c("#E41A1C", "#377EB8" ,"#4DAF4A", "#252525", "#FF7F00" ,"#FFFF33")
names(vc_cols) = c(  'Frame_Shift_Del','Missense_Mutation',
                     'Nonsense_Mutation',
                     'Multi_Hit',
                     'Frame_Shift_Ins',
                     'In_Frame_Del'
)
total = test_df$Total_protein[1]
# Plot
ggplot(test_df, aes(x=as.numeric(AA_mut), y=V1, label = paste(Amino_acids, AA_mut, sep = ": "))) +
  geom_segment( aes(x=as.numeric(AA_mut), xend=as.numeric(AA_mut), y=0, yend=V1), color ="grey") +
  geom_point(aes(colour = factor(Variant_Classification)), size=4) + 
  geom_text(nudge_x = 0, nudge_y = 0.4)+
  scale_x_continuous(limits = c(1,as.numeric(total)))+
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  )


plotProtein(gene = "OR8U1")


#maftools code 
geneID ="OR8U1"
#Protein domain source.
gff = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
gff = readRDS(file = gff)
data.table::setDT(x = gff)
prot = gff[HGNC %in% geneID]
txs = unique(prot$refseq.ID)
if(length(txs) > 1){
  prot = prot[which(prot$aa.length == max(prot$aa.length)),]
}

#len = as.numeric(max(prot$aa.length, na.rm = TRUE))
prot = prot[,domain_lenght := End - Start][order(domain_lenght, decreasing = TRUE)][,domain_lenght := NULL]

#xlimPos = pretty(0:max(prot$aa.length, na.rm = TRUE))
#xlimPos[length(xlimPos)] = max(prot$aa.length)
domains = unique(prot[,Label])
domain_cols = brewer.pal(n = 8, name = "Pastel1")

if(length(domains) > length(domain_cols)){
  domain_cols = sample(colours(), size = length(domains), replace = FALSE)
}
domain_cols = domain_cols[1:length(domains)]
#domain_cols = grDevices::adjustcolor(col = domain_cols)
names(domain_cols) = domains
prot[, domainCol := domain_cols[prot[, Label]]]

ggplot(test_df, aes(x=as.numeric(AA_mut), y=V1, label = paste(Amino_acids, AA_mut, sep = ": "))) +
  geom_segment( aes(x=as.numeric(AA_mut), xend=as.numeric(AA_mut), y=0, yend=V1), color ="grey") +
  geom_point(aes(colour = factor(Variant_Classification)), size=4) + 
  geom_text(nudge_x = 0, nudge_y = 0.4)+
  scale_x_continuous(limits = c(1,as.numeric(total)))+
  annotate("rect", xmin = 1, xmax = as.numeric(total), ymin = -0.25, ymax = 0.25, alpha = 0.75, fill = "grey")+
  annotate("rect", xmin = prot[,Start], xmax = prot[,End], ymin = -0.4, ymax = 0.4, fill=prot[,domainCol])+
  annotate("text", x = prot[,Start]+25, y = 0, label = prot[,Label])+
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  )


((prot[,End]-prot[,Start])/2) +prot[,Start]
#geom_rect(data = test_df,mapping=aes(xmin = 0, xmax = len, ymin = 0, ymax = 0.5), color = "grey" ,alpha=0.5, inherit.aes = FALSE)+
par(mar = c(0, 1, 2, 1))
plot(0, 0, pch = NA, ylim = c(0, 0.85), xlim = c(0, len), axes = FALSE, xlab = NA, ylab = NA)

rect(xleft = 0, ybottom = 0.5, xright = len, ytop = 0.7, col = "#95a5a6", border = "grey")

prot[, domainCol := domain_cols[prot[, Label]]]

rect(xleft = prot[,Start], ybottom = 0.4, xright = prot[,End], ytop = 0.8, col = prot[,domainCol], border = "grey")


text(x = xlimPos, y = 0.3, labels = xlimPos, col = "#34495e")
#rect(xleft = xlimPos, ybottom = 0.22, xright = xlimPos, ytop = 0.24, col = "gray90")

title(main = "OR8U1", adj = 0, font.main = 3, cex.main = 1, line = 0.8)
title(main = unique(prot[,refseq.ID]), adj = 0, font.main = 1, line = -0.5, cex.main = 0.8)


prot = prot[!duplicated(Label)]
prot$pos = rowMeans(x = prot[,.(Start, End)], na.rm = FALSE)
text(y = 0.6, x = prot$pos, labels = prot$Label, font = 3, cex = 1)

legend(x = "bottomleft", legend = names(domain_cols), col = domain_cols, pch = 15, ncol = legendNcol, bty = "n", xpd = TRUE, cex = legendTxtSize)

#end of test stuff ----------------------------------




# end of own lollipop plot ==============================



lollipopPlot(necrotic_mafset, gene = "TP53", AACol = "Amino_acids", proteinID = "NP_000537", refSeqID = "NM_000546")

necrotic_mafset@data$Amino_acids[108]
nec_genes = as.data.frame(necrotic_mafset@data$Hugo_Symbol)


getGeneSummary()
lollipopPlot(all_excBEN_mafset, gene = "OR8U1", showMutationRate = TRUE, AACol = "HGVSp")
plotProtein(gene = "OR8U1")
#Summary
ben_mafset
plotmafSummary(benPBMC_mafset)
all_mafset
necrotic_mafset

all = getGeneSummary(all_mafset)

#all_BP$Hugo_Symbol == "CHIT1"
CHIT1 = all[3,]
CHIT1BP = all_BP[887,]
#plot transition and transversion
all.titv = titv(maf = all_mafset, plot = FALSE, useSyn = TRUE)
plotTiTv(res = all.titv)

n.titv = titv(maf = necrotic_mafset, plot = FALSE, useSyn = TRUE)
plotTiTv(res = n.titv)

t1.titv = titv(maf = t1_mafset, plot = FALSE, useSyn = TRUE)
plotTiTv(res = t1.titv)

p.titv = titv(maf = peritumoral_mafset, plot = FALSE, useSyn = TRUE)
plotTiTv(res = p.titv)

#comapre mutation load against TCGA
#TO DO : mutation load schöner machen .... !!
par(mfrow=c(1,1))
gbm.mutload_ben = tcgaCompare(maf = benPBMC_mafset, cohortName = "BEN", logscale = TRUE, capture_size = 36.58)
gbm.mutload_nt1p = tcgaCompare(maf = all_excBEN_mafset, cohortName = 'GB', logscale = TRUE, capture_size =  36.58)
gbm.mutload_necrotic = tcgaCompare(maf = necrotic_mafset, cohortName = 'necrotic zone', logscale = TRUE, capture_size =  36.58)
gbm.mutload_t1 = tcgaCompare(maf = t1_mafset, cohortName = 'T1 zone', logscale = TRUE, capture_size =  36.58)
gbm.mutload_peritumoral = tcgaCompare(maf = peritumoral_mafset, cohortName = 'infiltartion zone', logscale = TRUE, capture_size = 36.58)

help("tcgaCompare")
#somatic interactions
#nachlesen was genau mir das sagen soll
somaticInteractions(maf = necrotic_mafset, top = 25, pvalue = c(0.05, 0.1))

somaticInteractions(maf = t1_mafset, top = 25, pvalue = c(0.05, 0.1))

somaticInteractions(maf = peritumoral_mafset, top = 25, pvalue = c(0.05, 0.1))

somaticInteractions(maf = all_mafset, top = 25, pvalue = c(0.05, 0.1))


#cancer driven genes
#vielleicht zu hohe R version 
necrotic.driver = oncodrive(maf = necrotic_mafset, AACol = NULL, pvalMethod = 'zscore')
plotOncodrive(res = necrotic.driver, fdrCutOff = 0.1, useFraction = T, labelSize = 0.6)

t1.driver = oncodrive(maf = t1_mafset, AACol = NULL, pvalMethod = 'zscore')
plotOncodrive(res = t1.driver, fdrCutOff = 0.1, useFraction = T, labelSize = 0.6)

peritumoral.driver = oncodrive(maf = peritumoral_mafset, AACol = NULL, pvalMethod = 'zscore')
plotOncodrive(res = peritumoral.driver, fdrCutOff = 0.1, useFraction = T, labelSize = 0.6)


#comparing 2 cohorts

n_vs_t1 <- mafCompare(m1 = necrotic_mafset, m2 = t1_mafset, m1Name = 'Necrotic', m2Name = 'T1')
n_vs_p <- mafCompare(m1 = necrotic_mafset, m2 = peritumoral_mafset, m1Name = 'Necrotic', m2Name = 'Peritumoral')
t1_vs_p <- mafCompare(m1 = t1_mafset, m2 = peritumoral_mafset, m1Name = 'T1', m2Name = 'Peritumoral')

p_vs_ben <- mafCompare(m1 = peritumoral_mafset, m2 = benPBMC_mafset, m1Name = "INF", m2Name = "BEn") 
n_vs_ben <- mafCompare(m1 = necrotic_mafset, m2 = benPBMC_mafset, m1Name = "NEC", m2Name = "BEn") 
t1_vs_ben <- mafCompare(m1 = t1_mafset, m2 = benPBMC_mafset, m1Name = "T1", m2Name = "BEn") 


forestPlot(mafCompareRes = n_vs_ben, pVal = 0.6)
forestPlot(mafCompareRes = t1_vs_ben, pVal = 0.3)
forestPlot(mafCompareRes = p_vs_ben, pVal = 0.6)
forestPlot(mafCompareRes = n_vs_t1, pVal = 0.1)
forestPlot(mafCompareRes = n_vs_p, pVal = 0.1)
forestPlot(mafCompareRes = t1_vs_p, pVal = 0.1)

help("forestPlot")

#oncogenic pathways
OncogenicPathways(necrotic_mafset)
OncogenicPathways(t1_mafset)
OncogenicPathways(peritumoral_mafset)
OncogenicPathways(all_mafset)
OncogenicPathways(benPBMC_mafset)


#TO DO: plot maf summary

plotmafSummary(maf = necrotic_mafset, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plotmafSummary(maf = t1_mafset, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plotmafSummary(maf = peritumoral_mafset, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plotmafSummary(maf = all_mafset, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#co oncoplot 

coOncoplot(m1 = benPBMC_mafset, m2 = all_mafset, m1Name = 'BEN', m2Name = 'all', removeNonMutated = TRUE)

#co barplot 
#sehen nice aus , mehr davon 

coBarplot(m1 = necrotic_mafset, m2 = t1_mafset, m1Name = "GBM necrotic Strelka", m2Name = "GBM t1 Strelka")
coBarplot(m1 = StrelkaBP_necrotic, m2 = StrelkaBP_t1, m1Name = "GBM necrotic Strelka", m2Name = "GBM t1 Strelka")


coBarplot(m1 = necrotic_mafset, m2 = peritumoral_mafset, m1Name = "GBM necrotic Strelka", m2Name = "GBM peritumoral Strelka")
coBarplot(m1 = StrelkaBP_necrotic, m2 = StrelkaBP_peritumoral, m1Name = "GBM necrotic Strelka", m2Name = "GBM peritumoral Strelka")

coBarplot(m1 = peritumoral_mafset, m2 = t1_mafset, m1Name = "GBM peritumoral Strelka", m2Name = "GBM t1 Strelka")
coBarplot(m1 = StrelkaBP_peritumoral, m2 = StrelkaBP_t1, m1Name = "GBM peritumoral Strelka", m2Name = "GBM t1 Strelka")

coBarplot(m1 = benPBMC_mafset, m2 = necrotic_mafset, m1Name = "BEN", m2Name = "NEC")

# drug gene interaction einfach mal so zum spa?
dgi = drugInteractions(maf = necrotic_mafset, fontSize = 0.75)




# ----- from Rike
# Calculate Varaint allele frequency
compute_vaf <- function(maf_obj){
  genes_vaf = subsetMaf(maf = maf_obj, query = "t_depth > 0", fields = c("t_depth","t_alt_count"), mafObj = FALSE)
  genes_vaf$VAF <- genes_vaf$t_alt_count/genes_vaf$t_depth
  genes_vaf <- genes_vaf[, mean(VAF), Hugo_Symbol]
  colnames(genes_vaf)[2] = "VAF"
  genes_vaf <- genes_vaf[which(genes_vaf$VAF > 0.05)]
  return(genes_vaf)
}
#VENN diagramms ----------------------------------------------

nec_df = getGeneSummary(necrotic_mafset)
t1_df = getGeneSummary(t1_mafset)
inf_df = getGeneSummary(peritumoral_mafset)
ben_df = getGeneSummary(benPBMC_mafset)
nec_genes <- nec_df$Hugo_Symbol
t1_genes <- t1_df$Hugo_Symbol
inf_genes <- inf_df$Hugo_Symbol
ben_genes <- ben_df$Hugo_Symbol

set_list = list( NEC = nec_genes, T1 = t1_genes, INF = inf_genes, BEN = ben_genes)

myCol = brewer.pal(4, "Pastel2")

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

drawVenn(set_list, "mutated genes")

#tumor heterogenity 
library("mclust")
help("inferHeterogeneity")

#GBM_necrotic_strelka.het = inferHeterogeneity(maf = necrotic_mafset, , vafCol = 'i_TumorVAF_WU')


#alles deinstallieren
#library(survtype)
#help(maf2matrix)

#mutmatrix = maf2matrix(Strelka_snvs_nt1p)

#install.packages("devtools")
#library(devtools) 
#install_github("hanjunwei-lab/SMDIC")
#library(SMDIC)

#mutmatrix = maf2matrix(Strelka_snvs_necrotic)

#mutmatrix_necrotic = mutCountMatrix(Strelka_snvs_necrotic)

#options(repos=c(getOption("repos"), "http://ccb.nki.nl/software/discover/repos/r"))
#install.packages("discover")
#library(discover)


#result.mutex <- pairwise.discover.test(mutmatrix_necrotic[mutmatrix_necrotic, ])
