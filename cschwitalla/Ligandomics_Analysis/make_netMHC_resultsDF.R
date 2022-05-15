library("readxl")
library("stringr")
library("tidyr")
library("dplyr")
library("tibble")
library("utils")

setwd("/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/")

## Dont curse and be so hard on yourself in your code

#HIER HABE ICH BULLSHIT GEMACHT ICH MUSS DIE CLASS II PEPTIDES AUF MHC II PREDICTEN , HIER HAB ICH DAS NUR AUF WAHRSCHEIBLICH HLA CLASS I GEMACHT ALSO 
#NETMHCIIPAN NOCH MACHEN 
#bei MHCII könnte ich eventuell auf den server laden aber mal sehen , weil ich hab a keine typisierung von den patienten 


makeResultsDF <- function(file){
  raw_file = readLines(file)
  column_names = c("Pos","MHC","Peptide","Core","Of","Gp","Gl","Ip","Il","Icore","Identity","Score_EL","%Rank_EL","Score_BA","%Rank_BA","Aff(nM)","nix","BindLevel")
  index_list = c()
  for(i in c(1:length(raw_file))){
    if(!startsWith(raw_file[i],"   1" )){
      index_list = append(index_list, i)
    }
  }
  rows = as.data.frame(as.data.frame(raw_file)[-index_list,])
  split_list = c()
  for(i in rows[,1]){
    str = strsplit(i, " ", fixed = TRUE)
    split_list = append(split_list, str)
  }
  cleaned_split_list = list()
  for( i in 1:length(split_list)){
    short_list = c()
    for (j in split_list[[i]] ){
      if ( j != ""){
        short_list = append(short_list, j)
      }
      cleaned_split_list[[i]] = short_list
    }
  }
  for(i in 1:length(cleaned_split_list)){
    if(length(cleaned_split_list[[i]]) != 18){
      cleaned_split_list[[i]] = append(cleaned_split_list[[i]], c("NO", "NO")) 
    }
  }
 
  final_df = do.call(rbind.data.frame, cleaned_split_list)
  names(final_df) =column_names
  final_df = final_df[-17]
  return(final_df)
}

#makeResults_I_small---------------------
raw_file_tumor_I = readLines(Tumor_I_file)
small_tumor_I = raw_file_tumor_I[grepl("<=", raw_file_tumor_I)]

Tumor_I_df = make_small_results(small_tumor_I)

make_small_results<- function(list){
  column_names = c("Pos","MHC","Peptide","Core","Of","Gp","Gl","Ip","Il","Icore","Identity","Score_EL","%Rank_EL","Score_BA","%Rank_BA","Aff(nM)","nix","BindLevel")
  
  rows = as.data.frame(list)
  split_list = c()
  for(i in rows[,1]){
    str = strsplit(i, " ", fixed = TRUE)
    split_list = append(split_list, str)
  }
  cleaned_split_list = list()
  for( i in 1:length(split_list)){
    short_list = c()
    for (j in split_list[[i]] ){
      if ( j != ""){
        short_list = append(short_list, j)
      }
      cleaned_split_list[[i]] = short_list
    }
  }
  for(i in 1:length(cleaned_split_list)){
    if(length(cleaned_split_list[[i]]) != 18){
      cleaned_split_list[[i]] = append(cleaned_split_list[[i]], c("NO", "NO")) 
    }
  }
  
  final_df = do.call(rbind.data.frame, cleaned_split_list)
  names(final_df) =column_names
  final_df = final_df[-17]
  return(final_df)
  
}

#----------------------------

makeResultsDF_II <- function(file){
  raw_file = readLines(file)
  column_names = c("Pos","MHC","Peptide","Of","Core","Core_rel","Identity","Score_EL","%Rank_EL","Exp_Bind ","Score_BA" , "Affinity(nM)", "%Rank_BA","BindLevel")
  
  
  index_list = c()
  for(i in c(1:length(raw_file))){
    if(str_detect(raw_file[i], "Sequence") == FALSE ){
      index_list = append(index_list, i)
    }
  }
  rows = as.data.frame(as.data.frame(raw_file)[-index_list,])
  
  split_list = c()
  for(i in rows){
    str = strsplit(i, " ", fixed = TRUE)
    split_list = append(split_list, str)
  }
  cleaned_split_list = list()
  for( i in 1:length(split_list)){
    short_list = c()
    for (j in split_list[[i]] ){
      if ( j != ""){
        short_list = append(short_list, j)
      }
      cleaned_split_list[[i]] = short_list
    }
  }
  for(i in 1:length(cleaned_split_list)){
    if(length(cleaned_split_list[[i]]) != max(lengths(cleaned_split_list))){
      difference = max(lengths(cleaned_split_list)) - length(cleaned_split_list[[i]])
      vec = rep("NO", times = difference)
      cleaned_split_list[[i]] = append(cleaned_split_list[[i]], vec) 
    }
  }
  
  final_df = do.call(rbind.data.frame, cleaned_split_list)
  names(final_df) =column_names
  return(final_df)
}

raw_tumor_II = readLines(Tumor_II_file)
small_tumor_II = raw_tumor_II[grepl("<=", raw_tumor_II)]

Tumor_II_df = small_II_df(small_tumor_II)


small_II_df <-function(list){
  column_names = c("Pos","MHC","Peptide","Of","Core","Core_rel","Identity","Score_EL","%Rank_EL","Exp_Bind ","Score_BA" , "Affinity(nM)", "%Rank_BA","BindLevel")
  rows = list
  
  split_list = c()
  for(i in rows){
    str = strsplit(i, " ", fixed = TRUE)
    split_list = append(split_list, str)
  }
  cleaned_split_list = list()
  for( i in 1:length(split_list)){
    short_list = c()
    for (j in split_list[[i]] ){
      if ( j != ""){
        short_list = append(short_list, j)
      }
      cleaned_split_list[[i]] = short_list
    }
  }
  for(i in 1:length(cleaned_split_list)){
    if(length(cleaned_split_list[[i]]) != max(lengths(cleaned_split_list))){
      difference = max(lengths(cleaned_split_list)) - length(cleaned_split_list[[i]])
      vec = rep("NO", times = difference)
      cleaned_split_list[[i]] = append(cleaned_split_list[[i]], vec) 
    }
  }
  
  final_df = do.call(rbind.data.frame, cleaned_split_list)
  names(final_df) =column_names
  return(final_df)
  
}



#files class I ---------------------------
#shared_I_file = file("shared_peptides_I_netMHCpan_res.txt")
NEC_I_file = file("NEC_I_netMHCpan_res.txt")
T1_I_file = file("T1_I_netMHCpan_res.txt")
INF_I_file = file("INF_I_netMHCpan_res.txt")
Tumor_I_file = file("Tumor_I_netMHCpan_res.txt")

#files class II  -----------------------------------
#shared_II_file = file("shared_peptides_II.txt")
NEC_II_file = file("NEC_II_netMHCIIpan_res.txt")
T1_II_file = file("T1_II_netMHCIIpan_res.txt")
INF_II_file = file("INF_II_netMHCIIpan_res.txt")
Tumor_II_file = file("Tumor_II_netMHCIIpan_res.txt")

#dataframes class I -----------------------------
shared_I_df = makeResultsDF(shared_I_file)
NEC_I_df = makeResultsDF(NEC_I_file)
T1_I_df = makeResultsDF(T1_I_file)
INF_I_df = makeResultsDF(INF_I_file)
#Tumor_I_df = makeResultsDF(Tumor_I_file)

#dataframes class II ---------------------------

shared_II_df = makeResultsDF_II(shared_II_file)
NEC_II_df = makeResultsDF_II(NEC_II_file)
T1_II_df = makeResultsDF_II(T1_II_file)
INF_II_df = makeResultsDF_II(INF_II_file)
Tumor_II_Df = makeResultsDF_II(Tumor_II_file)

#Results ------------------------------------------
# ich filter alle strog and weak binders raus 
# am besten darstellen --> peptide , binding level , HLA type 
# alle peptides die nicht da reinfallen auch --> zum filtern für die high frequency exclusive und shared peptides 


# filter Class I dfs =================
#shared_I_df_filterd = shared_I_df[shared_I_df$BindLevel != "NO",]
NEC_I_df_filterd = NEC_I_df[NEC_I_df$BindLevel != "NO",]
T1_I_df_filterd = T1_I_df[T1_I_df$BindLevel != "NO",]
INF_I_df_filterd = INF_I_df[INF_I_df$BindLevel != "NO",]

#filter Class II dfs ==================
shared_II_df_filterd = shared_II_df[shared_II_df$BindLevel != "NO",]
NEC_II_df_filterd = NEC_II_df[NEC_II_df$BindLevel != "NO",]
T1_II_df_filterd = T1_II_df[T1_II_df$BindLevel != "NO",]
INF_II_df_filterd = INF_II_df[INF_II_df$BindLevel != "NO",]

#Get unique peptides that are strong / weak binders =========================
# class I 
#shared_I_binders = unique(shared_I_df_filterd$Peptide)
NEC_I_binders = unique(NEC_I_df_filterd$Peptide)
T1_I_binders = unique(T1_I_df_filterd$Peptide)
INF_I_binders = unique(INF_I_df_filterd$Peptide)
Tumor_I_binder = unique(Tumor_I_df$Peptide)
# class II 
#shared_II_binders = unique(shared_II_df_filterd$Peptide)
NEC_II_binders = unique(NEC_II_df_filterd$Peptide)
T1_II_binders = unique(T1_II_df_filterd$Peptide)
INF_II_binders = unique(INF_II_df_filterd$Peptide)
Tumor_II_binders = unique(Tumor_II_df$Peptide)

# read in high frequent peptides tsv from bevore / that were input for NetMCH pan I and II -----------------------
#class I 
#shared_I_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/ shared_without_ben.tsv", sep = "\t")
NEC_I_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/NEC_excl_pepacc.tsv", sep = "\t")
T1_I_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/T1_excl_pepacc.tsv", sep = "\t")
INF_I_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/INF_excl_pepacc.tsv", sep = "\t")
Tumor_I_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/Tumor_I_pepacc.tsv", sep = "\t")


# class II 
NEC_II_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/NEC_excl_pepacc_II.tsv", sep = "\t")
T1_II_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/T1_excl_pepacc_II.tsv", sep = "\t")
INF_II_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/INF_excl_pepacc_II.tsv", sep = "\t")
Tumor_II_high_f = read.csv(file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Highfreq_peptides/Tumor_II_pepacc.tsv", sep = "\t")

# include predicted binders

is_binder <- function(df,list){
  binder = c()
  for(i in df$Sequence){
    if(i %in% list){
      binder = append(binder, "Binder")
    }
    else{
      binder = append(binder, "No binder")
    }
  }
  df$Binder_prediction = binder
  return(df)
}

shared_I_high_f = is_binder(shared_I_high_f,shared_I_binders)
NEC_I_high_f = is_binder(NEC_I_high_f,NEC_I_binders)
T1_I_high_f = is_binder(T1_I_high_f,T1_I_binders)
INF_I_high_f = is_binder(INF_I_high_f, INF_I_binders)
Tumor_I_high_f = is_binder(Tumor_I_high_f, Tumor_I_binder)

shared_II_high_f = is_binder(shared_II_high_f, shared_II_binders)
NEC_II_high_f = is_binder(NEC_II_high_f, NEC_II_binders)
T1_II_high_f = is_binder(T1_II_high_f, T1_II_binders)
INF_II_high_f = is_binder(INF_II_high_f, INF_II_binders)
Tumor_II_high_f = is_binder(Tumor_II_high_f, Tumor_II_binders)

#include hla type for binders 
NEC_I_high_f = include_HLA(NEC_I_df_filterd,NEC_I_high_f)
T1_I_high_f = include_HLA(T1_I_df_filterd,T1_I_high_f)
INF_I_high_f = include_HLA(INF_I_df_filterd, INF_I_high_f)
Tumor_I_high_f = include_HLA(Tumor_I_df, Tumor_I_high_f)

NEC_II_high_f = include_HLA_II(NEC_II_df_filterd, NEC_II_high_f)
T1_II_high_f = include_HLA_II(T1_II_df_filterd, T1_II_high_f)
INF_II_high_f = include_HLA_II(INF_II_df_filterd, INF_II_high_f)
Tumor_II_high_f = include_HLA_II(Tumor_II_df, Tumor_II_high_f)


# EXPORT RESULTS !
write.table(NEC_I_high_f, sep = "\t", row.names = FALSE, quote = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_I_high_f.tsv")
write.table(T1_I_high_f, sep = "\t", row.names = FALSE, quote = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_I_high_f.tsv")
write.table(INF_I_high_f, sep = "\t", row.names = FALSE, quote = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_I_high_f.tsv")
write.table(Tumor_I_high_f, sep = "\t", row.names = FALSE, quote = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_I_high_f.tsv")


write.table(NEC_II_high_f, sep = "\t", row.names = FALSE, quote = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/NEC_II_high_f.tsv")
write.table(T1_II_high_f, sep = "\t", row.names = FALSE, quote = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/T1_II_high_f.tsv")
write.table(INF_II_high_f, sep = "\t", row.names = FALSE, quote = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/INF_II_high_f.tsv")

Tumor_II_high_f[102,5]
too_long_string1 = Tumor_II_high_f[98,5]
too_long_string1 = str_replace(too_long_string1, "\n", "")
Tumor_II_high_f[98,5] = too_long_string

#Tumor_II_high_f[102,5]
too_long_string2 = Tumor_II_high_f[102,5]
too_long_string2 = str_replace(too_long_string2, "\n", "")
Tumor_II_high_f[102,5] = too_long_string

write.table(Tumor_II_high_f, sep = "\t", row.names = FALSE, quote = FALSE, file = "/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_II_high_f.tsv")

TUMOR_II_test = read.table("/Users/cschwitalla/Documents/Ligandomics_analysis/Peptide_selection_gioele/Tumor_II_high_f.tsv",sep = "\t" )


#### include hla 

include_HLA <- function(binding_prediction_df, high_freq_pep_df){
  binding_prediction_df_wb = binding_prediction_df[binding_prediction_df$BindLevel == "WB",][2:3]
  binding_prediction_df_sb = binding_prediction_df[binding_prediction_df$BindLevel == "SB",][2:3]
  
  binding_prediction_df_wb = binding_prediction_df_wb %>% group_by(Peptide) %>% summarise(MHC = unique(MHC))
  binding_prediction_df_wb = binding_prediction_df_wb %>% group_by(Peptide) %>% summarise(toString(c(MHC)))
  names(binding_prediction_df_wb) = c("Sequence", "Weak_binders")
  
  binding_prediction_df_sb = binding_prediction_df_sb %>% group_by(Peptide) %>% summarise(MHC = unique(MHC))
  binding_prediction_df_sb = binding_prediction_df_sb %>% group_by(Peptide) %>% summarise(toString(c(MHC)))
  names(binding_prediction_df_sb) = c("Sequence", "Strong_binders")
  
  high_freq_pep_df = merge(high_freq_pep_df,binding_prediction_df_wb, by = "Sequence", all.x = TRUE )
  high_freq_pep_df = merge(high_freq_pep_df,binding_prediction_df_sb, by = "Sequence", all.x = TRUE )
  
  return(high_freq_pep_df)
}




include_HLA_II <- function(binding_prediction_df, high_freq_pep_df){
  binding_prediction_df_wb = binding_prediction_df[binding_prediction_df$BindLevel == "<=WB",][2:3]
  binding_prediction_df_sb = binding_prediction_df[binding_prediction_df$BindLevel == "<=SB",][2:3]
  
  binding_prediction_df_wb = binding_prediction_df_wb %>% group_by(Peptide) %>% summarise(MHC = unique(MHC))
  binding_prediction_df_wb = binding_prediction_df_wb %>% group_by(Peptide) %>% summarise(paste(c(MHC), collapse = ","))
  
  names(binding_prediction_df_wb) = c("Sequence", "Weak_binders")
  
  binding_prediction_df_sb = binding_prediction_df_sb %>% group_by(Peptide) %>% summarise(MHC = unique(MHC))
  binding_prediction_df_sb = binding_prediction_df_sb %>% group_by(Peptide) %>% summarise(toString(c(MHC)))
  names(binding_prediction_df_sb) = c("Sequence", "Strong_binders")
  
  high_freq_pep_df = merge(high_freq_pep_df,binding_prediction_df_wb, by = "Sequence", all.x = TRUE )
  high_freq_pep_df = merge(high_freq_pep_df,binding_prediction_df_sb, by = "Sequence", all.x = TRUE )
  
  return(high_freq_pep_df)
}
NEC_I_high_f_hla = include_HLA(NEC_I_df_filterd,NEC_I_high_f)

NEC_II_high_f_hla = include_HLA_II(NEC_II_df_filterd,NEC_II_high_f)



T1_I_tes_wb2 = T1_I_df_filterd[T1_I_df_filterd$BindLevel == "WB",][2:3]
T1_I_test_sb = T1_I_df_filterd[T1_I_df_filterd$BindLevel == "SB",][2:3]
T1_I_test_wb_hla = T1_I_tes_wb2 %>% group_by(Peptide) %>% summarise(MHC = unique(MHC))
T1_I_test_wb_hla2 = T1_I_tes_wb2 %>% group_by(Peptide) %>% summarise(list(MHC))

names(T1_I_test_wb_hla2)[1] = "Sequence"

T1_test_merge = merge(T1_I_high_f,T1_I_test_wb_hla2, by = "Sequence", all.x = TRUE )
