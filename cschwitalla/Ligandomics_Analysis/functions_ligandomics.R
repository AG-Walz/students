################################################################################
###                                Purpose                                   ###
################################################################################
## Title: Functions for Immunopeptidomics Analysis
## Author: Carolin Schwitalla
##
## Description: This script contains all functions that are used for the
##              immunopeptidomics analysis for my master thesis
################################################################################
###                                                                          ###
################################################################################

## TITLE: multi waterfall data frames ------------------------------------------
## DESCRIPTION: makes dataframes for each tumor reagion with the frequency of
##              region specific peptides or source proteins as well as of the
##              combined tumor
##
## PARAMETERS: - df1(2,3): tumor region specific data frames with peptide +
##                         source protein information [data frame]
##             - with_seq: if true, frequencies will be computed for peptides
##                         if false, frequencies will be computed for source
##                         proteins [boolean]
##
## OUTPUT: List of dataframes with region specific peptide/protein frequencies
##         as well as the combinde tumor region frequency data frame
##
make_multi_waterfall_df <- function(df1, df2, df3, with_seq) {
  if (with_seq == TRUE) {
    uniq_seq <- unique(c(df1$Sequence, df2$Sequence, df3$Sequence))

    df1 <- df1 %>%
      group_by(Patient_ID) %>%
      summarise(Sequence = unique(Sequence))
    df2 <- df2 %>%
      group_by(Patient_ID) %>%
      summarise(Sequence = unique(Sequence))
    df3 <- df3 %>%
      group_by(Patient_ID) %>%
      summarise(Sequence = unique(Sequence))
    # neues df
    waterfall_df <- data.frame(uniq_seq)
    # ratio of peptide occurences für NEC und INF
    waterfall_df$set1 <- sapply(uniq_seq, function(seq) {
      mean(ifelse(is.na(tuple::matchAll(seq, df1$Sequence)), 0,
        (length(tuple::matchAll(seq, df1$Sequence)) / length(unique(df1$Patient_ID))) * 100
      ))
    })

    waterfall_df$set2 <- sapply(uniq_seq, function(seq) {
      mean(ifelse(is.na(tuple::matchAll(seq, df2$Sequence)), 0,
        (length(tuple::matchAll(seq, df2$Sequence)) / length(unique(df2$Patient_ID))) * 100
      ))
    })

    waterfall_df$set3 <- sapply(uniq_seq, function(seq) {
      mean(ifelse(is.na(tuple::matchAll(seq, df3$Sequence)), 0,
        (length(tuple::matchAll(seq, df3$Sequence)) / length(unique(df3$Patient_ID))) * 100
      ))
    })



    waterfall_df$Total <- sapply(uniq_seq, function(seq) {
      (length(tuple::matchAll(seq, c(df1$Sequence, df2$Sequence, df3$Sequence))) /
        length(c(unique(df1$Patient_ID), unique(df2$Patient_ID), unique(df3$Patient_ID)))) * 100
    })
  } else {
    uniq_acc <- unique(c(df1$Accessions, df2$Accessions, df3$Accessions))

    # unique sequence per sample weil wir sehen wollen in wie vielen samples kommt das pep vor und nicht wie oft im sample
    df1 <- df1 %>%
      group_by(Patient_ID) %>%
      summarise(Accessions = unique(Accessions))
    df2 <- df2 %>%
      group_by(Patient_ID) %>%
      summarise(Accessions = unique(Accessions))
    df3 <- df3 %>%
      group_by(Patient_ID) %>%
      summarise(Accessions = unique(Accessions))
    # neues df
    waterfall_df <- data.frame(uniq_acc)
    # ratio of peptide occurences für NEC und INF
    waterfall_df$set1 <- sapply(uniq_acc, function(acc) {
      mean(ifelse(is.na(tuple::matchAll(acc, df1$Accessions)), 0,
        (length(tuple::matchAll(acc, df1$Accessions)) / length(unique(df1$Patient_ID))) * 100
      ))
    })

    waterfall_df$set2 <- sapply(uniq_acc, function(acc) {
      mean(ifelse(is.na(tuple::matchAll(acc, df2$Accessions)), 0,
        (length(tuple::matchAll(acc, df2$Accessions)) / length(unique(df2$Patient_ID))) * 100
      ))
    })

    waterfall_df$set3 <- sapply(uniq_acc, function(acc) {
      mean(ifelse(is.na(tuple::matchAll(acc, df3$Accessions)), 0,
        (length(tuple::matchAll(acc, df3$Accessions)) / length(unique(df3$Patient_ID))) * 100
      ))
    })



    waterfall_df$Total <- sapply(uniq_acc, function(acc) {
      (length(tuple::matchAll(acc, c(df1$Accessions, df2$Accessions, df3$Accessions))) /
        length(c(unique(df1$Patient_ID), unique(df2$Patient_ID), unique(df3$Patient_ID)))) * 100
    })
  }

  # filtern des df
  exclusive_1 <- waterfall_df[which(waterfall_df$set2 == 0 & waterfall_df$set3 == 0), ]
  exclusive_2 <- waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set3 == 0), ]
  exclusive_3 <- waterfall_df[which(waterfall_df$set1 == 0 & waterfall_df$set2 == 0), ]

  shared <- waterfall_df[which(waterfall_df$set1 != 0 & waterfall_df$set2 != 0 & waterfall_df$set3 != 0), ]
  all <- waterfall_df

  # create new waterfall df waterfall_df = rbind(exclusive_1, exclusive_2, exclusive_3, exclusive_4, shared)
  waterfall_df2 <- list(Nec = exclusive_1, T1 = exclusive_2, INF = exclusive_3, Shared = shared, Tumor = all)
  return(waterfall_df2)
}




## TITLE: make Waterfall df ----------------------------------------------------
## DESCRIPTION: makes data frame for waterfall plots
##
## PARAMETERS: - df1 + df2: data frames of the two tumor regions peptides and
##                          source proteins that should be compared [data frame]
##             - with_seq: if true, frequencies will be computed for peptides
##                         if false, frequencies will be computed for source
##                         proteins [boolean]
##
## OUTPUT: data frame of peptide or protein frequencies for the tumor regions
##
make_waterfall_df <- function(df1, df2, with_seq) {
  if (with_seq == TRUE) {
    # union of unique sequencs --> x achse
    uniq_seq <- unique(c(df1$Sequence, df2$Sequence))
    df1 <- df1 %>%
      group_by(Patient_ID) %>%
      summarise(Sequence = unique(Sequence))
    df2 <- df2 %>%
      group_by(Patient_ID) %>%
      summarise(Sequence = unique(Sequence))
    # neues df
    waterfall_df <- data.frame(uniq_seq)
    # ratio of peptide occurences für NEC und INF
    waterfall_df$set1 <- sapply(uniq_seq, function(seq) {
      mean(ifelse(is.na(tuple::matchAll(seq, df1$Sequence)), 0,
        (length(tuple::matchAll(seq, df1$Sequence)) / length(unique(df1$Patient_ID))) * 100
      ))
    })

    waterfall_df$set2 <- sapply(uniq_seq, function(seq) {
      mean(ifelse(is.na(tuple::matchAll(seq, df2$Sequence)), 0,
        (length(tuple::matchAll(seq, df2$Sequence)) / length(unique(df2$Patient_ID))) * 100
      ))
    })

    waterfall_df$Total <- sapply(uniq_seq, function(seq) {
      (length(tuple::matchAll(seq, c(df1$Sequence, df2$Sequence))) /
        length(c(unique(df1$Patient_ID), unique(df2$Patient_ID)))) * 100
    })
  } else {
    (with_seq == FALSE)
  }
  df1 <- df1 %>%
    group_by(Patient_ID) %>%
    summarise(Accessions = unique(Accessions))
  df2 <- df2 %>%
    group_by(Patient_ID) %>%
    summarise(Accessions = unique(Accessions))

  # union of unique acc --> x achse
  uniq_acc <- unique(c(df1$Accessions, df2$Accessions))

  # neues df
  waterfall_df <- data.frame(uniq_acc)
  # ratio of peptide occurences für df1 und df2
  waterfall_df$set1 <- sapply(uniq_acc, function(acc) {
    mean(ifelse(is.na(tuple::matchAll(acc, df1$Accessions)), 0,
      (length(tuple::matchAll(acc, df1$Accessions)) / length(unique(df1$Patient_ID))) * 100
    ))
  })

  waterfall_df$set2 <- sapply(uniq_acc, function(acc) {
    mean(ifelse(is.na(tuple::matchAll(acc, df2$Accessions)), 0,
      (length(tuple::matchAll(acc, df2$Accessions)) / length(unique(df2$Patient_ID))) * 100
    ))
  })

  waterfall_df$Total <- sapply(uniq_acc, function(acc) {
    (length(tuple::matchAll(acc, c(df1$Accessions, df2$Accessions))) /
      length(c(unique(df1$Patient_ID), unique(df2$Patient_ID)))) * 100
  })


  # filtern des df
  exclusive_1 <- waterfall_df[which(waterfall_df$set2 == 0), ]
  exclusive_2 <- waterfall_df[which(waterfall_df$set1 == 0), ]
  # ordnen
  exclusive_1 <- exclusive_1[order(exclusive_1$set1, decreasing = TRUE), ]
  exclusive_2 <- exclusive_2[order(exclusive_2$set2), ]

  other <- waterfall_df[which(waterfall_df$set1 != 0 & waterfall_df$set2 != 0), ]
  other <- other[order(other$set1, other$set2), ]
  # create new waterfall df
  waterfall_df <- rbind(exclusive_1, other, exclusive_2)
  return(waterfall_df)
}



## TITLE: plot Waterfall plot---------------------------------------------------
## DESCRIPTION: function that plots a waterfall plot of a comparison
##
## PARAMETERS: - waterfall_df: data frame of two tumor regions and the respective
##                             peptide/protein frequencies [data frame]
##             - name1: name of the first tumor region [string]
##             - name2: name of the second tumor region [string]
##             - title: title of the plot [string]
##             - with_seq: if true, frequencies will be computed for peptides
##                         if false, frequencies will be computed for source
##                         proteins [boolean]
##
## OUTPUT: Waterfall plot of a pairwise comparison
##
plot_waterfall <- function(waterfall_df, name1, name2, title, with_seq) {
  if (with_seq == TRUE) {
    w <- waterfall_df %>% mutate(uniq_seq = factor(uniq_seq, levels = unique(uniq_seq)))
    x <- w$uniq_seq
    x_lab_text <- "Sequences"
  } else {
    (with_seq == FALSE)
  }
  w <- waterfall_df %>% mutate(uniq_acc = factor(uniq_acc, levels = unique(uniq_acc)))
  x <- w$uniq_acc
  x_lab_text <- "Source Proteins"


  # plot waterfall plot
  ggplot(w, aes(x = x)) +
    geom_bar(aes(y = set1, fill = name1), stat = "identity", width = 1) +
    geom_bar(aes(y = -set2, fill = name2), stat = "identity", width = 1) +
    scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 10)) +
    scale_fill_manual(values = c("salmon1", "steelblue3")) +
    ylab("Frequency of positive ligandomes (%)") +
    xlab(x_lab_text) +
    labs(title = title) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

## TITLE: plot Venn diagram for waterfall plot ---------------------------------
## DESCRIPTION: function that plots avenn diagram for the waterfallplot
##
## PARAMETERS: - listarea1: peptides/proteins of  first tumor region [vector]
##             - listarea2: peptides/proteins of the second tumor region [vector]
##             - title: title of the plot [string]
##
## OUTPUT: Venn diagram of a pairwise comparison
##
plot_venn_waterfall <- function(listarea1, listarea2, title) {
  VennDiagram::draw.pairwise.venn(
    area1 = length(unique(listarea1)), area2 = length(unique(listarea2)),
    cross.area = length(Reduce(dplyr::intersect, list(unique(listarea1), unique(listarea2)))),
    # category = category,
    fill = c("salmon1", "steelblue3"),
    cat.pos = c(0, 0), cat.dist = c(.05, .05),
    cex = 1.5,
    cat.cex = 1.5
  )
}




## TITLE: rewrite HLA types ----------------------------------------------------
## DESCRIPTION: function that rewrites HLA types to a format that is compatible
##              with netMHCpan HLA type input
## PARAMETERS:
##             - uniqe_HLA_types: vector with unique HLA types containing * [vector]
##             - as_string: set to FALSE will return new HLA types as vector
##                           set to TRUE will return a single string of new HLA
##                           types which can be used directly as netMHCpan HLA
##                           type input[boolean]
##
## OUTPUT: either vector or string of reformatted HLA types from A*01:01
##         to HLA-A01:01
##
rewrite_HLA_types <- function(uniqe_HLA_types, as_string) {
  new_HLA_types <- c()
  for (i in uniqe_HLA_types) {
    str2 <- strsplit(i, "*", fixed = TRUE)
    new_HLA_types <- append(
      new_HLA_types,
      paste0("HLA-", str2[[1]][1], str2[[1]][2])
    )
  }
  if (as_string == TRUE) {
    HLA_str <- ""
    for (i in unique(new_HLA_types)) {
      HLA_str <- paste0(HLA_str, i, ",")
    }
    return(HLA_str)
  } else {
    return(unique(new_HLA_types))
  }
}



## TITLE: create df from MHCquant results---------------------------------------
## DESCRIPTION: function to read in ligandomics data from tsv files
##
## PARAMETERS: - filelist: list of all files including the respective path [vector]
##
## OUTPUT: data frame with all detected peptides + source proteins + tumor region
##         information + patient ids
##
create_data_frame <- function(filelist) {
  # create empty classI dataframe
  file_count <- 0
  df <- data.frame(
    Patient_ID = character(),
    Tumor_region = character(),
    Sequence = character(),
    Accessions = character()
  )
  # loop over all files
  for (file in filelist) {
    file_count <- file_count + 1
    tempfile <- read.csv(file, header = TRUE, sep = "\t")
    file_list <- strsplit(file, split = "/")[[1]]
    patient_id <- strsplit(grep("ZH", file_list, value = TRUE), split = "_")[[1]][1]
    t_region <- strsplit(grep("ZH", file_list, value = TRUE), split = "_")[[1]][2]
    seq <- tempfile$sequence
    acc <- tempfile$accessions
    # temporaray df
    temp_df <- data.frame(
      Patient_ID = rep(patient_id, times = length(seq)),
      Sample_num = rep(paste("Sample", as.character(file_count), sep = "_"), times = length(seq)),
      Tumor_region = rep(t_region, times = length(seq)),
      Sequence = seq,
      Accessions = acc
    )
    # append to output df
    df <- rbind(df, temp_df)
  }
  df$Sequence <- str_replace_all(df$Sequence, "\\(Oxidation\\)", "") # get rid of the (Oxidation) string in sequences
  df$Tumor_region <- str_replace_all(df$Tumor_region, "\\.tsv", "")
  df$Sequence <- toupper(df$Sequence) # sequence all upper case
  return(df)
}


## TITLE: plot custom venn -----------------------------------------------------
## DESCRIPTION: function to plot venn diagrams
##
## PARAMETERS: - set: named list of data 
##             - title: title of the venn diagram [string ]
##
## OUTPUT: Venn diagram
##
plot_custom_venn <- function(set, title) {
  venn_data <- process_data(Venn(set))
  venn_plot <- ggVennDiagram(set,
    label_alpha = 0,
    label = "count",
    label_size = 4,
    edge_size = 0
  ) +
    scale_fill_gradient(low = "papayawhip", high = "paleturquoise4") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  plot(venn_plot)
  return(venn_data)
}



## TITLE: length distribution plot----------------------------------------------
## DESCRIPTION: function that plots the lenght distribution of identified peptides
##              in the raw data 
##
## PARAMETERS: - df: dataframe with peptide data and tumor region information
##             - as_bar: if true the len. distribution will be shown as bar plot
##                       if false the len. distribution will be shown as line plot
##             - title: title of the plot [string]
##
## OUTPUT:
##
plot_length_distribution <- function(df, as_bar = TRUE, title) {
  df$len <- apply(df, 2, nchar)[, 4]
  counts <- aggregate(df$len,
    by = list(Tumor_region = df$Tumor_region, len = df$len),
    FUN = table
  )

  if (as_bar == TRUE) {
    ggplot(counts, aes(x = len, y = x, group = Tumor_region, fill = Tumor_region)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = brewer.pal(4, "Set2")) +
      labs(x = "Peptide length", y = "number of HLA peptides") +
      theme_classic() +
      ggtitle(title) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.2))
  } else {
    ggplot(counts, aes(x = len, y = x, group = Tumor_region, color = Tumor_region)) +
      geom_line() +
      scale_color_manual(values = brewer.pal(4, "Set2")) +
      labs(x = "Peptide length", y = "number of HLA peptides") +
      theme_classic() +
      ggtitle(title)
  }
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.2))
}


## TITLE: check length of peptides----------------------------------------------
## DESCRIPTION: functions that takes a list of peptides as input and checks the
##              length of each peptide. Peptides that don't fullfill the length
##              criteria of netMHCpan restrictions will be excluded
##
## PARAMETERS: - list: list of peptide sequences [verctor of strings]
##
## OUTPUT: new list with all peptides that match the netMHCpan criteria
##
# for class II MHC peptides
check_len_II <- function(list) {
  new_list <- c()
  for (i in list) {
    if (nchar(i) >= 9 & nchar(i) < 20000) {
      new_list <- append(new_list, i)
    }
  }
  return(new_list)
}
# for class I MHC peptides
check_len_I <- function(list) {
  new_list <- c()
  for (i in list) {
    if (nchar(i) >= 8 & nchar(i) < 15) {
      new_list <- append(new_list, i)
    }
  }
  return(new_list)
}



## TITLE: saturation analysis plot----------------------------------------------
## DESCRIPTION:
##
## PARAMETERS: -
##
## OUTPUT:
##
# TODO: include saturation analysis 

## TITLE: get protein Accessions------------------------------------------------
## DESCRIPTION: function that rewrite protein Accessions to uniprot IDs
##
## PARAMETERS: - list: list of protein Accessions that should be rewritten to
##                     uniprot IDs [vector]
##
## OUTPUT: list of uniprot IDs
##
# TODO: dont need the for loop try out
get_protein_acc <- function(list) {
  acc_only <- c()
  for (i in list) {
    str <- strsplit(i, "|", fixed = TRUE)
    acc_only <- append(acc_only, str[[1]][2])
  }
  return(acc_only)
}


## TITLE: netMHCpan_to_df_I-----------------------------------------------------
## DESCRIPTION: function that creates a summarised dataframe from netMHCpan
##              binding prediction results for class I HLA peptides
##
## PARAMETERS: - file: text file with the netMHCpan class I results
##
## OUTPUT: - data frame with binding prediction results
##
netMHCpan_results_to_df_I <- function(file) {
  raw_file <- readLines(file)
  list <- raw_file[grepl("<=", raw_file)]

  column_names <- c(
    "Pos", "MHC", "Peptide", "Core", "Of", "Gp", "Gl", "Ip",
    "Il", "Icore", "Identity", "Score_EL", "%Rank_EL",
    "Score_BA", "%Rank_BA", "Aff(nM)", "nix", "BindLevel"
  )

  rows <- as.data.frame(list)
  split_list <- c()
  for (i in rows[, 1]) {
    str <- strsplit(i, " ", fixed = TRUE)
    split_list <- append(split_list, str)
  }
  cleaned_split_list <- list()
  for (i in 1:length(split_list)) {
    short_list <- c()
    for (j in split_list[[i]]) {
      if (j != "") {
        short_list <- append(short_list, j)
      }
      cleaned_split_list[[i]] <- short_list
    }
  }
  for (i in 1:length(cleaned_split_list)) {
    if (length(cleaned_split_list[[i]]) != 18) {
      cleaned_split_list[[i]] <- append(cleaned_split_list[[i]], c("NO", "NO"))
    }
  }

  final_df <- do.call(rbind.data.frame, cleaned_split_list)
  names(final_df) <- column_names
  final_df <- final_df[-17]
  return(final_df)
}


## TITLE: netMHCpan_to_df_II-----------------------------------------------------
## DESCRIPTION: function that creates a summarised dataframe from netMHCpan
##              binding prediction results for class II HLA peptides
##
## PARAMETERS: - file: text file with the netMHCpan class II results
##
## OUTPUT: - data frame with binding prediction results
##
netMHCpan_results_to_df_II <- function(file) {
  raw_file <- readLines(file)
  list <- raw_file[grepl("<=", raw_file)]

  column_names <- c(
    "Pos", "MHC", "Peptide", "Of", "Core", "Core_rel",
    "Identity", "Score_EL", "%Rank_EL", "Exp_Bind ",
    "Score_BA", "Affinity(nM)", "%Rank_BA", "BindLevel"
  )
  rows <- list

  split_list <- c()
  for (i in rows) {
    str <- strsplit(i, " ", fixed = TRUE)
    split_list <- append(split_list, str)
  }
  cleaned_split_list <- list()
  for (i in 1:length(split_list)) {
    short_list <- c()
    for (j in split_list[[i]]) {
      if (j != "") {
        short_list <- append(short_list, j)
      }
      cleaned_split_list[[i]] <- short_list
    }
  }
  for (i in 1:length(cleaned_split_list)) {
    if (length(cleaned_split_list[[i]]) != max(lengths(cleaned_split_list))) {
      difference <- max(lengths(cleaned_split_list)) - length(cleaned_split_list[[i]])
      vec <- rep("NO", times = difference)
      cleaned_split_list[[i]] <- append(cleaned_split_list[[i]], vec)
    }
  }

  final_df <- do.call(rbind.data.frame, cleaned_split_list)
  names(final_df) <- column_names
  return(final_df)
}

## TITLE: summarise_binder_prediction----------------------------------------------------
## DESCRIPTION:
##
## PARAMETERS: -
##
## OUTPUT:
##
summarise_binder_pred <- function(binding_prediction_df) {
  high_freq_pep_df <- data.frame(Peptide = unique(binding_prediction_df$Peptide))

  binding_prediction_df_wb <- binding_prediction_df[binding_prediction_df$BindLevel == "WB", ][2:3]
  binding_prediction_df_sb <- binding_prediction_df[binding_prediction_df$BindLevel == "SB", ][2:3]

  binding_prediction_df_wb <- binding_prediction_df_wb %>%
    group_by(Peptide) %>%
    summarise(MHC = unique(MHC))
  binding_prediction_df_wb <- binding_prediction_df_wb %>%
    group_by(Peptide) %>%
    summarise(toString(c(MHC)))
  names(binding_prediction_df_wb) <- c("Peptide", "Weak_binders")

  binding_prediction_df_sb <- binding_prediction_df_sb %>%
    group_by(Peptide) %>%
    summarise(MHC = unique(MHC))
  binding_prediction_df_sb <- binding_prediction_df_sb %>%
    group_by(Peptide) %>%
    summarise(toString(c(MHC)))
  names(binding_prediction_df_sb) <- c("Peptide", "Strong_binders")

  high_freq_pep_df <- merge(high_freq_pep_df,
    binding_prediction_df_wb,
    by.y = "Peptide",
    all.x = TRUE
  )
  high_freq_pep_df <- merge(high_freq_pep_df,
    binding_prediction_df_sb,
    by.y = "Peptide",
    all.x = TRUE
  )

  return(high_freq_pep_df)
}
