
# MET
pmp_preprocess <- function(pos_df, neg_df, metadata = NULL, samples_key = 'Sample', intens_cols = NULL, 
                           info_cols = NULL, metab_meta = NULL, pca_group = NULL, PCA_Title = NULL,
                           blank_name = NULL, qc_name = NULL,
                           blankFC = 5, max_perc_mv = 0.8, missingPeaksFraction = 0.8, max_rsd = 25, 
                           mv_imp_rowmax = 0.7, mv_imp_colmax = 0.7, mv_imp_method = 'knn'){
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'ggplot2', 'pmp', 'SummarizedExperiment', 'S4Vectors',
            'ggsci', 'stringr','dplyr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  metab_pos <- pos_df
  metab_neg <- neg_df
  
  # Get info columns and intensity columns
  if (is.null(info_cols)) {
    info_cols <- c("Alignment ID","Average Rt(min)","Average Mz","Metabolite name","Adduct type","Post curation result","Fill %","MS/MS assigned","Reference RT","Reference m/z",
                   "Formula","Ontology","INCHIKEY","SMILES","Annotation tag (VS1.0)","RT matched","m/z matched","MS/MS matched","Comment","Manually modified for quantification",
                   "Manually modified for annotation","Isotope tracking parent ID","Isotope tracking weight number","RT similarity","m/z similarity","Simple dot product",
                   "Weighted dot product","Reverse dot product","Matched peaks count","Matched peaks percentage","Total score","S/N average","Spectrum reference file name",
                   "MS1 isotopic spectrum","MS/MS spectrum")
    metab_pos_info <- metab_pos[, colnames(metab_pos) %in% info_cols]
    metab_neg_info <- metab_neg[, colnames(metab_neg) %in% info_cols]
  } else {
    metab_pos_info <- metab_pos[, info_cols]
    metab_neg_info <- metab_neg[, info_cols]
  }
  
  # Get QC, blanks, and sample columns
  if (is.null(intens_cols)) {
    intens_pattern <- paste('Blank_', 'QC_', paste0(samples_key, '_'), sep = '|')
    intens_cols <- str_detect(colnames(metab_pos), intens_pattern)
  }
  
  # Create dataframes for the SummarizedExperiment object
  if (length(intens_cols) != length(intersect(colnames(metab_pos),intens_cols)) | 
      length(intens_cols) != length(intersect(colnames(metab_neg),intens_cols))) {
    Mismatch1 <- setdiff(intens_cols,colnames(metab_pos))
    Mismatch2 <- setdiff(intens_cols,colnames(metab_neg))
    warning("The following 'intens_cols' names were not found in metab_pos: ", paste0(Mismatch1, sep = ' '))
    warning("The following 'intens_cols' names were not found in metab_neg: ", paste0(Mismatch2, sep = ' '))
    stop("Check names and run again")
  }
  
  metab_pos_counts <- metab_pos[, intens_cols] %>% sapply(.,as.numeric) %>% as.matrix()
  metab_neg_counts <- metab_neg[, intens_cols] %>% sapply(.,as.numeric) %>% as.matrix()
  
  # Remove MS/MS samples (not acquired in the same way)
  metab_pos_counts <- metab_pos_counts[, !(colnames(metab_pos_counts) %in% c('MSMS_pos', 'MSMS_neg'))]
  metab_neg_counts <- metab_neg_counts[, !(colnames(metab_neg_counts) %in% c('MSMS_pos', 'MSMS_neg'))]
  
  # Rename the data to indicate ionisation mode
  metab_pos_rownames <- paste0(metab_pos_info$`Alignment ID`, '_pos')
  metab_neg_rownames <- paste0(metab_neg_info$`Alignment ID`, '_neg')
  
  rownames(metab_pos_counts) <- metab_pos_rownames
  rownames(metab_neg_counts) <- metab_neg_rownames
  rownames(metab_pos_info) <- metab_pos_rownames
  rownames(metab_neg_info) <- metab_neg_rownames
  
  # Merge the positive and negative ionisation modes
  metab_counts <- rbind(metab_pos_counts, metab_neg_counts)
  metab_info <- rbind(metab_pos_info, metab_neg_info)
  
  # Create class and group vectors
  if (is.null(metab_meta)) {
    metab_class <- substr(colnames(metab_counts), start = 1, stop = 2)
  } else {
    if (is.null(metadata)) {
      stop("Please provide metadata dataframe")
    }
    metadata_temp <- metadata[colnames(metab_counts),]
    metadata_temp = metadata_temp[,metab_meta] %>% as.data.frame()
    if (length(metab_meta) > 1) {
      colnames(metadata_temp) = c("class",metab_meta[2:length(metab_meta)])
    } else {
      colnames(metadata_temp) = c("class")  
    }
  }
  
  # Alternate steps depending on whether metadata was provided
  if (is.null(metadata)) {
    # Create SummarizedExperiment object
    metab_SE <- SummarizedExperiment(assays = list(counts = metab_counts),
                                     rowData = list(info = metab_info),
                                     colData = DataFrame(class = metab_class))
  } else {
    # Check that the metadata matches the samples
    if (identical(rownames(metadata_temp), colnames(metab_counts)) == TRUE) {
      message("ColData + Counts Names Match")
    } else {
      stop("ColData + Counts Names DO NOT Match!")
    }
    
    # Create SummarizedExperiment object
    metab_SE <- SummarizedExperiment(assays = list(counts = metab_counts),
                                     metadata = list(metadata = metadata),
                                     rowData = list(info = metab_info),
                                     colData = metadata_temp)
  }
  
  message('SummarizedExperiment object created...')
  
  ###
  ### FILTERING AND NORMALISATION
  ###
  
  # Original number of features
  features0 <- dim(metab_SE)
  
  # Replace missing values with NA to be compatible with downstream filtering
  assay(metab_SE) <- replace(assay(metab_SE), assay(metab_SE) == 0, NA)
  
  message('Replaced missing values with NA...')
  
  if (is.null(blank_name)) {
    blank_name = "Bl"
    warning("blank_label = 'Bl' - ensure this is correct")
  }
  if (is.null(qc_name)) {
    qc_name = "QC"
    warning("qc_label = 'QC' - ensure this is correct")
  }
  
  # Filter peaks and samples based on blanks
  metab_filt <- filter_peaks_by_blank(df = metab_SE,
                                      fold_change = blankFC,
                                      classes = metab_SE$class,
                                      remove_samples = TRUE,
                                      remove_peaks = TRUE,
                                      blank_label = blank_name,
                                      qc_label = qc_name)
  
  message('Filtered peaks and samples based on blanks...')
  
  # Number of features
  features1 <- dim(metab_filt)
  
  # Filter samples based on missing values
  metab_filt <- filter_samples_by_mv(df = metab_filt,
                                     max_perc_mv = max_perc_mv)
  
  # Number of features
  features2 <- dim(metab_filt)
  
  # Filter peaks based on missing values
  metab_filt <- filter_peaks_by_fraction(df = metab_filt,
                                         min_frac = missingPeaksFraction,
                                         classes = metab_filt$class,
                                         method = 'across',
                                         qc_label = qc_name)
  
  # Number of features
  features3 <- dim(metab_filt)
  
  message('Filtered peaks and samples based on missing values...')
  
  # Filter peaks based on the % variation in the QC
  metab_filt <- filter_peaks_by_rsd(df = metab_filt,
                                    max_rsd = max_rsd,
                                    classes = metab_filt$class,
                                    qc_label = qc_name)
  
  # Number of features
  features4 <- dim(metab_filt)
  
  message('Filtered peaks based on the percentage variance in QC samples...')
  
  # Data normalisation
  metab_norm <- pqn_normalisation(df = metab_filt,
                                  classes = metab_filt$class,
                                  qc_label = qc_name)
  
  message('Normalised data using PQN...')
  
  if ( sum(colSums(assay(metab_norm) %>% is.na()) / nrow(assay(metab_norm)) > 0.60) > 0) {
    message('Removing samples with >60% missing values')
    Keep <- which(colSums(assay(metab_norm) %>% is.na()) / nrow(assay(metab_norm)) < 0.60) %>% names()
    metab_norm <- metab_norm[,Keep]
  }
  
  message('Imputing values now...')
  
  # Missing values imputation
  metab_imp <- mv_imputation(df = metab_norm,
                             rowmax = mv_imp_rowmax,
                             colmax = mv_imp_colmax,
                             method = mv_imp_method)
  
  message('Finished imputing values...')
  
  # Data scaling
  metab_glog <- glog_transformation(df = metab_imp,
                                    classes = metab_imp$class,
                                    qc_label = qc_name)
  
  message('Scaled data via glog algorithm...')
  
  opt_lambda <- processing_history(metab_glog)$glog_transformation$lambda_opt
  
  glog_plot <- glog_plot_optimised_lambda(df = metab_imp,
                                          optimised_lambda = opt_lambda,
                                          classes = metab_imp$class,
                                          qc_label = qc_name)
  
  # Perform PCA
  
  if (is.null(pca_group)) {
    pca_group = "class"
  }
  
  if (!(pca_group %in% names(metab_glog@colData@listData))) {
    stop("'pca_group' name not in colData - please specify using 'metab_meta' and 'pca_group' correctly")
  }
  
  PCA <- prcomp(t(assay(metab_glog)), center = TRUE)
  varexp <- c(summary(PCA)$importance[2,1]*100, summary(PCA)$importance[2,2]*100)
  
  # Create dataset
  data_PCA <- cbind(data.frame(Samples = rownames(PCA$x),
                               PC1 = PCA$x[,1],
                               PC2 = PCA$x[,2]),
                    plot_group = metab_glog@colData@listData[[paste0(pca_group)]])
  
  if (is.null(PCA_Title)) {
    PCA_Title <- "Metabolomics / Lipidomics"
  }
  
  # Plot results
  PCA_plot <- ggplot(data_PCA, aes(x = PC1, y = PC2)) +
    geom_point(aes_string(fill = paste0("plot_group"), color = paste0("plot_group"))) +
    stat_ellipse(aes_string(fill = paste0("plot_group")), geom = 'polygon', type = 't', level = 0.9, alpha = 0.2) +
    labs(title = paste0(PCA_Title),
         x = paste0('PC1 ', round(varexp[1], 2), '%'),
         y = paste0('PC2 ', round(varexp[2], 2), '%')) +
    scale_fill_jama(name = paste0(pca_group)) +
    scale_color_jama(name = paste0(pca_group))
  
  # Make filtering dimensions dataframe
  filtering_dims <- data.frame(rbind(features0, features1, features2, features3, features4))
  colnames(filtering_dims) <- c('Features', 'Samples')
  rownames(filtering_dims) <- c('Original', 'Blank_filtered', 'MV_sample_filtered', 
                                'MV_peak_filtered', 'QC_var_filtered')
  
  # Prepare elements in a list to return from function
  results <- list(imputed_results = metab_imp,
                  glog_results = metab_glog,
                  glog_plot = glog_plot,
                  PCA_plot = PCA_plot,
                  filtering_dimensions = filtering_dims)
  
  message('Done!')
  
  results
}

# LIP
pmp_preprocess <- function(pos_df, neg_df, metadata = NULL, samples_key = 'Sample', intens_cols = NULL, 
                           info_cols = NULL, metab_meta = NULL, pca_group = NULL, PCA_Title = NULL,
                           blank_name = NULL, qc_name = NULL,
                           blankFC = 5, max_perc_mv = 0.8, missingPeaksFraction = 0.8, max_rsd = 25, 
                           mv_imp_rowmax = 0.7, mv_imp_colmax = 0.7, mv_imp_method = 'knn'){
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'ggplot2', 'pmp', 'SummarizedExperiment', 'S4Vectors',
            'ggsci', 'stringr','dplyr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  metab_pos <- pos_df
  metab_neg <- neg_df
  
  # Get info columns and intensity columns
  if (is.null(info_cols)) {
    info_cols <- c("Alignment ID","Average Rt(min)","Average Mz","Metabolite name","Adduct type","Post curation result","Fill %","MS/MS assigned","Reference RT","Reference m/z",
                   "Formula","Ontology","INCHIKEY","SMILES","Annotation tag (VS1.0)","RT matched","m/z matched","MS/MS matched","Comment","Manually modified for quantification",
                   "Manually modified for annotation","Isotope tracking parent ID","Isotope tracking weight number","RT similarity","m/z similarity","Simple dot product",
                   "Weighted dot product","Reverse dot product","Matched peaks count","Matched peaks percentage","Total score","S/N average","Spectrum reference file name",
                   "MS1 isotopic spectrum","MS/MS spectrum")
    metab_pos_info <- metab_pos[, colnames(metab_pos) %in% info_cols]
    metab_neg_info <- metab_neg[, colnames(metab_neg) %in% info_cols]
  } else {
    metab_pos_info <- metab_pos[, info_cols]
    metab_neg_info <- metab_neg[, info_cols]
  }
  
  # Get QC, blanks, and sample columns
  if (is.null(intens_cols)) {
    intens_pattern <- paste('Blank_', 'QC_', paste0(samples_key, '_'), sep = '|')
    intens_cols <- str_detect(colnames(metab_pos), intens_pattern)
  }
  
  # Create dataframes for the SummarizedExperiment object
  if (length(intens_cols) != length(intersect(colnames(metab_pos),intens_cols)) | 
      length(intens_cols) != length(intersect(colnames(metab_neg),intens_cols))) {
    Mismatch1 <- setdiff(intens_cols,colnames(metab_pos))
    Mismatch2 <- setdiff(intens_cols,colnames(metab_neg))
    warning("The following 'intens_cols' names were not found in metab_pos: ", paste0(Mismatch1, sep = ' '))
    warning("The following 'intens_cols' names were not found in metab_neg: ", paste0(Mismatch2, sep = ' '))
    stop("Check names and run again")
  }
  
  metab_pos_counts <- metab_pos[, intens_cols] %>% sapply(.,as.numeric) %>% as.matrix()
  metab_neg_counts <- metab_neg[, intens_cols] %>% sapply(.,as.numeric) %>% as.matrix()
  
  # Remove MS/MS samples (not acquired in the same way)
  metab_pos_counts <- metab_pos_counts[, !(colnames(metab_pos_counts) %in% c('MSMS_pos', 'MSMS_neg'))]
  metab_neg_counts <- metab_neg_counts[, !(colnames(metab_neg_counts) %in% c('MSMS_pos', 'MSMS_neg'))]
  
  # Rename the data to indicate ionisation mode
  metab_pos_rownames <- paste0(metab_pos_info$`Alignment ID`, '_pos')
  metab_neg_rownames <- paste0(metab_neg_info$`Alignment ID`, '_neg')
  
  rownames(metab_pos_counts) <- metab_pos_rownames
  rownames(metab_neg_counts) <- metab_neg_rownames
  rownames(metab_pos_info) <- metab_pos_rownames
  rownames(metab_neg_info) <- metab_neg_rownames
  
  # Merge the positive and negative ionisation modes
  metab_counts <- rbind(metab_pos_counts, metab_neg_counts)
  metab_info <- rbind(metab_pos_info, metab_neg_info)
  
  # Create class and group vectors
  if (is.null(metab_meta)) {
    metab_class <- substr(colnames(metab_counts), start = 1, stop = 2)
  } else {
    if (is.null(metadata)) {
      stop("Please provide metadata dataframe")
    }
    metadata_temp <- metadata[colnames(metab_counts),]
    metadata_temp = metadata_temp[,metab_meta] %>% as.data.frame()
    if (length(metab_meta) > 1) {
      colnames(metadata_temp) = c("class",metab_meta[2:length(metab_meta)])
    } else {
      colnames(metadata_temp) = c("class")  
    }
  }
  
  # Alternate steps depending on whether metadata was provided
  if (is.null(metadata)) {
    # Create SummarizedExperiment object
    metab_SE <- SummarizedExperiment(assays = list(counts = metab_counts),
                                     rowData = list(info = metab_info),
                                     colData = DataFrame(class = metab_class))
  } else {
    # Check that the metadata matches the samples
    if (identical(rownames(metadata_temp), colnames(metab_counts)) == TRUE) {
      message("ColData + Counts Names Match")
    } else {
      stop("ColData + Counts Names DO NOT Match!")
    }
    
    # Create SummarizedExperiment object
    metab_SE <- SummarizedExperiment(assays = list(counts = metab_counts),
                                     metadata = list(metadata = metadata),
                                     rowData = list(info = metab_info),
                                     colData = metadata_temp)
  }
  
  message('SummarizedExperiment object created...')
  
  ###
  ### FILTERING AND NORMALISATION
  ###
  
  # Original number of features
  features0 <- dim(metab_SE)
  
  # Replace missing values with NA to be compatible with downstream filtering
  assay(metab_SE) <- replace(assay(metab_SE), assay(metab_SE) == 0, NA)
  
  message('Replaced missing values with NA...')
  
  if (is.null(blank_name)) {
    blank_name = "Bl"
    warning("blank_label = 'Bl' - ensure this is correct")
  }
  if (is.null(qc_name)) {
    qc_name = "QC"
    warning("qc_label = 'QC' - ensure this is correct")
  }
  
  # Filter peaks and samples based on blanks
  metab_filt <- filter_peaks_by_blank(df = metab_SE,
                                      fold_change = blankFC,
                                      classes = metab_SE$class,
                                      remove_samples = TRUE,
                                      remove_peaks = TRUE,
                                      blank_label = blank_name,
                                      qc_label = qc_name)
  
  message('Filtered peaks and samples based on blanks...')
  
  # Number of features
  features1 <- dim(metab_filt)
  
  # Filter samples based on missing values
  metab_filt <- filter_samples_by_mv(df = metab_filt,
                                     max_perc_mv = max_perc_mv)
  
  # Number of features
  features2 <- dim(metab_filt)
  
  # Filter peaks based on missing values
  metab_filt <- filter_peaks_by_fraction(df = metab_filt,
                                         min_frac = missingPeaksFraction,
                                         classes = metab_filt$class,
                                         method = 'across',
                                         qc_label = qc_name)
  
  # Number of features
  features3 <- dim(metab_filt)
  
  message('Filtered peaks and samples based on missing values...')
  
  # Filter peaks based on the % variation in the QC
  metab_filt <- filter_peaks_by_rsd(df = metab_filt,
                                    max_rsd = max_rsd,
                                    classes = metab_filt$class,
                                    qc_label = qc_name)
  
  # Number of features
  features4 <- dim(metab_filt)
  
  message('Filtered peaks based on the percentage variance in QC samples...')
  
  # Data normalisation
  metab_norm <- pqn_normalisation(df = metab_filt,
                                  classes = metab_filt$class,
                                  qc_label = qc_name)
  
  message('Normalised data using PQN...')
  
  if ( sum(colSums(assay(metab_norm) %>% is.na()) / nrow(assay(metab_norm)) > 0.60) > 0) {
    message('Removing samples with >60% missing values')
    Keep <- which(colSums(assay(metab_norm) %>% is.na()) / nrow(assay(metab_norm)) < 0.60) %>% names()
    metab_norm <- metab_norm[,Keep]
  }
  
  message('Imputing values now...')
  
  # Missing values imputation
  metab_imp <- mv_imputation(df = metab_norm,
                             rowmax = mv_imp_rowmax,
                             colmax = mv_imp_colmax,
                             method = mv_imp_method)
  
  message('Finished imputing values...')
  
  # Data scaling
  metab_glog <- glog_transformation(df = metab_imp,
                                    classes = metab_imp$class,
                                    qc_label = qc_name)
  
  message('Scaled data via glog algorithm...')
  
  opt_lambda <- processing_history(metab_glog)$glog_transformation$lambda_opt
  
  glog_plot <- glog_plot_optimised_lambda(df = metab_imp,
                                          optimised_lambda = opt_lambda,
                                          classes = metab_imp$class,
                                          qc_label = qc_name)
  
  # Perform PCA
  
  if (is.null(pca_group)) {
    pca_group = "class"
  }
  
  if (!(pca_group %in% names(metab_glog@colData@listData))) {
    stop("'pca_group' name not in colData - please specify using 'metab_meta' and 'pca_group' correctly")
  }
  
  PCA <- prcomp(t(assay(metab_glog)), center = TRUE)
  varexp <- c(summary(PCA)$importance[2,1]*100, summary(PCA)$importance[2,2]*100)
  
  # Create dataset
  data_PCA <- cbind(data.frame(Samples = rownames(PCA$x),
                               PC1 = PCA$x[,1],
                               PC2 = PCA$x[,2]),
                    plot_group = metab_glog@colData@listData[[paste0(pca_group)]])
  
  if (is.null(PCA_Title)) {
    PCA_Title <- "Metabolomics / Lipidomics"
  }
  
  # Plot results
  PCA_plot <- ggplot(data_PCA, aes(x = PC1, y = PC2)) +
    geom_point(aes_string(fill = paste0("plot_group"), color = paste0("plot_group"))) +
    stat_ellipse(aes_string(fill = paste0("plot_group")), geom = 'polygon', type = 't', level = 0.9, alpha = 0.2) +
    labs(title = paste0(PCA_Title),
         x = paste0('PC1 ', round(varexp[1], 2), '%'),
         y = paste0('PC2 ', round(varexp[2], 2), '%')) +
    scale_fill_jama(name = paste0(pca_group)) +
    scale_color_jama(name = paste0(pca_group))
  
  # Make filtering dimensions dataframe
  filtering_dims <- data.frame(rbind(features0, features1, features2, features3, features4))
  colnames(filtering_dims) <- c('Features', 'Samples')
  rownames(filtering_dims) <- c('Original', 'Blank_filtered', 'MV_sample_filtered', 
                                'MV_peak_filtered', 'QC_var_filtered')
  
  # Prepare elements in a list to return from function
  results <- list(imputed_results = metab_imp,
                  glog_results = metab_glog,
                  glog_plot = glog_plot,
                  PCA_plot = PCA_plot,
                  filtering_dimensions = filtering_dims)
  
  message('Done!')
  
  results
}


## GNPS
gnps_format_names_met <- function(gnps_pos_df, gnps_neg_df) {
  
  gnps_remove_tags <- function(vector) {
    vector <- gsub('Spectral Match to (.*) from NIST14', '\\1', vector) # fix "spectal match" tags
    vector <- gsub('ReSpect:[a-zA-Z0-9]* ([^|]*).*', '\\1', vector) # fix "ReSpect" tags
    vector <- gsub('Massbank:[a-zA-Z0-9]* ([^|]*).*', '\\1', vector) # fix "Massbank" tags
    vector <- gsub('MoNA:\\d* (.*)', '\\1', vector) # fix "MoNA" tags
    vector <- gsub('(.*) - \\d{1,3}.\\d{1,3} eV', '\\1', vector) # fix "eV" tags
  }
  
  # Change first column name from '#SCAN#'
  colnames(gnps_pos_df)[1] <- 'Alignment.ID'
  colnames(gnps_neg_df)[1] <- 'Alignment.ID'
  
  # Copy Compound Name column
  gnps_pos_df$compound_name_gnps <- gnps_pos_df$Compound_Name
  gnps_neg_df$compound_name_gnps <- gnps_neg_df$Compound_Name
  
  # Remove the extraneous text around metabolite feature names
  gnps_pos_df$compound_name_gnps <- gnps_remove_tags(gnps_pos_df$compound_name_gnps)
  gnps_neg_df$compound_name_gnps <- gnps_remove_tags(gnps_neg_df$compound_name_gnps)
  
  # Switch to title case
  gnps_pos_df$compound_name_gnps <- stringr::str_to_title(gnps_pos_df$compound_name_gnps)
  gnps_neg_df$compound_name_gnps <- stringr::str_to_title(gnps_neg_df$compound_name_gnps)
  
  # Create vector for new rownames to match pmp rownames
  new_rownames <- c(paste0(gnps_pos_df$Alignment.ID, '_pos'), paste0(gnps_neg_df$Alignment.ID, '_neg'))
  
  # rbind the data.frames
  gnps_df <- rbind(gnps_pos_df, gnps_neg_df) %>%
    dplyr::select(Alignment.ID, compound_name_gnps, everything()) %>%
    mutate(alignment_ionisation = new_rownames)
  rownames(gnps_df) <- new_rownames
  
  gnps_df
}

gnps_format_names_lip <- function(gnps_pos_df, gnps_neg_df) {
  
  gnps_remove_tags <- function(vector) {
    vector <- gsub('Spectral Match to (.*) from NIST14', '\\1', vector) # fix "spectal match" tags
    vector <- gsub('ReSpect:[a-zA-Z0-9]* ([^|]*).*', '\\1', vector) # fix "ReSpect" tags
    vector <- gsub('Massbank:[a-zA-Z0-9]* ([^|]*).*', '\\1', vector) # fix "Massbank" tags
    vector <- gsub('MoNA:\\d* (.*)', '\\1', vector) # fix "MoNA" tags
    vector <- gsub('(.*) - \\d{1,3}.\\d{1,3} eV', '\\1', vector) # fix "eV" tags
  }
  
  # Change first column name from '#SCAN#'
  colnames(gnps_pos_df)[21] <- 'Alignment.ID'
  colnames(gnps_neg_df)[21] <- 'Alignment.ID'
  
  # Copy Compound Name column
  gnps_pos_df$compound_name_gnps <- gnps_pos_df$Compound_Name
  gnps_neg_df$compound_name_gnps <- gnps_neg_df$Compound_Name
  
  # Remove the extraneous text around metabolite feature names
  gnps_pos_df$compound_name_gnps <- gnps_remove_tags(gnps_pos_df$compound_name_gnps)
  gnps_neg_df$compound_name_gnps <- gnps_remove_tags(gnps_neg_df$compound_name_gnps)
  
  # Switch to title case
  gnps_pos_df$compound_name_gnps <- stringr::str_to_title(gnps_pos_df$compound_name_gnps)
  gnps_neg_df$compound_name_gnps <- stringr::str_to_title(gnps_neg_df$compound_name_gnps)
  
  # Create vector for new rownames to match pmp rownames
  new_rownames <- c(paste0(gnps_pos_df$Alignment.ID, '_pos'), paste0(gnps_neg_df$Alignment.ID, '_neg'))
  
  # rbind the data.frames
  gnps_df <- rbind(gnps_pos_df, gnps_neg_df) %>%
    dplyr::select(Alignment.ID, compound_name_gnps, everything()) %>%
    mutate(alignment_ionisation = new_rownames)
  rownames(gnps_df) <- new_rownames 
  
  gnps_df
}

gnps_SE_names <- function(gnps_df, metab_SE) {
  
  # Get the metabolite feature info data from SE object
  #metab_info_temp <- as.data.frame(metab_SE@elementMetadata@listData)
  metab_info_temp <- as.data.frame(rowData(metab_SE))
  metab_info_temp$alignment_ionisation <- rownames(metab_info_temp)
  
  # Select just the GNPS name (and alignment_ionisation column for joining)
  gnps_to_match <- gnps_df %>% dplyr::select(alignment_ionisation, compound_name_gnps)
  
  # Left join the GNPS names to the main table in order to get the right order
  metab_info_temp <- metab_info_temp %>% 
    left_join(gnps_to_match, by = 'alignment_ionisation')
  
  # Obtain the vector with GNPS names
  compound_name_gnps <- metab_info_temp$compound_name_gnps
  
  # Return the vector
  compound_name_gnps
}


## Add HMDB 
## OLD 
add_hmdb <- function(metab_SE, hmdb, mass_tol = 0.002) {
  # Set ion mass
  ion <- 1.007276
  # Transform everything into a vector for faster looping
  hmdb$monisotopic_molecular_weight <- as.numeric(hmdb$monisotopic_molecular_weight)
  hmdb <- hmdb[!is.na(hmdb$monisotopic_molecular_weight),]
  
  Mass_db <- as.vector(as.numeric(hmdb$monisotopic_molecular_weight))
  KEGG_db <- as.vector(hmdb$kegg_id)
  Name_db <- as.vector(hmdb$name)
  HMDB_db <- as.vector(hmdb$accession)
  
  Mass_data <- as.vector(rowData(metab_SE)$`info.Average Mz`)
  Adduct_data <- rowData(metab_SE)$`info.Adduct type`
  
  HMDB_data <- rep(NA, length(Mass_data))
  KEGG_data <- rep(NA, length(Mass_data))
  HMDB_acc <- rep(NA, length(Mass_data)) 
  
  # Get masses corrected for ion precursors
  for (m in 1:length(Mass_data)){
    if (Adduct_data[m] %in% c('[M+2H]2+','[M+H]+','[M+Na]+','[M+NH4]+','[M+NH4]2+')){
      Mass_data[m] <- Mass_data[m]-ion}
    else {
      Mass_data[m] <- Mass_data[m]+ion}}
  # Run loop
  for (n in 1:length(Mass_db)){
    for (m in 1:length(Mass_data)){
      if (is.na(HMDB_data[m])==TRUE & between(Mass_db[n], Mass_data[m]-mass_tol, Mass_data[m]+mass_tol)==TRUE){
        HMDB_data[m] <- Name_db[n]
        KEGG_data[m] <- KEGG_db[n]
        HMDB_acc[m] <- HMDB_db[n]}
      else if (is.na(HMDB_data[m])==FALSE & between(Mass_db[n], Mass_data[m]-mass_tol, Mass_data[m]+mass_tol)==TRUE){
        HMDB_data[m] <- paste(HMDB_data[m],Name_db[n], sep=';')
        KEGG_data[m] <- paste(KEGG_data[m],KEGG_db[n], sep=';')
        HMDB_acc[m] <- paste(HMDB_acc[m],HMDB_db[n], sep=';')}}}
  # Add new information to SE experiment object
  rowData(metab_SE)$HMDB <- HMDB_data
  rowData(metab_SE)$KEGG <- KEGG_data
  rowData(metab_SE)$HMDB_ID <- HMDB_acc
  return(metab_SE)}

## NEW  
add_hmdb <- function(metab_SE, hmdb, mass_tol, cores) { #0.002
  pkgs <- c('foreach', 'doSNOW', 'itertools', 'dplyr', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # Set ion mass
  ion <- 1.007276
  # Transform everything into a vector for faster looping
  hmdb$monisotopic_molecular_weight <- as.numeric(hmdb$monisotopic_molecular_weight)
  hmdb <- hmdb[!is.na(hmdb$monisotopic_molecular_weight),]
  Mass_db <- as.vector(as.numeric(hmdb$monisotopic_molecular_weight))
  KEGG_db <- as.vector(hmdb$kegg)
  Name_db <- as.vector(hmdb$name)
  HMDB_id_db <- as.vector(hmdb$accession)
  chebi_db <- as.vector(hmdb$chebi_id)
  pubchem_db <- as.vector(hmdb$pubchem)
  bigg_db <- as.vector(hmdb$bigg)
  
  Mass_data <- as.vector(rowData(metab_SE)$`info.Average Mz`)
  Adduct_data <- rowData(metab_SE)$`info.Adduct type`
  HMDB_data <- rep(NA, length(Mass_data))
  HMDB_id <- rep(NA, length(Mass_data))
  KEGG_data <- rep(NA, length(Mass_data))
  chebi_data <- rep(NA, length(Mass_data))
  pubchem_data <- rep(NA, length(Mass_data))
  bigg_data <- rep(NA, length(Mass_data))
  
  # Get masses corrected for ion precursors
  for (m in 1:length(Mass_data)){
    if (Adduct_data[m] %in% c('[M+2H]2+','[M+H]+','[M+Na]+','[M+NH4]+','[M+NH4]2+')){
      Mass_data[m] <- Mass_data[m]-ion}
    else {
      Mass_data[m] <- Mass_data[m]+ion}}
  
  Matrix <- matrix(nrow = length(Mass_db), ncol = length(Mass_data))
  rownames(Matrix) <- Mass_db
  colnames(Matrix) <- Mass_data
  #Setup clusters
  cores=cores
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  message("HMDB Annotation Starting - Get a Coffee While You Wait :)")
  start_time <- Sys.time()
  TF_DF <- foreach(i = isplitCols(Matrix, chunks=cores), .combine = "cbind", .packages = c("dplyr")#, .options.snow=opts
  ) %dopar% {
    for (n in 1:nrow(i)) {
      for (m in 1:ncol(i)) {
        i[n,m] <- between(as.numeric(rownames(i)[n]), as.numeric(colnames(i)[m])-mass_tol, as.numeric(colnames(i)[m])+mass_tol)
      }
    }
    i
  }
  end_time <- Sys.time()
  message("Minutes Taken: ", round(end_time - start_time,2))
  
  stopCluster(cl)
  
  #Make names
  for (i in c(1:ncol(TF_DF))) {
    if (sum(TF_DF[,i]*1) > 0) {
      index = which(TF_DF[,i] == T) %>% as.numeric()
      HMDB_data[i] <- paste0(Name_db[index], collapse = ';')
      HMDB_id[i] <- paste0(HMDB_id_db[index], collapse = ';')
      KEGG_data[i] <- paste0(KEGG_db[index], collapse = ';')
      chebi_data[i] <- paste0(chebi_db[index], collapse = ';')
      pubchem_data[i] <- paste0(pubchem_db[index], collapse = ';')
      bigg_data[i] <- paste0(bigg_db[index], collapse = ';')
    }
  }
  # Add new information to SE experiment object
  rowData(metab_SE)$HMDB <- HMDB_data %>% dplyr::na_if(.,'')
  rowData(metab_SE)$KEGG <- KEGG_data %>% dplyr::na_if(.,'')
  rowData(metab_SE)$HMDB_accession <- HMDB_id %>% dplyr::na_if(.,'')
  rowData(metab_SE)$CHEBI_ID <- chebi_data %>% dplyr::na_if(.,'')
  rowData(metab_SE)$PubChem_ID <- pubchem_data %>% dplyr::na_if(.,'')
  rowData(metab_SE)$BiGG_ID <- bigg_data %>% dplyr::na_if(.,'')
  return(metab_SE)
}


## Add LMSD
# Function to list of all annotations from Lipid Maps structural database 
## OLD
add_lmsd <- function(metab_SE, lmsd, mass_tol = 0.002, cores = NA) {
  # Load required packages
  pkgs <- c('foreach', 'doParallel', 'data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # Set ion mass
  ion <- 1.007276
  lmsd$EXACT_MASS <- as.numeric(lmsd$EXACT_MASS)
  Mass_db <- lmsd$EXACT_MASS
  Mass_data <- data.frame(msdial_mz = rowData(metab_SE)$`info.Average Mz`)
  Mass_data$msdial_mz <- as.numeric(Mass_data$msdial_mz)
  rownames(Mass_data) <- rownames(rowData(metab_SE))
  Mass_data$Adduct_data <- rowData(metab_SE)$`info.Adduct type`
  Mass_data$corrected_mz <- NA
  # Get masses corrected for ion precursors
  for (m in 1:nrow(Mass_data)){
    if (Mass_data$Adduct_data[m] %in% c("[M+H]+","[M+NH4]+","[M+2H]2+","[M+H-H2O]+", "[M+Na]+")){
      Mass_data$corrected_mz[m] <- Mass_data$msdial_mz[m]-ion}
    else {
      Mass_data$corrected_mz[m] <- Mass_data$msdial_mz[m]+ion
    }
  }
  #setup parallel backend to use many processors
  if(is.na(cores) | cores > detectCores()){
    cores <- detectCores()
  }
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  ## Compare masses corrected with lmsd "exact masses" to find annotations that are with tolerance range
  full_lmsd_ann <- foreach(n=1:nrow(Mass_data), .combine=rbind, .packages = "dplyr") %dopar% {
    row_num <- c()
    for (m in 1:length(Mass_db)) {
      if(between(Mass_db[m], Mass_data$corrected_mz[n]-mass_tol, Mass_data$corrected_mz[n]+mass_tol)==TRUE)
        row_num <- c(row_num, m)
    }
    if(is.null(row_num)) { ## if no matches are found we want to add the data as NA
      temp <- data.frame(matrix(ncol = ncol(lmsd)))
      colnames(temp) <- colnames(lmsd)
      temp$LipidID <- rownames(Mass_data)[n]
      temp$corrected_mz <- Mass_data$corrected_mz[n]
      temp$delta <- 1 ## delta 1 = no match
    } else {
      temp <- lmsd[row_num,]
      temp$LipidID <- rownames(Mass_data)[n]
      temp$corrected_mz <- Mass_data$corrected_mz[n]
      temp$delta <- abs(Mass_data$corrected_mz[n] - lmsd$EXACT_MASS[row_num])
    }
    temp
  }
  stopCluster(cl)
  
  # arrange smallest delta for each lipidID
  full_lmsd_ann <- full_lmsd_ann %>%
    arrange(delta, by_group = T) %>%
    arrange(factor(LipidID, levels = rownames(Mass_data))) # maintain original order
  
  # Remove duplicate matches for lipidID
  distinct_lmsd_ann <- full_lmsd_ann %>%
    distinct(LipidID, NAME, SYSTEMATIC_NAME, .keep_all=T)
  
  
  # Select columns of interest
  lmsd_ann_sub <- distinct_lmsd_ann[, c("LipidID","corrected_mz","delta",
                                        "NAME","SYSTEMATIC_NAME", "CATEGORY",
                                        "MAIN_CLASS","EXACT_MASS", "ABBREVIATION",
                                        "SYNONYMS","KEGG_ID","HMDB_ID", "SUB_CLASS",
                                        "CLASS_LEVEL4")]
  # Change column names
  colnames(lmsd_ann_sub) <- c("LipidID","corrected_mz","delta",
                              "LMSD_NAME","LMSD_SYSTEMATIC_NAME", "LMSD_CATEGORY",
                              "LMSD_MAIN_CLASS","LMSD_EXACT_MASS", "LMSD_ABBREVIATION",
                              "LMSD_SYNONYMS","LMSD_KEGG_ID","LMSD_HMDB_ID", "LMSD_SUB_CLASS",
                              "LMSD_CLASS_LEVEL4")
  
  # Create aggregate string of all lipid annotations for each LipidID
  sub_lmsd <- lmsd_ann_sub[, c("LipidID", "delta", "LMSD_NAME", "LMSD_SYSTEMATIC_NAME",
                               "LMSD_ABBREVIATION", "LMSD_SYNONYMS")]
  sub_lmsd$delta <- round(sub_lmsd$delta,6)
  agg_lmsd <- aggregate(. ~ LipidID, data = sub_lmsd, function(x) paste(unique(x), collapse = " ; "), na.action = na.pass)
  agg_lmsd[agg_lmsd == "NA"] <- NA
  rownames(agg_lmsd) <- agg_lmsd$LipidID
  agg_lmsd <- agg_lmsd[rownames(metab_SE),]
  
  # Select 1st lipid for multiple lipid matches (lowest delta should be at the top of list)
  top_LMSD_match <- lmsd_ann_sub %>%
    group_by(LipidID) %>%
    slice_head() %>%
    arrange(factor(LipidID, levels = rownames(Mass_data)))
  
  # Get the metabolite feature info data from SE object
  metab_info_temp <- data.frame(metab_SE@elementMetadata@listData, check.names = F, stringsAsFactors = T) %>%
    mutate(LipidID = rownames(.))
  
  # Left join with to the main table using LipidID to get the correct ordering
  SE_metadata_added_lmsd <- metab_info_temp %>%
    left_join(top_LMSD_match, by = 'LipidID')
  rownames(SE_metadata_added_lmsd) <- SE_metadata_added_lmsd$LipidID
  
  lmsd_list <- list(
    "full_lmsd_ann" = distinct_lmsd_ann, ## full but with duplicates removed
    "agg_lmsd_df" = agg_lmsd,
    "metadata_lmsd_table" = SE_metadata_added_lmsd
  )
  
  return(lmsd_list)
}

## NEW
add_lmsd <- function(metab_SE, lmsd, mass_tol = 0.002, cores = NA) {
  # Load required packages
  pkgs <- c('foreach', 'doParallel', 'data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'doSNOW')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # Set ion mass
  ion <- 1.007276
  lmsd$EXACT_MASS <- as.numeric(lmsd$EXACT_MASS)
  Mass_db <- lmsd$EXACT_MASS
  Mass_data <- data.frame(msdial_mz = rowData(metab_SE)$`info.Average Mz`)
  Mass_data$msdial_mz <- as.numeric(Mass_data$msdial_mz)
  rownames(Mass_data) <- rownames(rowData(metab_SE))
  Mass_data$Adduct_data <- rowData(metab_SE)$`info.Adduct type`
  Mass_data$corrected_mz <- NA
  # Get masses corrected for ion precursors
  for (m in 1:nrow(Mass_data)){
    if (Mass_data$Adduct_data[m] %in% c("[M+H]+","[M+NH4]+","[M+2H]2+","[M+H-H2O]+", "[M+Na]+")){
      Mass_data$corrected_mz[m] <- Mass_data$msdial_mz[m]-ion}
    else {
      Mass_data$corrected_mz[m] <- Mass_data$msdial_mz[m]+ion
    }
  }
  #setup parallel backend to use many processors
  if(is.na(cores) | cores > detectCores()){
    cores <- detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
  } else {
    cl <- makeCluster(cores)
  }
  registerDoSNOW(cl)
  #Progress Bar
  iterations <- nrow(Mass_data)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  ## Compare masses corrected with lmsd "exact masses" to find annotations that are with tolerance range
  full_lmsd_ann <- foreach(n=1:nrow(Mass_data), .combine=rbind, .packages = "dplyr", .options.snow = opts) %dopar% {
    row_num <- c()
    for (m in 1:length(Mass_db)) {
      if(between(Mass_db[m], Mass_data$corrected_mz[n]-mass_tol, Mass_data$corrected_mz[n]+mass_tol)==TRUE)
        row_num <- c(row_num, m)
    }
    if(is.null(row_num)) { ## if no matches are found we want to add the data as NA
      temp <- data.frame(matrix(ncol = ncol(lmsd)))
      colnames(temp) <- colnames(lmsd)
      temp$LipidID <- rownames(Mass_data)[n]
      temp$corrected_mz <- Mass_data$corrected_mz[n]
      temp$delta <- 1 ## delta 1 = no match
    } else {
      temp <- lmsd[row_num,]
      temp$LipidID <- rownames(Mass_data)[n]
      temp$corrected_mz <- Mass_data$corrected_mz[n]
      temp$delta <- abs(Mass_data$corrected_mz[n] - lmsd$EXACT_MASS[row_num])
    }
    temp
  }
  close(pb)
  stopCluster(cl)
  
  # arrange smallest delta for each lipidID
  full_lmsd_ann <- full_lmsd_ann %>%
    arrange(delta, by_group = T) %>%
    arrange(factor(LipidID, levels = rownames(Mass_data))) # maintain original order
  
  # Remove duplicate matches for lipidID
  distinct_lmsd_ann <- full_lmsd_ann %>%
    distinct(LipidID, NAME, SYSTEMATIC_NAME, .keep_all=T)
  
  
  # Select columns of interest
  lmsd_ann_sub <- distinct_lmsd_ann[, c("LipidID","corrected_mz","delta",
                                        "NAME","SYSTEMATIC_NAME", "CATEGORY",
                                        "MAIN_CLASS","EXACT_MASS", "ABBREVIATION",
                                        "SYNONYMS","KEGG_ID","HMDB_ID", "SUB_CLASS",
                                        "CLASS_LEVEL4",
                                        "LIPIDMAPS_ID","CHEBI_ID","PUBCHEM_ID")] %>%
    rename(HMDB_accession = HMDB_ID, KEGG = KEGG_ID, PubChem_ID = PUBCHEM_ID, LMSD = NAME) # to make consistent with metabolomics dataset

  # Change column names
 # colnames(lmsd_ann_sub) <- c("LipidID","corrected_mz","delta",
 #                             "LMSD_NAME","LMSD_SYSTEMATIC_NAME", "LMSD_CATEGORY",
 #                             "LMSD_MAIN_CLASS","LMSD_EXACT_MASS", "LMSD_ABBREVIATION",
 #                             "LMSD_SYNONYMS","LMSD_KEGG_ID","LMSD_HMDB_ID", "LMSD_SUB_CLASS",
 #                             "LMSD_CLASS_LEVEL4",
 #                             "LIPIDMAPS_ID","CHEBI_ID","PUBCHEM_ID")
  
  # Create aggregate string of all lipid annotations for each LipidID
  sub_lmsd <- lmsd_ann_sub[, c("LipidID", "delta", "LMSD", "SYSTEMATIC_NAME",
                               "ABBREVIATION", "SYNONYMS")]
  sub_lmsd$delta <- round(sub_lmsd$delta,6)
  agg_lmsd <- aggregate(. ~ LipidID, data = sub_lmsd, function(x) paste(unique(x), collapse = " ; "), na.action = na.pass)
  agg_lmsd[agg_lmsd == "NA"] <- NA
  rownames(agg_lmsd) <- agg_lmsd$LipidID
  agg_lmsd <- agg_lmsd[rownames(metab_SE),]
  
  # Select 1st lipid for multiple lipid matches (lowest delta should be at the top of list)
  top_LMSD_match <- lmsd_ann_sub %>%
    group_by(LipidID) %>%
    slice_head() %>%
    arrange(factor(LipidID, levels = rownames(Mass_data)))
  
  # Get the metabolite feature info data from SE object
  metab_info_temp <- data.frame(metab_SE@elementMetadata@listData, check.names = F, stringsAsFactors = T) %>% rownames_to_column(var = "LipidID")
  
  # Left join with to the main table using LipidID to get the correct ordering
  SE_metadata_added_lmsd <- metab_info_temp %>%
    left_join(top_LMSD_match, by = 'LipidID')
  rownames(SE_metadata_added_lmsd) <- SE_metadata_added_lmsd$LipidID
  
  lmsd_list <- list(
    "full_lmsd_ann" = distinct_lmsd_ann, ## full but with duplicates removed
    "agg_lmsd_df" = agg_lmsd,
    "metadata_lmsd_table" = SE_metadata_added_lmsd
  )
  
  return(lmsd_list)
}

## Compare annotations 
## OLD
compare_annotations <- function(metab_SE) {
  
  # Prepare data.frame with alignment IDs and all four annotations and filter for at least one annotation
  msdial_gnps_hmdb <- data.frame('Alignment.ID' = rownames(metab_SE),
                                 'Retention.Time' = rowData(metab_SE)$`info.Average Rt(min)`,
                                 'HMDB_annotation' = rowData(metab_SE)$HMDB,
                                 'KEGG_annotation' = rowData(metab_SE)$KEGG,
                                 'HMDB_id' = rowData(metab_SE)$HMDB_ID,
                                 'MSDIAL_annotation' = rowData(metab_SE)$`info.Metabolite name`,
                                 'GNPS_annotation' = rowData(metab_SE)$compound_name_gnps
  ) %>%
    column_to_rownames(var = 'Alignment.ID') %>%
    mutate(MSDIAL_annotation = replace(MSDIAL_annotation, MSDIAL_annotation == 'Unknown', NA)) %>%
    filter(!is.na(MSDIAL_annotation) | !is.na(GNPS_annotation) | !is.na(HMDB_annotation) | !is.na(KEGG_annotation) | !is.na(HMDB_id))
  
  # Return the data.frame
  msdial_gnps_hmdb
}
# NEW
compare_annotations_met <- function(metab_SE) {
  # Load packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Prepare data.frame with alignment IDs and all four annotations and filter for at least one annotation
  msdial_gnps_hmdb <- data.frame('Alignment.ID' = rownames(metab_SE),
                                 'Retention.Time' = rowData(metab_SE)$`info.Average Rt(min)`,
                                 'Fill_percent' = rowData(metab_SE)$`info.Fill %`,
                                 'S/N_ratio' = rowData(metab_SE)$`info.S/N average`,
                                 'MSDIAL_annotation' = rowData(metab_SE)$`info.Metabolite name`,
                                 #'GNPS_annotation' = rowData(metab_SE)$compound_name_gnps,
                                 'HMDB_annotation' = rowData(metab_SE)$HMDB,
                                 'HMDB_accession' = rowData(metab_SE)$HMDB_accession,
                                 'KEGG_annotation' = rowData(metab_SE)$KEGG)
  msdial_gnps_hmdb <- msdial_gnps_hmdb %>%
    column_to_rownames(var = 'Alignment.ID') %>%
    mutate(MSDIAL_annotation = replace(MSDIAL_annotation, MSDIAL_annotation == 'Unknown', NA),
           KEGG_annotation= replace(KEGG_annotation, KEGG_annotation == '', NA)) %>%
    filter(!is.na(MSDIAL_annotation) | #!is.na(GNPS_annotation) | 
             !is.na(HMDB_annotation) | !is.na(KEGG_annotation))
  
  # Return the data.frame
  msdial_gnps_hmdb
}


## OLD
compare_annotations_lip <- function(metab_SE, agg_lmsd_ann) {
  # Load packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # LMSD_annotation will prioritise use of LMSD_NAME,then LMSD_ABBREVIATION and lastly LMSD_SYSTEMATIC_NAME
  LMSD_ann <- ifelse(is.na(agg_lmsd_ann$LMSD_NAME),
                     ifelse(is.na(agg_lmsd_ann$LMSD_ABBREVIATION),
                            ifelse(is.na(agg_lmsd_ann$LMSD_SYSTEMATIC_NAME), NA, agg_lmsd_ann$LMSD_SYSTEMATIC_NAME),
                            agg_lmsd_ann$LMSD_ABBREVIATION),
                     agg_lmsd_ann$LMSD_NAME)
  
  # Prepare data.frame with alignment IDs and all four annotations and filter for at least one annotation
  msdial_lmsd_gnps_hmdb <- data.frame('LipidID' = rowData(metab_SE)$LipidID,
                                      'Mz' = rowData(metab_SE)$`info.Average Mz`,
                                      'RT' = rowData(metab_SE)$`info.Average Rt(min)`,
                                      'MSDIAL_annotation' = rowData(metab_SE)$`info.Metabolite name`,
                                      'LMSD_annotation' = LMSD_ann,
                                      'GNPS_annotation' = rowData(metab_SE)$compound_name_gnps,
                                      'HMDB_annotation' = rowData(metab_SE)$HMDB,
                                      'KEGG_annotation' = rowData(metab_SE)$KEGG) %>%
    mutate(MSDIAL_annotation = replace(MSDIAL_annotation, MSDIAL_annotation == 'Unknown', NA),
           KEGG_annotation= replace(KEGG_annotation, KEGG_annotation == '', NA)) %>%
    filter(!is.na(MSDIAL_annotation) | !is.na(GNPS_annotation) | !is.na(HMDB_annotation) |
             !is.na(KEGG_annotation) | !is.na(LMSD_annotation))
  
  # Return the data.frame
  return(msdial_lmsd_gnps_hmdb)
}

# NEW
compare_annotations_lip <- function(metab_SE, agg_lmsd_ann) {
  # Load packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # LMSD_annotation will prioritise use of LMSD_NAME,then LMSD_ABBREVIATION and lastly LMSD_SYSTEMATIC_NAME
  LMSD_ann <- data.frame("LMSD_annotation" = ifelse(is.na(agg_lmsd_ann$LMSD),
                                        ifelse(is.na(agg_lmsd_ann$ABBREVIATION),
                                               ifelse(is.na(agg_lmsd_ann$SYSTEMATIC_NAME), NA, agg_lmsd_ann$SYSTEMATIC_NAME),
                                               agg_lmsd_ann$ABBREVIATION),
                                        agg_lmsd_ann$LMSD),
                         "LipidID" = agg_lmsd_ann$LipidID)

  # Prepare data.frame with alignment IDs and all four annotations and filter for at least one annotation
  msdial_lmsd_gnps_hmdb <- data.frame('LipidID' = rowData(metab_SE)$LipidID,
                                      'Mz' = rowData(metab_SE)$`info.Average Mz`,
                                      'RT' = rowData(metab_SE)$`info.Average Rt(min)`,
                                      'MSDIAL_annotation' = rowData(metab_SE)$`info.Metabolite name`) %>%
    left_join(LMSD_ann, by = join_by(LipidID) ) %>%
    mutate(MSDIAL_annotation = replace(MSDIAL_annotation, MSDIAL_annotation == 'Unknown', NA)# , 
           #KEGG_annotation= replace(KEGG_annotation, KEGG_annotation == '', NA)
    ) %>%
    filter(!is.na(MSDIAL_annotation) | !is.na(LMSD_annotation)  )
  
  # Return the data.frame
  return(msdial_lmsd_gnps_hmdb)
}


## Keep annotated
# Keep annotated only, and make a column with a consensus name 
## OLD
keep_annotated <-function(metab_SE) {
  # Keep only rows with at least one annotation
  metab_SE <- metab_SE[(!is.na(rowData(metab_SE)$HMDB) | !is.na(rowData(metab_SE)$compound_name_gnps) | !is.na(rowData(metab_SE)$HMDB_ID) |
                          rowData(metab_SE)$`info.Metabolite name` != 'Unknown'), ]
  # Create a rowData ionisation variable
  rowData(metab_SE)$ionisation <- gsub('\\d*_(pos|neg)', '\\1', rownames(metab_SE))
  
  # Add names to "shortname" by preference: HMDB > GNPS > MS-DIAL
  for (n in 1:nrow(metab_SE)) {
    if (!is.na(rowData(metab_SE)$HMDB)[n]) {
      rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$HMDB[n]
    } else if (!is.na(rowData(metab_SE)$compound_name_gnps)[n]) {
      rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$compound_name_gnps[n]
    } else {
      rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$`info.Metabolite name`[n]
    }
  }
  
  # Remove starting tags from the names
  rowData(metab_SE)$shortname <- gsub('w/o MS2:', '', 
                                      gsub('; LC-ESI-QQ; MS2; CE', '', 
                                           rowData(metab_SE)$shortname))
  # Keep only the first name out of a longer list
  rowData(metab_SE)$shortname <- gsub('([^;]*);.*', '\\1', rowData(metab_SE)$shortname)
  # Add the ionisation mode
  rowData(metab_SE)$shortname <- paste0(rowData(metab_SE)$shortname, 
                                        '_', rowData(metab_SE)$ionisation)
  # Return the SE object
  metab_SE
}
## NEW
keep_annotated_met <- function(metab_SE) {
  # Keep only rows with at least one annotation
  metab_SE <- metab_SE[(!is.na(rowData(metab_SE)$HMDB) | #!is.na(rowData(metab_SE)$compound_name_gnps) |
                          rowData(metab_SE)$`info.Metabolite name` != 'Unknown'), ]
  # Create a rowData ionisation variable
  rowData(metab_SE)$ionisation <- gsub('\\d*_(pos|neg)', '\\1', rownames(metab_SE))
  
  # Add names to "shortname" by preference: HMDB > GNPS > MS-DIAL
  for (n in 1:nrow(metab_SE)) {
    if (!is.na(rowData(metab_SE)$HMDB)[n]) {
      rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$HMDB[n]
    } #else if (!is.na(rowData(metab_SE)$compound_name_gnps)[n]) {
    #rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$compound_name_gnps[n]
    else {
      rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$`info.Metabolite name`[n]
    }
  }
  
  # Remove starting tags from the names
  rowData(metab_SE)$shortname <- gsub('w/o MS2:', '', 
                                      gsub('; LC-ESI-QQ; MS2; CE', '', 
                                           rowData(metab_SE)$shortname))
  # Remove leftover trailing whitespace
  rowData(metab_SE)$shortname <- trimws(rowData(metab_SE)$shortname)
  # Keep only the first name out of a longer list
  rowData(metab_SE)$shortname <- gsub('([^;]*);.*', '\\1', rowData(metab_SE)$shortname)
  # Add the ionisation mode, and then make unique
  rowData(metab_SE)$shortname <- make.unique(paste0(rowData(metab_SE)$shortname, 
                                                    '_', rowData(metab_SE)$ionisation))
  # Clean up kegg column
  rowData(metab_SE)$KEGG <- gsub('NA;|;NA|NA', '', rowData(metab_SE)$KEGG)
  rowData(metab_SE)$KEGG <- str_squish(rowData(metab_SE)$KEGG)

  # Return the SE object
  metab_SE
}


## OLD
keep_annotated_lip <- function(metab_SE) {
  # Keep only rows with at least one annotation
  metab_SE <- metab_SE[(!is.na(rowData(metab_SE)$HMDB) |
                          !is.na(rowData(metab_SE)$LMSD_NAME) | !is.na(rowData(metab_SE)$LMSD_SYSTEMATIC_NAME) |
                          !is.na(rowData(metab_SE)$LMSD_ABBREVIATION) | rowData(metab_SE)$`info.Metabolite name` != 'Unknown'), ]
  
  # Create a rowData ionisation variable
  rowData(metab_SE)$ionisation <- gsub('\\d*_(pos|neg)', '\\1', rownames(metab_SE))
  
  # Add names to "shortname" by preference: LMSD > HMDB > MS-DIAL > GNPS
  for (n in 1:nrow(metab_SE)) {
    
    if (!is.na(rowData(metab_SE)$LMSD_NAME[n]) | !is.na(rowData(metab_SE)$LMSD_ABBREVIATION[n])) {
      
      #rowData(metab_SE)$shortname[n] <- ifelse(is.na(rowData(metab_SE)$LMSD_NAME[n]),
      #  ifelse(is.na(rowData(metab_SE)$LMSD_ABBREVIATION[n]), rowData(metab_SE)$LMSD_SYSTEMATIC_NAME[n], rowData(metab_SE)$LMSD_ABBREVIATION[n]), rowData(metab_SE)$LMSD_NAME[n])
      
      rowData(metab_SE)$shortname[n] <- ifelse(is.na(rowData(metab_SE)$LMSD_NAME[n]), rowData(metab_SE)$LMSD_ABBREVIATION[n], rowData(metab_SE)$LMSD_NAME[n])
      
    } else if (!is.na(rowData(metab_SE)$HMDB[n]))  {
      
      rowData(metab_SE)$shortname[n] <- as.character(rowData(metab_SE)$HMDB[n])
      
    } else if  (rowData(metab_SE)$`info.Metabolite name`[n] != "Unknown" & !(str_detect(rowData(metab_SE)$`info.Metabolite name`[n], "RIKEN"))) {
      
      rowData(metab_SE)$shortname[n] <- as.character(rowData(metab_SE)$`info.Metabolite name`[n])
      
    } else if (!is.na(rowData(metab_SE)$LMSD_SYSTEMATIC_NAME[n])) {
      
      rowData(metab_SE)$shortname[n] <- as.character(rowData(metab_SE)$`LMSD_SYSTEMATIC_NAME`[n])
      
    } else { ## This is for the RIKEN match - which will be labeled NA and removed
      rowData(metab_SE)$shortname[n] <- NA
    }
  }
  
  #Remove any NA's
  metab_SE <- metab_SE[(!is.na(rowData(metab_SE)$shortname)), ]
  
  # Remove starting tags from the names
  rowData(metab_SE)$shortname <- gsub('w/o MS2:', '',
                                      gsub('; LC-ESI-QQ; MS2; CE', '',
                                           rowData(metab_SE)$shortname))
  # Keep only the first name out of a longer list
  rowData(metab_SE)$shortname <- gsub(' \\| .*', '', rowData(metab_SE)$shortname)
  rowData(metab_SE)$shortname <- gsub('\\|.*', '', rowData(metab_SE)$shortname)
  # rowData(metab_SE)$shortname <- gsub('([^;]*);[0-9]O', '\\1', rowData(metab_SE)$shortname) 
  # Add the ionisation mode, and then make unique
  rowData(metab_SE)$shortname <- make.unique(paste0(rowData(metab_SE)$shortname,
                                                    '|', rowData(metab_SE)$ionisation))
  
  # Clean up kegg column
  rowData(metab_SE)$KEGG <- gsub('[[:punct:]]', ' ', rowData(metab_SE)$KEGG) 
  rowData(metab_SE)$KEGG <- str_squish(rowData(metab_SE)$KEGG)
  rowData(metab_SE)$KEGG <- substring(rowData(metab_SE)$KEGG, 1, 6)
  
  rowData(metab_SE)$HMDB_ID <- gsub('[[:punct:]]', ' ', rowData(metab_SE)$HMDB_ID) 
  rowData(metab_SE)$HMDB_ID <- str_squish(rowData(metab_SE)$HMDB_ID)
  rowData(metab_SE)$HMDB_ID <- substring(rowData(metab_SE)$HMDB_ID, 1, 11)
  
  # Return the SE object
  return(metab_SE)
}

## NEW
keep_annotated_lip <- function(metab_SE) {
  
  # Keep only rows with at least one annotation
  metab_SE <- metab_SE[(#!is.na(rowData(metab_SE)$HMDB) | 
    rowData(metab_SE)$`info.Metabolite name` != 'Unknown' | 
      !is.na(rowData(metab_SE)$LMSD) | 
      !is.na(rowData(metab_SE)$SYSTEMATIC_NAME) |
      !is.na(rowData(metab_SE)$ABBREVIATION)), ]
  # Create a rowData ionisation variable
  rowData(metab_SE)$ionisation <- gsub('\\d*_(pos|neg)', '\\1', rownames(metab_SE))
  # Add names to "shortname" by preference: LMSD > MS-DIAL 
  for (n in 1:nrow(metab_SE)) {
    # LMSD 
    if  (!is.na(rowData(metab_SE)$LMSD[n]) | !is.na(rowData(metab_SE)$SYSTEMATIC_NAME[n]) | !is.na(rowData(metab_SE)$ABBREVIATION[n])) {
      rowData(metab_SE)$shortname[n] <- ifelse(is.na(rowData(metab_SE)$LMSD[n]),
                                               ifelse(is.na(rowData(metab_SE)$ABBREVIATION[n]),
                                                      rowData(metab_SE)$SYSTEMATIC_NAME[n],
                                                      rowData(metab_SE)$ABBREVIATION[n]),
                                               rowData(metab_SE)$LMSD[n])
     # # HMDB
     # } else if (!is.na(rowData(metab_SE)$HMDB[n]) )  {
     # rowData(metab_SE)$shortname[n] <- rowData(metab_SE)$HMDB[n]
      # MSDIAL
      } else if (rowData(metab_SE)$`info.Metabolite name`[n] != "Unknown" & !(str_detect(rowData(metab_SE)$`info.Metabolite name`[n], "RIKEN"))) {
        rowData(metab_SE)$shortname[n] <- as.character(rowData(metab_SE)$`info.Metabolite name`[n]) %>% gsub(";","#",.)
      } else { ## This is for the RIKEN match - which will be labeled NA and removed
      rowData(metab_SE)$shortname[n] <- NA 
      }
  }
  
  #Remove any NA's
  metab_SE <- metab_SE[(!is.na(rowData(metab_SE)$shortname)), ]
  
  # Remove starting tags from the names
  rowData(metab_SE)$shortname <- gsub('w/o MS2:', '',
                                      gsub('; LC-ESI-QQ; MS2; CE', '',
                                           rowData(metab_SE)$shortname))
  # Remove leftover trailing whitespace
  rowData(metab_SE)$shortname <- trimws(rowData(metab_SE)$shortname)
  # Keep only the first name out of a longer list
  rowData(metab_SE)$shortname <- gsub(' \\| .*', '', rowData(metab_SE)$shortname)
  rowData(metab_SE)$shortname <- gsub('([^;]*);.*', '\\1', rowData(metab_SE)$shortname)
  # Add the ionisation mode, and then make unique
  rowData(metab_SE)$shortname <- make.unique(paste0(rowData(metab_SE)$shortname,
                                                    '_', rowData(metab_SE)$ionisation))
  rowData(metab_SE)$shortname <- gsub("#",";", rowData(metab_SE)$shortname) # Fix MS-DIAL names back now
  # Return the SE object
  return(metab_SE)
}


# Function to save a curation table 
## OLD
save_curation_table <- function(metab_SE, filename) {
  
  df <- data.frame('Alignment_ID' = metab_SE@elementMetadata$`info.Alignment ID`,
                   'Metabolite_name' = metab_SE@elementMetadata$`info.Metabolite name`,
                   'Ionisation' = metab_SE@elementMetadata$ionisation,
                   'shortname' = metab_SE@elementMetadata$shortname,
                   'Fill%' = metab_SE@elementMetadata$`info.Fill %`,
                   'S/N' = metab_SE@elementMetadata$`info.S/N average`,
                   'ID_ion' = paste0(metab_SE@elementMetadata$`info.Alignment ID`, metab_SE@elementMetadata$ionisation))
  write_csv(df, file = filename)
}

## NEW
save_curation_table_met <- function(metab_SE, filename) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'readr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  df <- data.frame('Alignment_ID' = metab_SE@elementMetadata$`info.Alignment ID`,
                   'Metabolite_name' = metab_SE@elementMetadata$`info.Metabolite name`,
                   'Ionisation' = metab_SE@elementMetadata$ionisation,
                   'shortname' = metab_SE@elementMetadata$shortname,
                   'Fill' = metab_SE@elementMetadata$`info.Fill %`,
                   'S/N' = metab_SE@elementMetadata$`info.S/N average`,
                   'RT_avg' = metab_SE@elementMetadata$`info.Average Rt(min)`)
  write_csv(df, file = filename)
}

## NEW
save_curation_table_lip <- function(metab_SE, filename) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'readr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  df <- data.frame('Alignment_ID' = metab_SE@elementMetadata$`info.Alignment ID`,
                   'Lipid_name' = metab_SE@elementMetadata$`info.Metabolite name`,
                   'Ionisation' = metab_SE@elementMetadata$ionisation,
                   'shortname' = metab_SE@elementMetadata$shortname,
                   'Fill' = metab_SE@elementMetadata$`info.Fill %`,
                   'S/N' = metab_SE@elementMetadata$`info.S/N average`,
                   'RT_avg' = metab_SE@elementMetadata$`info.Average Rt(min)`)
  write_csv(df, file = filename)
}

metab_limma_categorical <- function(metab_SE, metadata_var, metadata_condition = NULL, model_matrix = NULL, 
                                    contrast_matrix = NULL, adjust_method = 'BH', rownames = NULL, adj_pval_threshold = 0.05, 
                                    logFC_threshold = 1, legend_metadata_string = NULL,
                                    volc_plot_title = NULL, volc_plot_subtitle = NULL,
                                    volc_plot_xlab = NULL, volc_plot_ylab = NULL, DoVenn = FALSE) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'limma', 'ggplot2', 'ggrepel', 'stringr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Filter assay data by the conditional statement if present
  if(!is.null(metadata_condition)) {
    metab_SE <- metab_SE[, metadata_condition]
  }
  
  # Prepare limma data.frame
  metab_limma_data <- data.frame(assay(metab_SE))
  
  # Use shortname as rownames unless rownames provided
  if (is.null(rownames)) {
    rownames(metab_limma_data) <- make.unique(rowData(metab_SE)$shortname)
  } else {
    rownames(metab_limma_data) <- rownames
  }
  
  # Ensure that NA values are actually NA, and not character type 'NA'
  if (!is.null(metadata_condition)) {
    test_var <- metab_SE@metadata$metadata[[metadata_var]][metadata_condition]
  } else {
    test_var <- as.vector(metab_SE@metadata$metadata[[metadata_var]])
  }
  
  ensure_NA <- function(vector) {
    v <- replace(vector, vector == 'NA', NA)
    v <- replace(v, v == 'TRUE', 'Yes')
    v <- replace(v, v == 'FALSE', 'No')
    return(v)
  }
  
  test_var <- factor(ensure_NA(test_var))
  
  # Filter limma data to remove NA values
  metab_limma_data <- metab_limma_data[,!is.na(test_var)]
  
  # Create design matrix if none provided
  if (is.null(model_matrix)) {
    metab_limma_design <- model.matrix(~ 0 + test_var)
    colnames(metab_limma_design) <- levels(test_var)
  } else {
    metab_limma_design <- model_matrix
  }
  
  # Fit the expression matrix to a linear model
  fit <- lmFit(metab_limma_data, metab_limma_design)
  
  # Define function to handle creation of the contrasts matrix
  make_contrasts_vector <- function(design_matrix) {
    # Define levels and make new lists
    levels <- colnames(design_matrix)
    num_levels <- length(levels)
    ci <- list()
    cj <- list()
    
    # For loops for generate contrast statements
    for (i in 2:num_levels) {
      ci[length(ci) + 1] <- paste0(levels[1], '-', levels[i])
      for (j in (i + 1):num_levels) {
        cj[length(cj) + 1] <- paste0(levels[i], '-', levels[j])
      }
    }
    
    # Unlist elements and remove NA values, then remove the last element (contrasts itself)
    cx <- c(unlist(ci), unlist(cj)); cx <- cx[!str_detect(cx, 'NA')]
    cx <- cx[1:length(cx) - 1]
    
    # Generate contrasts matrix and return it
    cont_matrix <- makeContrasts(contrasts = cx, levels = levels)
    cont_matrix
  }
  
  # Create the contrasts matrix if none provided
  if (is.null(contrast_matrix)) {
    cont_matrix <- make_contrasts_vector(metab_limma_design)
  } else {
    cont_matrix <- contrast_matrix
  }
  
  # Bayes statistics of differential expression
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2, robust = TRUE, trend = FALSE)
  
  # Generate a list of the differentially abundant metabolites
  get_all_topTables <- function(fit2, top = TRUE) {
    # Prepare an empty list object
    all_topTables <- list()
    
    # Run through each of the coefficients and add topTable to the list
    for (coef in 1:dim(fit2$contrasts)[2]) {
      coef_name <- colnames(fit2$contrasts)[coef]
      metab_topTable <- topTable(fit2, number = dim(fit2)[1], adjust.method = adjust_method, coef = coef)
      if (top == TRUE) {
        metab_topTable <- metab_topTable %>%
          filter(adj.P.Val < adj_pval_threshold) %>%
          filter(abs(logFC) >= logFC_threshold)
      }
      all_topTables[[coef_name]] <- metab_topTable
    }
    
    # Return list object
    all_topTables
  }
  
  # Get both significant results and all results lists
  metab_limma_signif <- get_all_topTables(fit2, top = TRUE)
  metab_limma_all <- get_all_topTables(fit2, top = FALSE)
  
  # Generate a venn diagram of the results
  results <- decideTests(fit2)
  if (DoVenn) {
    venn_diagram <- vennDiagram(results)
  } else {
    venn_diagram <- "DoVenn = FALSE"
  }
  
  
  # Add a direction column to the 'metab_limma_all' topTables
  limma_add_categorical_direction <- function(metab_limma_all) {
    for (i in 1:length(metab_limma_all)) {
      metab_limma_all[[i]] <- data.frame(metab_limma_all[[i]]) %>%
        mutate(direction = case_when(
          adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ 'Decreased',
          adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ 'Increased',
          adj.P.Val >= adj_pval_threshold | abs(logFC) < logFC_threshold ~ 'NS'
        ))
    }
    metab_limma_all
  }
  
  metab_limma_all <- limma_add_categorical_direction(metab_limma_all)
  
  # Assign test_string if provided
  if (!is.null(legend_metadata_string)) {
    test_string <- legend_metadata_string
  } else {
    test_string <- 'Test Variable'
  }
  
  # Create a vector of strings for color legend
  plot_color_names <- function(test_string, metab_limma_all) {
    vector = ''
    for (i in 1:length(metab_limma_all)) {
      contrast_name <- names(metab_limma_all[i])
      contrast_name <- gsub('-', ' vs. ', contrast_name)
      contrast_name <- paste0(test_string, ':\n', contrast_name)
      vector <- c(vector, contrast_name)
    }
    vector <- vector[2:length(vector)]
    vector
  }
  
  plot_color_labels <- plot_color_names(test_string, metab_limma_all)
  
  # Define function to plot multiple volcano plots
  limma_volcano_plots <- function(topTable, plot_title, plot_subtitle, plot_xlab, plot_ylab, plot_color_labels) {
    # Create a blank list to hold plots
    plots <- list()
    
    # Make volcano plot for each topTable
    for (i in 1:length(metab_limma_all)) {
      plot <- ggplot(metab_limma_all[[i]], aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = direction)) +
        geom_hline(yintercept = -log10(adj_pval_threshold), linetype = 2) +
        geom_vline(xintercept = -logFC_threshold, linetype = 2) +
        geom_vline(xintercept = logFC_threshold, linetype = 2) +
        scale_color_manual(values = c('Decreased' = 'blue', 'Increased' = 'red', 'NS' = 'grey70'), name = plot_color_labels[i]) +
        geom_text_repel(data = metab_limma_all[[i]][metab_limma_all[[i]]$adj.P.Val < adj_pval_threshold & 
                                                      abs(metab_limma_all[[i]]$logFC) >= logFC_threshold, ],
                        label = rownames(metab_limma_all[[i]][metab_limma_all[[i]]$adj.P.Val < adj_pval_threshold & 
                                                                abs(metab_limma_all[[i]]$logFC) >= logFC_threshold, ]),
                        size = 2) +
        labs(title = plot_title,
             subtitle = plot_subtitle,
             x = plot_xlab,
             y = plot_ylab)
      
      # Add plot to list
      var_name <- names(metab_limma_all[i])
      plots[[var_name]] <- plot
    }
    
    # Return plots
    plots
  }
  
  # Prepare labelling
  if (is.null(volc_plot_title)) {
    plot_title <- 'Differential Intensity Metabolites'
  } else {
    plot_title <- volc_plot_title
  }
  
  if (is.null(volc_plot_subtitle)) {
    plot_subtitle <- NULL
  } else {
    plot_subtitle <- volc_plot_subtitle
  }
  
  if (is.null(volc_plot_xlab)) {
    plot_xlab <- 'log2FC'
  } else {
    plot_xlab <- volc_plot_xlab
  }
  
  if (is.null(volc_plot_ylab)) {
    plot_ylab <- 'Significance\n-log10(adjusted p-value)'
  } else {
    plot_ylab <- volc_plot_ylab
  }
  
  # Plot all volcano plots
  volcano_plots <- limma_volcano_plots(metab_limma_all, plot_title, plot_subtitle, plot_xlab, plot_ylab, plot_color_labels)
  
  # Make a list of all components above to return from function
  return_list <- list(input_data = metab_limma_data,
                      input_metadata = metab_SE@metadata$metadata[colnames(metab_limma_data),],
                      test_variable = test_var,
                      model_matrix = metab_limma_design,
                      contrast_matrix = cont_matrix,
                      limma_significant = metab_limma_signif,
                      limma_all = metab_limma_all,
                      volcano_plots = volcano_plots,
                      venn_diagram = venn_diagram,
                      limma_type = 'categorical')
  
  # Return the list
  return_list
}

# Metabolomics/Lipidomics
ag.fun <- function(x) {
  if (is.numeric(x))
    max(x)
  else if (is.character(x))
    cbind(paste0(x, collapse = ", "))
}

ag.fun <- function(x) {
  if (is.numeric(x))
    cbind(x[1]) #max(x)
  else if (is.character(x))
    cbind(x[1])} # keep first one only

mean_bind <- function(x) {
  if (is.numeric(x))
    mean(x)
  else if (is.character(x) | is.factor(x))
    cbind(x[1])
}

# Merge
dmetab.agg <- dmetab.agg %>%
  aggregate(., by = list(.$shortname), FUN = ag.fun) %>%
  column_to_rownames(., var = "Group.1") 

dmetab <- dmetab.agg %>%
  dplyr::select(rownames(met.metadata))

dmetab.info <- dmetab.agg %>%
  dplyr::select(colnames(dmetab.info))


