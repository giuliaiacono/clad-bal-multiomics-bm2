plot_fill <- function(Group = FALSE,
                      Transplant_indication = FALSE,
                      Time_interval_detailed = FALSE,
                      Time_interval_combined = FALSE,
                      Time_interval = FALSE,
                      Days_interval = FALSE,
                      Days_interval2 = FALSE,
                      Months_grouped = FALSE) {
  list(

    if (Group)
      scale_fill_manual(values = c("CLAD-free" = "#5C6378", "CLAD" = "#EA9010")),

    if (Transplant_indication)
      scale_fill_manual(values = c("CF/NCFB" = "#EA9010", "COPD" = "#150578", "ILD" = "#BEE5BF", "Other" = "#5C6378", "CLAD" = "#197278")),

    if (Time_interval_detailed)
      scale_fill_manual(values = c("<1.5 months" = "#EEEEEE", "1.5-3 months" = "#D9E1E8",
                                   "3-6 months" = "#B2CADC", "6-12 months" = "#82ADCE", ">12 months" = "#6099C4")),
    if (Time_interval)
      scale_fill_manual(values = c("<6 months" = "#D9E1E8", "6-12 months" = "#82ADCE", ">12 months" = "#78aafa")),

    if (Time_interval_combined)
      scale_fill_manual(values = c("<3 months" = "#EEEEEE", "3-6 months" = "#D9E1E8", "6-12 months" = "#82ADCE", ">12 months" = "#6099C4")),

    if (Days_interval)
      scale_fill_manual(values = c("<6 months" = "lightgrey", ">6 months" = "#3a86ff")),

    if (Days_interval2)
      scale_fill_manual(values = c("<3" = "lightgrey", ">3" = "#389ce8")),

    if (Months_grouped)
      scale_fill_manual(values = c("<1.5" = "#EEEEEE", "1.5-3" = "darkgrey", "4-6" = "#389ce8", "7-9" = "#266fa6", "10-12" = "#005da3", "13-18" = "#014b82", ">18" = "#01223b") )
  )
}

plot_colours <- function(Group = FALSE,
                      Transplant_indication = FALSE,
                      Time_interval_detailed = FALSE,
                      Time_interval_combined = FALSE,
                      Time_interval = FALSE,
                      Days_interval = FALSE,
                      Days_interval2 = FALSE,
                      Months_grouped = FALSE) {
  list(

    if (Group)
      scale_colour_manual(values = c("CLAD-free" = "#5C6378", "CLAD" = "#EA9010")),

    if (Transplant_indication)
      scale_colour_manual(values = c("CF/NCFB" = "#EA9010", "COPD" = "#150578", "ILD" = "#BEE5BF", "Other" = "#5C6378", "CLAD" = "#197278")),

    if (Time_interval_detailed)
      scale_colour_manual(values = c("<1.5 months" = "#EEEEEE", "1.5-3 months" = "#D9E1E8",
                                   "3-6 months" = "#B2CADC", "6-12 months" = "#82ADCE", ">12 months" = "#6099C4")),
    if (Time_interval)
      scale_colour_manual(values = c("<6 months" = "#D9E1E8", "6-12 months" = "#82ADCE", ">12 months" = "#78aafa")),

    if (Time_interval_combined)
      scale_colour_manual(values = c("<3 months" = "#EEEEEE", "3-6 months" = "#D9E1E8", "6-12 months" = "#82ADCE", ">12 months" = "#6099C4")),

    if (Days_interval)
      scale_colour_manual(values = c("<6 months" = "lightgrey", ">6 months" = "#3a86ff")),

    if (Days_interval2)
      scale_colour_manual(values = c("<3" = "lightgrey", ">3" = "#389ce8")),

    if (Months_grouped)
      scale_colour_manual(values = c("<1.5" = "#EEEEEE", "1.5-3" = "darkgrey", "4-6" = "#389ce8", "7-9" = "#266fa6", "10-12" = "#005da3", "13-18" = "#014b82", ">18" = "#01223b") )
  )
}

# Select the metadata to add as as column annotation and choose colors
annotate_heatmap <- function(metadata, annotations){

  #metadata <- bact.metadata
  #annotations <- c("Donor_age", "Code")

  # Select metadata columns
  list_annotation <- list()
  for (i in annotations){
    list_annotation[i] <- metadata %>% dplyr::select(i)}
  list_annotation_df <- data.frame(dplyr::bind_cols(list_annotation), check.names = FALSE)
  rownames(list_annotation_df) <- metadata$Patient_days

  if (any(grepl("_", annotations))){
    colnames(list_annotation_df)[which(grepl("_", colnames(list_annotation_df)))] <- gsub("_", " ", colnames(list_annotation_df))
    annotations[which(grepl("_", annotations))] <- gsub("_", " ", annotations)
  }

  # Select annotation colours
  colour_annotation <- list()

  if (isTRUE(any(annotations == "Group"))){
    colour_annotation[["Group"]] <- c("CLAD-free" = "#5C6378", "CLAD" = "#EA9010")
  }

  if (isTRUE(any(annotations == "Transplant_indication"))){
    colour_annotation[["Transplant_indication"]] <- c("CF/NCFB" = "#EA9010", "COPD" = "#150578", "ILD" = "#BEE5BF",
                                                      "Other" = "#5C6378", "CLAD" = "#197278")}


  if (isTRUE(any(annotations == "Status"))){
    colour_annotation[["Status"]] <- c("Before" = "#5C6378", "After" = "#3a86ff")
  }

  if (isTRUE(any(annotations == "Days interval"))){
    colour_annotation[["Days interval"]] <- c("<6 months" = "lightgrey", ">6 months" = "#3a86ff")
  }

  if (isTRUE(any(annotations == "Days"))){
    col_fun <- colorRamp2(c(min(list_annotation_df$Days), mean(list_annotation_df$Days), max(list_annotation_df$Days)), c("lightgrey", "darkgrey", "black"))
    colour_annotation[["Days"]] <- col_fun}

  if (isTRUE(any(annotations == "Months"))){
    col_fun <- colorRamp2(c(min(metadata$Months), mean(metadata$Months), max(metadata$Months)), c("lightgrey", "darkgrey", "black"))
    colour_annotation[["Months"]] <- col_fun}

  if (isTRUE(any(annotations == "Interval"))){
    colour_annotation[["Interval"]] <- c("<1.5 months" = "#EEEEEE", "1.5-3 months" = "#D9E1E8", "3-6 months" = "#B2CADC", "6-12 months" = "#82ADCE", ">12 months" = "#6099C4")
  }

  if (isTRUE(any(annotations == "Immunosuppression"))){
    colour_annotation[["Immunosuppression"]] <- c("High" = "#626262", "Mid" = "#A4A4A4", "Low" = "#E1E1E1")
  }

  if (isTRUE(any(annotations == "Time post-transplant (months)"))){
    colour_annotation[["Time post-transplant (months)"]] <- c("<1.5" = "#EEEEEE", "1.5-3" = "darkgrey", "4-6" = "#389ce8", "7-9" = "#266fa6",
                                                              "10-12" = "#005da3", "13-18" = "#014b82", ">18" = "#01223b")
  }

  if (isTRUE(any(annotations == "Cohort"))){
    colour_annotation[["Cohort"]] <- c("Stable" = "lightgrey", "Biomarker" = "#626262")
  }

  colAnn <- HeatmapAnnotation(df = list_annotation_df,
                              col = colour_annotation,
                              annotation_name_side = "right",
                              height = unit(0.1, "cm"),
                              annotation_legend_param = list(direction = "vertical")) #annotation_label = c(" ", " ", " ") title_position = "leftcenter-rot"

  return(colAnn)
}

annotate_heatmap_row <- function(metadata, row_id = row_id, row_or = NULL, annotations, path_palette = NULL, phylum_palette = NULL){

  # Select metadata columns
  list_annotation <- list()
  list_colour_annotation <- list()

  for (i in annotations){

    list_annotation[i] <- metadata %>% dplyr::select(all_of(i)) }

  annotation_df <- data.frame(dplyr::bind_cols(list_annotation))
  rownames(annotation_df) <- metadata %>% pull({{row_id}})
  if (!is.null(row_or)) {
    annotation_df <- annotation_df %>% arrange(match(rownames(.), row_or))
  }

  # Select annotation colours
  if (isTRUE(any(annotations == "Pathway"))){
    col_fun <- c(colorRampPalette(path_palette, space = 'rgb')(length(levels(metadata$Pathway))-1), "#EAEAEA")
    list_colour_annotation[["Pathway"]] <- data.frame(levels(metadata$Pathway), col_fun) %>% tibble::deframe()}

  if (isTRUE(any(annotations == "Class"))){
    list_colour_annotation[["Class"]] <- path_palette}

  if (isTRUE(any(annotations == "Cluster"))){
    col_fun <- c("#f78436", "#6099c4")
    list_colour_annotation[["Cluster"]] <- data.frame(levels(metadata$Cluster), col_fun) %>% tibble::deframe()}

  if (isTRUE(any(annotations == "Phylum"))){
    #col_fun <- c("#f26a8d", "#3F3D89", "#8339B8", "#CC6AD8", "#EA9010", "#1B8293", "#80ed99", "#38B2FB", "#5C6378")
    list_colour_annotation[["Phylum"]] <- phylum_palette}


  if (isTRUE(any(annotations %in% c("Tacrolimus", "Prednisolone", "Cyclosporine", "Mycophenolate", "Azathioprine")))){

    min <- corr_metadata %>% summarise_if(is.numeric, ~min) %>% min()
    max <- corr_metadata %>% summarise_if(is.numeric, ~max) %>% max()
    col_fun <- colorRamp2(c(min, 0, max), c("#6099c4" , "white", "#f78436"))

    for (i in 2:length(annotations)){
      list_colour_annotation[[annotations[i]]] <- col_fun}
  }

  rowAnn <- rowAnnotation(df = annotation_df,
                          #annotation_legend_param = list(Tacrolimus = list(title = "Pearson correlation")),
                          annotation_label = gpar(fontsize = 6),
                          #border = TRUE,
                          col = list_colour_annotation)

  return(rowAnn)
}
