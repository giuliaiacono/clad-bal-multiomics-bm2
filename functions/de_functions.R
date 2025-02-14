
# DE functions
## Deseq2 function
deseq2 <- function(countdata, metadata, genes_info, test_variable, 
                   design, control_contrast = NULL, treat_contrast = NULL,
                   logFC_threshold, adj_pval_threshold, top_pval_threshold,
                   categorical = FALSE, continuous = FALSE){
  
  de_list <- list()
  de_list[["Design"]] <- design
  de_list[["Variable"]] <- test_variable
  de_list[["logFC_threshold"]] <- as.numeric(logFC_threshold)
  
  # Create dds object
  des_dataset <- DESeq2::DESeqDataSetFromMatrix(countData = countdata, 
                                                colData = metadata, 
                                                rowData = genes_info,
                                                design = as.formula(design)) 
  rownames(des_dataset) <- rowData(des_dataset)$Gene.Name
  de_list[["des_dataset"]] <- des_dataset
  
  # Run Deseq2
  des_run <- DESeq(des_dataset, minReplicatesForReplace = Inf) 
  de_list[["des_run"]] <- des_run
  
  # CATEGORICAL
  if (isTRUE(categorical)){
    #de_results <- DESeq2::results(des_run, contrast = c(as.character(test_variable), as.character(treat_contrast), as.character(control_contrast)), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = qf(.99, as.numeric(length(resultsNames(des_run))), as.numeric(nrow(metadata)) - as.numeric(length(resultsNames(des_run)))))
    
    # Apply fold change shrinkage
    de_results <- lfcShrink(des_run, coef = paste(as.character(test_variable), as.character(treat_contrast), "vs", as.character(control_contrast), sep = "_"), type = "apeglm")}
  
  # CONTINUOUS
  else if (isTRUE(continuous)){
    de_results <- DESeq2::results(des_run, alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = qf(.99, as.numeric(length(resultsNames(des_run))), as.numeric(nrow(metadata)) - as.numeric(length(resultsNames(des_run)))) )}
  
  # Check results
  summary(de_results)
  print(data.frame(cbind("max" = max(de_results$log2FoldChange), 
                         "min" = min(de_results$log2FoldChange),
                         "hits_below" = data.frame(de_results) %>% filter(padj < 0.05 & abs(log2FoldChange) > logFC_threshold) %>% count() %>% pull(n))))
  
  # Add extra gene info
  de_results$SYMBOL <- rowData(des_run)$Gene.Name
  de_results$ENSEMBL <- rowData(des_run)$Gene.ID
  de_results$GENENAME <- mapIds(org.Hs.eg.db, keys = de_results$ENSEMBL, keytype = "ENSEMBL",  column = "GENENAME", multiVals="first")
  de_list[["de_results"]] <- de_results
  
  # Extract results table
  de_results_df <- data.frame(de_results)
  de_results_df <- de_results_df %>%
    mutate(Direction = case_when(
      padj < adj_pval_threshold & log2FoldChange <= -logFC_threshold ~ 'Decreased',
      padj < adj_pval_threshold & log2FoldChange >= logFC_threshold ~ 'Increased',
      padj >= adj_pval_threshold | abs(log2FoldChange) < logFC_threshold ~ 'NS',
      is.na(de_results_df$padj) == TRUE | is.na(de_results_df$pvalue) == TRUE ~ "NA")) %>% 
    filter(Direction != "NA")
  de_list[["de_results_df"]] <- de_results_df
  
  # Show if no significant genes detected
  if (isTRUE(all(de_results_df$padj >= 0.05))){
    print("No DE hits - all padj values above 0.05")
    break}
  
  # Get sign hits
  sig_hits_df <- de_results_df %>% filter(Direction != "NS") 
  de_list[["sig_hits"]] <- sig_hits_df
  
  # Get top significant genes based on logFC expression
  top_hits_df <- de_results_df %>% filter(Direction != "NS") %>% filter(log2FoldChange > abs(logFC_threshold) & padj < top_pval_threshold)
  de_list[["top_hits"]] <- top_hits_df
  
  # Extract normalised dataset for visualisation
  vsd <- assay(vst(des_run, blind = FALSE))
  #vsd <- data.frame(normalizeVSN(assay(des_dataset)), check.names = FALSE)
  de_list[["normalised_counts"]] <- vsd
  
  de_list[["heat_all"]] <- t(scale(t(vsd))) #need samples in columns
  
  # Subset for any significant genes
  da_genes_int <- vsd[rownames(vsd) %in% rownames(sig_hits_df),]
  de_list[["heat_sig_hits"]] <- t(scale(t(da_genes_int)))
  
  # Subset da top genes with logfc > thr and adjpval < thr
  da_top_genes_int <- vsd[rownames(vsd) %in% rownames(top_hits_df),]
  de_list[["heat_top_hits"]] <- t(scale(t(da_top_genes_int)))
  
  return(de_list)}


## Limma process functions
### Met
# IN USE
limma_met_process <- function(list, data, metadata, 
                              data_info, coef_name, test_var1,
                              title, efit,  
                              hmdb_subset, folder_id, 
                              logFC_threshold, adj_pval_threshold = 0.05, top_pval_threshold = 0.01,
                              fella, plimit, threshold){
  
  list[[title]][[coef_name]]$metadata <- metadata
  list[[title]][[coef_name]]$efit <- efit
  
  de_results <- topTable(efit, coef = as.character(coef_name), number = Inf, adjust = 'BH') %>%
    mutate(direction = case_when(
      adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ as.character(levels(test_var1)[1]), # control
      adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ as.character(levels(test_var1)[2]), # treatment
      adj.P.Val >= adj_pval_threshold | abs(logFC) < logFC_threshold ~ 'n.s.',
      is.na(adj.P.Val) == TRUE | is.na(P.Value) == TRUE ~ "NA")) %>% 
    filter(direction != "NA")
  list[[title]][[coef_name]]$de_results <- as.data.frame(de_results)
  
  sig_hits <- de_results %>%
    filter(adj.P.Val < adj_pval_threshold) %>%
    filter(abs(logFC) >= logFC_threshold) %>% 
    rownames_to_column("shortname") %>% 
    left_join(data_info[, c("shortname", "KEGG")])
  list[[title]][[coef_name]]$sig_hits <- sig_hits 
  
  top_hits <- de_results %>% 
    filter(abs(logFC) >= logFC_threshold & adj.P.Val < top_pval_threshold) %>%
    arrange(logFC) %>% 
    rownames_to_column("shortname") %>%  
    left_join(data_info[, c("shortname", "KEGG")])
  list[[title]][[coef_name]]$top_hits <- top_hits
  
  list[[title]][[coef_name]]$heat_all <- t(scale(t(data)))
  list[[title]][[coef_name]]$heat_sig_hits <- t(scale(t(as.data.frame(data[rownames(data) %in% sig_hits$shortname, ]))))
  list[[title]][[coef_name]]$heat_top_hits <- t(scale(t(as.data.frame(data[rownames(data) %in% top_hits$shortname, ]))))
  
  # Pathway analysis
  path.info <- data_info %>% filter(rownames(.) %in% sig_hits$shortname) %>% dplyr::select(`info.Alignment ID`, `info.Average Rt(min)`, `info.Average Mz`, `info.Metabolite name`, `info.Adduct type`, `info.Fill %`, HMDB, "accession" = HMDB_Accession, KEGG, shortname) %>%
    left_join(hmdb_subset, by = "accession")
  
  list[[title]][[coef_name]]$path.info <- path.info
  list[[title]][[coef_name]]$path.ids <- path.info %>% dplyr::select("HMDB_ID" = accession) %>% distinct() %>% pull(HMDB_ID)
  
  write.csv(path.info, paste0(folder_id, Sys.Date(), "-", title, "-",  coef_name, "-path.csv"))
  write.csv(sig_hits, paste0(folder_id, Sys.Date(), "-", title, "-",  coef_name,  "-sig-hits-table.csv"))
  
  met.ontology <- sig_hits %>%
    left_join(path.info, by = "shortname") 
  rownames(met.ontology) <- met.ontology$shortname
  met.ontology <- met.ontology[rownames(list[[title]][[coef_name]]$heat_sig_hits),]
  list[[title]][[coef_name]]$ontology <- met.ontology
  
  if (any(str_detect(path.info$KEGG, "C")) & !is.null(fella)){
    
    # Introduce compounds
    s.path <- defineCompounds(compounds = path.info$KEGG, data = fella)
    s.path.analysis <- enrich(compounds = s.path@userinput@metabolites, data = fella, methods =  listMethods(), approx = "normality")
    
    # Generate results table and filter
    s.hy.table <- generateResultsTable(object = s.path.analysis, method = "hypergeom", data = fella, plimit = plimit, threshold = threshold)
    
    if (!is.null(s.hy.table)){
      
      colnames(s.hy.table) <- c("Pathway", "Pathway_name", "CompoundHits", "CompoundsInPathway", "pvalue")
      
      # Generate hypergeom graph and merge with table so compounds are matched to pathways
      s.path.hy <- generateResultsGraph(object = s.path.analysis, method = "hypergeom", data = fella, plimit = plimit, threshold = threshold)
      st.hy <- data.frame(get.edgelist(s.path.hy))
      colnames(st.hy) <- c("name", "Pathway")
      st.hy <- st.hy %>% inner_join(data.frame(vertex_attr(s.path.hy)), by = "name")
      
      s.hy.table <- st.hy %>%
        full_join(s.hy.table, by = "Pathway") %>%
        na.omit() %>%
        dplyr::select("Compound" = name, Pathway, "Metabolite" = label, Pathway_name, CompoundHits, CompoundsInPathway, pvalue)
      s.hy.table$Pathway_name <- gsub(" - Homo sapiens (human)", "", s.hy.table$Pathway_name, fixed = TRUE)
      
      list[[title]][[coef_name]]$go_res <- s.hy.table
    }
  }
  
  return(list)
  
}


### Lip
# IN USE (previously limma_lip_process)
limma_molecules_process <- function(list, data, metadata, coef_name, test_var1,
                              title, efit, folder_id, logFC_threshold, adj_pval_threshold = 0.05, top_pval_threshold = 0.01){
  
  list[[title]]$efit <- efit
  
  for (i in 1:length(coef_name)) {
    
    coef <- coef_name[i]
    
    de_results <- topTable(efit, coef = as.character(coef), number = Inf, adjust = 'BH') %>%
      mutate(direction = case_when(
        adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ as.character(levels(test_var1)[1]), # control
        adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ as.character(levels(test_var1)[2]), # treatment
        adj.P.Val >= adj_pval_threshold | abs(logFC) < logFC_threshold ~ 'n.s.',
        is.na(adj.P.Val) == TRUE | is.na(P.Value) == TRUE ~ "NA")) %>%
      rownames_to_column("Molecule") 
    list[[title]][[coef]]$de_results <- as.data.frame(de_results)
    
    sig_hits <- de_results %>%
      filter(adj.P.Val < adj_pval_threshold) %>%
      filter(abs(logFC) >= logFC_threshold)
    list[[title]][[coef]]$sig_hits <- sig_hits
    
    top_hits <- de_results %>% 
      filter(abs(logFC) >= logFC_threshold & adj.P.Val < top_pval_threshold) 
    list[[title]][[coef]]$top_hits <- top_hits
    
    list[[title]][[coef]]$heat_all <- t(scale(t(data)))
    list[[title]][[coef]]$heat_sig_hits <- t(scale(t(as.data.frame(data[rownames(data) %in% sig_hits$Molecule, ]))))
    list[[title]][[coef]]$heat_top_hits <- t(scale(t(as.data.frame(data[rownames(data) %in% top_hits$Molecule, ]))))
    
    write.csv(sig_hits, paste0(folder_id, title, "-",  coef, "-sig-hits-table.csv"))
    
  }
  
  return(list)

}

molecules_pathway_analysis <- function(list, data_info, title, molecules, folder_id, fella) {
  
  molecules_info <- data_info %>% 
  filter(shortname %in% molecules) %>% 
  dplyr::select("Molecule" = shortname, CATEGORY, MAIN_CLASS, SUB_CLASS, info.Ontology, ABBREVIATION, SYNONYMS, HMDB_Accession, KEGG, LIPIDMAPS_ID)
  list[[title]]$molecules_info <- molecules_info
  write.csv(molecules_info, paste0(folder_id, title, "-path.csv"), row.names = FALSE)
  
  if (any(str_detect(molecules_info$KEGG, "C"))){
    
    # Introduce compounds
    compounds <- molecules_info %>% pull(KEGG) %>% gsub(";.*", "", .) %>% na.omit()
    s.path <- defineCompounds(compounds = molecules_info$KEGG, data = fella)
    s.path.analysis <- enrich(compounds = s.path@userinput@metabolites, data = fella, methods =  listMethods(), approx = "normality")
    
    # Generate results table and filter
    s.hy.table <- generateResultsTable(object = s.path.analysis, method = "hypergeom", data = fella, plimit = 49, threshold = 1)
    
    if (!is.null(s.hy.table)){
      
      colnames(s.hy.table) <- c("Pathway", "Pathway_name", "CompoundHits", "CompoundsInPathway", "pvalue")
      
      # Generate hypergeom graph and merge with table so compounds are matched to pathways
      s.path.hy <- generateResultsGraph(object = s.path.analysis, method = "hypergeom", data = fella, plimit = 49, threshold = 1)
      st.hy <- data.frame(get.edgelist(s.path.hy))
      colnames(st.hy) <- c("name", "Pathway")
      st.hy <- st.hy %>% inner_join(data.frame(vertex_attr(s.path.hy)), by = "name")
      
      s.hy.table <- st.hy %>%
        full_join(s.hy.table, by = "Pathway") %>%
        na.omit() %>%
        dplyr::select("KEGG" = name, Pathway, "Molecule" = label, Pathway_name, CompoundHits, CompoundsInPathway, pvalue)
      s.hy.table$Pathway_name <- gsub(" - Homo sapiens (human)", "", s.hy.table$Pathway_name, fixed = TRUE)
      
      list[[title]]$fella_output <- s.hy.table
    }
  }
  
  return(list)
  
}



### Mic


# IN USE
limma_mic_process <- function(list, data, metadata, coef_name, test_var1,
                              title, efit, folder_id, 
                              logFC_threshold, adj_pval_threshold = 0.05, top_pval_threshold = 0.01){
  
  list[[title]]$efit <- efit
  
  de_results <- topTable(efit, coef = as.character(coef_name), number = Inf, adjust = 'BH') %>%
    mutate(direction = case_when(
      adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ as.character(levels(test_var1)[1]),
      adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ as.character(levels(test_var1)[2]),
      adj.P.Val >= adj_pval_threshold | abs(logFC) < logFC_threshold ~ 'n.s.',
      is.na(adj.P.Val) == TRUE | is.na(P.Value) == TRUE ~ "NA")) %>% 
    filter(direction != "NA")
  list[[title]]$de_results <- as.data.frame(de_results)
  
  sig_hits <- de_results %>%
    filter(adj.P.Val < adj_pval_threshold) %>%
    filter(abs(logFC) >= logFC_threshold)
  list[[title]]$sig_hits <- sig_hits
  
  top_hits <- de_results %>% 
    filter(abs(logFC) >= logFC_threshold & adj.P.Val < top_pval_threshold) %>%
    arrange(logFC)
  list[[title]]$top_hits <- top_hits
  
  tbi.da <- data %>% psmelt() %>% filter(OTU %in% rownames(list[[paste0(title)]]$sig_hits))
  tbi.da.heat <- reshape2::dcast(tbi.da, Sample ~ OTU, value.var = "Abundance") %>% 
    column_to_rownames(., var = "Sample")
  tbi.da.heat <- t(tbi.da.heat[rownames(metadata),])
  list[[paste0(title)]]$heat_sig_hits <- tbi.da.heat
  #list[[title]]$heat_top_hits <- tbi.da.heat[rownames(list[[paste0(title)]]$top_hits), ]
  
  write.csv(sig_hits, paste0(folder_id, "mic/", Sys.Date(), "-", title, "-sig-hits-table.csv"))
  
  return(list)
  
}



### RNA


## RNA limma process
limma_rna_process <- function(list, folder_id, test_variable, countdata, coef_name, logFC_threshold, test_var, pval_threshold = 0.05, adj_pval_threshold = 0.05, top_pval_threshold = 0.01, cpm = FALSE){
  
  list[[test_variable]]$efit <- efit
  de_results <- topTable(efit, coef = coef_name, number = Inf, adjust = 'BH') %>%
    mutate(direction = case_when(
      P.Value < pval_threshold & adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ as.character(levels(test_var)[1]),
      P.Value < pval_threshold & adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ as.character(levels(test_var)[2]),
      P.Value >= pval_threshold | adj.P.Val >= adj_pval_threshold | abs(logFC) < logFC_threshold ~ 'n.s.',
      is.na(adj.P.Val) == TRUE | is.na(P.Value) == TRUE ~ "NA")) %>% 
    filter(direction != "NA") %>%
    mutate(SYMBOL = rownames(.))
  list[[test_variable]]$de_results <- as.data.frame(de_results)
  
  sig_hits <- de_results %>%
    filter(P.Value < pval_threshold, adj.P.Val < adj_pval_threshold, abs(logFC) >= logFC_threshold) 
  list[[test_variable]]$sig_hits <- sig_hits
  
  top_hits <- de_results %>% 
    filter(P.Value < pval_threshold, adj.P.Val < top_pval_threshold, abs(logFC) >= logFC_threshold) 
  list[[test_variable]]$top_hits <- top_hits
  
  if (isTRUE(cpm)){
    heat <- as.data.frame(cpm(countdata, log = TRUE))
  }
  
  if (isFALSE(cpm)){
    heat <- countdata 
  }
  
  list[[test_variable]]$vsd <- heat
  list[[test_variable]]$heat_all <- t(scale(t(heat)))
  list[[test_variable]]$heat_sig_hits <- t(scale(t(as.data.frame(heat[rownames(heat) %in% rownames(sig_hits), ]))))
  list[[test_variable]]$heat_top_hits <- t(scale(t(as.data.frame(heat[rownames(heat) %in% rownames(top_hits), ]))))
  
  #write.csv(list[[test_variable]]$de_results, paste0(folder_id, "all_hits.csv"))
  #write.csv(list[[test_variable]]$sig_hits, paste0(folder_id, "sig_hits.csv"))
  #write.csv(list[[test_variable]]$top_hits, paste0(folder_id, "top_hits.csv"))
  
  return(list)
}




### Heatmaps


# HEATMAP
limma_heat_plots <- function(list, title, metadata, annotation, folder_id, order = NULL,
                             row_title, fontsize = 10, width, height, res = 600, palette = NULL){
  
  ## HEAT
  if (row_title == "Metabolite"){
    
    if (isTRUE(order)){
      
      png(paste0(folder_id, "met/", Sys.Date(), "-", title, "-heat-me.png"), width = width, height = height, res = res)
      
      draw(Heatmap(list[[paste0(title)]]$heat_sig_hits, 
                   row_title = row_title, 
                   heatmap_legend_param = list(title = "Intensity"),
                   show_column_names = FALSE,
                   column_order = colnames(list[[paste0(title)]]$heat_sig_hits[, metadata %>% group_by(Timepoint, AgeGroup) %>% arrange(AgeGroup, .by_group = TRUE) %>% pull(SeqID)]),
                   top_annotation = annotate_tbi(metadata, annotations = annotation),
                   row_names_gp = gpar(fontsize = fontsize)), 
           padding = unit(c(1, 1, 1, 1.5), "cm"))
      
      invisible(dev.off()) 
      
      p <- grid.grabExpr(draw(Heatmap(list[[paste0(title)]]$heat_sig_hits, 
                                      row_title = row_title, 
                                      heatmap_legend_param = list(title = "Intensity"),
                                      show_column_names = FALSE,
                                      column_order = colnames(list[[paste0(title)]]$heat_sig_hits[, metadata %>% group_by(Timepoint, AgeGroup) %>% arrange(AgeGroup, .by_group = TRUE) %>% pull(SeqID)]),
                                      top_annotation = annotate_tbi(metadata, annotations = annotation),
                                      row_names_gp = gpar(fontsize = fontsize)), 
                              padding = unit(c(1, 1, 1, 1.5), "cm")))
      
    }
    
    if (is.null(order)){
      
      png(paste0(folder_id, "met/", Sys.Date(), "-", title, "-heat-me.png"), width = width, height = height, res = res)
      
      draw(Heatmap(list[[paste0(title)]]$heat_sig_hits, 
                   row_title = row_title, 
                   heatmap_legend_param = list(title = "Intensity"),
                   show_column_names = FALSE,
                   top_annotation = annotate_tbi(metadata, annotations = annotation),
                   row_names_gp = gpar(fontsize = fontsize)), 
           padding = unit(c(1, 1, 1, 1.5), "cm"))
      
      invisible(dev.off()) 
      
      p <- grid.grabExpr(draw(Heatmap(list[[paste0(title)]]$heat_sig_hits, 
                                      row_title = row_title, 
                                      heatmap_legend_param = list(title = "Intensity"),
                                      show_column_names = FALSE,
                                      top_annotation = annotate_tbi(metadata, annotations = annotation),
                                      row_names_gp = gpar(fontsize = fontsize)), 
                              padding = unit(c(1, 1, 1, 1.5), "cm")))
    }
    
    return(p)
    
  }
  
  if (row_title == "ASV"){
    
    if (isTRUE(order)){
      
      png(paste0(folder_id, "mic/", Sys.Date(), "-", title, "-heat-mi.png"),  width = width, height = height, res = res)
      
      draw(Heatmap(list[[paste0(title)]]$heat_sig_hits, 
                   row_title = row_title, 
                   heatmap_legend_param = list(title = "Abundance"),
                   show_column_names = FALSE,
                   column_order = colnames(list[[paste0(title)]]$heat_sig_hits[, metadata %>% group_by(Timepoint, AgeGroup) %>% arrange(AgeGroup, .by_group = TRUE) %>% pull(SeqID)]),
                   top_annotation = annotate_tbi(metadata, annotations = annotation),
                   col = palette,
                   row_names_gp = gpar(fontsize = fontsize)), 
           padding = unit(c(1, 1, 1, 1.5), "cm"))
      
      invisible(dev.off()) 
      
      p <- grid.grabExpr(draw(Heatmap(list[[paste0(title)]]$heat_sig_hits, 
                                      row_title = row_title, 
                                      heatmap_legend_param = list(title = "Abundance"),
                                      show_column_names = FALSE,
                                      column_order = colnames(list[[paste0(title)]]$heat_sig_hits[, metadata %>% group_by(Timepoint, AgeGroup) %>% arrange(AgeGroup, .by_group = TRUE) %>% pull(SeqID)]),
                                      top_annotation = annotate_tbi(metadata, annotations = annotation),
                                      col = palette,
                                      row_names_gp = gpar(fontsize = fontsize)), 
                              padding = unit(c(1, 1, 1, 1.5), "cm")))
    }
    
    if (is.null(order)){
      
      png(paste0(folder_id, "mic/", Sys.Date(), "-", title, "-heat-mi.png"),  width = width, height = height, res = res)
      
      draw(Heatmap(list[[paste0(title)]]$heat_sig_hits, 
                   row_title = row_title, 
                   heatmap_legend_param = list(title = "Abundance"),
                   show_column_names = FALSE,
                   top_annotation = annotate_tbi(metadata, annotations = annotation),
                   col = palette,
                   row_names_gp = gpar(fontsize = fontsize)), 
           padding = unit(c(1, 1, 1, 1.5), "cm"))
      
      invisible(dev.off()) 
      
      p <- grid.grabExpr(draw(Heatmap(list[[paste0(title)]]$heat_sig_hits, 
                                      row_title = row_title, 
                                      heatmap_legend_param = list(title = "Abundance"),
                                      show_column_names = FALSE,
                                      top_annotation = annotate_tbi(metadata, annotations = annotation),
                                      col = palette,
                                      row_names_gp = gpar(fontsize = fontsize)), 
                              padding = unit(c(1, 1, 1, 1.5), "cm")))
    }
    
    return(p)
  }
  
}



### Hits


# HITS
limma_hits_plot <- function(list, data, metadata, title, type, folder_id,
                            x = x){
  
  if (type == "Metabolite"){
    
    top_hits <- list[[paste0(title)]][["sig_hits"]] %>% arrange(desc(logFC)) %>% rownames_to_column(var = "Metabolite")
    hits_figure <- list()
    
    for (i in 1:dim(top_hits)[1]){
      
      var_select <- top_hits$Metabolite[i] #i
      tab <- data[var_select,] %>%
        t() %>%
        data.frame() %>%
        rownames_to_column("SeqID") 
      colnames(tab) <- c("SeqID", "Intensity")
      tab <- tab %>%
        inner_join(metadata, by = "SeqID")
      
      p <-  ggplot(tab, aes(x = {{x}}, y = Intensity)) +
        geom_boxplot(aes(fill = AgeGroup), alpha = 0.2, width = 0.75, outlier.shape = NA) +
        geom_point(aes(colour = AgeGroup, shape = Sex), position = position_dodge(width=0.75), size = 1.5, stroke = 1) +
        labs(title = var_select) + 
        theme(plot.title = 6) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) +
        theme_light() +
        tbi_plot_fill(AgeGroup = TRUE) +
        tbi_plot_colours(AgeGroup = TRUE) +
        theme(legend.position = "bottom", legend.justification = c("left"))
      
      ggsave(paste0(folder_id, "met/", Sys.Date(), "-", paste0(gsub("/", "", var_select)), ".png"), p, width = 12, height = 8, units = 'cm')
      
      hits_figure[[var_select]] <- p}
    
    # Report 4 top hits in one figure
    top_hits <- rbind(head(top_hits, n = 2), tail(top_hits, n = 2)) %>% distinct()
    
    figure[[4]] <- ggarrange(hits_figure[[top_hits$Metabolite[1]]], hits_figure[[top_hits$Metabolite[2]]], hits_figure[[top_hits$Metabolite[3]]], hits_figure[[top_hits$Metabolite[4]]],
                             ncol = 2, nrow = 2, 
                             align = "v",
                             common.legend = TRUE, legend = "bottom")
    
    ggsave(paste0(folder_id, "met/", Sys.Date(), "-", title, "-hits-me.png"), figure[[4]], width = 22, height = 12, units = 'cm')
    
  }
  
  if (type == "ASV"){
    
    top_hits <- list[[paste0(title)]][["sig_hits"]] %>% arrange(desc(logFC)) %>% rownames_to_column(var = "ASV")
    hits_figure <- list()
    
    for (i in 1:dim(top_hits)[1]){
      
      var_select <- top_hits$ASV[i] #i
      
      tab <- otu_table(data)[var_select,] %>%
        t() %>%
        data.frame() %>%
        rownames_to_column("SeqID") 
      colnames(tab) <- c("SeqID", "Abundance")
      tab <- tab %>%
        inner_join(metadata, by = "SeqID")
      
      p <-  ggplot(tab, aes(x = {{x}}, y = Abundance)) +
        geom_boxplot(aes(fill = AgeGroup), alpha = 0.2, outlier.shape = NA) +
        geom_jitter(aes(colour = AgeGroup, shape = Sex), size = 1.5, stroke = 1) +
        labs(title = var_select) + 
        theme(plot.title = 6) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) +
        theme_light() +
        tbi_plot_fill(AgeGroup = TRUE) +
        tbi_plot_colours(AgeGroup = TRUE) +
        facet_wrap(~Treatment)
      theme(legend.position = "bottom", legend.justification = c("left"))
      
      ggsave(paste0(folder_id, "mic/", Sys.Date(), "-", paste0(gsub("/", "", var_select)), ".png"), p, width = 14, height = 8, units = 'cm')
      
      hits_figure[[var_select]] <- p}
    
    # Report 4 top hits in one figure
    top_hits <- rbind(head(top_hits, n = 2), tail(top_hits, n = 2)) %>% distinct()
    
    figure[[4]] <- ggarrange(hits_figure[[top_hits$ASV[1]]], hits_figure[[top_hits$ASV[2]]], hits_figure[[top_hits$ASV[3]]], hits_figure[[top_hits$ASV[4]]],
                             ncol = 2, nrow = 2, 
                             align = "v",
                             common.legend = TRUE, legend = "bottom")
    
    ggsave(paste0(folder_id, "mic/", Sys.Date(), "-", title, "-hits-mi.png"), figure[[4]], width = 22, height = 12, units = 'cm')
  }
  
}


### Analyse GO


limma_rna_go <- function(list, test_variable, folder_id, logFC_threshold, top = FALSE, sig = FALSE,
                         org.Hs.eg.db, genes_bg, toptable_rows = 500,
                         adj_pval_threshold = 0.05, top_pval_threshold = 0.01){
  
  gene.list <- list()
  
  if (isTRUE(sig)){
    dset <- list[[test_variable]]$sig_hits
  }
  
  if (isTRUE(top)){
    dset <- list[[test_variable]]$top_hits
  }
  
  for (i in 1:length(dset$SYMBOL)){
    gene <- geneinfo_2_html(dset$SYMBOL[i])
    gene <- sub(".*GeneCards database: <a href = ", "", gene)
    gene <- sub(" target.*", "", gene)
    gene <- gsub('"', '', gene)
    gene.list[[i]] <- gene
  }
  gene.df <- dplyr::bind_cols(gene.list) %>% t()
  gene.df <- cbind(gene.df, dset$SYMBOL, dset$GENENAME, dset$Group)
  write.csv(gene.df, paste0(folder_id, Sys.Date(), "-", test_variable, "-gene-links.csv"))
  
  ## GOANA
  # Annotation for each go
  anno.go <- AnnotationDbi::select(org.Hs.eg.db, keytype = "ENTREZID", keys = dset$ENTREZID, columns = "GO") %>% filter(ONTOLOGY == "BP") %>% left_join(dset, by = "ENTREZID") %>%
    dplyr::select(ENTREZID, GO, SYMBOL) %>% distinct()
  list[[test_variable]]$anno_goana <- anno.go
  goana <- goana(dset$ENTREZID, species="Hs", FDR = 0.01) %>% filter(Ont == "BP")
  list[[test_variable]]$goana_res <- goana
  
  # Merge annotations
  list[[test_variable]]$go_res <- goana %>%
    rownames_to_column(., var = "GO") %>%
    right_join(anno.go, by = "GO") %>%
    dplyr::select(GO, Term, N, DE, P.DE, SYMBOL) %>% 
    filter(P.DE < 0.05) %>%
    arrange(desc(DE))
  
  ## GENE TONIC
  # Vector of de gene names with fdr threshold < 0.01
  de_results_sub <- dset %>%
    mutate(log2FoldChange = logFC)
  de_symb <- dset$SYMBOL
  # Background genes
  bg_ids <- genes_bg$Gene.Name
  
  # Enrich
  suppressMessages(top_go_de <- pcaExplorer::topGOtable(de_symb,
                                                        bg_ids,
                                                        ontology = "BP",
                                                        mapping = "org.Hs.eg.db",
                                                        geneID = "symbol",
                                                        topTablerows = toptable_rows))
  
  # Convert for gene tonic
  res_enrich <- shake_topGOtableResult(top_go_de)
  
  # Annotation df
  anno_df <- data.frame(
    gene_id =  genes_bg$Geneid,
    gene_name =  genes_bg$Gene.Name,
    stringsAsFactors = FALSE)
  rownames(anno_df) <- anno_df$gene_id
  list[[test_variable]]$anno_df <- anno_df
  
  # Uses rownames, need to be the same as annotation object
  rownames(de_results_sub) <- de_results_sub$ENSEMBL
  
  # Aggregate (adds two columns)
  res_enrich <- get_aggrscores(res_enrich = res_enrich,
                               res_de = de_results_sub,
                               annotation_obj = anno_df,
                               aggrfun = mean)
  list[[test_variable]]$res_enrich <- res_enrich
  
  res_enrich_sig <- res_enrich %>% filter(gs_pvalue < 0.05 & gs_de_count > 2) %>% arrange(desc(gs_de_count))
  
  list[[test_variable]]$res_enrich_sig <- res_enrich_sig
  
  #GeneTonic(dds = des_dataset,
  #res_de = deseq2_results,
  #res_enrich = res_enrich,
  #annotation_obj = anno_df,
  #project_id = "RNA-seq")
  
  # de_list[["gt_object"]] <- GeneTonic(dds = de_list$dataset, res_de = de_list$de_results, res_enrich = res_enrich_sig, annotation_obj = anno_df, 
  #project_id = project_id)
  
  return(list)
  
  #go.wgcna <- enrichPathway(gene = rna.list[[test_variable]]$sig_hits$ENTREZID, pvalueCutoff = 1, qvalueCutoff = 1)
  #rna.list[[test_variable]]$go.wgcna.df <- data.frame(go.wgcna)
  #rna.list[[test_variable]]$netplot <- cnetplot(go.wgcna, categorySize = 'pvalue', showCategory = 6, color_gene = 'indianred', color_category = 'grey', layout = 'fr')
}



select_top_pathways <- function(net_hits, list, test_variable){
  
  #test_variable <- test_variable
  #net_hits <- net_hits
  #list <- rna.list[c(1,10)]
  
  # Add columns
  net_hits$description <-  as.character(paste0("[", abbreviate(net_hits$gs_description, 4, named = FALSE, method = "both", dot = FALSE), "]")) 
  net_hits$Colour_id <- colorRampPalette(c("#046C9A","#ECCBAE", "#692895", "#08b6b6", "#ff8a5b","#faa613"), space = 'rgb')(length(net_hits$gs_description))
  net_hits$log_gs_pvalue <- abs(log10(net_hits$gs_pvalue))
  
  # Split
  s <- strsplit(net_hits$gs_genes, split = ",")
  nedf <- data.frame("SYMBOL" = unlist(s),
                     "Pathway_description" = rep(net_hits$gs_description,sapply(s, length)),
                     "Colour_id" = rep(net_hits$Colour_id,sapply(s, length)),
                     "Pathway_id" = rep(net_hits$gs_id,sapply(s, length)), 
                     "de_count" = rep(net_hits$gs_de_count,sapply(s, length)), 
                     "z_score" = rep(net_hits$z_score,sapply(s, length)),
                     "gs_pvalue" = rep(net_hits$gs_pvalue,sapply(s, length)),
                     "log_gs_pvalue" = rep(net_hits$log_gs_pvalue,sapply(s, length)))
  
  ne <- nedf %>% 
    left_join(list[[test_variable]]$sig_hits[,c("SYMBOL", "direction")], by = "SYMBOL") %>%
    rename(Gene = SYMBOL)
  
  list[[test_variable]]$go_summary <- ne
  
  return(list)
}




