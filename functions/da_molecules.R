
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


