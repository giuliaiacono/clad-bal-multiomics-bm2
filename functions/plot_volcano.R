plot_volcano <- function(list, folder_id, width, height){ 
  
  g <- ggplot(list[[test_variable]]$de_results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(colour = Group)) + 
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_vline(xintercept = -logFC_threshold, linetype = 2) +
    geom_vline(xintercept = logFC_threshold, linetype = 2) +
    labs(colour = "Upregulated in") +
    geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val)), max.overlaps = 25, label = ifelse(list[[test_variable]]$de_results$Group != 'n.s.', as.character(rownames(list[[test_variable]]$de_results)),''), cex=1.5) +
    plot_colours(Group = TRUE) +
    theme(legend.position = "bottom", legend.justification = c("left")) 
  
  ggsave(paste0(folder, Sys.Date(), "-vol-", test_variable, ".png"), g, width = width, height = height, units = 'cm')
  
  return(g)}