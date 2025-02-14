pca_process <- function(data, metadata, k_max = FALSE){
  
  pca_list <- list()
  
  # Determine best clustering
  k <- fviz_nbclust(t(data), kmeans, method = "silhouette")$data
  
  if (isFALSE(k_max)){
    k_max <- as.numeric(k$clusters[which.max(k$y)])}
  
  # Create clusters using partitioning around medoid
  pam <- pam(t(data), k_max, metric = "euclidean", stand = FALSE)
  pam.cluster <- t(data.frame(as.list(pam$clustering), check.names = FALSE))
  pam.cluster <- cbind(rownames(pam.cluster), pam.cluster)
  colnames(pam.cluster) <- c("Patient_days", "Cluster")
  
  PCA <- prcomp(t(data), center = TRUE, scale = TRUE)
  varexp <- c(summary(PCA)$importance[2,1]*100,summary(PCA)$importance[2,2]*100)
  
  data.PCA <- cbind(data.frame(Patient_days = rownames(PCA$x),
                               PC1 = PCA$x[,1],
                               PC2 = PCA$x[,2])) %>% 
    left_join(metadata, by = "Patient_days") %>% 
    arrange(Record.ID, Days) %>% 
    left_join(data.frame(pam.cluster), by = "Patient_days") %>%
    mutate(Cluster = factor(Cluster, levels = seq(1:k_max)))
  
  metadata <- metadata %>%
    left_join(data.PCA[, c("Patient_days", "Cluster")], by = "Patient_days")
  
  pca_list$PCA <- PCA
  pca_list$varexp <- varexp
  pca_list$data.PCA <- data.PCA
  pca_list$metadata <- metadata
  
  return(pca_list)
}
