#meta <- mofa_imputed@samples_metadata
#rownames(meta) <- meta$Patient_days
# Extract latent factors
#latent_factors <- as.data.frame(get_factors(mofa_imputed))
#colnames(latent_factors) <- gsub("group1.","",colnames(latent_factors))

# Fit LMM for each latent factor
fit_lmm_and_extract_pvals <- function(latent_factors, metadata, variables_of_interest) {
  results_list <- list()
  for (var in variables_of_interest) {
    message(var)
    pb <- progress::progress_bar$new(total = ncol(latent_factors))
    results <- vapply(colnames(latent_factors), FUN = function(factor_name) {
      feat <- latent_factors[, factor_name]
      df.test <- cbind(metadata, feat = feat)
      formula <- as.formula(paste0("feat ~ ", var, "+ (1|Record.ID)"))
      fit <- tryCatch({
        suppressMessages(suppressWarnings(
          lmerTest::lmer(formula, data = df.test, REML = FALSE)
        ))
      }, error = function(e) {
        return(NULL)
      })
      if (is.null(fit)) {
        pb$tick()
        return(c('estimate' = NA_real_, 'p.value' = NA_real_))
      }
      res <- tryCatch({
        coef(summary(fit))
      }, error = function(e) {
        return(NULL)
      })
      if (is.null(res) || nrow(res) < 2) {
        pb$tick()
        return(c('estimate' = NA_real_, 'p.value' = NA_real_))
      }
      pb$tick()
      return(c('estimate' = res[2, "Estimate"], 'p.value' = res[2, "Pr(>|t|)"]))
    }, FUN.VALUE = double(2))
    results <- t(results) %>% 
      as_tibble(rownames = 'Factor') %>% 
      mutate(variable = var)
    results_list[[(length(results_list) + 1)]] <- results
  }
  final_results <- bind_rows(results_list) %>% 
    mutate(Adjusted_p_value = p.adjust(p.value, method = "BH"))
  return(final_results)
}

# Set these to whatever variables you would like to test for association with latent factor
#variables_of_interest <- c("SampleLabel", "Group", "Batch", "CLAD_type", "Timepoint", "days_to_clad", "WCC_blood")
# Apply function
#df.res.assoc <- fit_lmm_and_extract_pvals(latent_factors, meta, variables_of_interest)
