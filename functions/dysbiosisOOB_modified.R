dysbiosisOBB_modified <- function(x= NULL,
                         group_col = NULL,
                         case_label = NULL,
                         seed_value = 1235,
                         add_tuneRF_params = list(ntreeTry=1000,
                                                  stepFactor=1.5,
                                                  improve=0.01,
                                                  trace=TRUE,
                                                  dobest=FALSE),
                         ntree = 1000,
                         plot_roc = TRUE,
                         ...){
  
  prop <- group.var <- arg.list.tuneRF.all <- bestmtry <- min.oob <- NULL
  OOBError <- forest.res <- case.prob <- roc2 <- mtry <- NULL
  if(is.null(x) ||
     is.null(group_col) ||
     is.null(case_label)){
    stop("All arguments must be specified")
  }
  message("The random seed of your session for reproducibility is: ", seed_value)
  
  # Get abundance
  prop <- abundances(x)
  
  if(is.factor(meta(x)[,group_col])){
    group.var <- meta(x)[,group_col]
  } else {
    group.var <- as.factor(meta(x)[,group_col])
  }
  
  arg.list.tuneRF.all <- c(list(x=t(prop),
                                y=group.var,
                                sampsize=table(group.var)), add_tuneRF_params)
  set.seed(seed_value)
  bestmtry <- as.data.frame(do.call("tuneRF", arg.list.tuneRF.all))
  
  min.oob <- bestmtry |>
    dplyr::filter(OOBError ==  min(bestmtry$OOBError)) |>
    dplyr::pull(mtry)
  
  message("The mtry used is : ", min.oob)
  
  forest.res <- randomForest(x=t(prop),
                             y=group.var,
                             sampsize=table(group.var),
                             importance = TRUE,
                             ntree = ntree,
                             mtry = min.oob,
                             ...)
  
  case.prob <- stats::predict(forest.res,  type='prob')[, case_label]
  if(plot_roc){
    roc2 <- pROC::roc(group.var,
                      case.prob ,
                      plot=TRUE,
                      ci = TRUE,
                      auc.polygon=TRUE,
                      max.auc.polygon=TRUE,
                      grid=TRUE,
                      print.auc=TRUE)
    #message(ci(roc2))
  }
  
  res.df <- data.frame(oob.score = case.prob)
  
  #importance <- data.frame(varImp(forest.res, scale=FALSE)$importance)
  
  l <- list()
  l[["dataframe"]] <- cbind(res.df, meta(x))
  l[["model"]] <- forest.res
  
  return(l)
  
}
