# install.packages("tidyverse")
# install.packages("caret")
# install.packages("leaps")
# install.packages("MASS")
# install.packages("UBL")
# install.packages("kernlab")
# install.packages("doParallel")
# install.packages("elasticnet")
# install.packages("obliqueRF")
# install.packages("quantregForest")

library(rcdk)
library(caret) 
library(leaps) 
library(MASS) 
library(caTools)
library(tidyverse) 
library(reshape2)
library(RRF)
library(kernlab)

# LOAD DATA AND REMOVE NA COLUMNS
load_data = function(old_data_path,
                     pcb_data_path, 
                     test_data_path) {

  df_old <- read_delim(
    old_data_path,
    delim = ";",
    col_names = TRUE,
    trim_ws = TRUE
  ) %>% as.data.frame()
  
  df_pcb <- read_delim(
    pcb_data_path,
    delim = ";",
    col_names = TRUE,
    trim_ws = TRUE
  ) %>% as.data.frame()
  
  df_pcb_test <- read_delim(
    test_data_path,
    delim = ";",
    col_names = TRUE,
    trim_ws = TRUE
  ) %>% as.data.frame()
  
  nanCols <-
    c(unique(which(is.na(
      rbind(df_old, df_pcb, df_pcb_test)
    ), arr.ind = TRUE)[,2]))
  
  df_old <- dplyr::select(df_old,-nanCols)
  df_pcb <- dplyr::select(df_pcb,-nanCols)
  df_pcb_test <- dplyr::select(df_pcb_test,-nanCols)
  
  return (list(df_old=df_old, df_pcb=df_pcb, df_pcb_test=df_pcb_test))
}


# SPLIT DATA FUNCTION 
split_data = function(df) {
  
  unique_data <- df %>%
    dplyr::select(name) %>%
    unique()
  
  set.seed(123) 
  split <- sample.split(unique_data$name, SplitRatio = 0.8)
  training_set <- unique_data %>%
    filter(split == TRUE) %>%
    mutate(split = "TRUE") %>%
    left_join(df) %>%
    na.omit()
  
  test_set <- unique_data %>%
    filter(split == FALSE) %>%
    mutate(split = "FALSE") %>%
    left_join(df) %>%
    na.omit()
  
  training_set <- training_set %>%
    dplyr::select(-split)
  
  test_set <- test_set %>%
    dplyr::select(-split)
  
  return (list(train=training_set, test=test_set))
}

# DATA PREPROCESSING FUNCTION
preprocess_data = function(df) {
  name <- df$name
  logIE <- df$logIE
  
  # REMOVE NON DESCRIPTORS
  df <- df %>%
    dplyr::select(
      -c(
        logIE,
        logIE_pred,
        SMILES,
        organic_modifier,
        additive,
        instrument,
        name,
        source,
      )
    )
  
  # REMOVE NZV & HIGHLY CORRELATED DESCRIPTORS
  df <- df %>% dplyr::select(-all_of(nearZeroVar(df)))
  
  descr_cor <- cor(df, use = "complete.obs")
  highly_cor_descr <- findCorrelation(descr_cor, cutoff = .90)
  df <- df %>% dplyr::select(-all_of(highly_cor_descr))
  
  df$name <- name
  df$logIE <- logIE
  return(df)
}

# FOLDS IN MODEL TRAINING 
get_folds = function(df, k) {
  set.seed(123)
  folds = groupKFold(df$name, k = k)
  return (folds)
}

# TRAIN MODEL FUNCTION
train_models = function(models,
                        df,
                        folds,
                        control_method = "boot",
                        multi_proccessing = FALSE) {
  
  fitControl <-
    trainControl(method = control_method,
                 index = folds)
  
  if (multi_proccessing) {
    library(doParallel)
    # Use nrCores -1 to not overload CPU
    cl <- makePSOCKcluster(detectCores()[1] - 1)
    registerDoParallel(cl)
  }
  fitted_models = list()
  set.seed(123) 
  for (m in models) {
    fit <- train(
      logIE ~ .,
      data <- df %>% dplyr::select(-c(name)),
      method = m,
      ntree = 100,
      trControl = fitControl,
      importance = "impurity",
      verbose = FALSE
    )
    cat("## DONE TRAINING", m, "MODEL ## \n")
    fitted_models[[m]] = fit
  }
  if (multi_proccessing) {
    stopCluster(cl)
  }
  return(fitted_models)
}
