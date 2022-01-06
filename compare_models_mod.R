library("plotly")
library("RRF")
library("stringr")
library("Metrics")
source("sara_functions_mod.R")

########################################
# TRAIN MODELS AND COMPARE TO BASELINE #
########################################

# Load data and remove NA columns
data <- load_data("Data/descs_old_data.csv", "Data/descs_pcb.csv", "Data/descs_pcb_test.csv")

# Unique: Split both and concatenate train/train and test/test
old_splits <- split_data(data$df_old) 
train <- rbind(old_splits$train, data$df_pcb) 

old_splits$test$pcb <- FALSE
data$df_pcb_test$pcb <- TRUE

test <- rbind(old_splits$test, data$df_pcb_test)

# Preprocess
train <- preprocess_data(train)
folds <- get_folds(train, 10)

# Training models
models <-
  list("rf") #(models: xgbTree, rf, RRF, xgbDART, xgbLinear, qrf)

fitted_models <- train_models(models,
                              train,
                              folds,
                              multi_proccessing = TRUE)

# Predicting test data and evaluating the model
preds <-
  data.frame(logIE = test$logIE, pcb = test$pcb) 

for (m in fitted_models) {
  preds[[m$method]] = predict(m, newdata = test, predict.all = TRUE)
}

annots = preds[preds$pcb == TRUE,]

ymax = max(max(preds$logIE)) + 0.1
ymin = min(min(preds$logIE)) - 0.1

xmax = max(max(apply(preds[3:length(colnames(preds))], 2, max))) + 0.1
xmin = min(min(apply(preds[3:length(colnames(preds))], 2, min))) - 0.1

figs <- list()
for (c in colnames(preds)[3:length(colnames(preds))]) {
  fig <- plot_ly(
    x = preds[, c],
    y = preds$logIE,
    type = 'scatter',
    mode = 'markers'
  ) %>% add_annotations(
    x = xmax - 0.7,
    y = ymin + 0.3,
    text = fitted_models[[c]]$method,
    showarrow = F,
    xref = 'x',
    yref = "y",
    font = list(
      color = '#264E86',
      family = 'sans serif',
      size = 12
    ),
    showlegend = FALSE
  ) %>% add_annotations(
    x = xmin + 0.4,
    y = ymax - 0.3,
    text = str_c("RMSE", round(rmse(
      exp(preds$logIE), exp(preds[, c])
    ), digits = 4), sep = ": "),
    showarrow = F,
    xref = 'x',
    yref = "y",
    font = list(
      color = '#264E86',
      family = 'sans serif',
      size = 12
    ),
    showlegend = FALSE
  ) %>% add_annotations(
    x = annots[, c],
    y = annots$logIE,
    text = "pcb",
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 5,
    arrowsize = .8,
    ax = 20,
    ay = -40
  ) %>% layout(
    xaxis = list(title = "pred", range = c(xmin, xmax)),
    yaxis = list(title = "true", range = c(ymin, ymax))
  )
  figs[[c]] <- fig
}

subplot(figs,
        titleX = TRUE,
        titleY = TRUE,
        nrows = 1) %>% layout(showlegend = FALSE)


# TABLE DESCRIBING IMPORTANT DESCRIPTORS
varian <- varImp(fitted_models$rf)
variance_table <- varian$importance %>%  
  rownames_to_column() %>%
  filter(Overall > 0)

saveRDS(fitted_models$rf, "Model/0.90rfbest.rds")
