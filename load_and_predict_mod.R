source("sara_functions_mod.R")
library("Metrics")
library("tidyverse")
library("readr")

########################################
#       LOAD MODEL AND PREDICT         #
########################################

# Load all data
data <- load_data("Data/descs_old_data.csv", "Data/descs_pcb.csv", "Data/descs_pcb_test.csv")
old_splits <- split_data(data$df_old)

# Split both and concatenate train/train and test/test
train <- rbind(old_splits$train, data$df_pcb)
test <- rbind(old_splits$test, data$df_pcb_test)

# Load models
model_1 <- readRDS("Model/0.90rfbest.rds")
model_0 <- readRDS("Model/regressor_neg_new.rds")

# Make predictions with old model
prediction_old = predict(model_0, newdata = test, predict.all = TRUE)
prediction_old = prediction_old$aggregate

# Make predictions with new model 
predictions <- test %>%
  dplyr::select(logIE, name, organic) %>%
  mutate(prediction_new = predict(model_1, newdata = test, predict.all = TRUE)) %>%
  mutate(prediction_old = prediction_old)

# Evaluation
rmse_new <- round(10^rmse((predictions$logIE), (predictions$prediction_new)), digits = 4)
rmse_old <- round(10^rmse((predictions$logIE), (predictions$prediction_old)), digits = 4)

r2_new <- cor(test$logIE, predictions$prediction_new) ^ 2
r2_old <- cor(test$logIE, predictions$prediction_old) ^ 2

# Save results
write.csv(predictions,'predictions.csv')