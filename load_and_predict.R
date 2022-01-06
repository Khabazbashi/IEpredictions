library("tidyverse")
library("readr")

############################################
# LOAD MODEL AND DATA AND MAKE PREDICTIONS #
############################################
data <- read_delim("Data/descs_temp_pred.csv",
                   delim = ";",
                   col_names = TRUE)

regressor_new = readRDS("Model/model_seed_1337_rf_cv_55.rds")
regressor_old = readRDS("Model/regressor_neg_new.rds")

prediction_new = predict(regressor_new, newdata = data, predict.all = TRUE)
prediction_old = predict(regressor_old, newdata = data, predict.all = TRUE) 
prediction_old = prediction_old$aggregate 

result <- data %>%
  dplyr::select(SMILES, name, logIE) %>%
  mutate(logIE_pred_old = prediction_old) %>%
  mutate(logIE_pred_new = prediction_new) %>%
  mutate(error_new = abs(exp(logIE)-exp(prediction_new))) %>%
  mutate(error_old = abs(exp(logIE)-exp(prediction_old))) %>%
  mutate(old_new_diff_abs = abs(prediction_old - prediction_new))

############################################
#             MODEL EVALUATION             #
############################################
mean_new <- mean(result$error_new)
mean_old <- mean(result$error_old)
r2_new <- cor(result$logIE, result$logIE_pred_new) ^ 2
r2_old <- cor(result$logIE, result$logIE_pred_old) ^ 2

ggplot() +
  geom_point(
    data = result,
    aes(x = logIE, y = logIE_pred_old)
  ) +
  labs(x = "Calculated", y = "Predicted old model")

ggplot() +
  geom_point(
    data = result,
    aes(x = logIE, y = logIE_pred_new)
  ) +
  labs(x = "Calculated", y = "Predicted new model")

write_delim(result,
            "Result/IEpred_temp_pred.csv",
           delim = ";")