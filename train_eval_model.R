#install.packages("tidyverse")
#install.packages("caret")
#install.packages("leaps")
#install.packages("MASS")

library(rcdk)
library(caret) 
library(leaps) 
library(MASS) 
library(caTools)
library(tidyverse) 
library(reshape2)

# LOAD DATA
df <- read_delim(
  'Data/descs_cal_data.csv',
  delim = ";",
  col_names = TRUE,
  trim_ws = TRUE
) %>%
  drop_na()

# REMOVE NON DESCRIPTORS
temp_df <- df %>%
  dplyr::select(
    -c(
      logIE_pred,
      SMILES,
      name,
      organic_modifier,
      additive,
      instrument,
      source,
      logIE,
    )
  )

# REMOVE NEAR-ZERO VARIANCE AND HIGHLY CORRELATED DESCRIPTORS
nzv <- nearZeroVar(temp_df)
temp_df <- temp_df %>%
  dplyr::select(-nearZeroVar(temp_df))  

descr_cor <- cor(temp_df, use = "complete.obs")
temp_df <- temp_df %>%
  dplyr::select(-findCorrelation(descr_cor, cutoff = 0.75))  

temp_df <- temp_df %>%
  dplyr::mutate(logIE = df$logIE, name = df$name)


# COMPUTE INDICES FOR TRAIN PARTITION
unique_data <- temp_df %>%
  dplyr::select(name) %>%
  unique()

split <- sample.split(unique_data$name, SplitRatio = 0.8)
training_set <- unique_data %>%
  filter(split == TRUE) %>%
  mutate(split = "TRUE") %>%
  left_join(temp_df) %>% 
  na.omit()

test_set <- unique_data %>%
  filter(split == FALSE) %>%
  mutate(split = "FALSE") %>%
  left_join(temp_df) %>% 
  na.omit()

temp_df <- rbind(training_set, test_set) %>%
  na.omit()

training_set <- training_set %>%
  dplyr::select(-name, -split) 

test_set <- test_set %>%
  dplyr::select(-name, -split) 
  

# CROSS-VALIDATION OR BOOT
fitControl <- trainControl(
  method = "repeatedcv",
  repeats = 5,
  number = 5,
  verboseIter = TRUE,
)

set.seed(1337)
rfFit <- train(
  logIE ~ .,
  data = training_set,
  method = "xgbTree", #qrf #rf #xgbTree
  #nTrees = 100,
  trControl = fitControl,
  objective = "reg:squarederror",
  verbose = FALSE
)

# ASSESSING THE MODEL 
training_set <- training_set %>%
  mutate(logIE_pred = predict(rfFit, newdata = training_set))

test_set <- test_set %>%
  mutate(logIE_pred = predict(rfFit, newdata = test_set))

ggplot() +
  geom_point(data = training_set,
             aes(x = logIE, y = logIE_pred, color = "Train")) +
  geom_point(data = test_set,
             aes(x = logIE, y = logIE_pred, color = "Test")) +
  labs(x = "Calculated", y = "Predicted")

write_delim(training_set,"Result/training_1337_pred_rf_cv_55_.csv", delim = ";")
write_delim(test_set,"Result/test_1337_pred_rf_cv_55_.csv", delim = ";")

print(rfFit)
plot(rfFit)
r2_cv <- cor(test_set$logIE, test_set$logIE_pred) ^ 2

saveRDS(rfFit, "Model/model_seed_1337_rf_cv_55.rds")
