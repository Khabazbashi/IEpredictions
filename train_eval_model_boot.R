#install.packages("tidyverse")
#install.packages("caret")
#install.packages("leaps")
#install.packages("MASS")

library(rcdk)
library(caret) #machine learning work-flow
library(leaps) #library for step-wise regression
library(MASS) #contains some important linear regression tools
library(caTools)
library(tidyverse) # helps us to write concise code
library(reshape2)

# LOAD DATA
df <- read_delim(
  'Data/data.csv',
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
      solvent,
      error_abs,
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

set.seed(1337)
uniquerows <- sample(nrow(unique_data))
unique_data <- unique_data[uniquerows,]

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
  dplyr::select(-split) 

test_set <- test_set %>%
  dplyr::select(-split) 
  

# CROSS-VALIDATION
set.seed(1337)
  folds = groupKFold(training_set$name, k = 5)  
  fitControl <- trainControl(method = "boot", index = folds)

regressor_fish_boot <- train(logIE ~ .,
  data = training_set %>%
  dplyr::select(-name),
  objective = "reg:squarederror",
  method = "xgbTree",  #"RRF", "xgbDART", "xgbLinear", "svmLinear"
  trControl = fitControl,
  importance = "impurity") 


# ASSESSING THE MODEL 
varian <- varImp(regressor_fish_boot)
variance_table <- varian$importance %>%  #table describing which descriptors are more important
  rownames_to_column() %>%
  filter(Overall > 0)


training_set <- training_set %>%
  mutate(logIE_pred = predict(regressor_fish_boot, newdata = training_set))

test_set <- test_set %>%
  mutate(logIE_pred = predict(regressor_fish_boot, newdata = test_set))

ggplot() +
  geom_point(data = training_set,
             aes(x = logIE, y = logIE_pred, color = "Train")) +
  geom_point(data = test_set,
             aes(x = logIE, y = logIE_pred, color = "Test")) +
  labs(x = "Calculated", y = "Predicted")

write_delim(training_set,"Result/training_1337_pred_xgbTree_boot.csv", delim = ";")
write_delim(test_set,"Result/test_1337_pred_xgbTree_boot.csv", delim = ";")

print(rfFit)
plot(rfFit)
r2_model <- cor(test_set$logIE, test_set$logIE_pred) ^ 2

saveRDS(regressor_fish_boot, "Model/model_seed_1337_xgbTree_boot.rds")
