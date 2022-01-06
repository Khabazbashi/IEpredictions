# IEpredictions

This code was used in my thesis project when working on predictions of ionization efficiency for OH-PCBs.
The code is based on previous work by Anneli Kruve and colleagues at Stockholms University.

## Data folder:
Contains files data and holdout_test. Data includes all previous data + 24 OH-PCBs (LC and FIA).
Holdout_test includes 12 OH-PCBs (LC and FIA) which model has never been introduced too.

## Scripts:
 1. `load_and_predict`: Load model from Model-folder and predict IE-values. Saves result in Result- folder.

 1. `train_eval_model`: Script for training a new model. Separating test and train so that no compound appears in both. Saves completed model in Model-folder. Saves result of predictions in Result-folder.

 1. `descriptor_calculator`: Only calculates descriptors. Saves result in Result-folder.
