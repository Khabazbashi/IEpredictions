library("tidyverse")
source("PaDEL_descs_calculator.R")
source("compound_eluent.R")

############################################
# CALCULATE DESCRIPTORS FOR YOUR COMPOUNDS #
############################################

# READ AND LOAD .CSV FILE TO STANDARDS 
standards = read_delim("Data/old_data.csv",
                       delim = ";",
                       col_names = TRUE)

# SELECT COLUMNS FROM .CSV FILE 
standards <- standards %>% 
  dplyr::select(
    c(
      SMILES,
      organic_modifier,
      organic,
      pH.aq.,
      additive,
      instrument,
      source,
      name,
      NH4,
      logIE,
      logIE_pred,
      additive_concentration_mM
    )
  )

# CALCULATE PaDEL DESCRIPTORS 
descriptor_calculated = PaDEL_original(standards)

# ADD VISCOSITY, POLARITY INDEX AND SURFACE TENSTION TO DESCRIPTORS 
descriptor_calculated <- descriptor_calculated %>%
  dplyr::select(SMILES, organic, organic_modifier, additive, additive_concentration_mM, name, logIE, instrument, source, logIE_pred, pH.aq., NH4, everything()) %>%
  mutate(
    viscosity = viscosity(organic, organic_modifier),
    surface_tension = surfacetension(organic, organic_modifier),
    polarity_index = polarityindex(organic, organic_modifier))

# PRINT RESULTS IN RESULT FOLDER
write_delim(
  descriptor_calculated,
  "Data/descs_old_data.csv",
  delim = ";"
)
