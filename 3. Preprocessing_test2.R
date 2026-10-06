# Hannah-Marie Lamle, Victor Ruano-Rodriguez
# Preprocessing of the spectral data for analysis 
# more methods and comparing the spread of TP1 to other TP's 
# with these other methods: 
# 10/6/26


library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(tidyverse)
library(patchwork)


# Load metadata: 
# READ ME: You don't need to do this if you are already working straight from
# script 1. If you're picking up with the preprocessing then you will need to 
# reload the metadata sheets

tp1_metadata_final <- read.csv("TP1_reflectance_lookup.csv")
tp2_metadata_final <- read.csv("TP2_reflectance_lookup.csv")
tp3_metadata_final <- read.csv("TP3_reflectance_lookup.csv")
tp4_metadata_final <- read.csv("TP4_reflectance_lookup.csv")
tp5_metadata_final <- read.csv("TP5_reflectance_lookup.csv")





# Load Data using MAP: -----------------------------------------

## function to load the data from the text files: 

read_spectrum <- function(file) {
  
  # If there is no file path, return an empty spectrum
  if (is.na(file) || file == "") {
    return(
      tibble(
        wavelength = numeric(),
        reflectance = numeric()
      )
    )
  }
  
  lines <- readLines(file)
  start <- grep("Begin Spectral Data", lines)
  
  spectrum <- read.table(
    file,
    skip = start,
    header = FALSE
  )
  
  spectrum %>%
    select(
      wavelength = V1,
      reflectance = V2
    )
}


# Make Dataframes: 

tp1_spectra <- tp1_metadata_final %>%
  mutate(
    spectrum = map(full_path, read_spectrum)
  )

tp2_spectra <- tp2_metadata_final %>%
  mutate(
    spectrum = map(full_path, read_spectrum)
  )

tp3_spectra <- tp3_metadata_final %>%
  mutate(
    spectrum = map(full_path, read_spectrum)
  )

tp4_spectra <- tp4_metadata_final %>%
  mutate(
    spectrum = map(full_path, read_spectrum)
  )

tp5_spectra <- tp5_metadata_final %>%
  mutate(
    spectrum = map(full_path, read_spectrum)
  )


# --------- combine TP df's -----------------------------

library(prospectr)

tp1_spectra <- tp1_spectra %>%
  mutate(timepoint = 1)

tp2_spectra <- tp2_spectra %>%
  mutate(timepoint = 2)

tp3_spectra <- tp3_spectra %>%
  mutate(timepoint = 3)

tp4_spectra <- tp4_spectra %>%
  mutate(timepoint = 4)

tp5_spectra <- tp5_spectra %>%
  mutate(timepoint = 5)

all_spectra <- bind_rows(
  tp1_spectra,
  tp2_spectra,
  tp3_spectra,
  tp4_spectra,
  tp5_spectra
)

spectra_long <- all_spectra %>%
  unnest(spectrum)


spectra_wide <- spectra_long %>%
  select(fragment_ID, site, tank, treatment, timepoint, scan, wavelength, reflectance) %>%
  pivot_wider(
    names_from = wavelength,
    values_from = reflectance
  )


# ----------------- Mean Centering Normalization -----------------

# average triplicate scans per fragment: 

spectra_avg <- spectra_long %>%
  select(fragment_ID, site, tank, treatment, timepoint, scan, wavelength, reflectance) %>%
  group_by(fragment_ID, timepoint, wavelength) %>% 
  summarise(
    frag_avg = mean(reflectance), 
    frag_sd = sd(reflectance)
  )


# First remove standards and skeleton scans: 

spectra_avg_clean <- spectra_avg %>%
  ungroup() %>%
  filter(
    !str_detect(
      str_to_lower(fragment_ID),
      "spectralon|skeleton"
    )
  )



## mean centering normalization ---------------------------------

avg_wide <- spectra_avg_clean %>%
  pivot_wider(
    id_cols = c(fragment_ID, timepoint),
    names_from = wavelength,
    values_from = frag_avg
  )

mcn <- scale(avg_wide[, 3:ncol(avg_wide)], center = TRUE, scale = FALSE)


# reattach matrix: 
mcn_wide <- cbind(
  avg_wide[, 1:2],
  as.data.frame(mcn)
)


### plot results: ----------------------------------------------------

mcn_long <- mcn_wide %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(
    wavelength = as.numeric(wavelength)
  )


mcn_long %>%
  filter(timepoint == 1) %>%
  ggplot(
    aes(
      x = wavelength,
      y = reflectance,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean-centered reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


mcn_long %>%
  filter(timepoint == 2) %>%
  ggplot(
    aes(
      x = wavelength,
      y = reflectance,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean-centered reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


mcn_long %>%
  filter(timepoint == 3) %>%
  ggplot(
    aes(
      x = wavelength,
      y = reflectance,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean-centered reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


mcn_long %>%
  filter(timepoint == 4) %>%
  ggplot(
    aes(
      x = wavelength,
      y = reflectance,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean-centered reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


# next steps: -----------------------

# line 135 how to group by and average but keep other columns too.

# try SNV and see how that data looks 
# try Rolo&Ryan index thing? 

# redo the pca's and compare (and separate out TP1 specifically?)

