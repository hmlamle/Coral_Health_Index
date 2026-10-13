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
library(prospectr)


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



# ------------------- Test for temporal autocorrelation ------------------------
# pause on this for now -------------------------------------------
# According to Otto 1999 (Chemometrics book 3rd edition) for analyzing 
# time-series data, should first check for correlation in the data. 

# IF NO CORRELATION: can proceed to do statistical tests without transformation. 

# IF CORRELATION: 



# ----------------- Average scans & Clean df:  ---------------------------------

# average triplicate scans per fragment: 
spectra_avg <- spectra_long %>%
  group_by(
    fragment_ID, site, tank, treatment, timepoint, wavelength
  ) %>%
  summarise(
    frag_avg = mean(reflectance, na.rm = TRUE),
    frag_sd = sd(reflectance, na.rm = TRUE),
    n_scans = sum(!is.na(reflectance)),
    .groups = "drop"
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

tp1_unaverage <- spectra_long %>%
  filter(wavelength >= 350, wavelength <= 750, 
         timepoint == 1) %>%
  ggplot(aes(wavelength, reflectance, group = interaction(fragment_ID, scan), color = fragment_ID)) +
  geom_line() +
  scale_x_continuous(
    breaks = seq(350, 750, by = 50)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  facet_wrap(~ tank) +
  theme_classic() +
  theme(legend.position = "none")


tp1_avg <- spectra_avg %>%
  filter(wavelength >= 350, wavelength <= 750, 
         timepoint == 1) %>%
  ggplot(aes(wavelength, frag_avg, color = fragment_ID)) +
  geom_line() +
  scale_x_continuous(
    breaks = seq(350, 750, by = 50)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  facet_wrap(~ tank) +
  theme_classic() +
  theme(legend.position = "none")


tp1_unaverage | tp1_avg
# filters out some of the bad geometry, but still a lot of spread in the data
# which is due to bad geometry 

# Mean Centering Normalization:  ---------------------------------

avg_wide <- spectra_avg_clean %>%
  pivot_wider(
    id_cols = c(fragment_ID, site, tank, treatment, timepoint),
    names_from = wavelength,
    values_from = frag_avg
  )

mcn <- scale(avg_wide[, 6:1112], center = TRUE, scale = FALSE)


# reattach matrix: 
mcn_wide <- cbind(
  avg_wide[, 1:5],
  as.data.frame(mcn)
)


## plot results: ----------------------------------------------------

mcn_long <- mcn_wide %>%
  pivot_longer(
    cols = 6:ncol(.),
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
  coord_cartesian(ylim = c(-50, 150)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 1 Global Mean Center Normalization",
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
  coord_cartesian(ylim = c(-50, 150)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 2 Global Mean Center Normalization",
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
  coord_cartesian(ylim = c(-50, 150)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 3 Global Mean Center Normalization",
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
  coord_cartesian(ylim = c(-50, 150)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 4 Global Mean Center Normalization",
    x = "Wavelength (nm)",
    y = "Mean-centered reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


mcn_long %>%
  filter(timepoint == 5) %>%
  ggplot(
    aes(
      x = wavelength,
      y = reflectance,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-50, 150)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 5 Global Mean Center Normalization",
    x = "Wavelength (nm)",
    y = "Mean-centered reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")

# time point 1 still has an insane amount of variation in the data due to 
# bad geometry not due to biology. this is not a good preprocessing method 
# (in my eyes) to be able to reliably predict traits in the field 
 
# Standard Normal Variate: ----------------------------------------------------

avg_spectral_matrix <- avg_wide %>%
  select(-fragment_ID, -site, -tank, -treatment, -timepoint) %>%
  as.matrix()
dim(avg_spectral_matrix) # should be 899 fragments x 2048 wavelengths 

spectra_SNV <- standardNormalVariate(avg_spectral_matrix)  # only SNV to control scatter 

svn_df <- cbind(
  avg_wide[, 1:5],
  as.data.frame(spectra_SNV)
)

svn_df_long <- svn_df %>%
  pivot_longer(
    cols = 6:ncol(.),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(
    wavelength = as.numeric(wavelength)
  )

## plot results: ----------------------------------------------------


svn_df_long %>%
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
  coord_cartesian(ylim = c(-2, 4)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 1 Standard Normal Variate",
    x = "Wavelength (nm)",
    y = "SVN reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


svn_df_long %>%
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
  coord_cartesian(ylim = c(-2, 4)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 2 Standard Normal Variate",
    x = "Wavelength (nm)",
    y = "SVN reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


svn_df_long %>%
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
  coord_cartesian(ylim = c(-2, 4)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 3 Standard Normal Variate",
    x = "Wavelength (nm)",
    y = "SVN reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


svn_df_long %>%
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
  coord_cartesian(ylim = c(-2, 4)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 4 Standard Normal Variate",
    x = "Wavelength (nm)",
    y = "SVN reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


svn_df_long %>%
  filter(timepoint == 5) %>%
  ggplot(
    aes(
      x = wavelength,
      y = reflectance,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-2, 4)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 5 Standard Normal Variate",
    x = "Wavelength (nm)",
    y = "SVN reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")

# from the SVN transformation, no change in the amplitude of time points in 
# 32C tanks, but changes in the spectral shape. HOWEVER, as Rolo & Ryan 
# pointed out, we cannot use this method individually and be able to compare 
# relative changes in repeated measures 


# Z-score Normalization (global SVN) ------------------------------------------

zscore_spectra <- scale(avg_wide[, 6:1112], center = TRUE, scale = TRUE) # rescale data, creates a matrix


# reattach matrix: 
zscore_wide <- cbind(
  avg_wide[, 1:5],
  as.data.frame(zscore_spectra)
)



## plot results: ----------------------------------------------------

# change the df format: 
zscore_long <- zscore_wide %>%
  pivot_longer(
    cols = 6:ncol(.),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(
    wavelength = as.numeric(wavelength)
  )



zscore_long %>%
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
  coord_cartesian(ylim = c(-2, 6)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 1 Z-score normalization",
    x = "Wavelength (nm)",
    y = "Z-scored reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


zscore_long %>%
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
  coord_cartesian(ylim = c(-2, 6)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 2 Z-score normalization",
    x = "Wavelength (nm)",
    y = "Z-scored reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


zscore_long %>%
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
  coord_cartesian(ylim = c(-2, 6)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 3 Z-score normalization",
    x = "Wavelength (nm)",
    y = "Z-scored reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


zscore_long %>%
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
  coord_cartesian(ylim = c(-2, 6)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 4 Z-score normalization",
    x = "Wavelength (nm)",
    y = "Z-scored reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


zscore_long %>%
  filter(timepoint == 5) %>%
  ggplot(
    aes(
      x = wavelength,
      y = reflectance,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-2, 6)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 5 Z-score normalization",
    x = "Wavelength (nm)",
    y = "Z-scored reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")



# R&R Index: ------------------------------------------------------------------
# lastly trying Rolo and Ryan's index: 

## Step 1: average the triplicate scans for frags and skeletons ---------------
# we already did this above (lines 145-155)

# but we are redoing this because we need to ensure that fragment scans are 
# matchd to the right skeleton scans. 






## Step 2: z-score normalization for the skeleton scans as well? ---------------

# remove spectralon: 
spectra_avg2 <- spectra_avg %>%
  ungroup() %>%
  filter(
    !str_detect(
      str_to_lower(fragment_ID),
      "spectralon"
    )
  )

# Keep only wavelengths from 350 to 750 nm
spectra_avg2_350_750 <- spectra_avg2 %>%
  filter(wavelength >= 350, wavelength <= 750)

# Wide format
avg_wide2 <- spectra_avg2_350_750 %>%
  pivot_wider(
    id_cols = c(fragment_ID, site, tank, treatment, timepoint),
    names_from = wavelength,
    values_from = frag_avg
  )

# Identify wavelength columns automatically
metadata_cols <- c(
  "fragment_ID", "site", "tank", "treatment", "timepoint"
)

wavelength_cols <- setdiff(names(avg_wide2), metadata_cols)

# Z-score each wavelength across samples
zscore_spectra2 <- scale(
  as.matrix(avg_wide2[, wavelength_cols]),
  center = TRUE,
  scale = TRUE
)

# Reattach metadata
zscore_wide2 <- bind_cols(
  avg_wide2 %>% select(all_of(metadata_cols)),
  as.data.frame(zscore_spectra2)
)

# Return to long format
zscore_long2 <- zscore_wide2 %>%
  pivot_longer(
    cols = all_of(wavelength_cols),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(wavelength))


## Step 3: create index (divide by skeleton and log transform) -----------------

# isolate skeleton results 
skeleton_ref <- zscore_long2 %>%
  filter(str_detect(str_to_lower(fragment_ID), "^skeleton_tank")) %>%
  select(
    tank, timepoint, wavelength,
    skeleton_reflectance = reflectance
  )

RR_index <- zscore_long2 %>%
  filter(!str_detect(str_to_lower(fragment_ID), "^skeleton_tank")) %>%
  left_join(
    skeleton_ref,
    by = c("tank", "timepoint", "wavelength")
  ) %>%
  mutate(
    ratio = reflectance / skeleton_reflectance
  )

# log ratio isn't going to work because there's too many negative values because 
# we z-scored the data... 

### plot z-score ratio: ------------------------------------------------

RR_index %>%
  filter(timepoint == 1) %>%
  ggplot(
    aes(
      x = wavelength,
      y = ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 1 RR Index",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")



RR_index %>%
  filter(timepoint == 2) %>%
  ggplot(
    aes(
      x = wavelength,
      y = ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 2 RR Index",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


RR_index %>%
  filter(timepoint == 3) %>%
  ggplot(
    aes(
      x = wavelength,
      y = ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 3 RR Index",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")


RR_index %>%
  filter(timepoint == 4) %>%
  ggplot(
    aes(
      x = wavelength,
      y = ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 4 RR Index",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")



RR_index %>%
  filter(timepoint == 5) %>%
  ggplot(
    aes(
      x = wavelength,
      y = ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 5 RR Index",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")

## non z-scored data log transform: -------------------------------------------

RR_index_nolog <- spectra_avg2

# Extract the skeleton reference spectra
skeleton_ref <- RR_index_nolog %>%
  filter(grepl("^skeleton_tank", fragment_ID)) %>%
  select(
    tank, timepoint, wavelength,
    skeleton_reflectance = frag_avg
  )

# Match each coral measurement to its tank/timepoint skeleton
RR_index_nolog2 <- RR_index_nolog %>%
  filter(!grepl("^skeleton_tank", fragment_ID)) %>%
  left_join(
    skeleton_ref,
    by = c("tank", "timepoint", "wavelength")
  ) %>%
  mutate(
    log_ratio = log(frag_avg / skeleton_reflectance)
  )


### plot raw log transform: ----------------------------------------------------


RR_index_nolog2 %>%
  filter(timepoint == 1) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
   coord_cartesian(ylim = c(-4, 3)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 1 Raw log transform",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")




RR_index_nolog2 %>%
  filter(timepoint == 2) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
   coord_cartesian(ylim = c(-4, 3)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 2 Raw log transform",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")



RR_index_nolog2 %>%
  filter(timepoint == 3) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
   coord_cartesian(ylim = c(-4, 3)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 3 Raw log transform",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")



RR_index_nolog2 %>%
  filter(timepoint == 4) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
   coord_cartesian(ylim = c(-4, 3)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 4 Raw log transform",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")




RR_index_nolog2 %>%
  filter(timepoint == 5) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_ratio,
      group = fragment_ID,
      color = fragment_ID
    )
  ) +
  geom_line(alpha = 0.5) +
  coord_cartesian(ylim = c(-4, 3)) +
  facet_wrap(~ treatment) +
  labs(
    title = "Time Point 5 Raw log transform",
    x = "Wavelength (nm)",
    y = "Biologically filtered Reflectance"
  ) +
  theme_classic() +
  theme(legend.position = "none")




# PCA for preprocessing methods: -----------------------------------------------

# How are the time point 1 fragments clustering within each treatment using 
# the different preprocessing methods? 


## Mean centering normalization --------------------------------------------
# Step 1: Filter wavelength range and transform to wide format 
pca_wide <- spectra_avg_clean %>%
  filter(wavelength >= 350, wavelength <= 750) %>%
  pivot_wider(
    id_cols = c(
      fragment_ID, site, tank, treatment, timepoint
    ),
    names_from = wavelength,
    values_from = frag_avg
  )

# Step 2: Identify wavelength columns
metadata_cols <- c(
  "fragment_ID", "site", "tank", "treatment", "timepoint"
)

wavelength_cols <- setdiff(names(pca_wide), metadata_cols)

# Step 3: Extract numeric spectral data
spectral_matrix <- as.matrix(pca_wide[, wavelength_cols])

# Remove rows with missing or non-finite spectral values
keep <- complete.cases(spectral_matrix) &
  apply(is.finite(spectral_matrix), 1, all)

pca_metadata <- pca_wide[keep, metadata_cols]
spectral_matrix <- spectral_matrix[keep, , drop = FALSE]

# Step 4: Run PCA
pca_raw <- prcomp(
  spectral_matrix,
  center = TRUE,
  scale. = FALSE
)

# Step 5: Combine PCA scores with metadata
pca_scores <- bind_cols(
  pca_metadata,
  as.data.frame(pca_raw$x[, 1:2])
) %>%
  mutate(
    timepoint = factor(timepoint),
    treatment = factor(treatment)
  )

# Calculate variance explained by each PC
var_explained <- summary(pca_raw)$importance[2, ] * 100

ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = timepoint
  )
) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ treatment) +
  labs(
    title = "PCA of Coral Reflectance Spectra",
    subtitle = "Raw reflectance, 350–750 nm",
    x = paste0(
      "PC1 (", round(var_explained[1], 1), "% variance)"
    ),
    y = paste0(
      "PC2 (", round(var_explained[2], 1), "% variance)"
    ),
    color = "Timepoint"
  ) +
  theme_classic()


## Z-score normalization: ------------------------------------------------------


pca_zscore <- prcomp(
  spectral_matrix,
  center = TRUE,
  scale. = TRUE
)

pca_z_scores <- bind_cols(
  pca_metadata,
  as.data.frame(pca_zscore$x[, 1:2])
) %>%
  mutate(
    timepoint = factor(timepoint),
    treatment = factor(treatment)
  )

z_var_explained <- summary(pca_zscore)$importance[2, ] * 100

ggplot(
  pca_z_scores,
  aes(x = PC1, y = PC2, color = timepoint)
) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ treatment) +
  labs(
    title = "PCA of Wavelength-wise Z-scored Coral Spectra",
    subtitle = "350–750 nm",
    x = paste0(
      "PC1 (", round(z_var_explained[1], 1), "% variance)"
    ),
    y = paste0(
      "PC2 (", round(z_var_explained[2], 1), "% variance)"
    ),
    color = "Timepoint"
  ) +
  theme_classic()


## z-scored RR ratio: -----------------------------------------------------

# Prepare the data for PCA
pca_ratio_wide <- RR_index %>%
  filter(
    !str_detect(str_to_lower(fragment_ID), "spectralon|skeleton"),
    wavelength >= 350,
    wavelength <= 750
  ) %>%
  select(
    fragment_ID, site, tank, treatment, timepoint,
    wavelength, ratio
  ) %>%
  pivot_wider(
    id_cols = c(fragment_ID, site, tank, treatment, timepoint),
    names_from = wavelength,
    values_from = ratio
  )

metadata_cols <- c(
  "fragment_ID", "site", "tank", "treatment", "timepoint"
)

# Exclude metadata columns first
candidate_wavelengths <- setdiff(
  names(pca_ratio_wide),
  metadata_cols
)

# Keep wavelengths between 350 and 750 nm
wavelength_numbers <- as.numeric(candidate_wavelengths)

wavelength_cols <- candidate_wavelengths[
  !is.na(wavelength_numbers) &
    wavelength_numbers >= 350 &
    wavelength_numbers <= 750
]

# Check
length(wavelength_cols)
head(wavelength_cols)
tail(wavelength_cols)

# Create the matrix directly
ratio_matrix <- as.matrix(
  pca_ratio_wide[, wavelength_cols, drop = FALSE]
)

# Confirm it is a numeric matrix
class(ratio_matrix)
typeof(ratio_matrix)

# Count missing or infinite values
sum(!is.finite(ratio_matrix))



# wavelength_cols <- setdiff(names(pca_ratio_wide), metadata_cols)
# 
# 
# # Convert each wavelength column to numeric
# pca_ratio_wide <- pca_ratio_wide %>%
#   mutate(
#     across(
#       all_of(wavelength_cols),
#       ~ as.numeric(as.character(.x))
#     )
#   )
# 
# # Create numeric matrix
# ratio_matrix <- as.matrix(
#   pca_ratio_wide[, wavelength_cols, drop = FALSE]
# )
# 
# # Check the result
# class(ratio_matrix)
# typeof(ratio_matrix)

# Number of coral samples before filtering
nrow(ratio_matrix)

# Number of samples with complete, finite spectra
sum(apply(ratio_matrix, 1, function(x) all(is.finite(x))))

# Number of samples being removed
sum(!apply(ratio_matrix, 1, function(x) all(is.finite(x))))

# Number of invalid values per sample
invalid_per_sample <- rowSums(!is.finite(ratio_matrix))

summary(invalid_per_sample)

# Remove samples with missing or non-finite values
keep <- apply(ratio_matrix, 1, function(x) all(is.finite(x)))

ratio_matrix_clean <- ratio_matrix[keep, , drop = FALSE]
ratio_metadata <- pca_ratio_wide[keep, metadata_cols]

# Run PCA
pca_ratio <- prcomp(
  ratio_matrix_clean,
  center = TRUE,
  scale. = FALSE
)

# Calculate variance explained
variance_ratio <- 100 * pca_ratio$sdev^2 /
  sum(pca_ratio$sdev^2)

# Combine PCA scores with sample metadata
pca_ratio_scores <- bind_cols(
  ratio_metadata,
  as.data.frame(pca_ratio$x)
)

# Plot PC1 vs PC2
ggplot(
  pca_ratio_scores,
  aes(
    x = PC1,
    y = PC2,
    color = factor(timepoint)
  )
) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~ treatment) +
  labs(
    title = "PCA of Z-Scored Coral-to-Skeleton Ratios",
    x = paste0("PC1 (", round(variance_ratio[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_ratio[2], 1), "%)"),
    color = "Timepoint"
  ) +
  theme_classic()



## Raw spectra log(RR index): ---------------------------------------------------

# TP2 Tank 7 has multiple skeletons that i needed to redo so its all messed up. 
# Need to load "tp2_metadata_final" from script 1 into enviornment to run this
# code:

library(stringr)

# Create the skeleton reference lookup
skeleton_ref <- spectra_avg2 %>%
  filter(str_detect(str_to_lower(fragment_ID), "^skeleton_tank")) %>%
  mutate(
    calibration = str_extract(
      fragment_ID,
      "calibration\\d+"
    ),
    calibration = coalesce(calibration, "original")
  ) %>%
  select(
    tank,
    timepoint,
    wavelength,
    calibration,
    skeleton_reflectance = frag_avg
  )

# Prepare coral spectra and assign their calibration
coral_spectra <- spectra_avg2 %>%
  filter(!str_detect(str_to_lower(fragment_ID), "^skeleton_tank|^spectralon_tank")) %>%
  left_join(
    tp2_metadata_final %>%
      select(fragment_ID, tank, reflectance_file, notes) %>%
      distinct(),
    by = c("fragment_ID", "tank")
  ) %>%
  mutate(
    calibration = str_extract(
      notes,
      "calibration\\d+"
    ),
    calibration = coalesce(calibration, "original")
  )

# Join each coral spectrum to its matching skeleton calibration
RR_index_nolog2 <- coral_spectra %>%
  left_join(
    skeleton_ref,
    by = c("tank", "timepoint", "wavelength", "calibration")
  ) %>%
  mutate(
    log_ratio = log(frag_avg / skeleton_reflectance)
  )







# Prepare the data for PCA
pca_logratio_wide <- RR_index_nolog2 %>%
  filter(
    !str_detect(str_to_lower(fragment_ID), "spectralon|skeleton"),
    wavelength >= 350,
    wavelength <= 750,
    is.finite(log_ratio)
  ) %>%
  select(
    fragment_ID, site, tank, treatment, timepoint,
    wavelength, log_ratio
  ) %>%
  pivot_wider(
    id_cols = c(fragment_ID, site, tank, treatment, timepoint),
    names_from = wavelength,
    values_from = log_ratio
  )

# Identify wavelength columns
metadata_cols <- c(
  "fragment_ID", "site", "tank", "treatment", "timepoint"
)

wavelength_cols <- setdiff(names(pca_logratio_wide), metadata_cols)

# Convert to a numeric matrix
logratio_matrix <- as.matrix(
  pca_logratio_wide[, wavelength_cols]
)
storage.mode(logratio_matrix) <- "numeric"

# Remove samples with missing or non-finite values
keep <- apply(logratio_matrix, 1, function(x) all(is.finite(x)))

logratio_matrix_clean <- logratio_matrix[keep, , drop = FALSE]
logratio_metadata <- pca_logratio_wide[keep, metadata_cols]

# Run PCA
pca_logratio <- prcomp(
  logratio_matrix_clean,
  center = TRUE,
  scale. = FALSE
)

# Calculate variance explained
variance_logratio <- 100 * pca_logratio$sdev^2 /
  sum(pca_logratio$sdev^2)

# Combine scores with metadata
pca_logratio_scores <- bind_cols(
  logratio_metadata,
  as.data.frame(pca_logratio$x)
)

# Plot PC1 vs PC2
ggplot(
  pca_logratio_scores,
  aes(
    x = PC1,
    y = PC2,
    color = factor(timepoint)
  )
) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~ treatment) +
  labs(
    title = "PCA of Log-Transformed Raw Coral-to-Skeleton Ratios",
    x = paste0("PC1 (", round(variance_logratio[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_logratio[2], 1), "%)"),
    color = "Timepoint"
  ) +
  theme_classic()
