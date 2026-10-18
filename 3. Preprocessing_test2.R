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
  select(fragment_ID, site, tank, treatment, timepoint, scan, skeleton_standard, wavelength, reflectance) %>%
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
    fragment_ID, site, tank, treatment, timepoint, skeleton_standard, wavelength
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
# we already did this above (lines 145-155), so can use the spectra_avg df


# first remove the spectralon only, we want to keep skeleton in for this 
spectra_avg_RR <- spectra_avg %>%
  ungroup() %>%
  filter(
    !str_detect(
      str_to_lower(fragment_ID),
      "spectralon"
    )
  ) %>%
  filter(wavelength >= 350, wavelength <= 750) # also filter to only visible range



# transform the df 
avg_RR_wide <- spectra_avg_RR %>%
  pivot_wider(
    id_cols = c(fragment_ID, site, tank, treatment, timepoint, skeleton_standard),
    names_from = wavelength,
    values_from = frag_avg
  )



## Step 2: z-score normalization for the skeleton scans as well? ---------------


# Identify wavelength columns automatically
metadata_cols <- c(
  "fragment_ID", "site", "tank", "treatment", "timepoint", "skeleton_standard"
)

wavelength_cols <- setdiff(names(avg_RR_wide), metadata_cols)

# Z-score each wavelength across samples
RR_zscore <- scale(
  as.matrix(avg_RR_wide[, wavelength_cols]),
  center = TRUE,
  scale = TRUE
)

# Reattach metadata
zscore_RR <- bind_cols(
  avg_RR_wide %>% select(all_of(metadata_cols)),
  as.data.frame(RR_zscore)
)

# Return to long format
zscore_RR_long <- zscore_RR %>%
  pivot_longer(
    cols = all_of(wavelength_cols),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(wavelength = as.numeric(wavelength))


## Step 3: create index (divide by skeleton and log transform) -----------------

# Extract the skeleton spectra to use as a lookup table
skeleton_lookup <- zscore_RR_long %>%
  filter(str_detect(fragment_ID, "^skeleton_tank")) %>%
  select(
    fragment_ID,
    timepoint,
    wavelength,
    skeleton_reflectance = reflectance
  )

# Match each fragment to its assigned skeleton standard
zscore_RR_long <- zscore_RR_long %>%
  left_join(
    skeleton_lookup,
    by = c(
      "skeleton_standard" = "fragment_ID",
      "timepoint",
      "wavelength"
    )
  ) %>%
  mutate(
    fragment_skeleton_ratio = if_else(
      !str_detect(fragment_ID, "^skeleton_tank") &
        !is.na(skeleton_reflectance) &
        skeleton_reflectance != 0,
      reflectance / skeleton_reflectance,
      NA_real_
    )
  )

# log ratio isn't going to work because there's too many negative values because 
# we z-scored the data... 

### plot z-score ratio: ------------------------------------------------

zscore_RR_long %>%
  filter(timepoint == 1) %>%
  ggplot(
    aes(
      x = wavelength,
      y = fragment_skeleton_ratio,
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
      y = fragment_skeleton_ratio,
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
      y = fragment_skeleton_ratio,
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
      y = fragment_skeleton_ratio,
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
      y = fragment_skeleton_ratio,
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

# 1. Create a lookup table containing the averaged skeleton spectra
skeleton_lookup2 <- spectra_avg_RR %>%
  filter(str_detect(fragment_ID, "^skeleton_tank")) %>%
  select(
    fragment_ID,
    timepoint,
    wavelength,
    skeleton_avg = frag_avg
  )


# Extract the skeleton reference spectra
skeleton_ref <- spectra_avg_RR %>%
  filter(grepl("^skeleton_tank", fragment_ID)) %>%
  select(
    tank, timepoint, wavelength,
    skeleton_avg = frag_avg
  )

# 2. Match each fragment to its correct skeleton spectrum
spectra_avg_RR <- spectra_avg_RR %>%
  left_join(
    skeleton_lookup2,
    by = c(
      "skeleton_standard" = "fragment_ID",
      "timepoint",
      "wavelength"
    )
  ) %>%
  mutate(
    log_frag_skeleton_ratio = if_else(
      !str_detect(fragment_ID, "^skeleton_tank") &
        !is.na(skeleton_avg) &
        !is.na(frag_avg) &
        frag_avg > 0 &
        skeleton_avg > 0,
      log(frag_avg / skeleton_avg),
      NA_real_
    )
  )

### plot raw log transform: ----------------------------------------------------


spectra_avg_RR %>%
  filter(timepoint == 1) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_frag_skeleton_ratio,
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




spectra_avg_RR %>%
  filter(timepoint == 2) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_frag_skeleton_ratio,
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



spectra_avg_RR %>%
  filter(timepoint == 3) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_frag_skeleton_ratio,
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



spectra_avg_RR %>%
  filter(timepoint == 4) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_frag_skeleton_ratio,
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




spectra_avg_RR %>%
  filter(timepoint == 5) %>%
  ggplot(
    aes(
      x = wavelength,
      y = log_frag_skeleton_ratio,
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
    subtitle = "Mean Centered reflectance, 350–750 nm",
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


# Helper function: reshape spectra, run PCA, and plot scores

run_spectral_pca <- function(data, value_col, dataset_name) {
  
  # Keep coral fragments only; exclude skeleton reference scans
  coral_data <- data %>%
    filter(!str_detect(fragment_ID, "^skeleton_tank")) %>%
    select(
      fragment_ID, site, tank, treatment, timepoint,
      wavelength, value = all_of(value_col)
    )
  
  # Reshape: one row per fragment/timepoint, one column per wavelength
  spectral_wide <- coral_data %>%
    mutate(wavelength = paste0("wl_", wavelength)) %>%
    pivot_wider(
      id_cols = c(
        fragment_ID, site, tank, treatment, timepoint
      ),
      names_from = wavelength,
      values_from = value
    )
  
  # Separate metadata from wavelength variables
  metadata <- spectral_wide %>%
    select(fragment_ID, site, tank, treatment, timepoint)
  
  spectral_matrix <- spectral_wide %>%
    select(starts_with("wl_")) %>%
    as.data.frame()
  
  # Retain observations with complete spectral data
  keep <- complete.cases(spectral_matrix) &
    apply(spectral_matrix, 1, function(x) all(is.finite(x)))
  
  metadata <- metadata[keep, , drop = FALSE]
  spectral_matrix <- spectral_matrix[keep, , drop = FALSE]
  
  # Run PCA: center wavelengths but do not scale each wavelength
  pca <- prcomp(
    spectral_matrix,
    center = TRUE,
    scale. = FALSE
  )
  
  # Calculate variance explained
  variance_explained <- 100 *
    pca$sdev^2 / sum(pca$sdev^2)
  
  # Store PCA scores alongside coral metadata
  scores <- bind_cols(
    metadata,
    as_tibble(pca$x)
  )
  
  # Plot the first two principal components
  # Plot PC1 vs PC2, faceted by treatment
  p <- ggplot(
    scores,
    aes(
      x = PC1,
      y = PC2,
      color = factor(timepoint),
      shape = factor(timepoint)
    )
  ) +
    geom_point(size = 2.5, alpha = 0.8) +
    facet_wrap(~ treatment) +
    labs(
      title = paste("PCA:", dataset_name),
      x = paste0(
        "PC1 (", round(variance_explained[1], 1), "%)"
      ),
      y = paste0(
        "PC2 (", round(variance_explained[2], 1), "%)"
      ),
      color = "Timepoint",
      shape = "Timepoint"
    ) +
    theme_classic(base_size = 13)
  
  print(p)
  
  return(list(
    pca = pca,
    scores = scores,
    plot = p,
    variance_explained = variance_explained,
    n_excluded = sum(!keep)
  ))
}

# run PCAs: 

# A. Raw averaged reflectance
pca_raw <- run_spectral_pca(
  data = spectra_avg_RR,
  value_col = "frag_avg",
  dataset_name = "Raw averaged reflectance"
)

# B. Log-transformed raw averaged reflectance
spectra_avg_RR <- spectra_avg_RR %>%
  mutate(
    log_frag_avg = if_else(
      !is.na(frag_avg) & frag_avg > 0,
      log(frag_avg),
      NA_real_
    )
  )

pca_log_raw <- run_spectral_pca(
  data = spectra_avg_RR,
  value_col = "log_frag_avg",
  dataset_name = "Log-transformed raw reflectance"
)

# C. Z-scored spectra
pca_zscore <- run_spectral_pca(
  data = zscore_RR_long,
  value_col = "reflectance",
  dataset_name = "Z-scored reflectance"
)

# D. Ratio of z-scored fragment and skeleton spectra
pca_zscore_ratio <- run_spectral_pca(
  data = zscore_RR_long,
  value_col = "fragment_skeleton_ratio",
  dataset_name = "Ratio of z-scored spectra"
)

# E. Log-transformed averaged fragment-to-skeleton ratio
pca_log_ratio <- run_spectral_pca(
  data = spectra_avg_RR,
  value_col = "log_frag_skeleton_ratio",
  dataset_name = "Log fragment-to-skeleton ratio"
)


(pca_raw$plot | pca_log_raw$plot) /
  (pca_zscore$plot | pca_zscore_ratio$plot) /
  (pca_log_ratio$plot)





# -------------- Compare PC2 results through time -----------------------------


## Extract PC2 scores from each PCA ---------------------------
pc2_mean_centered <- pca_raw$scores %>%
  mutate(method = "Mean-centered")

pc2_zscore <- pca_zscore$scores %>%
  mutate(method = "Z-score")

pc2_zscore_ratio <- pca_zscore_ratio$scores %>%
  mutate(method = "Z-score ratio")

pc2_log_ratio <- pca_log_ratio$scores %>%
  mutate(method = "Raw log ratio")

## Combine into one data frame 
pc2_all <- bind_rows(
  pc2_mean_centered,
  pc2_zscore,
  pc2_zscore_ratio,
  pc2_log_ratio
) %>%
  mutate(
    timepoint = as.numeric(as.character(timepoint)),
    method = factor(
      method,
      levels = c(
        "Mean-centered",
        "Z-score",
        "Z-score ratio",
        "Raw log ratio"
      )
    )
  )

# Check that the columns are present
names(pc2_all)
glimpse(pc2_all)


## Make plots: ---------------------------------------
ggplot(
  pc2_all,
  aes(x = timepoint, y = PC2)
) +
  geom_point(
    aes(color = treatment),
    alpha = 0.5,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    linewidth = 0.9
  ) +
  facet_wrap(~method, scales = "free_y") +
  scale_x_continuous(
    breaks = 1:5,
    labels = paste0("TP", 1:5)
  ) +
  labs(
    x = "Timepoint",
    y = "PC2 score",
    color = "Treatment",
    title = "Temporal trends in PC2 across spectral preprocessing methods",
    subtitle = "Black lines show linear regressions; shaded bands show 95% confidence intervals"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

pc2_models <- pc2_all %>%
  group_by(method) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(PC2 ~ timepoint, data = .x)),
    results = map(model, broom::tidy, conf.int = TRUE),
    fit = map(model, broom::glance)
  )

# Regression coefficients, including slope and confidence interval
pc2_regression_results <- pc2_models %>%
  select(method, results) %>%
  unnest(results) %>%
  filter(term == "timepoint") %>%
  select(
    method,
    slope = estimate,
    std_error = std.error,
    statistic = statistic,
    p_value = p.value,
    conf_low = conf.low,
    conf_high = conf.high
  )

# Model fit statistics
pc2_fit_results <- pc2_models %>%
  select(method, fit) %>%
  unnest(fit) %>%
  select(
    method,
    r_squared = r.squared,
    adj_r_squared = adj.r.squared,
    residual_sd = sigma,
    n = nobs
  )

# Combine the results
pc2_comparison <- pc2_regression_results %>%
  left_join(pc2_fit_results, by = "method")

pc2_comparison
