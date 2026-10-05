# Hannah-Marie Lamle 
# Preprocessing of the spectral data for analysis 


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



# Quality Control: ----------------------------------


## TP1 Standards:  --------------------

tp1_qc <- tp1_spectra %>%
  filter(
    str_detect(
      fragment_ID,
      regex("^(spectralon|skeleton)_", ignore_case = TRUE)
    )
  ) %>%
  mutate(
    qc_type = case_when(
      str_detect(fragment_ID, regex("^spectralon_", ignore_case = TRUE)) ~ "Spectralon",
      str_detect(fragment_ID, regex("^skeleton_", ignore_case = TRUE)) ~ "Skeleton"
    )
  ) %>%
  unnest(spectrum)



tp1_qc_plot <- ggplot(
  tp1_qc,
  aes(
    x = wavelength,
    y = reflectance,
    group = interaction(fragment_ID, scan),
    color = factor(tank)
  )
) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ tank) +
  scale_color_viridis_d(option = "turbo") +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    color = "Tank"
  )

tp1_qc_plot
# tanks 10, 11, and 12 the spectralon scans are super high.... we may have to do 
# some transformation? Let's continue and see how noisy the rest of the time points are
# but here's a list of things to check in on: 
# tank 2 super low spectralon? 
# tank 4 confirm that these scans are dead coral? they look alive
# tank 5 skeleton looks alive
# tank 6 two scans definitely live coral
# tank 8 confirm skeleton
# tank 9 confirm all skeleton scans and remove the live coral spectrum
# tank 10 all scans high, confirm skeleton
# tank 11 all scans high, confirm skeleton
# tank 12 all scans high, confirm skeleton



### Correct scan to frag lookup: --------------------------


### Plot indivdiual standards and remove coral:

plot_tp1_tank <- function(tank_number) {
  
  tp1_qc %>%
    filter(tank == tank_number) %>%
    ggplot(aes(
      x = wavelength,
      y = reflectance,
      group = reflectance_file,
      color = factor(reflectance_file)
    )) +
    geom_line(linewidth = 0.7) +
    labs(
      title = paste("TP1 Spectralon & Skeleton Scans — Tank", tank_number),
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "Scan ID"
    ) +
    theme_classic()
}

plot_tp1_tank(1)
# this looks pretty good. no need to recalibrate anything I don't think. 

plot_tp1_tank(2)
# spectralon is bad. way low. bad geometry? 

plot_tp1_tank(3)
# no skeletons, oops but spectralon looks pretty good

plot_tp1_tank(4)
# good. checked with the data sheets and confirmed those scans are skeleton

plot_tp1_tank(5)
# confirmed skeleton scans 

plot_tp1_tank(6)
# fixed skeleton scans (and coral scans that were incorrect)

plot_tp1_tank(7)
# good. 

plot_tp1_tank(8)
# good. 

plot_tp1_tank(9)
# fixed skeleton scans (and coral scans that were incorrect)

plot_tp1_tank(10)
# confirmed those are skeleton scans 

plot_tp1_tank(11)
# confirmed that is true spectralon and skeleton scans 

plot_tp1_tank(12)
# good. 


### Summary stats of standards: --------------------

tp1_qc_summary <- tp1_qc %>%
  filter(wavelength >= 450, wavelength <= 900) %>%
  group_by(tank, qc_type, scan, fragment_ID, file_time) %>%
  summarise(
    mean_reflectance = mean(reflectance, na.rm = TRUE),
    median_reflectance = median(reflectance, na.rm = TRUE),
    sd_reflectance = sd(reflectance, na.rm = TRUE),
    max_reflectance = max(reflectance, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(
  tp1_qc_summary,
  aes(
    x = file_time,
    y = mean_reflectance,
    color = factor(tank)
  )
) +
  geom_point(size = 3) +
  facet_wrap(~qc_type, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Acquisition time",
    y = "Mean reflectance (450–900 nm)",
    color = "Tank"
  )


ggplot(tp1_qc, aes(wavelength, reflectance, color = factor(tank), group = interaction(tank, scan))) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ qc_type) +
  theme_classic()



ggplot(
  tp1_skeleton,
  aes(
    x = wavelength,
    y = reflectance,
    group = scan
  )
) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ tank, ncol = 4) +
  theme_classic() +
  labs(
    title = "TP1 Skeleton scans by tank",
    x = "Wavelength (nm)",
    y = "Reflectance (%)"
  ) +
  scale_x_continuous(
    limits = c(350, 1050)
  )


## TP2 Standards: -------------------------


tp2_qc <- tp2_spectra %>%
  filter(
    str_detect(
      fragment_ID,
      regex("^(spectralon|skeleton)_", ignore_case = TRUE)
    )
  ) %>%
  mutate(
    qc_type = case_when(
      str_detect(fragment_ID, regex("^spectralon_", ignore_case = TRUE)) ~ "Spectralon",
      str_detect(fragment_ID, regex("^skeleton_", ignore_case = TRUE)) ~ "Skeleton"
    )
  ) %>%
  unnest(spectrum)



tp2_qc_plot <- ggplot(
  tp2_qc,
  aes(
    x = wavelength,
    y = reflectance,
    group = interaction(fragment_ID, scan),
    color = factor(tank)
  )
) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ tank) +
  scale_color_viridis_d(option = "turbo") +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    color = "Tank"
  )

tp2_qc_plot
# TP2 looks muchhhhh better with the newly implemented methodology which is good. 
# no coral fragments in standards scans


## TP3 Standards: -----------------------------------------


tp3_qc <- tp3_spectra %>%
  filter(
    str_detect(
      fragment_ID,
      regex("^(spectralon|skeleton)_", ignore_case = TRUE)
    )
  ) %>%
  mutate(
    qc_type = case_when(
      str_detect(fragment_ID, regex("^spectralon_", ignore_case = TRUE)) ~ "Spectralon",
      str_detect(fragment_ID, regex("^skeleton_", ignore_case = TRUE)) ~ "Skeleton"
    )
  ) %>%
  unnest(spectrum)

tp3_qc_no12 <- tp3_qc %>%
  filter(!tank == 12)

tp3_qc_plot <- ggplot(
  tp3_qc,
  aes(
    x = wavelength,
    y = reflectance,
    group = interaction(fragment_ID, scan),
    color = factor(tank)
  )
) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ tank) +
  scale_color_viridis_d(option = "turbo") +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    color = "Tank"
  )

tp3_qc_plot
# Woah, what is going on with tank 12???? Let's go through each tank. 


### Correct scan to frag lookup: --------------------------


### Plot indivdiual standards and remove coral:

plot_tp3_tank <- function(tank_number) {
  
  tp3_qc %>%
    filter(tank == tank_number) %>%
    ggplot(aes(
      x = wavelength,
      y = reflectance,
      group = reflectance_file,
      color = factor(reflectance_file)
    )) +
    geom_line(linewidth = 0.7) +
    labs(
      title = paste("TP3 Spectralon & Skeleton Scans — Tank", tank_number),
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "Scan ID"
    ) +
    theme_classic()
}

plot_tp3_tank(1)
# looks good, yes there was 5 skeleton scans as well that is correct. 

plot_tp3_tank(2)
# good. 

plot_tp3_tank(3)
# good. 

plot_tp3_tank(4)
# good. 

plot_tp3_tank(5)
# good. 

plot_tp3_tank(6)
# good. 

plot_tp3_tank(7)
# good. 

plot_tp3_tank(8)
# good. 

plot_tp3_tank(9)
# good.

plot_tp3_tank(10)
# the naming scheme from oceanview got messed up for this tank.
# plot all the scans for tank 10: 

tp3_tank10 <- tp3_spectra %>%
  filter(tank == 10) %>%
  unnest(spectrum)

ggplot(
  tp3_tank10,
  aes(
    x = wavelength,
    y = reflectance,
    color = reflectance_file)
) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ fragment_ID) +
  scale_color_viridis_d(option = "turbo") +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    color = "Tank"
  )
# went back through and removed the extra scans with weird names for tank 10. 
# still showing a weird curve for some of the skeleton scans though. not sure why. 


plot_tp1_tank(11)
# lot of baseline variation but looks good

plot_tp3_tank(12)
# crazy high. this tank also has weird naming scheme. removed the wrong names
# and reloaded the dataframe, have a look: 

tp3_tank12 <- tp3_spectra %>%
  filter(tank == 12) %>%
  unnest(spectrum)

ggplot(
  tp3_tank12,
  aes(
    x = wavelength,
    y = reflectance,
    color = reflectance_file)
) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ fragment_ID) +
  scale_color_viridis_d(option = "turbo") +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    color = "Tank"
  )
# yay, looks good now!! 




## TP4 Standards: ---------------------------------------------------


tp4_qc <- tp4_spectra %>%
  filter(
    str_detect(
      fragment_ID,
      regex("^(spectralon|skeleton)_", ignore_case = TRUE)
    )
  ) %>%
  mutate(
    qc_type = case_when(
      str_detect(fragment_ID, regex("^spectralon_", ignore_case = TRUE)) ~ "Spectralon",
      str_detect(fragment_ID, regex("^skeleton_", ignore_case = TRUE)) ~ "Skeleton"
    )
  ) %>%
  unnest(spectrum)



tp4_qc_plot <- ggplot(
  tp4_qc,
  aes(
    x = wavelength,
    y = reflectance,
    group = interaction(fragment_ID, scan),
    color = factor(tank)
  )
) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ tank) +
  scale_color_viridis_d(option = "turbo") +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    color = "Tank"
  )

tp4_qc_plot
# amazing, all look really good. no need to transform




## TP5 Standards: ----------------------------------------------------


tp5_qc <- tp5_spectra %>%
  filter(
    str_detect(
      fragment_ID,
      regex("^(spectralon|skeleton)_", ignore_case = TRUE)
    )
  ) %>%
  mutate(
    qc_type = case_when(
      str_detect(fragment_ID, regex("^spectralon_", ignore_case = TRUE)) ~ "Spectralon",
      str_detect(fragment_ID, regex("^skeleton_", ignore_case = TRUE)) ~ "Skeleton"
    )
  ) %>%
  unnest(spectrum)



tp5_qc_plot <- ggplot(tp5_qc, aes(wavelength, reflectance, group = interaction(fragment_ID, scan), color = factor(tank))) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ tank) +
  scale_color_viridis_d(option = "turbo") +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    color = "Tank"
  )

tp5_qc_plot
# tank 7 has a live coral scan in there
# tank 4 look at the scan that's super flat? weird 


### Correct scan to frag lookup: --------------------------


### Plot indivdiual standards and remove coral:

plot_tp5_tank <- function(tank_number) {
  
  tp5_qc %>%
    filter(tank == tank_number) %>%
    ggplot(aes(
      x = wavelength,
      y = reflectance,
      group = reflectance_file,
      color = factor(reflectance_file)
    )) +
    geom_line(linewidth = 0.7) +
    labs(
      title = paste("TP5 Spectralon & Skeleton Scans — Tank", tank_number),
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "Scan ID"
    ) +
    theme_classic()
}

plot_tp5_tank(1)
# good.

plot_tp5_tank(2)
# good. 

plot_tp5_tank(3)
# good. 

plot_tp5_tank(4)
# good. 

plot_tp5_tank(5)
# good. 

plot_tp5_tank(6)
# good. 

plot_tp5_tank(7)
# coral scan in there. 
# but fixed it. 

plot_tp5_tank(8)
# good. 

plot_tp5_tank(9)
# good.

plot_tp5_tank(10)
# good.

plot_tp5_tank(11)
# good.

plot_tp5_tank(12)
# good. 


# ---------------------- Preprocessing -----------------------------

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

# not sure which of these worked or if I need both so keeping both:
spectra_wide <- spectra_wide %>%
  mutate(
    across(
      where(is.list),
      ~ as.numeric(unlist(.x))
    )
  )

# see above ^
spectra_wide <- spectra_wide %>%
  mutate(
    across(
      where(is.list),
      ~ map_dbl(.x, ~ .x[[1]])
    )
  )

head(spectra_wide) # looks good.


spectral_matrix <- spectra_wide %>%
  select(-fragment_ID, -site, -tank, -treatment, -timepoint, -scan) %>%
  as.matrix()
dim(spectral_matrix)


# testing how well PCA separates the matrices based on the different pre-processing modes: 
spectra_raw <- spectral_matrix # no preprocessing 
spectra_SNV <- standardNormalVariate(spectral_matrix)  # only SNV to control scatter 
spectra_SG <- savitzkyGolay(spectral_matrix, 1, 2, 15) # only savitzky-golay with first order derivative, second binomial and window size of 15
spectra_both <- spectral_matrix %>%                    #  
  standardNormalVariate() %>%
  savitzkyGolay(1, 2, 15)



# PCA on all four different processed spectras: -------------------------------



# 1. Remove spectralon and skeleton scans: 

# Identify rows that are NOT Spectralon or skeleton scans
keep_rows <- !str_detect(
  str_to_lower(spectra_wide$fragment_ID),
  "spectralon|skeleton"
)

# Check how many scans are being retained/removed
sum(keep_rows)
sum(!keep_rows)

# See exactly what was removed
spectra_wide %>%
  filter(!keep_rows) %>%
  count(fragment_ID, sort = TRUE)


# 2. APPLY THE SAME FILTER TO ALL FOUR SPECTRAL MATRICES

spectra_list <- list(
  Raw = spectra_raw[keep_rows, ],
  SNV = spectra_SNV[keep_rows, ],
  SG = spectra_SG[keep_rows, ],
  SG_SNV = spectra_both[keep_rows, ]
)

# Check that all matrices have the same number of rows
lapply(spectra_list, nrow)


# 3. CREATE METADATA FOR THE REMAINING CORAL SCANS

metadata_clean <- spectra_wide %>%
  filter(keep_rows) %>%
  select(
    fragment_ID,
    tank,
    timepoint,
    treatment
  ) %>%
  mutate(
    spectrum_group = paste(
      fragment_ID,
      timepoint,
      sep = "_"
    )
  )


# Check number of scans per fragment × timepoint
metadata_clean %>%
  count(fragment_ID, timepoint) %>%
  count(n)


# 4. AVERAGE REPLICATE SCANS FOR EACH FRAGMENT × TIMEPOINT

spectra_avg <- list()

for (method in names(spectra_list)) {
  
  # Convert spectral matrix to dataframe
  x <- as.data.frame(spectra_list[[method]])
  
  # Add fragment × timepoint grouping variable
  x <- x %>%
    mutate(
      spectrum_group = metadata_clean$spectrum_group
    )
  
  # Average scans within each fragment × timepoint
  x_avg <- x %>%
    group_by(spectrum_group) %>%
    summarise(
      across(
        everything(),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  # Save averaged spectra
  spectra_avg[[method]] <- x_avg
}


# Check dimensions after averaging
lapply(spectra_avg, dim)


# 5. REATTACH METADATA TO THE AVERAGED SPECTRA

avg_metadata <- metadata_clean %>%
  distinct(
    spectrum_group,
    fragment_ID,
    tank,
    timepoint,
    treatment
  )


pca_data <- list()

for (method in names(spectra_avg)) {
  
  pca_data[[method]] <- avg_metadata %>%
    left_join(
      spectra_avg[[method]],
      by = "spectrum_group"
    )
}


# 6. RUN PCA FOR EACH PREPROCESSING METHOD

pca_results <- list()
pca_scores <- list()

for (method in names(pca_data)) {
  
  # Identify spectral columns
  spectral_cols <- setdiff(
    names(pca_data[[method]]),
    c(
      "fragment_ID",
      "tank",
      "timepoint",
      "treatment",
      "spectrum_group"
    )
  )
  
  # Create spectral matrix
  x <- as.matrix(
    pca_data[[method]][, spectral_cols]
  )
  
  # Run PCA
  pca <- prcomp(
    x,
    center = TRUE,
    scale. = FALSE
  )
  
  # Save PCA object
  pca_results[[method]] <- pca
  
  # Calculate variance explained
  variance_explained <- summary(pca)$importance[2, 1:2] * 100
  
  # Extract PC scores and metadata
  scores <- as.data.frame(
    pca$x[, 1:2]
  ) %>%
    bind_cols(
      pca_data[[method]] %>%
        select(
          fragment_ID,
          tank,
          timepoint,
          treatment
        )
    ) %>%
    mutate(
      method = method,
      PC1_variance = variance_explained[1],
      PC2_variance = variance_explained[2]
    )
  
  # Save scores
  pca_scores[[method]] <- scores
}


# 7. MAKE PCA PLOTS

p_raw <- ggplot(
  pca_scores$Raw,
  aes(
    x = PC1,
    y = PC2,
    color = factor(timepoint),
    shape = factor(treatment)
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  labs(
    title = "Raw Spectra",
    x = paste0(
      "PC1 (",
      round(unique(pca_scores$Raw$PC1_variance), 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(unique(pca_scores$Raw$PC2_variance), 1),
      "%)"
    ),
    color = "Timepoint",
    shape = "Treatment"
  ) +
  theme_classic()


p_snv <- ggplot(
  pca_scores$SNV,
  aes(
    x = PC1,
    y = PC2,
    color = factor(timepoint),
    shape = factor(treatment)
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  labs(
    title = "SNV",
    x = paste0(
      "PC1 (",
      round(unique(pca_scores$SNV$PC1_variance), 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(unique(pca_scores$SNV$PC2_variance), 1),
      "%)"
    ),
    color = "Timepoint",
    shape = "Treatment"
  ) +
  theme_classic()


p_sg <- ggplot(
  pca_scores$SG,
  aes(
    x = PC1,
    y = PC2,
    color = factor(timepoint),
    shape = factor(treatment)
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  labs(
    title = "Savitzky-Golay",
    x = paste0(
      "PC1 (",
      round(unique(pca_scores$SG$PC1_variance), 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(unique(pca_scores$SG$PC2_variance), 1),
      "%)"
    ),
    color = "Timepoint",
    shape = "Treatment"
  ) +
  theme_classic()


p_both <- ggplot(
  pca_scores$SG_SNV,
  aes(
    x = PC1,
    y = PC2,
    color = factor(timepoint),
    shape = factor(treatment)
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  labs(
    title = "SNV + Savitzky-Golay",
    x = paste0(
      "PC1 (",
      round(unique(pca_scores$SG_SNV$PC1_variance), 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(unique(pca_scores$SG_SNV$PC2_variance), 1),
      "%)"
    ),
    color = "Timepoint",
    shape = "Treatment"
  ) +
  theme_classic()




# 8. COMBINE ALL FOUR PCA PLOTS



all_pca <- p_raw +
  p_snv +
  p_sg +
  p_both
all_pca

ggsave("figs/pca_preprocessing_results.png", all_pca, width = 11.5, height = 5.5, units = "in")
