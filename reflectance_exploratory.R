# Hannah-Marie Lamle
# graphing reflectance data
# 7/8/2026


library(tidyverse)
library(readxl)

# # --------------------- Load Data: ----------------------------
# metadata <- read.csv("TP1_reflectance_lookup.csv")
# 
# 
# reflectance_folder <- "C:/Users/hanna/OneDrive - Florida International University/1.Research/2.Bleaching_spectral_traits/1.data/Reflectance/TP1"
# 
# metadata <- metadata %>%
#   mutate(
#     file_path = file.path(reflectance_folder, reflectance_file)
#   )
# 
# 
# 
# # ------------ Function to import txt files: --------------------
# 
# read_oceanview <- function(file) {
#   
#   lines <- readLines(file)
#   
#   start_line <- which(str_detect(lines, "Begin Spectral Data")) + 1
#   
#   read_tsv(
#     file,
#     skip = start_line - 1,
#     col_names = c("wavelength", "reflectance"),
#     show_col_types = FALSE
#   )
# }
# 
# 
# # ----------------------- read files: --------------------------
# 
# 
# spectra_long <- metadata %>%
#   filter(!is.na(reflectance_file)) %>%
#   mutate(
#     file_path = file.path(reflectance_folder, reflectance_file),
#     data = map(file_path, read_oceanview)
#   ) %>%
#   unnest(data)
# 
# 
# # -------------------- plot raw data: -------------------------
# 
# 
# ggplot(
#   spectra_long,
#   aes(x = wavelength, y = reflectance, group = reflectance_file, color = tank.y)
# ) +
#   geom_line(alpha = 0.25) +
#   theme_classic() +
#   labs(
#     x = "Wavelength (nm)",
#     y = "Reflectance (%)",
#     title = "Raw reflectance spectra"
#   )




# ---------------- Function to import OceanView txt files ----------------

read_oceanview <- function(file) {
  lines <- readLines(file)
  start_line <- which(str_detect(lines, "Begin Spectral Data")) + 1
  
  read_tsv(
    file,
    skip = start_line - 1,
    col_names = c("wavelength", "reflectance"),
    show_col_types = FALSE
  )
}

# ---------------- Function to read one timepoint ----------------

read_timepoint_reflectance <- function(timepoint, lookup_file, reflectance_folder) {
  
  metadata <- read.csv(lookup_file) %>%
    mutate(
      timepoint = timepoint,
      file_path = file.path(reflectance_folder, reflectance_file)
    )
  
  spectra <- metadata %>%
    filter(!is.na(reflectance_file), reflectance_file != "") %>%
    mutate(data = map(file_path, read_oceanview)) %>%
    unnest(data)
  
  return(spectra)
}


# ---------------- define files: --------------------------------

tp_info <- tribble(
  ~timepoint, ~lookup_file,                  ~reflectance_folder,
  "TP1",      "TP1_reflectance_lookup.csv",   "C:/Users/hanna/OneDrive - Florida International University/1.Research/2.Bleaching_spectral_traits/1.data/Reflectance/TP1",
  "TP5",      "TP5_reflectance_lookup.csv",   "C:/Users/hanna/OneDrive - Florida International University/1.Research/2.Bleaching_spectral_traits/1.data/Reflectance/TP5"
)


# ------------------ Read everything: ------------------------

spectra_all <- tp_info %>%
  mutate(
    spectra = pmap(
      list(timepoint, lookup_file, reflectance_folder),
      read_timepoint_reflectance
    )
  ) %>%
  select(spectra) %>%
  unnest(spectra)


# tp5_test <- read.csv("TP5_reflectance_lookup.csv") %>%
#   mutate(
#     reflectance_file = str_trim(reflectance_file),
#     file_path = file.path(
#       "C:/Users/hanna/OneDrive - Florida International University/1.Research/2.Bleaching_spectral_traits/1.data/Reflectance/TP5",
#       reflectance_file
#     )
#   )
# 
# tp5_test %>%
#   filter(!file.exists(file_path)) %>%
#   select(fragment_ID, tank, scan, reflectance_file, file_path)


ggplot(
  spectra_all,
  aes(
    x = wavelength,
    y = reflectance,
    group = interaction(timepoint, reflectance_file),
    color = timepoint
  )
) +
  geom_line(alpha = 0.25) +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    title = "Raw reflectance spectra across timepoints"
  )



ggplot(
  spectra_all,
  aes(
    x = wavelength,
    y = reflectance,
    group = reflectance_file,
    color = treatment
  )
) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ timepoint) +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    title = "Raw reflectance spectra by timepoint"
  )

head(spectra_all)


spectra_test <- spectra_all %>%
  # filter(tank == 11) %>%
  drop_na(site) %>%
  group_by(timepoint, wavelength) %>%
  summarise(mean = mean(reflectance, na.rm = TRUE))


ggplot(spectra_test, aes(wavelength, mean, color = timepoint)) +
  geom_line() +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance (%)",
    title = "Raw reflectance spectra by timepoint"
  )
