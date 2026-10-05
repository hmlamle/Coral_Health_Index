# Hannah-Marie Lamle 
# making the reflectance metadata excel sheet for interns to data entry
# 6/30/26


library(tidyverse)
library(readxl)
library(writexl)


# ---------------------- Just time point 1 -----------------------------


# time point 1: 
df1 <- read_excel("C:\\Users\\hanna\\Florida International University\\Coral Reef Fisheries - 1. EPA_SEAGRANT Coral Spectral\\data\\raw\\CORAL COLOR EXPERIMENT\\1. CORALS\\3. Reflectance\\TP1\\reflectance_metadata.xlsx")


# 1. Duplicate each fragment row 3 times and assign scan 1, 2, 3
frag_scans <- df %>%
  select(-scan, -reflectance_file) %>%   # removes empty columns if they already exist
  crossing(scan = 1:3) %>%
  mutate(
    reflectance_file = NA_character_
  )

spectralon <- df %>%
  distinct(tank, treatment) %>%
  mutate(
    fragment_ID = "spectralon",
    site = NA_character_,
    scan = 1,
    reflectance_file = NA_character_
  )


# 2. Skeleton (3 scans per tank)

skeleton <- df %>%
  distinct(tank, treatment) %>%
  crossing(scan = 1:3) %>%
  mutate(
    fragment_ID = "skeleton",
    site = NA_character_,
    reflectance_file = NA_character_
  )

#---------------------------------
# Combine everything
#---------------------------------
reflectance_metadata <- bind_rows(
  frag_scans,
  spectralon,
  skeleton
) %>%
  arrange(tank, scan)

write_xlsx(reflectance_metadata, "C:\\Users\\hanna\\Florida International University\\Coral Reef Fisheries - 1. EPA_SEAGRANT Coral Spectral\\data\\raw\\CORAL COLOR EXPERIMENT\\1. CORALS\\3. Reflectance\\TP1\\reflectance_metadata_long_TP1.xlsx")




# from master sheet for all time points -------------------



master <- read.csv("C:\\Users\\hanna\\OneDrive - Florida International University\\1.Research\\2.Bleaching_spectral_traits\\3.analysis\\frag_placement\\experiment_randomization\\final_ssid_assignments_06062026.csv")


make_reflectance_sheet <- function(master, timepoint_num) {
  
  removed_tps <- paste0("T", seq_len(timepoint_num - 1))
  
  df_tp <- master %>%
    filter(
      if (timepoint_num == 1) TRUE else !sample_time %in% removed_tps
    ) %>%
    select(fragment_ID, site, tank, treatment)
  
  frag_scans <- df_tp %>%
    crossing(scan = 1:3) %>%
    mutate(reflectance_file = NA_character_)
  
  spectralon <- df_tp %>%
    distinct(tank, treatment) %>%
    mutate(
      fragment_ID = paste0("spectralon_tank", tank),
      site = NA_character_,
      scan = 1,
      reflectance_file = NA_character_
    )
  
  skeleton <- df_tp %>%
    distinct(tank, treatment) %>%
    crossing(scan = 1:3) %>%
    mutate(
      fragment_ID = paste0("skeleton_tank", tank),
      site = NA_character_,
      reflectance_file = NA_character_
    )
  
  bind_rows(frag_scans, spectralon, skeleton) %>%
    arrange(tank, fragment_ID, scan)
}


reflectance_sheets <- list(
  TP1 = make_reflectance_sheet(master, 1),
  TP2 = make_reflectance_sheet(master, 2),
  TP3 = make_reflectance_sheet(master, 3),
  TP4 = make_reflectance_sheet(master, 4),
  TP5 = make_reflectance_sheet(master, 5)
)

write_xlsx(
  reflectance_sheets[1],
  "C:\\Users\\hanna\\Florida International University\\Coral Reef Fisheries - 1. EPA_SEAGRANT Coral Spectral\\data\\raw\\CORAL COLOR EXPERIMENT\\1. CORALS\\3. Reflectance\\TP1\\reflectance_metadata_long_TP1.xlsx")

write_xlsx(
  reflectance_sheets[2],
  "C:\\Users\\hanna\\Florida International University\\Coral Reef Fisheries - 1. EPA_SEAGRANT Coral Spectral\\data\\raw\\CORAL COLOR EXPERIMENT\\1. CORALS\\3. Reflectance\\TP2\\reflectance_metadata_long_TP2.xlsx")

write_xlsx(
  reflectance_sheets[3],
  "C:\\Users\\hanna\\Florida International University\\Coral Reef Fisheries - 1. EPA_SEAGRANT Coral Spectral\\data\\raw\\CORAL COLOR EXPERIMENT\\1. CORALS\\3. Reflectance\\TP3\\reflectance_metadata_long_TP3.xlsx")

write_xlsx(
  reflectance_sheets[4],
  "C:\\Users\\hanna\\Florida International University\\Coral Reef Fisheries - 1. EPA_SEAGRANT Coral Spectral\\data\\raw\\CORAL COLOR EXPERIMENT\\1. CORALS\\3. Reflectance\\TP4\\reflectance_metadata_long_TP4.xlsx")

write_xlsx(
  reflectance_sheets[5],
  "C:\\Users\\hanna\\Florida International University\\Coral Reef Fisheries - 1. EPA_SEAGRANT Coral Spectral\\data\\raw\\CORAL COLOR EXPERIMENT\\1. CORALS\\3. Reflectance\\TP5\\reflectance_metadata_long_TP5.xlsx")

