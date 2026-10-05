# Hannah-Marie Lamle
# Loading reflectance data 
# 7/8/2026


library(readxl)
library(tidyverse)

# ---------------------- Time Point 1 Lookup Table -------------------------
# Separate lookup because this tp the data is messy (swapped names for tanks, etc)
## ------------------------ Load data: ----------------------------------

tp1_meta <- readxl::read_xlsx("data/reflectance_metadata_long_TP1.xlsx")

tp1_files <- list.files(
  path = "data/TP1",
  pattern = "\\.txt$",
  full.names = TRUE
)


## ------------------ Make lookup table for files: --------------------

file_lookup <- tibble(
  full_path = tp1_files, 
  reflectance_file = basename(tp1_files),
  file_time = file.info(tp1_files)$mtime
) %>%
  mutate(
    tank = str_extract(reflectance_file, regex("TANK\\d+", ignore_case = TRUE)),
    tank = str_remove(tank, regex("TANK", ignore_case = TRUE)),
    tank = as.integer(tank)
  ) %>%
  arrange(tank, file_time) %>%
  group_by(tank) %>%
  mutate(order_in_tank = row_number()) %>%
  ungroup()

# correct the tanks that i accidentally swapped: 

file_lookup <- file_lookup %>%
  mutate(
    tank_from_filename = tank,
    tank_corrected = case_when(
      tank_from_filename == 1  ~ 12,
      tank_from_filename == 12 ~ 1,
      tank_from_filename == 2  ~ 11,
      tank_from_filename == 11 ~ 2,
      TRUE ~ tank_from_filename
    )
  )





## --------------------- bind data  -------------------------------


tp1_metadata_final <- tp1_meta %>%
  left_join(
    file_lookup,
    by = c("tank" = "tank_corrected", "order_in_tank")
  ) %>%
  filter(is.na(notes) | !str_detect(notes, regex("circled", ignore_case = TRUE))) # remove circled rows, these were bad scans 





write.csv(tp1_metadata_final, "TP1_reflectance_lookup.csv")




# ------------------------- All other TP's: ----------------------------
## ------------------------ Load data: ----------------------------------

tp2_meta <- read_xlsx("data/reflectance_metadata_long_TP2.xlsx")

tp2_files <- list.files(
  path = "data/TP2",
  pattern = "\\.txt$",
  full.names = TRUE
)

tp3_meta <- read_xlsx("data/reflectance_metadata_long_TP3.xlsx")

tp3_files <- list.files(
  path = "data/TP3",
  pattern = "\\.txt$",
  full.names = TRUE
)


tp4_meta <- read_xlsx("data/reflectance_metadata_long_TP4.xlsx")

tp4_files <- list.files(
  path = "data/TP4",
  pattern = "\\.txt$",
  full.names = TRUE
)


tp5_meta <- read_xlsx("data/reflectance_metadata_long_TP5.xlsx")

tp5_files <- list.files(
  path = "data/TP5",
  pattern = "\\.txt$",
  full.names = TRUE
)


## ------------------ Make lookup table for files: --------------------

# TP2
TP2_lookup <- tibble(
  full_path = tp2_files, 
  reflectance_file = basename(tp2_files),
  file_time = file.info(tp2_files)$mtime
) %>%
  mutate(
    tank = str_extract(reflectance_file, regex("TANK\\d+", ignore_case = TRUE)),
    tank = str_remove(tank, regex("TANK", ignore_case = TRUE)),
    tank = as.integer(tank)
  ) %>%
  arrange(tank, file_time) %>%
  group_by(tank) %>%
  mutate(order_in_tank = row_number()) %>%
  ungroup()



#TP3
TP3_lookup <- tibble(
  full_path = tp3_files, 
  reflectance_file = basename(tp3_files),
  file_time = file.info(tp3_files)$mtime
) %>%
  mutate(
    tank = str_extract(reflectance_file, regex("TANK\\d+", ignore_case = TRUE)),
    tank = str_remove(tank, regex("TANK", ignore_case = TRUE)),
    tank = as.integer(tank)
  ) %>%
  arrange(tank, file_time) %>%
  group_by(tank) %>%
  mutate(order_in_tank = row_number()) %>%
  ungroup()




#TP4
TP4_lookup <- tibble(
  full_path = tp4_files, 
  reflectance_file = basename(tp4_files),
  file_time = file.info(tp4_files)$mtime
) %>%
  mutate(
    tank = str_extract(reflectance_file, regex("TANK\\d+", ignore_case = TRUE)),
    tank = str_remove(tank, regex("TANK", ignore_case = TRUE)),
    tank = as.integer(tank)
  ) %>%
  arrange(tank, file_time) %>%
  group_by(tank) %>%
  mutate(order_in_tank = row_number()) %>%
  ungroup()



#TP5
TP5_lookup <- tibble(
  full_path = tp5_files, 
  reflectance_file = basename(tp5_files),
  file_time = file.info(tp5_files)$mtime
) %>%
  mutate(
    tank = str_extract(reflectance_file, regex("TANK\\d+", ignore_case = TRUE)),
    tank = str_remove(tank, regex("TANK", ignore_case = TRUE)),
    tank = as.integer(tank)
  ) %>%
  arrange(tank, file_time) %>%
  group_by(tank) %>%
  mutate(order_in_tank = row_number()) %>%
  ungroup()




## --------------------- bind data ------------------------------------

#TP2
tp2_metadata_final <- tp2_meta %>%
  left_join(
    TP2_lookup,
    by = c("tank", "order_in_tank")
  ) %>%
  filter(is.na(notes) | !str_detect(notes, regex("circled", ignore_case = TRUE))) # remove circled rows, these were bad scans 


write.csv(tp2_metadata_final, "TP2_reflectance_lookup.csv")




#TP3
tp3_metadata_final <- tp3_meta %>%
  left_join(
    TP3_lookup,
    by = c("tank", "order_in_tank")
  ) %>%
  filter(is.na(notes) | !str_detect(notes, regex("circled", ignore_case = TRUE))) # remove circled rows, these were bad scans 


write.csv(tp3_metadata_final, "TP3_reflectance_lookup.csv")



#TP4
tp4_metadata_final <- tp4_meta %>%
  left_join(
    TP4_lookup,
    by = c("tank", "order_in_tank")
  ) %>%
  filter(is.na(notes) | !str_detect(notes, regex("circled", ignore_case = TRUE))) # remove circled rows, these were bad scans 


write.csv(tp4_metadata_final, "TP4_reflectance_lookup.csv")




#TP5
tp5_metadata_final <- tp5_meta %>%
  left_join(
    TP5_lookup,
    by = c("tank", "order_in_tank")
  ) %>%
  filter(is.na(notes) | !str_detect(notes, regex("circled", ignore_case = TRUE))) # remove circled rows, these were bad scans 


write.csv(tp5_metadata_final, "TP5_reflectance_lookup.csv")


# Yay, all done!!! 