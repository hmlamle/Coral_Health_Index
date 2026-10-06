library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(tidyverse)
library(patchwork)
library(prospectr)



# time point 1 fragment plotting to look at scatter effects: ---------


tank1_tp1 <- tp1_spectra %>%
  filter(
    tank == 1,
    !str_detect(fragment_ID, regex("^spectralon_", ignore_case = TRUE))
  ) %>%
  unnest(spectrum)

ggplot(tank1_tp1, aes(x = wavelength, y = reflectance, group = interaction(fragment_ID, scan), color = fragment_ID)) +
  geom_line(alpha = 0.7) +
  theme_classic() +
  facet_wrap(~ fragment_ID) +
  labs(
    title = "TP1 Reflectance — Tank 1",
    x = "Wavelength (nm)",
    y = "Reflectance"
  ) +
  theme(legend.position = "none")



# standard normal variate preprocess: 

tank1_wide <- tank1_tp1 %>%
  select(fragment_ID, site, tank, treatment, scan, wavelength, reflectance) %>%
  pivot_wider(
    names_from = wavelength,
    values_from = reflectance
  )

head(tank1_wide)

tank1_matrix <- tank1_wide %>%
  select(-fragment_ID, -site, -tank, -treatment, -scan) %>%
  as.matrix()
dim(tank1_matrix)

tank1_SNV <- standardNormalVariate(tank1_matrix)  # only SNV to control scatter 

# bind back to metadata: 

tank1_metadata <- tank1_wide %>%
  select(fragment_ID, site, tank, treatment, scan)

tank1_SNV_df <- bind_cols(
  tank1_metadata,
  as.data.frame(tank1_SNV)
)

tank1_SNV_long <- tank1_SNV_df %>%
  pivot_longer(
    cols = -c(fragment_ID, site, tank, treatment, scan),
    names_to = "wavelength",
    values_to = "reflectance"
  ) %>%
  mutate(
    wavelength = as.numeric(wavelength)
  )

ggplot(tank1_SNV_long, aes(x = wavelength, y = reflectance, group = interaction(fragment_ID, scan), color = fragment_ID)) +
  geom_line(alpha = 0.7) +
  theme_classic() +
  labs(
    title = "TP1 Reflectance — Tank 1",
    x = "Wavelength (nm)",
    y = "Reflectance"
  ) +
  theme(legend.position = "none")


tank1_SNV_long %>%
  filter(fragment_ID == "cave-ssid26-04162026", 
         scan == 1,
         wavelength >= 350,
         wavelength <= 750) %>%
  ggplot(aes(x = wavelength, y = reflectance)) +
  geom_line() +
  theme_classic() +
  labs(
    x = "Wavelength (nm)",
    y = "Normalized Reflectance"
  )
