# Calculate the likely share of missing deaths due to late registration
#
# We are interested in the likely share of missing deaths for 2021
# due to late registration by April 26, 2022. In order to calculate
# this index, we check the proportion of 2020 deaths that have been
# reported by April 26, 2021, assuming that the current (April 30, 2022)
# data for 2020 is complete. We perform this check based on historical
# github archives of the World Mortality Dataset. We perform this check
# separately for various countries.

# Init ------------------------------------------------------------

library(yaml)
library(httr)
library(readr)
library(dplyr)
library(ggplot2)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  wmd_snapshot_2021_04_26 = 'https://raw.githubusercontent.com/akarlinsky/world_mortality/4a08b8d449c43de797f4aca1229901f432845073/world_mortality.csv',
  wmd_snapshot_2022_04_30 = 'https://raw.githubusercontent.com/akarlinsky/world_mortality/main/world_mortality.csv'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  late_registration_proportions =
    './out/93-late_registration_proportions.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# constants specific to this analysis
cnst <- within(list(), {

})

# list containers for analysis artifacts
dat <- list()

# Download WMD github archive -------------------------------------

# download April 2021 and April 2022 snapshots of
# World Mortality Dataset from github
dat$wmd_2021 <-
  GET(url = paths$input$wmd_snapshot_2021_04_26, progress()) %>%
  content(as = 'parsed', type = 'text/csv', encoding = 'UTF-8') %>%
  rename(deaths_as_of_2021_04_26 = deaths)
dat$wmd_2022 <-
  GET(url = paths$input$wmd_snapshot_2022_04_30, progress()) %>%
  content(as = 'parsed', type = 'text/csv', encoding = 'UTF-8') %>%
  rename(deaths_as_of_2022_04_30 = deaths)

dat$registration_completeness_2020 <-
  dat$wmd_2022 %>%
  left_join(dat$wmd_2021) %>%
  group_by(iso3c, year) %>%
  summarise(
    deaths_as_of_2021_04_26 = sum(deaths_as_of_2021_04_26),
    deaths_as_of_2022_04_30 = sum(deaths_as_of_2022_04_30)
  ) %>%
  mutate(
    p =
      (deaths_as_of_2021_04_26 / deaths_as_of_2022_04_30)*100
  ) %>%
  filter(year == 2020) %>%
  ungroup()

dat$registration_completeness_2020 %>%
  filter(iso3c == 'USA') %>%
  pull(p)

dat$registration_completeness_2020 %>%
  filter(iso3c == 'HUN') %>%
  pull(p)

dat$registration_completeness_2020 %>%
  filter(!is.na(p) & substr(iso3c, 1,2) %in% config$skeleton$region) %>%
  mutate(
    iso3c = forcats::fct_reorder(iso3c, p)
  ) %>%
  ggplot() +
  geom_point(aes(x = p, y = iso3c))

# Export ----------------------------------------------------------

# export results of analysis
