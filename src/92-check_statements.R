# Check specific statements made in paper

# Init ------------------------------------------------------------

library(yaml); library(tidyverse); library(purrr)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  region_metadata = './cfg/region_metadata.csv',
  figspecs = './cfg/figure_specification.R',
  lifetables = './out/40-lifetables.rds',
  cntfc_lt_debug = './tmp/24-counterfactual_lt_debug.rds',
  skeleton = './tmp/10-harmonized_skeleton.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  fig_e0forecast = './out',
  rds_e0forecastT = './out/90-e0forecastT.rds',
  rds_e0forecastF = './out/90-e0forecastF.rds',
  rds_e0forecastM = './out/90-e0forecastM.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  regions_for_analysis = config$regions_for_all_cause_analysis
})

# meta data on regions
region_meta <- read_csv(paths$input$region_meta, na = '.')

# list containers for analysis artifacts
dat <- list()
fig <- list()

# Functions -------------------------------------------------------

# figure specifications
source(paths$input$figspec)

# Load data -------------------------------------------------------

dat$lifetables <-
  readRDS(paths$input$lifetables) %>%
  left_join(region_meta, by = c(region_iso = 'region_code_iso3166_2')) %>%
  filter(region_iso %in% cnst$regions_for_analysis, quarter == 'annual')
dat$cntfc_lt_debug <- readRDS(paths$input$cntfc_lt_debug)
dat$skeleton <- readRDS(paths$input$skeleton)

# 2 subsequent losses
dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual') %>%
  select(region_name, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_diff_mean_2020 < 0, ex_diff_mean_2021 < 0) %>%
  pull(region_name)

# incomplete bounce-back
dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual') %>%
  select(region_name, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_diff_mean_2020 < 0, ex_diff_mean_2021 > 0, ex_mean_2019 > ex_mean_2021)

# exceeding pre-pandemic LE in 2021
dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual') %>%
  select(region_name, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_mean_2019 < ex_mean_2021)

# still loss in 2021 compared with 2019
dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual') %>%
  select(region_name, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_mean_2019 > ex_mean_2021)

# first drop below 2019 in 2021
dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual') %>%
  select(region_name, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_diff_mean_2020 >= 0, ex_diff_mean_2021 < 0, ex_mean_2019 > ex_mean_2021)

# lower LE in 2021 than expected
dat$lifetables %>%
  filter(age == 0, sex == 'T') %>%
  select(region_name, year, projected, ex_mean) %>%
  pivot_wider(names_from = c(projected, year), values_from = c(ex_mean)) %>%
  filter(projected_2021 > actual_2021)

# any bounce back
dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual') %>%
  select(region_name, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_diff_mean_2020 < 0, ex_diff_mean_2021 > 0) %>%
  pull(region_name)

# younger age group contributed more than older age group to
# e0 losses in 2021 given losses in 2021
dat$lifetables %>%
  filter(sex == 'T', projected == 'actual') %>%
  mutate(
    age_group = cut(age, c(0, 60, Inf), right = FALSE,
                    labels = c('<60', '60+'))
  ) %>%
  select(
    region_iso, age, age_group, year, e0_cntrb_t_mean
  ) %>%
  group_by(region_iso, year, age_group) %>%
  summarise(e0_cntrb_t_mean = sum(e0_cntrb_t_mean)) %>%
  ungroup() %>%
  filter(year > 2019) %>%
  pivot_wider(names_from = c(year, age_group), values_from = e0_cntrb_t_mean) %>%
  filter((`2021_<60` + `2021_60+` < 0), `2021_<60` < 0, `2021_<60` < `2020_<60`)
