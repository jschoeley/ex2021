# Harmonize external life table estimates

# Init ------------------------------------------------------------

library(here); library(glue)
library(yaml)
library(readr); library(dplyr); library(tidyr)

# Constants -------------------------------------------------------

wd <- here()

config <- read_yaml(glue('{wd}/cfg/config.yaml'))
region_meta <- read_csv(glue('{wd}/cfg/region_metadata.csv'), na = '.')

cnst <- list()
cnst <- within(cnst, {
  # path to global objects
  path_global = glue('{wd}/src/00-global.R')
  # where to put the harmonized data
  path_harmonized = glue('{wd}/tmp')
  # where to find the life expectancy data
  path_wpp = glue('{wd}/dat/wpp/wpp_ex.rds')
  path_hmd_females = glue('{wd}/dat/hmdhfd/fltper_1x1.rds')
  path_hmd_males = glue('{wd}/dat/hmdhfd/mltper_1x1.rds')
  path_hmd_total = glue('{wd}/dat/hmdhfd/bltper_1x1.rds')
  # skeleton path
  path_skeleton = glue('{wd}/tmp/10-harmonized_skeleton.rds')
  # translation of ex sex code to harmonized sex code
  code_sex_wpp =
    c(`Male` = config$skeleton$sex$Male,
      `Female` = config$skeleton$sex$Female,
      `Total` = config$skeleton$sex$Total)
  # lookup table for region codes
  # only countries defined in skeleton
  region_lookup_wpp = 
    region_meta %>%
    filter(region_code_iso3166_2 %in% config$skeleton$region) %>%
    select(region_code_iso3166_2, region_code_wpp) %>%
    drop_na()
  # lookup table for region codes
  # only countries defined in skeleton
  region_lookup_hmd = 
    region_meta %>%
    filter(region_code_iso3166_2 %in% config$skeleton$region) %>%
    select(region_code_iso3166_2, region_code_hmd) %>%
    drop_na()
})

# Functions -------------------------------------------------------

source(cnst$path_global)

# Data ------------------------------------------------------------

dat <- list()

dat$skeleton <- readRDS(cnst$path_skeleton)

# wpp estimates
dat$wpp_ex <- readRDS(cnst$path_wpp)
# hmd estimates
dat$hmd_lt_females <- readRDS(cnst$path_hmd_females)
dat$hmd_lt_males <- readRDS(cnst$path_hmd_males)
dat$hmd_lt_total <- readRDS(cnst$path_hmd_total)

# Harmonize WPP ---------------------------------------------------

dat$wpp_clean <-
  dat$wpp_ex %>%
  # select columns of interest
  select(
    ex, sex = Sex,
    region_code_wpp = LocID, iso_year = MidPeriod, 
    age_start = AgeGrpStart, age_width = AgeGrpSpan
  ) %>%
  mutate(
    sex =
      factor(sex, levels = names(cnst$code_sex_wpp),
             labels = cnst$code_sex_wpp) %>%
      as.character(),
    region_iso = factor(
      region_code_wpp,
      levels = cnst$region_lookup_wpp$region_code_wpp,
      labels = cnst$region_lookup_wpp$region_code_iso3166_2
    ) %>% as.character()
  ) %>%
  # add row id
  mutate(id = GenerateRowID(region_iso, sex, iso_year, age_start)) %>%
  select(id, ex_wpp = ex)

# Harmonize HMD ---------------------------------------------------

# merge female, male, and total, prepare for LCA forecast
dat$hmd_clean <-
  bind_rows(
    '{config$skeleton$sex$Male}' := dat$hmd_lt_males,
    '{config$skeleton$sex$Female}' := dat$hmd_lt_females,
    '{config$skeleton$sex$Total}' := dat$hmd_lt_total,
    .id = 'sex'
  ) %>%
  as_tibble() %>% 
  select(
    nmx_hmd = mx, lx_hmd = lx, Tx_hmd = Tx, ex_hmd = ex,
    sex, region_code_hmd, year = Year,
    age_start = Age
  ) %>%
  mutate(
    region_iso = factor(
      region_code_hmd,
      levels = cnst$region_lookup_hmd$region_code_hmd,
      labels = cnst$region_lookup_hmd$region_code_iso3166_2
    ) %>% as.character()
  ) %>%
  filter(age_start <= 100) %>%
  mutate(nmx_hmd = ifelse(age_start == 100, lx_hmd/Tx_hmd, nmx_hmd)) %>%
  # add row id
  mutate(id = GenerateRowID(region_iso, sex, year, age_start)) %>%
  select(id, ex_hmd = ex_hmd, nmx_hmd = nmx_hmd)
  
# Join with skeleton ----------------------------------------------

# Please note that the WPP life expectancy estimates use 
# 0, 1-4, & 5-year age groups instead of single-year age groups

# join the different sources of population count data
# with the skeleton
dat$joined <-
  dat$skeleton %>% 
  left_join(
    dat$wpp_clean,
    by = 'id'
  ) %>%
  left_join(
    dat$hmd_clean,
    by = 'id'
  ) %>%
  select(id, ex_wpp_estimate = ex_wpp, ex_hmd_estimate = ex_hmd,
         nmx_hmd_estimate = nmx_hmd)

# Export ----------------------------------------------------------

saveRDS(
  dat$hmd_clean,
  file = glue('{cnst$path_harmonized}/23-harmonized_hmd_lifetables.rds')
)

saveRDS(
  dat$joined,
  file = glue('{cnst$path_harmonized}/23-harmonized_external_lifetables.rds')
)
