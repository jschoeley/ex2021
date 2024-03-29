# Harmonize mid-year population counts
#
# (1) Mid-year population estimates for all nation states from
#     World Population Prospects
#     - years < 2020 from estimates
#     - year 2020 from projections
# (2) Mid-year population estimates for England & Wales,
#     Northern Ireland and Scotland from HMD/HFD
#     - years < 2019 from HMD midyear pop estimates
#     - years 2019,2020 projected based on HMD and HFD population,
#      mortality and fertility data

# Init ------------------------------------------------------------

library(glue)
library(readr); library(dplyr); library(tidyr); library(purrr); library(yaml)
library(ggplot2)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  global = './src/00-global.R',
  region_meta = './cfg/region_metadata.csv',
  figspecs = './cfg/figure_specification.R',
  skeleton = './tmp/10-harmonized_skeleton.rds',
  wpp_population = './dat/wpp/wpp_population.rds',
  hmdhfd_gb_projection = './dat/hmdhfd/hmdhfd_gb_projection.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  harmonized_population = './tmp/20-harmonized_population.rds',
  out = './out'
)

# global configuration
config <- read_yaml(paths$input$config)

# meta data on regions
region_meta <- read_csv(paths$input$region_meta, na = '.')

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  # translation of wpp sex code to harmonized sex code
  code_sex_wpp =
    c(`male` = config$skeleton$sex$Male,
      `female` = config$skeleton$sex$Female)
  # translation of hmdhfd sex code to harmonized sex code
  code_sex_hmdhfd =
    c(`Male` = config$skeleton$sex$Male,
      `Female` = config$skeleton$sex$Female)
  # lookup table for region codes
  # only countries defined in skeleton
  region_lookup_wpp = 
    region_meta %>%
    filter(region_code_iso3166_2 %in% config$skeleton$region) %>%
    select(region_code_iso3166_2, region_code_wpp) %>%
    drop_na()
  region_lookup_hmd =
    region_meta %>%
    filter(region_code_iso3166_2 %in% config$skeleton$region) %>%
    select(region_code_iso3166_2, region_code_hmd)
  # first year in harmonized data set
  skeleton_first_year = config$skeleton$year$min
  # last year in harmonized data set
  skeleton_final_year = config$skeleton$year$max
  
  # number of years to project GB population data from 2018
  uk_projection_n_years = 3
})

# list containers for analysis artifacts
dat <- list()
fig <- list()

# Functions -------------------------------------------------------

# figure specifications
source(paths$input$figspecs)

# global functions
source(paths$input$global)

#' Project a Population via Stable Assumption
StableProjection1x1 <- function (
  pop_m, pop_f, mx_m, mx_f, fx_f, n = 1, srb = 1.04
) {
  
  nage = length(pop_m)
  
  # first project the female population for n steps
  # Leslie projection matrix
  A_f <- matrix(0, nrow = nage, ncol = nage)
  diag(A_f[-1,-nage]) <- head(exp(-mx_f), -1)
  A_f[1,] <- fx_f * 1/(1+srb)
  # population matrix, first column is jump off population
  P_f <- matrix(0, nrow = nage, ncol = n+1)
  P_f[,1] <- pop_f
  # project population
  for (i in 1:n+1) {
    P_f[,i] <- A_f%*%P_f[,i-1]
  }
  
  # now project the male population with male births derived from
  # projected female births via sex ratio
  A_m <- matrix(0, nage, ncol = nage)
  diag(A_m[-1,-nage]) <- head(exp(-mx_m), -1)
  P_m <- matrix(0, nrow = nage, ncol = n+1)
  P_m[,1] <- pop_m
  for (i in 1:n+1) {
    P_m[,i] <- A_m%*%P_m[,i-1]
    P_m[1,i] <- P_f[1,i]*srb
  }
  
  projection <-
    expand_grid(n = 1:n, age_start = 1:nage-1) %>%
    mutate(Female = c(P_f[,-1]), Male = c(P_m[,-1])) %>%
    pivot_longer(
      cols = c(Female, Male),
      names_to = 'sex',
      values_to = 'population_midyear'
    )
  
  return(projection)
  
}

#' Make a Grid of Population Pyramids
PopPyramids <- function(
  dat, population, age, sex, year, highlight, facet, title
) {
  require(ggplot2); require(dplyr)
  dat %>%
    transmute(
      pop = ifelse(sex == 'Male', -{{population}}, {{population}}),
      age = {{age}}, sex = {{sex}}, year = {{year}}, hl = {{highlight}},
      fct = {{facet}}
    ) %>%
    ggplot(
      aes(x = age, y = pop, color = hl,
          group = interaction(sex, year))
    ) +
    geom_line(show.legend = FALSE) +
    annotate('text', x = 90, y = -Inf, label = '\u2642',
             hjust = -1, size = 10) +
    annotate('text', x = 90, y = Inf, label = '\u2640',
             hjust = 1, size = 10) +
    geom_hline(yintercept = 0) +
    scale_colour_manual(values = c('grey60', 'red')) +
    scale_x_continuous(
      breaks = seq(0, 100, 20),
      labels = function (x) ifelse(x == 100, '100+', x)
    ) +
    scale_y_continuous(
      labels = function(x) {formatC(abs(x*1e-3), format = 'd')}
    ) +
    coord_flip() +
    facet_wrap(~fct, scales = 'free_x') +
    labs(x = 'Age', y = 'Midyear population in 1000s')
}

# Data ------------------------------------------------------------

dat$skeleton <- readRDS(paths$input$skeleton)

dat$wpp <- readRDS(paths$input$wpp_population)
dat$hmdhfd <- readRDS(paths$input$hmdhfd_gb_projection)

# Harmonize WPP population ----------------------------------------

dat$wpp_clean <-
  dat$wpp %>%
  # select columns of interest
  select(
    population_source,
    region_code_wpp = LocID,
    iso_year = Time, age = AgeGrpStart,
    male = PopMale, female = PopFemale
  ) %>%
  # sex to long format
  pivot_longer(
    cols = c(female, male),
    names_to = 'sex',
    values_to = 'population_july1st'
  ) %>%
  # population is in factor 1000, convert to actual number
  mutate(
    population_midyear = population_july1st*1000
  ) %>%
  # ensure proper names of factor variables
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
  mutate(id = GenerateRowID(region_iso, sex, iso_year, age))

# check if data makes sense
fig$wpp_population <-
  PopPyramids(
    dat = dat$wpp_clean,
    population = population_midyear,
    age = age, sex = sex, year = iso_year, highlight = population_source,
    facet = region_iso
  ) +
  labs(title = 'WPP population estimates 2015-18 (grey) and population projections 2019-21 (red)') +
  fig_spec$MyGGplotTheme(scaler = 0.8)

# Harmonize pop data for UK regions -------------------------------

# for years not yet in the data we get population estimates via a
# Leslie-Matrix projection of midyear population assuming a stable
# population, i.e. no migration and constant fertility and mortality
# rates. we use the 'female dominant' projection of the male population.

# prepare a data set with required variables for
# Leslie matrix population projection for GB sub-regions
dat$gb_pop_estimates <-
  # skeleton
  expand_grid(
    year = 2015:2018,
    region_code_hmd = c('GBRTENW', 'GBR_NIR', 'GBR_SCO'),
    sex = c('Female', 'Male'),
    age_start = 0:110
  ) %>%
  # midyear population
  left_join(
    dat$hmdhfd$hmd_pop %>%
      pivot_longer(
        c(Female, Male),
        names_to = 'sex', values_to = 'population_midyear'
      ) %>%
      select(region_code_hmd, year = Year, age_start = Age,
             sex, population_midyear)
  ) %>%
  # death rates
  left_join(
    dat$hmdhfd$hmd_mort %>%
      pivot_longer(
        c(Female, Male),
        names_to = 'sex', values_to = 'death_rate'
      ) %>%
      select(region_code_hmd, year = Year, age_start = Age, sex, death_rate)
  ) %>%
  # female fertility rates
  left_join(
    dat$hmdhfd$hfd_fert %>%
      select(region_code_hmd, year = Year, age_start = Age,
             female_fertility_rate = ASFR)
  ) %>%
  mutate(
    # fertility rates are NA outside the age range [12, 55],
    # replace with 0
    female_fertility_rate =
      ifelse(is.na(female_fertility_rate), 0, female_fertility_rate),
    # same with death rates in highest ages
    death_rate =
      ifelse(is.na(death_rate), 0, death_rate)
  )

# perform the projections by sub-region
# start 2018 and project 2019, 20, 21
cnst$jumpoff_year <- 2018
dat$gb_pop_projections <-
  dat$gb_pop_estimates %>%
  filter(year == cnst$jumpoff_year) %>%
  group_by(region_code_hmd) %>%
  group_modify(~{
    
    male <- filter(.x, sex == 'Male')
    female <- filter(.x, sex == 'Female')
    
    projection <- StableProjection1x1(
      pop_m = male$population_midyear, pop_f = female$population_midyear,
      mx_m = male$death_rate, mx_f = female$death_rate,
      fx_f = female$female_fertility_rate, n = cnst$uk_projection_n_years, srb = 1.04
    )
    
    return(projection)
    
  }) %>%
  ungroup() %>%
  mutate(year = cnst$jumpoff_year+n)

# bind and harmonize estimates and projections
dat$gb_clean <-
  bind_rows(
    hmd_estimates = select(
      dat$gb_pop_estimates,
      region_code_hmd, year, sex, age_start, population_midyear
    ),
    hmd_projections = select(
      dat$gb_pop_projections,
      region_code_hmd, year, age_start, sex, population_midyear
    ),
    .id = 'population_source'
  ) %>%
  # make age 100 an open age group
  mutate(
    age_start = ifelse(age_start > 100, 100, age_start)
  ) %>%
  group_by(population_source, region_code_hmd, year, sex, age_start) %>%
  summarise(population_midyear = sum(population_midyear)) %>%
  ungroup() %>%
  # ensure proper names of factor variables
  mutate(
    sex =
      factor(sex, levels = names(cnst$code_sex_hmdhfd),
             labels = cnst$code_sex_hmdhfd) %>%
      as.character(),
    region_iso = factor(
      region_code_hmd,
      levels = cnst$region_lookup_hmd$region_code_hmd,
      labels = cnst$region_lookup_hmd$region_code_iso3166_2
    ) %>% as.character()
  ) %>%
  # add row id
  mutate(id = GenerateRowID(region_iso, sex, year, age_start))

# check if projection went well
fig$gb_population <-
  PopPyramids(
    dat = dat$gb_clean,
    population = population_midyear,
    age = age_start, sex = sex, year = year,
    highlight = population_source, facet = region_iso
  ) +
  labs(title = 'HMD population estimates 2015-18 (grey) and stable population projection 2019-21 (red) for U.K. regions') +
  fig_spec$MyGGplotTheme(scaler = 0.8)

# Join population with skeleton -----------------------------------

# join the different sources of population count data
# with the skeleton
dat$pop_joined <-
  dat$skeleton %>% 
  left_join(
    select(
      dat$wpp_clean,
      id,
      population_midyear_wpp = population_midyear,
      population_source_wpp = population_source
    ),
    by = 'id'
  ) %>%
  left_join(
    select(
      dat$gb_clean,
      id, population_midyear_hmd = population_midyear,
      population_source_hmd = population_source
    ),
    by = 'id'
  ) %>%
  mutate(
    # use wpp estimates unless they are missing (as for GB regions)
    population_midyear = case_when(
      !is.na(population_midyear_wpp) ~ round(population_midyear_wpp, 1),
      !is.na(population_midyear_hmd) ~ round(population_midyear_hmd, 1),
      TRUE ~ as.numeric(NA)
    ),
    population_source = case_when(
      !is.na(population_midyear_wpp) ~ population_source_wpp,
      !is.na(population_midyear_hmd) ~ population_source_hmd,
      TRUE ~ as.character(NA)
    )
  )

# populate sex category "total"
dat$row_female <- which(dat$pop_joined$sex == 'Female')
dat$row_male <- which(dat$pop_joined$sex == 'Male')
dat$row_total <- which(dat$pop_joined$sex == 'Total')

dat$pop_joined$population_midyear[dat$row_total] <-
  dat$pop_joined$population_midyear[dat$row_female] +
  dat$pop_joined$population_midyear[dat$row_male]

# select variables of interest for joined data
dat$joined <-
  dat$pop_joined %>%
  select(id, population_midyear, population_source)

# Export ----------------------------------------------------------

saveRDS(dat$joined, file = paths$output$harmonized_population)

fig_spec$ExportFigure(
  fig$wpp_population, path = paths$output$out, device = 'pdf',
  filename = '20-wpp_population', scale = 2
)
fig_spec$ExportFigure(
  fig$gb_population, path = paths$output$out, device = 'pdf',
  filename = '20-gb_population'
)
