# Harmonize the cumulative COVID-19 death counts for STMF countries

# Init ------------------------------------------------------------

library(here); library(glue)
library(yaml); library(readr)
library(dplyr); library(tidyr); library(lubridate)
library(ggplot2)
library(ungroup)

# Constants -------------------------------------------------------

wd <- here()

config <- read_yaml(glue('{wd}/cfg/config.yaml'))
region_meta <- read_csv(glue('{wd}/cfg/region_metadata.csv'), na = '.')
source(glue('{wd}/cfg/figure_specification.R'))

cnst <- list()
cnst <- within(cnst, {
  # path to global objects
  path_global = glue('{wd}/src/00-global.R')
  # where to put the harmonized data
  path_harmonized = glue('{wd}/tmp')
  # where to put the figures
  path_fig = glue('{wd}/tmp')
  # where to find the coverage covid data
  path_coverage = glue('{wd}/dat/coverage')
  # skeleton path
  path_skeleton = glue('{wd}/tmp/10-harmonized_skeleton.rds')
  # translation of COVerAGE-DB sex code to harmonized sex code
  code_sex_coverage =
    c(`m` = config$skeleton$sex$Male,
      `f` = config$skeleton$sex$Female)
  # lookup table for region codes
  # only countries defined in skeleton
  region_lookup = 
    region_meta %>%
    filter(region_code_iso3166_2 %in% config$skeleton$region) %>%
    select(region_code_iso3166_2, region_code_coverage) %>%
    drop_na()
  # pclm life-table closeout age
  pclm_highest_age = 110
})

dat <- list()
fig <- list()

# Functions -------------------------------------------------------

source(cnst$path_global)

# Data ------------------------------------------------------------

dat$skeleton <- readRDS(cnst$path_skeleton)

dat$coverage_raw <- readRDS(glue('{cnst$path_coverage}/coverage.rds'))

# Prepare data for age harmonization ------------------------------

dat$coverage_cleaned <-
  dat$coverage_raw %>%
  # select columns of interest
  select(
    region_code_coverage = Country, date = Date, sex = Sex,
    age_start = Age, age_width = AgeInt, death_covid = Deaths
  ) %>%
  # ensure proper names of factor variables
  mutate(
    sex =
      factor(sex, levels = names(cnst$code_sex_coverage),
             labels = cnst$code_sex_coverage) %>%
      as.character(),
    region_iso = factor(
      region_code_coverage,
      levels = cnst$region_lookup$region_code_coverage,
      labels = cnst$region_lookup$region_code_iso3166_2
    ) %>% as.character(),
    date =
      as_date(date, format = '%d.%m.%Y'),
    year =
      year(date),
    age_width =
      as.numeric(ifelse(age_start == 100, 'Inf', age_width))
  ) %>%
  # we're only interested in 2020 and later
  filter(
    year >= 2020
  ) %>%
  select(
    region_iso, year, date, sex, age_start, age_width, death_covid
  ) %>% # calculate covid deaths per year (de-cumulate)
  drop_na(death_covid) %>%
  group_by(region_iso, year, sex) %>%
  mutate(date_range = case_when(date == min(date) ~ 'min', date == max(date) ~ 'max', TRUE ~ 'NA')) %>%
  filter(date_range %in% c('min', 'max')) %>%
  pivot_wider(names_from = date_range, values_from = c(death_covid, date)) %>%
  ungroup() %>%
  mutate(
    death_covid = death_covid_max - death_covid_min
  ) %>%
  select(region_iso, year, date = date_max, sex, age_start, age_width, death_covid) %>%
  # truncate eventual negative artifacts to 0
  mutate(death_covid = ifelse(death_covid < 0, 0, death_covid))

# Age harmonization -----------------------------------------------

# for each country and sex in a year
# - extract some data quality information indicating the completeness
# of the covid death information
# - extract the cumulative number of covid deaths at the latest
# available date by age group
# - ungroup into single year age groups ages 0 to 110
# - sum up deaths across ages 100 to 110
dat$coverage_ungrouped <-
  dat$coverage_cleaned %>%
  group_by(region_iso, year, sex) %>% 
  group_modify(~{
    cat('PCLM on', .y$region_iso, .y$sex, '\n')
    
    # slice to latest date where data is
    # available for all age groups
    most_recent_complete_data <-
      .x %>%
      group_by(date) %>%
      mutate(any_na = anyNA(death_covid)) %>%
      filter(!any_na) %>%
      ungroup() %>%
      filter(date == max(date))
    # flag indicating that there is not a single complete
    # age profile of cumulative covid deaths for this country-sex
    no_data = nrow(most_recent_complete_data) == 0
    
    if (!no_data) {

      # number of age groups
      n_agegroups_raw <- length(most_recent_complete_data$death_covid)
      # width of the last age group
      nlast <- cnst$pclm_highest_age -
        tail(most_recent_complete_data$age_start, 1)+1
            
      fit_pclm <- pclm(
        x = most_recent_complete_data$age_start,
        y = most_recent_complete_data$death_covid,
        nlast = nlast, out.step = 1
        #control = list(lambda = 0.1)
      )
      single_age_death <- round(fitted.values(fit_pclm), 1)
      pclm_lambda <- fit_pclm$smoothPar[1]
      the_date <- most_recent_complete_data$date[1]
      
    } else {
      n_agegroups_raw <- NA
      single_age_death <- NA
      pclm_lambda <- NA
      the_date <- NA
    }
    
    ungrouped_deaths <- tibble(
      age_start = 0:cnst$pclm_highest_age,
      death_covid = single_age_death,
      death_covid_lambda = pclm_lambda,
      death_covid_date = the_date,
      death_covid_nageraw = n_agegroups_raw
    )
    
    return(ungrouped_deaths)
  }) %>%
  ungroup()

# PCLM returns NaN if there are convergence problems. replace with NA
dat$coverage_ungrouped <-
  dat$coverage_ungrouped %>%
  mutate(
    death_covid = ifelse(is.nan(death_covid), NA, death_covid)
  )

# Summarize all ages past 100 in single age group -----------------

dat$coverage_on_skeleton <-
  dat$coverage_ungrouped %>%
  mutate(age_start = ifelse(age_start >= 100, 100, age_start)) %>%
  group_by(region_iso, year, sex, age_start) %>%
  summarise(
    death_covid = sum(death_covid),
    death_covid_date = death_covid_date[1],
    death_covid_nageraw = death_covid_nageraw[1]
  ) %>%
  ungroup() %>%
  # add id
  mutate(
    id = GenerateRowID(region_iso, sex, year, age_start)
  ) %>%
  select(id, death_covid, death_covid_date, death_covid_nageraw)
dat$coverage_on_skeleton <-
  dat$skeleton %>%
  left_join(dat$coverage_on_skeleton, by = 'id')

# Populate sex category "Total" -----------------------------------

dat$row_female <- which(dat$coverage_on_skeleton$sex == 'Female')
dat$row_male <- which(dat$coverage_on_skeleton$sex == 'Male')
dat$row_total <- which(dat$coverage_on_skeleton$sex == 'Total')

dat$coverage_on_skeleton$death_covid[dat$row_total] <-
  dat$coverage_on_skeleton$death_covid[dat$row_female] +
  dat$coverage_on_skeleton$death_covid[dat$row_male]

# Finalize --------------------------------------------------------

dat$covid <-
  dat$coverage_on_skeleton %>%
  # set COVID deaths pre 2020 to 0
  mutate(
    death_covid = ifelse(year >= 2020, death_covid, 0)
  ) %>%
  select(id, death_covid, death_covid_date, death_covid_nageraw)
  
# Diagnostic plots ------------------------------------------------

fig$covid_raw_data <-
  dat$coverage_cleaned %>%
  group_by(region_iso, date, sex) %>%
  summarise(death_covid = sum(death_covid)) %>%
  ggplot(aes(x = date, y = death_covid, color = sex)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = as_date('2020-12-31')) +
  scale_x_date(date_labels = '%b') +
  scale_color_manual(values = fig_spec$sex_colors) +
  facet_wrap(~region_iso, scales = 'free_y') +
  fig_spec$MyGGplotTheme(panel_border = TRUE, grid = 'xy') +
  labs(
    title = 'Cumulative COVID deaths, 2020',
    x = 'Age', y = 'Cumulative COVID deaths'
  )

dat$coverage_most_recent_date <-
  dat$coverage_ungrouped %>%
  group_by(region_iso) %>%
  summarise(min_date = max(death_covid_date))

fig$covid_pclm <-
  dat$coverage_ungrouped %>%
  ggplot() +
  geom_line(aes(x = age_start, y = death_covid,
                color = sex, lty = as.factor(year))) +
  geom_label(
    aes(
      x = 0, y = 0,
      label = min_date
    ),
    hjust = 0, vjust = 0,
    label.size = 0, alpha = 0.6,
    data = dat$coverage_most_recent_date
  ) +
  scale_color_manual(values = fig_spec$sex_colors) +
  facet_wrap(~region_iso, scales = 'free_y') +
  fig_spec$MyGGplotTheme(panel_border = TRUE, grid = 'xy') +
  labs(
    title = 'Ungrouped covid deaths as of <date> by country and age',
    subtitle = 'Starting date of data collection unknown in all cases',
    x = 'Age', y = 'COVID deaths'
  )

# Export ----------------------------------------------------------

saveRDS(
  dat$covid,
  file = glue('{cnst$path_harmonized}/22-harmonized_covid.rds')
)

fig_spec$ExportFigure(
  fig$covid_pclm, path = cnst$path_fig, device = 'pdf',
  filename = '22-covid_pclm', scale = 2
)
fig_spec$ExportFigure(
  fig$covid_raw_data, path = cnst$path_fig, device = 'pdf',
  filename = '22-covid_raw_data', scale = 2
)
