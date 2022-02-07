# Calculate vaccination uptake statistics by age

# Init ------------------------------------------------------------

library(glue); library(yaml)
library(dplyr); library(tidyr); library(readr)
library(ggplot2); library(ggrepel); library(ggpmisc)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  region_meta = './cfg/region_metadata.csv',
  vaccination = './dat/coverage/rates_v2.rds',
  arriaga_cntf = './out/arriaga_cntfc.rds',
  figspecs = './cfg/figure_specification.R'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  vaxe0_rds = './out/vaxe0.rds',
  vaxe0_fig = './out/vaxe0.pdf'
)

# global configuration
config <- read_yaml(paths$input$config)

# meta data on regions
region_meta <- read_csv(paths$input$region_meta, na = '.')

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  regions_for_vax_analysis = config$regions_for_vax_analysis}
)

tmp <- list()
dat <- list()

# Function --------------------------------------------------------

source(paths$input$figspecs)

# Load data -------------------------------------------------------

dat$arriaga_cntf <- readRDS(paths$input$arriaga_cntf)
dat$vaccination <- readRDS(paths$input$vaccination)

# Prepare vaccination data ----------------------------------------

dat$vaccination_sub <-
  dat$vaccination %>%
  mutate(
    date = as.Date(Date),
    region_iso = factor(
      Country,
      region_meta$region_code_coverage[region_meta$region_name_short != 'RUS'],
      region_meta$region_code_iso3166_2[region_meta$region_name_short != 'RUS']
    ) %>% as.character()
  ) %>%
  filter(Measure == 'Vaccination2', lubridate::year(date) == 2021) %>%
  select(region_iso, date = Date, rate) %>%
  complete(region_iso, date = seq.Date(as.Date('2021-01-01'), as.Date('2021-12-31'), by = 'day')) %>%
  group_by(region_iso) %>%
  mutate(rate = zoo::na.approx(rate, rule = 2)) %>%
  ungroup() %>%
  arrange(region_iso, date)

dat$vaccination_sub %>%
  ggplot() +
  geom_point(aes(x = date, y = rate), size = 0.1) +
  facet_wrap(~region_iso)

# Add vaccination uptake measures ---------------------------------

dat$vaccination_summary <-
  dat$vaccination_sub %>%
  group_by(region_iso) %>%
  summarise(
    # share fully vaccinated by july 1st
    vax_jul = rate[date == '2021-07-01'],
    # share fully vaccinated by oct 1st
    vax_oct = rate[date == '2021-10-01'],
    # integral vaccination measure
    vax_int = mean(rate)
  )

# Join vaccination and e0 data ------------------------------------

dat$arriaga_cntf_sub <-
  dat$arriaga_cntf %>%
  filter(age == 0, year == 2021, sex == 'T',
         region_iso %in% cnst$regions_for_vax_analysis) %>%
  select(region_iso, ex_actual_minus_expected_mean)

# Plot vaccination efficacy ---------------------------------------

vaxe0 <- list()
vaxe0$data <-
  left_join(
    dat$vaccination_summary,
    dat$arriaga_cntf_sub
  ) %>%
  mutate(
    #vax_measure = vax_int,
    vax_measure = vax_oct,
    ex_measure = ex_actual_minus_expected_mean
  ) %>%
  filter(is.finite(ex_actual_minus_expected_mean))


vaxe0$data %>%
  ggplot(aes(x = vax_measure, y = ex_measure)) +
  geom_point() +
  geom_text_repel(aes(label = region_iso)) +
  geom_smooth(method = 'lm', se = FALSE) +
  stat_poly_eq() +
  fig_spec$MyGGplotTheme()
