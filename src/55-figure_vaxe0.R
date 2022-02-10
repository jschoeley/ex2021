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
  vaccination_old = './tmp/rates_v2.rds',
  vaccination_young = './tmp/rates_v2_under60.rds',
  arriaga_cntf = './out/arriaga_cntfc.rds',
  figspecs = './cfg/figure_specification.R'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  vaxe0_rds = './out/vaxe0.rds',
  vaxe0_fig = './out'
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
fig <- list()

# Function --------------------------------------------------------

source(paths$input$figspecs)

# Load data -------------------------------------------------------

dat$arriaga_cntf <- readRDS(paths$input$arriaga_cntf)
dat$vaccination_old <- readRDS(paths$input$vaccination_old)
dat$vaccination_young <- readRDS(paths$input$vaccination_young)

# Prepare vaccination data ----------------------------------------

dat$vaccination_sub <-
  bind_rows(
    `60+` = dat$vaccination_old,
    `<60` = dat$vaccination_young,
    .id = 'age'
  ) %>%
  mutate(
    date = as.Date(Date),
    region_iso = factor(
      Country,
      region_meta$region_code_coverage[region_meta$region_name_short != 'RUS'],
      region_meta$region_code_iso3166_2[region_meta$region_name_short != 'RUS']
    ) %>% as.character()
  ) %>%
  select(age, region_iso, date, Measure, rate) %>%
  pivot_wider(names_from = Measure, values_from = rate) %>%
  filter(lubridate::year(date) == 2021) %>%
  complete(age, region_iso, date = seq.Date(as.Date('2021-01-01'), as.Date('2021-12-31'), by = 'day')) %>%
  group_by(region_iso) %>%
  mutate(
    across(
      .cols = c(Vaccination1, Vaccination2, Vaccination3, VaccinationBooster),
      .fns = zoo::na.approx, rule = 2, maxgap = 100
    )
  ) %>%
  ungroup() %>%
  arrange(region_iso, date)

dat$vaccination_sub %>%
  filter(age == '60+') %>%
  ggplot() +
  geom_point(aes(x = date, y = Vaccination2), size = 0.1) +
  facet_wrap(~region_iso)

# Prepare e0 data -------------------------------------------------

dat$e0_diff_summary <-
  dat$arriaga_cntf %>%
  mutate(age = cut(as.integer(age), c(0, 60, Inf),
                   right = FALSE, include.lowest = TRUE,
                   labels = c('<60', '60+'))) %>%
  group_by(region_iso, sex, age, year) %>%
  summarise(
    e0_cntrb_t_mean = sum(e0_cntrb_t_mean)
  ) %>%
  ungroup()

# Add vaccination uptake measures ---------------------------------

dat$vaccination_summary <-
  dat$vaccination_sub %>%
  group_by(region_iso, age) %>%
  summarise(
    # share fully vaccinated by july 1st
    vax2_jul = Vaccination2[date == '2021-07-01'],
    # share fully vaccinated by oct 1st
    vax2_oct = Vaccination2[date == '2021-10-01'],
    # integral vaccination measure
    vax2_int = mean(Vaccination2)
  ) %>%
  ungroup()

# Join vaccination and e0 data ------------------------------------

dat$vaxe0 <-
  dat$e0_diff_summary %>%
  select(region_iso, year, sex, age, e0_cntrb_t_mean) %>%
  left_join(dat$vaccination_summary) %>%
  filter(year == 2021, sex == 'T',
         region_iso %in% cnst$regions_for_vax_analysis)

# Plot vaccination efficacy ---------------------------------------

fig$vaxe0 <- list()
fig$vaxe0$data <-
  dat$vaxe0 %>%
  mutate(
    #vax_measure = vax2_int,
    vax_measure = vax2_oct,
    ex_measure = e0_cntrb_t_mean
  ) %>%
  filter(is.finite(e0_cntrb_t_mean)) %>%
  left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2'))

fig$vaxe0$plot <-
  fig$vaxe0$data %>%
  ggplot(aes(x = vax_measure, y = -ex_measure, group = age)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = 'lm', se = FALSE, color = 'grey', size = 0.75) +
  geom_text_repel(aes(label = region_name_short),
                  family = 'robotocondensed', size = 3, color = 'grey50') +
  geom_point(aes(fill = age), shape = 21, color = 'black') +
  #stat_poly_eq() +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(`<60` = 'white', `60+` = 'black')) +
  coord_cartesian(xlim = c(0.10, 1), ylim = c(-1.1, 3), expand = c(0,0)) +
  fig_spec$MyGGplotTheme(grid = 'xy', axis = '') +
  labs(
    y = 'Years of Life expectancy deficit\nin 2021 contributed by age group',
    x = '% fully vaccinated in age group by Oct 1st 2021'
  )

fig$vaxe0$plot

fig_spec$ExportFigure(
  fig$vaxe0$plot, device = 'pdf',
  filename = 'vaxe0',
  path = paths$output$vaxe0_fig,
  width = fig_spec$width, height = 0.7*fig_spec$width,
  scale = 1
)