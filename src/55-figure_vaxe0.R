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
  arriaga_cntf = './out/40-arriaga_cntfc.rds',
  figspecs = './cfg/figure_specification.R'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  vaxe0_rds = './out/55-vaxe0.rds',
  vaxe0_csv = './tmp/55-vaxe0.csv',
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
      .fns = zoo::na.approx, rule = 2, yleft = 0, maxgap = 100
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
  group_by(region_iso, sex, age, year, quarter) %>%
  summarise(
    e0_cntrb_t_mean = sum(e0_cntrb_t_mean)
  ) %>%
  ungroup()

# Add vaccination uptake measures ---------------------------------

dat$vaccination_summary <-
  dat$vaccination_sub %>%
  group_by(region_iso, age) %>%
  summarise(
    # share fully vaccinated by jan 1st
    Q1 = Vaccination2[date == '2021-01-01'],
    # share fully vaccinated by apr 1st
    Q2 = Vaccination2[date == '2021-04-01'],
    # share fully vaccinated by jul 1st
    Q3 = Vaccination2[date == '2021-07-01'],
    # share fully vaccinated by oct 1st
    Q4 = Vaccination2[date == '2021-10-01'],
    # integral vaccination measure
    annual = mean(Vaccination2)
  ) %>%
  ungroup() %>%
  pivot_longer(
    c(Q1, Q2, Q3, Q4, annual),
    names_to = 'quarter',
    values_to = 'p_vax'
  )

# Join vaccination and e0 data ------------------------------------

dat$vaxe0 <-
  dat$e0_diff_summary %>%
  select(region_iso, year, sex, age, quarter, e0_cntrb_t_mean) %>%
  left_join(dat$vaccination_summary) %>%
  filter(year == 2021, sex == 'T',
         region_iso %in% cnst$regions_for_vax_analysis)

# Plot vaccination efficacy ---------------------------------------

fig$vaxe0 <- list()
fig$vaxe0$data <-
  dat$vaxe0 %>%
  mutate(
    vax_measure = p_vax,
    ex_measure = e0_cntrb_t_mean
  ) %>%
  filter(is.finite(e0_cntrb_t_mean), quarter == 'Q4') %>%
  left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2'))

fig$vaxe0$plot <-
  fig$vaxe0$data %>%
  ggplot(aes(x = vax_measure, y = -ex_measure, group = age)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = 'lm', se = FALSE, color = 'grey', size = 0.75) +
  geom_text_repel(aes(label = region_name_short),
                  family = 'robotocondensed', size = 2, color = 'grey50') +
  geom_point(color = 'black', size = 1) +
  stat_poly_eq() +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.2)) +
  facet_wrap(~age) +
  fig_spec$MyGGplotTheme(grid = 'xy', axis = '') +
  labs(
    y = 'Years of life expectancy deficit\nin Q4 2021 contributed by age group',
    x = '% fully vaccinated in age group by Oct 1st 2021'
  )
fig$vaxe0$plot

# Export ----------------------------------------------------------

fig_spec$ExportFigure(
  fig$vaxe0$plot, device = 'pdf',
  filename = '55-vaxe0',
  path = paths$output$vaxe0_fig,
  width = fig_spec$width, height = 0.5*fig_spec$width,
  scale = 1
)
saveRDS(fig$vaxe0, file = paths$output$vaxe0_rds)
fig$vaxe0$data %>%
  select(region_code_iso3166_2_alpha3, age, vax_measure, ex_measure) %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$vaxe0_csv)
