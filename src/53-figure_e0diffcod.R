# Generate Figure of annual e0 changes

# Init ------------------------------------------------------------

library(yaml); library(tidyverse); library(readr)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  region_metadata = './cfg/region_metadata.csv',
  figspecs = './cfg/figure_specification.R',
  lifetables = './out/40-codecomp.rds',
  e0avgdiff = './out/40-e0avgdiff.rds',
  arriaga_cntfc = './out/40-arriaga_cntfc.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  data = './dat/output_data.rds',
  fig_e0diff = './out',
  rds_e0diffT = './out/53-e0diffcodT.rds',
  csv_e0diffT = './tmp/53-e0diffcodT.csv',
  rds_e0diffF = './out/53-e0diffcodF.rds',
  csv_e0diffF = './tmp/53-e0diffcodF.csv',
  rds_e0diffM = './out/53-e0diffcodM.rds',
  csv_e0diffM = './tmp/53-e0diffcodM.csv'
)

# global configuration
config <- read_yaml(paths$input$config)

# meta data on regions
region_meta <- read_csv(paths$input$region_meta, na = '.')

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  regions_for_analysis = config$regions_for_covid_cause_analysis
})

# list containers for analysis artifacts
dat <- list()
fig <- list()

# Functions -------------------------------------------------------

# figure specifications
source(paths$input$figspec)

# Load decomposition results --------------------------------------

dat$arriaga_cntfc <- readRDS(paths$input$arriaga_cntfc)

# Create figure ---------------------------------------------------

name <- 'e0diffcod'
strata <- c('T', 'F', 'M')

fig <- map(strata, ~{
  
  const <- list(); const <- within(const, {
    age_breaks = c(0, 20, 40, 60, 80, Inf)
    age_names = c('0-19', '20-39', '40-59', '60-79', '80+')
    color_vline = 'grey50'
    font_plot = 'roboto'
    n_age = length(age_breaks)
    vertical_gap = 0.3
    fill_covid = '#f2c84b'
    fill_noncovid = 'grey40'
  })
  
  data <-
    dat$arriaga_cntfc %>%
    filter(
      sex == .x,
      year %in% 2020:2021,
      region_iso %in% cnst$regions_for_analysis
    ) %>%
    mutate(
      age = as.integer(age),
      age = as.character(cut(age, breaks = const$age_breaks,
                             labels = const$age_names,
                             right = TRUE, include.lowest = TRUE))
    ) %>%
    bind_rows(mutate(.,age = 'Total')) %>%
    group_by(region_iso, sex, age, year) %>%
    summarise(
      e0_cntrb_t_mean = sum(e0_cntrb_t_mean),
      e0_cntrb_t_covid_mean = sum(e0_cntrb_t_covid_mean),
      e0_cntrb_t_noncovid_mean = sum(e0_cntrb_t_noncovid_mean)
    ) %>%
    select(
      region_iso, sex, age, year, e0_cntrb_t_mean,
      e0_cntrb_t_covid_mean, e0_cntrb_t_noncovid_mean
    ) %>%
    pivot_wider(
      names_from = year,
      values_from =
        c(e0_cntrb_t_mean, e0_cntrb_t_covid_mean, e0_cntrb_t_noncovid_mean)
    ) %>%
    transmute(
      region_iso, sex, total = ifelse(age == 'Total', TRUE, FALSE),
      e0_measure_covid = e0_cntrb_t_covid_mean_2021,
      e0_measure_noncovid = e0_cntrb_t_noncovid_mean_2021
    ) %>%
    pivot_longer(
      cols = c(e0_measure_covid, e0_measure_noncovid)
    ) %>%
    ungroup() %>%
    mutate(
      name =
        case_when(
          name == 'e0_measure_covid' ~ 'COVID-19',
          name == 'e0_measure_noncovid' ~ 'Other'
        )
    ) %>%
    left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2'))
  
  plot <-
    data %>%
    ggplot(aes(y = age)) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 5.4, ymax = 6.6, fill = 'grey90') +
    geom_vline(
      xintercept = seq(-2.5, 0.5, 0.5)*12,
      color = '#FFFFFF', size = 0.2
    ) +
    geom_col(
      aes(fill = name, x = value*12, group = age),
      position = position_stack(), data = . %>% filter(!total),
      width = 0.4
    ) +
    geom_col(
      aes(fill = name, x = value*12, group = age),
      position = position_stack(), data = . %>% filter(total)
    ) +
    geom_vline(xintercept = 0, color = 'grey50') +
    facet_wrap(~region_name) +
    scale_fill_manual(values = c(
      'COVID-19' = const$fill_covid, 'Other' = const$fill_noncovid
    )) +
    scale_x_continuous(breaks = seq(-2.5, 0.5, 0.5)*12,
                       labels = c('', '-24', '', '-12', '', '0', '+6')) +
    fig_spec$MyGGplotTheme(axis = '', grid = '', panel_border = F) +
    theme(
      panel.background = element_rect(fill = 'grey95', color = 'grey95'),
      #strip.background = element_rect(fill = 'grey90', color = 'grey90'),
      legend.position = c(0.08,0.1)
    ) +
    labs(
      x = 'Contributions by age and cause of death to months of life expectancy deficit in 2021',
      y = 'Age group',
      fill = 'Cause of death'
    )
  
  list(cnst = const, data = data, plot = plot)
  
})
names(fig) <- paste0(name, strata)


# Just the US -----------------------------------------------------


const <- list(); const <- within(const, {
  age_breaks = c(0, 20, 40, 60, 80, Inf)
  age_names = c('0-19', '20-39', '40-59', '60-79', '80+')
  color_vline = 'grey50'
  font_plot = 'roboto'
  n_age = length(age_breaks)
  vertical_gap = 0.3
  fill_covid = '#f2c84b'
  fill_noncovid = 'grey40'
})
lt <- readRDS('./out/40-lifetables.rds')
lt_plot <-
  lt %>%
  filter(
    sex == 'T',
    year %in% 2020:2021,
    region_iso %in% 'US',
    projected == 'actual'
  ) %>%
  mutate(
    age = as.integer(age),
    age = as.character(cut(age, breaks = const$age_breaks,
                           labels = const$age_names,
                           right = TRUE, include.lowest = TRUE))
  ) %>%
  bind_rows(mutate(.,age = 'Total')) %>%
  group_by(region_iso, sex, age, year) %>%
  summarise(
    e0_cntrb_t_mean = sum(e0_cntrb_t_mean),
    e0_cntrb_t_covid_mean = sum(e0_cntrb_t_covid_mean),
    e0_cntrb_t_noncovid_mean = sum(e0_cntrb_t_noncovid_mean)
  ) %>%
  select(
    region_iso, sex, age, year, e0_cntrb_t_mean,
    e0_cntrb_t_covid_mean, e0_cntrb_t_noncovid_mean
  ) %>%
  pivot_wider(
    names_from = year,
    values_from =
      c(e0_cntrb_t_mean, e0_cntrb_t_covid_mean, e0_cntrb_t_noncovid_mean)
  ) %>%
  transmute(
    region_iso, sex, total = ifelse(age == 'Total', TRUE, FALSE),
    e0_measure_covid = e0_cntrb_t_covid_mean_2021 + e0_cntrb_t_covid_mean_2020,
    e0_measure_noncovid = e0_cntrb_t_noncovid_mean_2021 + e0_cntrb_t_noncovid_mean_2020
  ) %>%
  pivot_longer(
    cols = c(e0_measure_covid, e0_measure_noncovid)
  ) %>%
  ungroup() %>%
  mutate(
    name =
      case_when(
        name == 'e0_measure_covid' ~ 'COVID-19',
        name == 'e0_measure_noncovid' ~ 'Other'
      )
  ) %>%
  left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2'))

fig_us <-
  lt_plot %>%
  ggplot(aes(y = age)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 5.4, ymax = 6.6, fill = 'grey90') +
  geom_vline(
    xintercept = seq(-2.5, 0.5, 0.5)*12,
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = name, x = value*12, group = age),
    position = position_stack(), data = . %>% filter(!total),
    width = 0.4
  ) +
  geom_col(
    aes(fill = name, x = value*12, group = age),
    position = position_stack(), data = . %>% filter(total)
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  facet_wrap(~region_name) +
  scale_fill_manual(values = c(
    'COVID-19' = const$fill_covid, 'Other' = const$fill_noncovid
  )) +
  scale_x_continuous(breaks = seq(-2.5, 0.5, 0.5)*12,
                     labels = c('', '-24', '', '-12', '', '0', '+6')) +
  fig_spec$MyGGplotTheme(axis = '', grid = '', panel_border = F) +
  theme(
    panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    #strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = c(0.08,0.1)
  ) +
  labs(
    x = 'Months of LE change since 2019',
    y = 'Age group',
    fill = 'Cause of death',
    title = 'U.S. life expectancy change since 2019 contributed by age and cause of death',
    #subtitle = 'Life expectancy changes are explained by mortality increases or decreases in certain age groups and causes of death',
    caption = 'Data: HMD-STMF, WHO-WPP, COVerAGE-DB, own calculations'
  )
fig_spec$ExportFigure(
  fig_us, device = 'pdf',
  filename = '53-us',
  path = paths$output$fig_e0diff,
  width = fig_spec$width*0.8, height = fig_spec$width*0.7
)

# Export ----------------------------------------------------------

fig_spec$ExportFigure(
  fig$e0diffcodT$plot, device = 'pdf',
  filename = '53-e0diffcodT',
  path = paths$output$fig_e0diff,
  width = fig_spec$width, height = fig_spec$width
)
saveRDS(fig$e0diffcodT, file = paths$output$rds_e0diffT)
fig$e0diffcodT$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0diffT)

fig_spec$ExportFigure(
  fig$e0diffcodF$plot, device = 'pdf',
  filename = '53-e0diffcodF',
  path = paths$output$fig_e0diff,
  width = fig_spec$width, height = fig_spec$width
)
saveRDS(fig$e0diffcodF, file = paths$output$rds_e0diffF)
fig$e0diffcodF$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0diffF)

fig_spec$ExportFigure(
  fig$e0diffcodM$plot, device = 'pdf',
  filename = '53-e0diffcodM',
  path = paths$output$fig_e0diff,
  width = fig_spec$width, height = fig_spec$width
)
saveRDS(fig$e0diffcodM, file = paths$output$rds_e0diffM)
fig$e0diffcodM$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0diffM)
