# Generate Figure of annual e0 forecasts

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
  csv_e0forecastT = './tmp/90-e0forecastT.csv',
  rds_e0forecastF = './out/90-e0forecastF.rds',
  csv_e0forecastF = './tmp/90-e0forecastF.csv',
  rds_e0forecastM = './out/90-e0forecastM.rds',
  csv_e0forecastM = './tmp/90-e0forecastM.csv'
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

dat$lifetables <- readRDS(paths$input$lifetables)
dat$cntfc_lt_debug <- readRDS(paths$input$cntfc_lt_debug)
dat$skeleton <- readRDS(paths$input$skeleton)

# Plot ------------------------------------------------------------

name <- 'e0forecast'
strata <- c('T', 'F', 'M')

fig$e0forecast <- map(strata, ~{
  data <- 
    dat$lifetables %>%
    filter(age == 0, sex == .x, region_iso %in% cnst$regions_for_analysis,
           quarter == 'annual') %>%
    left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2')) %>%
    select(year, region_name, sex, age, projected, ex_mean) %>%
    group_by(region_name) %>%
    mutate(ex_diff = ex_mean-ifelse(is.na(ex_mean[1]), ex_mean[2], ex_mean[1])) %>%
    ungroup()
  
  plot <-
    data %>%
    filter(
      (projected == 'projected' & year > 2019) | (projected == 'actual')
    ) %>%
    ggplot(aes(x = year)) +
    geom_hline(aes(yintercept = 0), color = '#666666') +
    geom_point(aes(y = ex_diff, fill = projected), shape = 21) +
    scale_fill_manual(values = c(projected = 'white', actual = 'black')) +
    facet_wrap(~region_name) +
    fig_spec$MyGGplotTheme() +
    theme(
      panel.background = element_rect(fill = 'grey95', color = 'grey95'),
      legend.position = c(0.92,0.1)
    ) +
    labs(y = 'Projected and observed change in life expectancy since 2015', x = NULL,
         fill = NULL)
  
  list(data = data, plot = plot)
})
names(fig$e0forecast) <- paste0(name, strata)

fig$e0forecast$e0forecastT$plot

# Export ----------------------------------------------------------

fig_spec$ExportFigure(
  fig$e0forecast$e0forecastT$plot, device = 'pdf',
  filename = '90-e0forecastT',
  path = paths$output$fig_e0forecast,
  width = fig_spec$width, height = 0.9*fig_spec$width
)
saveRDS(fig$e0forecast$e0forecastT, file = paths$output$rds_e0forecastT)
fig$e0forecast$e0forecastT$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0forecastT)

fig_spec$ExportFigure(
  fig$e0forecast$e0forecastF$plot, device = 'pdf',
  filename = '90-e0forecastF',
  path = paths$output$fig_e0forecast,
  width = fig_spec$width, height = 0.9*fig_spec$width
)
saveRDS(fig$e0forecast$e0forecastF, file = paths$output$rds_e0forecastF)
fig$e0forecast$e0forecastF$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0forecastF)

fig_spec$ExportFigure(
  fig$e0forecast$e0forecastM$plot, device = 'pdf',
  filename = '90-e0forecastM',
  path = paths$output$fig_e0forecast,
  width = fig_spec$width, height = 0.9*fig_spec$width
)
saveRDS(fig$e0forecast$e0forecastM, file = paths$output$rds_e0forecastM)
fig$e0forecast$e0forecastM$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0forecastM)
