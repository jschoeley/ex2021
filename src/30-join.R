# Assemble final data set

# Init ------------------------------------------------------------

library(glue)
library(dplyr); library(readr); library(openxlsx)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  skeleton = './tmp/harmonized_skeleton.rds',
  global = './src/00-global.R',
  population = './tmp/harmonized_population.rds',
  death = './tmp/harmonized_death.rds',
  covid = './tmp/harmonized_covid.rds',
  external_lifetables = './tmp/harmonized_external_lifetables.rds',
  counterfactual_mortality = './tmp/harmonized_counterfactual_mortality.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  out_rds = './out/lt_input.rds',
  out_csv = './out/lt_input.csv',
  out_xlsx = './out/lt_input.xlsx'
)

# list containers for analysis artifacts
dat <- list()

# Functions -------------------------------------------------------

source(paths$input$global)

# Data ------------------------------------------------------------

dat$skeleton <- readRDS(paths$input$skeleton)
dat$population <- readRDS(paths$input$population)
dat$death <- readRDS(paths$input$death)
dat$covid <- readRDS(paths$input$covid)
dat$external_lifetables <- readRDS(paths$input$external_lifetables)
dat$counterfactual_mortality <- readRDS(paths$input$counterfactual_mortality)

# Join ------------------------------------------------------------

dat$lt_input <-
  dat$skeleton %>% 
  mutate(nweeks_year = ifelse(YearHasIsoWeek53(year), 53L, 52L)) %>%
  arrange(sex, region_iso, year, age_start) %>%
  left_join(dat$death, by = 'id') %>%
  left_join(dat$population, by = 'id') %>%
  left_join(dat$covid, by = 'id') %>%
  left_join(dat$external_lifetables, by = 'id') %>%
  left_join(dat$counterfactual_mortality, by = 'id')

# Export ----------------------------------------------------------

saveRDS(dat$lt_input, file = paths$output$out_rds)

write_csv(dat$lt_input, file = paths$output$out_csv)

write.xlsx(dat$lt_input, file = paths$output$out_xlsx,
           keepNA = TRUE, na.string = '.',
           firstRow = TRUE, firstCol = TRUE, overwrite = TRUE)
