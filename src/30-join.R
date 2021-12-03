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
  ex = './tmp/harmonized_ex.rds'
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
dat$ex <- readRDS(paths$input$ex)

# Join ------------------------------------------------------------

dat$lt_input <-
  dat$skeleton %>% 
  mutate(nweeks_year = ifelse(YearHasIsoWeek53(year), 53L, 52L)) %>%
  arrange(sex, region_iso, year, age_start) %>%
  left_join(dat$death, by = 'id') %>%
  left_join(dat$population, by = 'id') %>%
  left_join(dat$covid, by = 'id') %>%
  left_join(dat$ex, by = 'id')

# Create sex category "total" -------------------------------------

dat$row_female <- which(dat$lt_input$sex == 'Female')
dat$row_male <- which(dat$lt_input$sex == 'Male')
dat$row_total <- which(dat$lt_input$sex == 'Total')

dat$lt_input$death_total[dat$row_total] <-
  dat$lt_input$death_total[dat$row_female] +
  dat$lt_input$death_total[dat$row_male]
dat$lt_input$population_py[dat$row_total] <-
  dat$lt_input$population_py[dat$row_female] +
  dat$lt_input$population_py[dat$row_male]
dat$lt_input$population_midyear[dat$row_total] <-
  dat$lt_input$population_midyear[dat$row_female] +
  dat$lt_input$population_midyear[dat$row_male]
dat$lt_input$death_covid[dat$row_total] <-
  dat$lt_input$death_covid[dat$row_female] +
  dat$lt_input$death_covid[dat$row_male]

# Export ----------------------------------------------------------

saveRDS(dat$lt_input, file = paths$output$out_rds)

write_csv(dat$lt_input, file = paths$output$out_csv)

write.xlsx(dat$lt_input, file = paths$output$out_xlsx,
           keepNA = TRUE, na.string = '.',
           firstRow = TRUE, firstCol = TRUE, overwrite = TRUE)
