# Download data on COVID-19 death counts from COVerAGE-DB

# Init ------------------------------------------------------------

library(glue)
library(readr); library(yaml)
library(httr); library(dplyr); library(tidyr)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  region_meta = './cfg/region_metadata.csv',
  url_covid = 'https://osf.io/7tnfh/download',
  coverage = './dat/coverage'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  coverage = './dat/coverage/coverage.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# meta data on regions
region_meta <- read_csv(paths$input$region_meta, na = '.')

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  # lookup table for region codes
  # only countries defined in skeleton
  region_lookup = 
    region_meta %>%
    filter(region_code_iso3166_2 %in% config$skeleton$region) %>%
    select(region_code_iso3166_2, region_code_coverage) %>%
    drop_na()
})

# list containers for analysis artifacts
dat <- list()

# Download --------------------------------------------------------

# download international COVID-19 death by age and sex counts
dat$coverage_zip <- GET(paths$input$url_covid, progress())

# Unzip and bind to single table ----------------------------------

# save downloaded zip to file
writeBin(
  object = content(dat$coverage_zip, 'raw'),
  con = glue('{paths$input$tmpdir}/coverage.zip')
)

# unzip coverage data and read table
dat$coverage <-
  unz(
    glue('{paths$input$tmpdir}/coverage.zip'),
    filename = 'Data/Output_5.csv'
  ) %>%
  read_csv(
    col_types = 'ccccciiddd',
    skip = 3
  )

# subset to data of interest
dat$coverage_sub <-
  dat$coverage %>%
  filter(
    Sex %in% c('f', 'm'),
    # only country level
    Country %in% cnst$region_lookup$region_code_coverage,
    Region == 'All'
  )

# Export ----------------------------------------------------------

saveRDS(dat$coverage_sub, file = paths$output$coverage)
