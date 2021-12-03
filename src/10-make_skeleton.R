# Create data base skeleton
#
# Here we define the "skeleton" of the data base used
# for analysis. It's a definition of years, ages, sexes, and
# regions that we wish to acquire data for.

# Init ------------------------------------------------------------

library(yaml)
library(dplyr); library(tidyr); library(glue)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  global = './src/00-global.R',
  config = './cfg/config.yaml'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  harmonized_skeleton = glue('{paths$input$tmp}/harmonized_skeleton.rds')
)

# global configuration
config <- read_yaml(paths$input$config)

# list containers for analysis artifacts
dat <- list()

# Functions -------------------------------------------------------

source(paths$input$global)

# Generate skeleton -----------------------------------------------

dat$skeleton <-
  expand_grid(
    region_iso =
      config$skeleton$region,
    sex =
      unlist(config$skeleton$sex),
    year =
      seq(config$skeleton$year$min, config$skeleton$year$max, 1) %>%
      as.integer(),
    tibble(
      age_start = seq(config$skeleton$age$min, config$skeleton$age$max, 1),
      age_width = c(diff(age_start), Inf)
    )
  )

# Add unique row id -----------------------------------------------

dat$skeleton <-
  dat$skeleton %>%
  mutate(
    id = GenerateRowID(region_iso, sex, year, age_start)
  )

# Define order of rows and columns --------------------------------

col_order <- quos(id, region_iso, sex, year, age_start, age_width)
dat$skeleton <-
  dat$skeleton %>%
  arrange(id) %>%
  select(!!!col_order)

# Export ---------------------------------------------------------

saveRDS(dat$skeleton, file = paths$output$harmonized_skeleton)
