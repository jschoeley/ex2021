# Add counterfactual mortality estimates

# Init ------------------------------------------------------------

library(dplyr); library(tidyr); library(yaml)
library(demography)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  skeleton = './tmp/harmonized_skeleton.rds',
  global = './src/00-global.R',
  death = './tmp/harmonized_death.rds',
  population = './tmp/harmonized_population.rds',
  hmd_lifetables = './tmp/harmonized_hmd_lifetables.rds'
  #covid = './tmp/harmonized_covid.rds',
  #ex = './tmp/harmonized_ex.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  harmonized_counterfactual_mortality = './tmp/harmonized_counterfactual_mortality.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# list containers for analysis artifacts
dat <- list()

# Functions -------------------------------------------------------

# global functions
source(paths$input$global)

# Load data -------------------------------------------------------

dat$skeleton <- readRDS(paths$input$skeleton)
dat$hmd_lifetables <- readRDS(paths$input$hmd_lifetables)
dat$death <- readRDS(paths$input$death)
dat$population <- readRDS(paths$input$population)

# Generate dataset for mortality forecasts ------------------------

dat$counterfactual_mortality_skeleton <-
  dat$skeleton %>%
  select(-id, -age_width) %>%
  complete(region_iso, sex, year = 1990:2021, age_start) %>%
  mutate(id = GenerateRowID(region_iso, sex, year, age_start)) %>%
  left_join(dat$hmd_lifetables) %>%
  left_join(dat$death) %>%
  left_join(dat$population)

# Add nmx forecasts -----------------------------------------------

# Lee-Carter forecast of HMD death rates adjusted to jump off
# at our STMF nmx estimates in 2019

dat$counterfactual_mortality <-
  dat$counterfactual_mortality_skeleton %>%
  group_by(region_iso, sex) %>%
  group_modify(~{
    
    cat(.y$region_iso, ' ', .y$sex, '\n')
    
    years = 1990:config$skeleton$year$max
    jumpoff_year = 2019
    h = config$skeleton$year$max-jumpoff_year
    n_years = length(years)
    age = config$skeleton$age$min:config$skeleton$age$max
    n_age = length(age)
    fudge = 1e-6
    
    dat <- .x
    
    M <- matrix(
      dat$nmx_hmd,
      nrow = n_age,
      ncol = n_years,
      dimnames = list(age, years)
    )
    M[M == 0] <- fudge
    P <- matrix(NA, nrow = n_age, ncol = n_years)
    DD <-
      demogdata(
        M, pop = P, type = 'mortality', ages = age, years = years,
        label = '', name = 'nmx'
      )
    LC <-
      lca(DD, series = 'nmx', interpolate = TRUE, adjust = 'e0',
          years = 1990:jumpoff_year)
    ax <- LC$ax
    bx <- LC$bx
    kt <- LC$kt
    ax_adjusted_jumpoff <-
      log(.x$death_total/.x$population_py)[which(.x$year==jumpoff_year)]
    fcst <- forecast::rwf(kt, h = h, drift = TRUE)
    drift <- fcst$model$par$drift
    x <- 1:h
    nmx_fcst <- exp(
      apply(t(x*drift), 2, function (kt_fcst) kt_fcst*bx) +
        ax_adjusted_jumpoff      
    )

    dat$nmx_cntfc <- 
      c(
        c(LC$fitted$y),
        c(nmx_fcst)
      )
    
    return(dat)
    
  }) %>%
  ungroup() %>%
  # add row id
  mutate(id = GenerateRowID(region_iso, sex, year, age_start)) %>%
  right_join(dat$skeleton, by = 'id') %>%
  select(id, nmx_cntfc)

# Export ----------------------------------------------------------

saveRDS(dat$counterfactual_mortality,
        paths$output$harmonized_counterfactual_mortality)
