# Add counterfactual mortality estimates

# Init ------------------------------------------------------------

library(dplyr); library(tidyr); library(yaml)
library(demography); library(ggplot2)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  skeleton = './tmp/10-harmonized_skeleton.rds',
  global = './src/00-global.R',
  death = './tmp/21-harmonized_death.rds',
  population = './tmp/20-harmonized_population.rds',
  hmd_lifetables = './tmp/23-harmonized_hmd_lifetables.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  harmonized_counterfactual_mortality = './tmp/24-harmonized_counterfactual_mortality.rds'
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
  left_join(dat$population) %>%
  mutate(
    nmx_pclm = death_total / population_py,
    nmx_to_forecast = ifelse(year < 2015, nmx_hmd, nmx_pclm)
  )

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
    fudge = 0
    
    dat <- .x
    
    M <- matrix(
      dat$nmx_to_forecast,
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
          years = 2015:jumpoff_year)

    nmx_fcst <- forecast::forecast(LC, h = h, jumpchoice = 'actual')$rate[['nmx']]

    l2 <- length(c(exp(c(LC$fitted$y)), c(nmx_fcst)))
    l1 <- nrow(dat)
    dat$nmx_cntfc <- 
      c(
        rep(NA, l1-l2),
        exp(c(LC$fitted$y)),
        c(nmx_fcst)
      )
    
    return(dat)
    
  }) %>%
  ungroup()

# Check forecast --------------------------------------------------

dat$counterfactual_mortality %>%
  #select(year, region_iso, sex, age_start, projected, ex_mean) %>%
  filter(age_start %in% c(40, 60, 80, 90), sex == 'Total') %>%
  mutate(age_start = as.factor(age_start)) %>%
  #pivot_wider(names_from = projected, values_from = ex_mean) %>%
  ggplot(aes(x = year, group = age_start)) +
  geom_point(aes(y = nmx_pclm, color = age_start)) +
  geom_vline(xintercept = 2014.5) +
  geom_vline(xintercept = 2019.5) +
  geom_line(aes(y = nmx_cntfc)) +
  coord_cartesian(xlim = c(2015, 2021)) +
  scale_y_log10() +
  facet_wrap(~region_iso)

dat$counterfactual_lt_debug <-
  dat$counterfactual_mortality %>%
  group_by(region_iso, sex, year) %>%
  mutate(
    npx_cntfc = exp(-nmx_cntfc),
    lx_cntfc = head(cumprod(c(1, npx_cntfc)), -1),
    ndx_cntfc = c(-diff(lx_cntfc), tail(lx_cntfc, 1)),
    nLx_cntfc = ifelse(nmx_cntfc == 0, lx_cntfc, ndx_cntfc/nmx_cntfc),
    Tx_cntfc = rev(cumsum(rev(nLx_cntfc))),
    ex_cntfc = Tx_cntfc/lx_cntfc,
    npx_pclm = exp(-nmx_pclm),
    lx_pclm = head(cumprod(c(1, npx_pclm)), -1),
    ndx_pclm = c(-diff(lx_pclm), tail(lx_pclm, 1)),
    nLx_pclm = ifelse(nmx_pclm == 0, lx_pclm, ndx_pclm/nmx_pclm),
    Tx_pclm = rev(cumsum(rev(nLx_pclm))),
    ex_pclm = Tx_pclm/lx_pclm
  )

dat$counterfactual_lt_debug %>%
  filter(age_start == 0, sex == 'Total', year >= 2015) %>%
  ggplot(aes(x = year)) +
  geom_point(aes(y = ex_hmd), color = 'red') +
  geom_point(aes(y = ex_pclm), color = 'blue') +
  geom_point(aes(y = ex_cntfc), shape = 21, data = . %>% filter(year > 2019)) +
  facet_wrap(~region_iso, scales = 'free_y')

# Prepare for export ----------------------------------------------

# add row id
dat$counterfactual_mortality <-
  dat$counterfactual_mortality %>%
  mutate(id = GenerateRowID(region_iso, sex, year, age_start)) %>%
  right_join(dat$skeleton, by = 'id') %>%
  select(id, nmx_cntfc)

# Export ----------------------------------------------------------

saveRDS(dat$counterfactual_mortality,
        paths$output$harmonized_counterfactual_mortality)
saveRDS(
  dat$counterfactual_lt_debug,
  paste0(paths$output$tmpdir, '/24-counterfactual_lt_debug.rds')
)
