# Estimate life tables and associated statistics with CIs
#
# Simulation based inference based on Poisson samples of
# observed death counts

# Init ------------------------------------------------------------

library(glue); library(yaml)
library(dplyr); library(tidyr)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  lt_input = './out/lt_input.rds',
  figspecs = './cfg/figure_specification.R'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  lifetables = './out/lifetables.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  regions_for_analysis = config$regions_for_all_cause_analysis
  # number of Poisson life-table replicates
  n_sim = 500
  # quantiles for CI's
  quantiles = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)
})

# Function --------------------------------------------------------

source(paths$input$figspecs)

# This function returns TRUE wherever elements are the same, including NA's,
# and FALSE everywhere else.
compareNA <- function(v1, v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

# Data ------------------------------------------------------------

lt_input <- list()

# input data for life-table calculation
# harmonized death counts and population exposures with open age group 100+
lt_input$openage_100 <- readRDS(paths$input$lt_input)

# Create open age group 85+ ---------------------------------------

lt_input$openage_85 <-
  lt_input$openage_100 %>%
  # for each life-table input stratum create age group 85+
  group_by(region_iso, sex, year) %>%
  group_modify(~{
    input_sorted <- arrange(.x, age_start)
    lt85 <- filter(input_sorted, age_start <= 85)
    lt85p <- filter(input_sorted, age_start > 85)
    
    lt85[nrow(lt85), 'age_width'] <- Inf
    lt85[nrow(lt85), 'death_total'] <- sum(lt85p$death_total)
    lt85[nrow(lt85), 'population_midyear'] <- sum(lt85p$population_midyear)
    lt85[nrow(lt85), 'population_py'] <- sum(lt85p$population_py)
    lt85[nrow(lt85), 'death_covid'] <- sum(lt85p$death_covid)
    
    return(lt85)
  }) %>%
  ungroup() %>%
  # harmonize order of variables
  select(all_of(names(lt_input$openage_100))) %>%
  arrange(sex, region_iso, year, age_start)

# Create Poisson replicates of counts -----------------------------

# create life table replicates by region, sex, and year
# based on repeatedly sampling death counts from a Poisson
# distribution with mean equal to estimated mean from PCLM

lifetables <- list()
lifetables$cnst <- list(
  nage = 86,
  nyears = config$skeleton$year$max-config$skeleton$year$min+1,
  nregions = length(config$skeleton$region)
)

lifetables$input <-
  lt_input$openage_85 %>%
  select(id, age_start, age_width, population_py,
         death_total, death_covid)
# vector ordered by sex, region, year, age distributed into
# 6D array [age x year, region_id, sex, sim_id, var_id]
# var_id:
# (1) person years exposure
# (2) total death count
# (3) covid death count
# (4) nmx
# (5) npx
# (6) nqx
# (7) lx
# (8) ndx
# (9) nLx
# (10) Tx
# (11) ex
# (12) ex_lag
# (13) ex_diff
# (14) ex_diff_lag
# (15) bbi
lifetables$simulation <-
  array(
    dim = c(
      age = lifetables$cnst$nage,
      year = lifetables$cnst$nyears,
      region_iso = lifetables$cnst$nregions,
      sex = 3,
      sim_id = cnst$n_sim+1,
      var_id = 15
    ),
    dimnames = list(
      0:85,
      config$skeleton$year$min:config$skeleton$year$max,
      sort(config$skeleton$region),
      c('F', 'M', 'T'),
      1:(cnst$n_sim+1),
      c('population_py', 'death_total', 'death_covid',
        'nmx', 'npx', 'nqx', 'lx', 'ndx', 'nLx', 'Tx', 'ex',
        'ex_lag', 'ex_diff', 'ex_diff_lag', 'bbi')
    )
  )
lifetables$simulation[,,,,,'population_py'] <-
  lifetables$input$population_py
lifetables$simulation[,,,,1,'death_total'] <-
  lifetables$input$death_total
lifetables$simulation[,,,,1,'death_covid'] <-
  lifetables$input$death_covid
# simulate total death counts
lifetables$simulation[,,,,-1,'death_total'] <-
  apply(lifetables$simulation[,,,,1,'death_total'], MARGIN = 1:4,
        function (lambda) rpois(n = cnst$n_sim, lambda),
        simplify = TRUE) %>%
  aperm(c(2,3,4,5,1))
# simulate covid death counts
lifetables$simulation[,,,,-1,'death_covid'] <-
  apply(lifetables$simulation[,,,,1,'death_covid'], MARGIN = 1:4,
        function (lambda) rpois(n = cnst$n_sim, lambda),
        simplify = TRUE) %>%
  aperm(c(2,3,4,5,1))

# Calculate lifetables over simulated counts ----------------------

# nmx
lifetables$simulation[,,,,,'nmx'] <-
  lifetables$simulation[,,,,,'death_total'] /
  lifetables$simulation[,,,,,'population_py']
# npx
lifetables$simulation[,,,,,'npx'] <-
  exp(-lifetables$simulation[,,,,,'nmx'])
# nqx
lifetables$simulation[,,,,,'nqx'] <-
  1-lifetables$simulation[,,,,,'npx']
# lx
lifetables$simulation[,,,,,'lx'] <-
  apply(
    lifetables$simulation[,,,,,'npx'], MARGIN = 2:5,
    function (npx) head(cumprod(c(1, npx)), -1)
  )
# ndx
lifetables$simulation[,,,,,'ndx'] <-
  apply(
    lifetables$simulation[,,,,,'lx'], MARGIN = 2:5,
    function (lx) c(-diff(lx), tail(lx, 1))
  )
# nLx = ifelse(mx==0, lx*nx, ndx/nmx)
lifetables$simulation[,,,,,'nLx'] <-
  lifetables$simulation[,,,,,'ndx']/lifetables$simulation[,,,,,'nmx']
I <- compareNA(lifetables$simulation[,,,,,'nmx'],0)
lifetables$simulation[,,,,,'nLx'][I] <-
  lifetables$simulation[,,,,,'lx'][I]
# Tx = rev(cumsum(rev(nLx)))
lifetables$simulation[,,,,,'Tx'] <-
  apply(
    lifetables$simulation[,,,,,'nLx'], MARGIN = 2:5,
    function (nLx) rev(cumsum(rev(nLx)))
  )
# ex = Tx/lx
lifetables$simulation[,,,,,'ex'] <-
  lifetables$simulation[,,,,,'Tx']/lifetables$simulation[,,,,,'lx']

# Calculate annual ex change --------------------------------------

lifetables$simulation[,,,,,'ex_lag'] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'ex']

lifetables$simulation[,,,,,'ex_diff'] <-
  lifetables$simulation[,,,,,'ex']-lifetables$simulation[,,,,,'ex_lag']

lifetables$simulation[,,,,,'ex_diff_lag'] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'ex_diff']

# bounce back indicator
lifetables$simulation[,,,,,'bbi'] <-
  1 - (
    lifetables$simulation[,,,,,'ex_diff_lag'] +
    lifetables$simulation[,,,,,'ex_diff']
  ) /
  lifetables$simulation[,,,,,'ex_diff_lag']
I <- compareNA(sign(lifetables$simulation[,,,,,'ex_diff_lag']), 1)
lifetables$simulation[,,,,,'bbi'][I] <-
  -lifetables$simulation[,,,,,'bbi'][I]

# Calculate sex differences ---------------------------------------

sexdiff <- list()
sexdiff$simulation <-
  lifetables$simulation[,,,'F',,]-
  lifetables$simulation[,,,'M',,]

# Calculates CI over simulations ----------------------------------

lifetables$ci <-
  apply(
    lifetables$simulation[
      ,,,,-1,
      c('nmx', 'npx', 'nqx', 'lx', 'ex', 'ex_diff', 'bbi')
    ],
    c(1:4, 6),
    quantile, prob = cnst$quantiles, names = FALSE, na.rm = TRUE,
    simplify = TRUE
  ) %>%
  aperm(c(2:6,1))
V <- dimnames(lifetables$ci); V[[6]] <- paste0('q', cnst$quantiles)
names(attr(lifetables$ci, 'dim'))[6] <- 'quantile'
dimnames(lifetables$ci) <- V

# Test ------------------------------------------------------------

lifetables$simulation['0','2021',,'T',1,'population_py'] ==
  lt_input$openage_85 %>%
  filter(age_start == 0, year == 2021, sex == 'Total') %>%
  select(population_py)

# Export ----------------------------------------------------------

# export results of analysis
lifetables$ci_df <-
  as.data.frame.table(lifetables$ci)
names(lifetables$ci_df) <-
  c(names(attr(lifetables$ci, 'dim')), 'value')

lifetables$ci_df <-
  lifetables$ci_df %>%
  as_tibble() %>%
  pivot_wider(id_cols = c(age, year, region_iso, sex),
              names_from = c(var_id, quantile),
              values_from = value)

saveRDS(lifetables$ci_df, paths$output$lifetables)