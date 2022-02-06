# Estimate life tables and associated statistics with CIs
#
# Simulation based life table inference based on Poisson samples of
# observed death counts. Implemented in a huge array.

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
  lifetables = './out/lifetables.rds',
  sexdiff = './out/sexdiff.rds',
  e0avgdiff = './out/e0avgdiff.rds',
  codecomp = './out/codecomp.rds',
  arriaga_cntfc = './out/arriaga_cntfc.rds'
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

tmp <- list()

# Function --------------------------------------------------------

source(paths$input$figspecs)

# this function returns TRUE wherever elements are the same,
# including NA's, and FALSE everywhere else
compareNA <- function(v1, v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

# return a vector of mean and quantiles
QuantileWithMean <- function (x, prob = cnst$quantiles) {
  q <- quantile(x, prob = prob, names = FALSE, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  result <- c(m, q)
  names(result) <- c('mean', paste0('q', prob))
  return(result)
}

# interpolate a range of values with its mean
InterpolateWithMean <- function (x) {
  rep(mean(x), length(x))
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
# 7D array [age, year, region_id, sex, sim_id, var_id, projected]
# DESCRIPTION OF VAR ID'S
# LIFE TABLE
# (1)  <population_py> person years exposure
# (2)  <death_total>   total death count
# (3)  <death_covid>   covid death count
# (4)  <nmx>           death rate
# (5)  <npx>           conditional probability of surviving age x
# (6)  <nqx>           conditional probability of dying within age x
# (7)  <lx>            probability of surviving to age x
# (8)  <ndx>           probability of dying within age x
# (9)  <nLx>           life table person years of exposure in age x
# (10) <Tx>            life table person years of exposure above age x            
# (11) <ex>            life expectancy at age x
# ARRIAGA DECOMPOSITION
# (12) <ex_lag>        1 year lag in ex
# (13) <ex_diff>       1 year difference in ex         
# (14) <ex_diff_lag>   1 year lag in ex difference
# (15) <lx_lag>        1 year lag in lx
# (16) <Lx_lag>        1 year lag in Lx
# (17) <Tx_lag>        1 year lag in Tx
# (18) <e0_cntrb_d>    direct contribution of nmx changes to e0 changes
# (19) <e0_cntrb_i>    indirect contribution of nmx changes to e0 changes
# (20) <e0_cntrb_t>    total contrib. of nmx changes to e0 changes
# CAUSE OF DEATH DECOMPOSITION
# (21) <nmx_lag>       1 year lag in nmx
# (22) <R_covid>       share of covid deaths
# (23) <R_covid_lag>   1 year lag in covid death
# (24) <e0_cntrb_t_covid>    total contribution of covid mortality
#                            changes to e0 changes
# (25) <e0_cntrb_t_noncovid> total contrib. of noncovid mortality
#                            changes to e0 changes
# ADDITIONAL INDICATORS
# (26) <bbi>                Bounce-back indicator
# (27) <ex_diff_2_year>     Life-expectancy change from 2 years prior
lifetables$simulation <-
  array(
    dim = c(
      age = lifetables$cnst$nage,
      year = lifetables$cnst$nyears,
      region_iso = lifetables$cnst$nregions,
      sex = 3,
      sim_id = cnst$n_sim+1,
      var_id = 27,
      projected = 2
    ),
    dimnames = list(
      0:85,
      config$skeleton$year$min:config$skeleton$year$max,
      sort(config$skeleton$region),
      c('F', 'M', 'T'),
      1:(cnst$n_sim+1),
      c('population_py', 'death_total', 'death_covid',
        'nmx', 'npx', 'nqx', 'lx', 'ndx', 'nLx', 'Tx', 'ex',
        'ex_lag', 'ex_diff', 'ex_diff_lag', 'lx_lag',
        'nLx_lag', 'Tx_lag',
        'e0_cntrb_d', 'e0_cntrb_i', 'e0_cntrb_t',
        'nmx_lag', 'R_covid', 'R_covid_lag',
        'e0_cntrb_t_covid', 'e0_cntrb_t_noncovid', 'bbi', 'ex_diff_2_year'),
      c('actual', 'projected')
    )
  )

# actual observables
lifetables$simulation[,,,,,'population_py','actual'] <-
  lifetables$input$population_py
lifetables$simulation[,,,,1,'death_total','actual'] <-
  lifetables$input$death_total
lifetables$simulation[,,,,1,'death_covid','actual'] <-
  lifetables$input$death_covid

# projected observables based on 5 year average nmx change
# here we store all-cause death counts how we would
# expect them under continuing pre-pandemic trends
lifetables$simulation[,,,,,'population_py','projected'] <-
  lifetables$input$population_py
tmp$nmx <-
  lifetables$simulation[,,,,1,'death_total','actual'] /
  lifetables$simulation[,,,,1,'population_py','actual']
tmp$nmx_lag <-
  tmp$nmx[,c(NA,1:(lifetables$cnst$nyears-1)),,]
tmp$nmx_diff <- tmp$nmx-tmp$nmx_lag
tmp$nmx_diff_avg_pre2020 <-
  apply(
    tmp$nmx_diff[,2:5,,],
    # apply function to vector of data by year
    MARGIN = c(1,3,4), function (x) rep(mean(x), lifetables$cnst$nyears)
  ) %>%
  aperm(c(2,1,3,4))
tmp$nmx_cntf <-
  apply(
    tmp$nmx_diff_avg_pre2020,
    # apply function to vector of data by year
    c(1,3,4), function (x) x*seq(-4,2,1)
  ) %>%
  aperm(c(2,1,3,4)) +
  tmp$nmx[,rep(5, 7),,]
lifetables$simulation[,,,,1,'death_total','projected'] <-
  tmp$nmx_cntf*lifetables$simulation[,,,,1,'population_py','projected']
lifetables$simulation[,,,,,'death_covid','projected'] <- 0

# simulate total death counts
lifetables$simulation[,,,,-1,'death_total',] <-
  apply(lifetables$simulation[,,,,1,'death_total',],
        MARGIN = 1:5, function (lambda) rpois(n = cnst$n_sim, lambda),
        simplify = TRUE) %>%
  aperm(c(2,3,4,5,1,6))
# simulate covid death counts
lifetables$simulation[,,,,-1,'death_covid','actual'] <-
  apply(lifetables$simulation[,,,,1,'death_covid','actual'],
        MARGIN = 1:4, function (lambda) rpois(n = cnst$n_sim, lambda),
        simplify = TRUE) %>%
  aperm(c(2,3,4,5,1))

# Calculate lifetables over simulated counts ----------------------

# nmx
lifetables$simulation[,,,,,'nmx',] <-
  lifetables$simulation[,,,,,'death_total',] /
  lifetables$simulation[,,,,,'population_py',]
# npx, using constant hazard assumption,
# i.e. npx = exp(-nmx)) for single year age groups
lifetables$simulation[,,,,,'npx',] <-
  exp(-lifetables$simulation[,,,,,'nmx',])
lifetables$simulation[86,,,,,'npx',] <- 0
# nqx
lifetables$simulation[,,,,,'nqx',] <-
  1-lifetables$simulation[,,,,,'npx',]
# lx
lifetables$simulation[,,,,,'lx',] <-
  apply(
    lifetables$simulation[,,,,,'npx',],
    # apply function to vector of data by age
    MARGIN = 2:6, function (npx) head(cumprod(c(1, npx)), -1)
  )
# ndx
lifetables$simulation[,,,,,'ndx',] <-
  apply(
    lifetables$simulation[,,,,,'lx',],
    # apply function to vector of data by age
    MARGIN = 2:6, function (lx) c(-diff(lx), tail(lx, 1))
  )
# nLx = ifelse(mx==0, lx*nx, ndx/nmx)
lifetables$simulation[,,,,,'nLx',] <-
  lifetables$simulation[,,,,,'ndx',]/lifetables$simulation[,,,,,'nmx',]
tmp$I <- compareNA(lifetables$simulation[,,,,,'nmx',],0)
lifetables$simulation[,,,,,'nLx',][tmp$I] <-
  lifetables$simulation[,,,,,'lx',][tmp$I]
# Tx = rev(cumsum(rev(nLx)))
lifetables$simulation[,,,,,'Tx',] <-
  apply(
    lifetables$simulation[,,,,,'nLx',],
    # apply function to vector of data by age
    MARGIN = 2:6, function (nLx) rev(cumsum(rev(nLx)))
  )
# ex = Tx/lx
lifetables$simulation[,,,,,'ex',] <-
  lifetables$simulation[,,,,,'Tx',] /
  lifetables$simulation[,,,,,'lx',]

# Calculate annual ex change --------------------------------------

lifetables$simulation[,,,,,'ex_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'ex',]
lifetables$simulation[,,,,,'ex_diff',] <-
  lifetables$simulation[,,,,,'ex',]-lifetables$simulation[,,,,,'ex_lag',]
lifetables$simulation[,,,,,'ex_diff_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'ex_diff',]

# Calculate Arriaga decomposition ---------------------------------

# decompose annual changes in e0 into age specific mortality changes

# see Arriaga (1984)
# Measuring and explaining the change in life expectancies
# DOI 10.2307/2061029

lifetables$simulation[,,,,,'lx_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'lx',]
lifetables$simulation[,,,,,'nLx_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'nLx',]
lifetables$simulation[,,,,,'Tx_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'Tx',]

lifetables$simulation[,,,,,'e0_cntrb_d',] <-
  (
    lifetables$simulation[,,,,,'nLx',]/lifetables$simulation[,,,,,'lx',]-
      lifetables$simulation[,,,,,'nLx_lag',]/lifetables$simulation[,,,,,'lx_lag',]
  ) * lifetables$simulation[,,,,,'lx',]

lifetables$simulation[,,,,,'e0_cntrb_i',] <-
  (
    lifetables$simulation[,,,,,'lx_lag',]/
      lifetables$simulation[,,,,,'lx',]-
      apply(lifetables$simulation[,,,,,'lx_lag',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))/
      apply(lifetables$simulation[,,,,,'lx',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))
  ) * apply(lifetables$simulation[,,,,,'Tx',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))
lifetables$simulation[86,,,,,'e0_cntrb_i',] <- 0

lifetables$simulation[,,,,,'e0_cntrb_t',] <-
  lifetables$simulation[,,,,,'e0_cntrb_d',] +
  lifetables$simulation[,,,,,'e0_cntrb_i',]

# Calculate cause of death decomposition --------------------------

lifetables$simulation[,,,,,'nmx_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'nmx',]

# R_covid = D_covid / D_total
lifetables$simulation[,,,,,'R_covid',] <-
  lifetables$simulation[,,,,,'death_covid',]/lifetables$simulation[,,,,,'death_total',]
# R_noncovid = 1 - R_covid
lifetables$simulation[,,,,,'R_covid_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'R_covid',]

lifetables$simulation[,,,,,'e0_cntrb_t_covid',] <-
  lifetables$simulation[,,,,,'e0_cntrb_t',] * (
    (
      lifetables$simulation[,,,,,'R_covid',] *
        lifetables$simulation[,,,,,'nmx',] -
        lifetables$simulation[,,,,,'R_covid_lag',] *
        lifetables$simulation[,,,,,'nmx_lag',]
    ) /
      (lifetables$simulation[,,,,,'nmx',] -
         lifetables$simulation[,,,,,'nmx_lag',])
  )

lifetables$simulation[,,,,,'e0_cntrb_t_noncovid',] <-
  lifetables$simulation[,,,,,'e0_cntrb_t',] -
  lifetables$simulation[,,,,,'e0_cntrb_t_covid',]

# Calculate additional indicators ---------------------------------

# bounce back indicator
lifetables$simulation[,,,,,'bbi',] <-
  1 - (
    lifetables$simulation[,,,,,'ex_diff_lag',] +
      lifetables$simulation[,,,,,'ex_diff',]
  ) /
  lifetables$simulation[,,,,,'ex_diff_lag',]
tmp$I <- compareNA(sign(lifetables$simulation[,,,,,'ex_diff_lag',]), 1)
lifetables$simulation[,,,,,'bbi',][tmp$I] <-
  -lifetables$simulation[,,,,,'bbi',][tmp$I]

# 2 year life expectancy difference
lifetables$simulation[,,,,,'ex_diff_2_year',] <-
  lifetables$simulation[,,,,,'ex_diff',] +
  lifetables$simulation[,,,,,'ex_diff_lag',]

# Calculate sex differences ---------------------------------------

sexdiff <- list()
sexdiff$simulation <-
  lifetables$simulation[,,,'F',,,]-
  lifetables$simulation[,,,'M',,,]

# add additional variables to array
# var_id:
# (27) <ex_diff_sign>                sign of sex difference in ex change
# (28) <ex_diff_change_from_2019>    change in sex difference from 2019
# (29) <ex_diff_drop_from_2019_flag> drop in sex difference from 2019
# (30) <ex_diff_rise_from_2019_flag> rise in sex difference from 2019
D <- dim(sexdiff$simulation)
D[5] <- D[5]+4
Dn <- dimnames(sexdiff$simulation)
Dn[[5]] <- c(
  Dn[[5]], 'ex_diff_sign', 'ex_diff_change_from_2019',
  'ex_diff_drop_from_2019_flag',
  'ex_diff_rise_from_2019_flag'
)
sexdiff$temp <- array(NA, D, Dn)
sexdiff$temp[,,,,-(27:30),] <- sexdiff$simulation
sexdiff$simulation <- sexdiff$temp

# +2: ex increase for both sexes
# -2: ex decrease for both sexes
# 0: mixed increase-decrease by sex
sexdiff$simulation[,,,,'ex_diff_sign',] <-
  sign(lifetables$simulation[,,,'F',,'ex_diff',])+
  sign(lifetables$simulation[,,,'M',,'ex_diff',])

sexdiff$simulation[,,,,'ex_diff_change_from_2019',] <-
  sexdiff$simulation[,,,,'ex_diff',] -
  sexdiff$simulation[,rep('2019', dim(sexdiff$simulation)[2]),,,'ex_diff',]

sexdiff$simulation[,,,,'ex_diff_drop_from_2019_flag',] <-
  sexdiff$simulation[,,,,'ex_diff_change_from_2019',] < 0

sexdiff$simulation[,,,,'ex_diff_rise_from_2019_flag',] <-
  sexdiff$simulation[,,,,'ex_diff_change_from_2019',] >= 0

# Calculate average ex change 2016 to 2019 ------------------------

e0avgdiff <- list()

e0avgdiff$simulation <- apply(
  lifetables$simulation[,2:5,,,,'ex_diff',],
  # apply function to vector of data by year
  MARGIN = c(1, 3:6), mean
) %>%
  aperm(c(1,3,2,4,5))

# Calculate cause contribution to annual e0 change ----------------

# simply sum the age specific cause contributions

codecomp <- list()

codecomp$simulation <-
  apply(
    lifetables$simulation[
      ,,,,,
      c('e0_cntrb_t_covid', 'e0_cntrb_t_noncovid', 'e0_cntrb_t'),
    ],
    # apply function to vector of data by age
    MARGIN = c(2:7), FUN = sum
  )


# Counterfactual Arriaga decomposition ----------------------------

# decompose annual changes in e0 into age specific mortality changes

# see Arriaga (1984)
# Measuring and explaining the change in life expectancies
# DOI 10.2307/2061029

arriaga_cntfc <- list()
arriaga_cntfc$simulation <-
  array(
    dim = c(
      age = lifetables$cnst$nage,
      year = lifetables$cnst$nyears,
      region_iso = lifetables$cnst$nregions,
      sex = 3,
      sim_id = cnst$n_sim+1,
      var_id = 22
    ),
    dimnames = list(
      0:85,
      config$skeleton$year$min:config$skeleton$year$max,
      sort(config$skeleton$region),
      c('F', 'M', 'T'),
      1:(cnst$n_sim+1),
      c('death_total_actual', 'death_total_expected',
        'death_covid_actual', 'death_covid_expected',
        'nmx_actual', 'nmx_expected',
        'ex_actual_minus_expected',
        'ex_actual', 'ex_expected',
        'lx_actual', 'lx_expected',
        'nLx_actual', 'nLx_expected',
        'Tx_actual', 'Tx_expected',
        'e0_cntrb_d', 'e0_cntrb_i', 'e0_cntrb_t',
        'R_covid_actual', 'R_covid_expected',
        'e0_cntrb_t_covid', 'e0_cntrb_t_noncovid')
    )
  )

arriaga_cntfc$simulation[,,,,,'death_total_actual'] <-
  lifetables$simulation[,,,,,'death_total','actual']
arriaga_cntfc$simulation[,,,,,'death_total_expected'] <-
  lifetables$simulation[,,,,,'death_total','projected']
arriaga_cntfc$simulation[,,,,,'death_covid_actual'] <-
  lifetables$simulation[,,,,,'death_covid','actual']
arriaga_cntfc$simulation[,,,,,'death_covid_expected'] <-
  lifetables$simulation[,,,,,'death_covid','projected']

arriaga_cntfc$simulation[,,,,,'nmx_actual'] <-
  lifetables$simulation[,,,,,'nmx','actual']
arriaga_cntfc$simulation[,,,,,'nmx_expected'] <-
  lifetables$simulation[,,,,,'nmx','projected']

arriaga_cntfc$simulation[,,,,,'lx_actual'] <-
  lifetables$simulation[,,,,,'lx','actual']
arriaga_cntfc$simulation[,,,,,'lx_expected'] <-
  lifetables$simulation[,,,,,'lx','projected']

arriaga_cntfc$simulation[,,,,,'nLx_actual'] <-
  lifetables$simulation[,,,,,'nLx','actual']
arriaga_cntfc$simulation[,,,,,'nLx_expected'] <-
  lifetables$simulation[,,,,,'nLx','projected']

arriaga_cntfc$simulation[,,,,,'Tx_actual'] <-
  lifetables$simulation[,,,,,'Tx','actual']
arriaga_cntfc$simulation[,,,,,'Tx_expected'] <-
  lifetables$simulation[,,,,,'Tx','projected']

arriaga_cntfc$simulation[,,,,,'ex_actual'] <-
  lifetables$simulation[,,,,,'ex','actual']
arriaga_cntfc$simulation[,,,,,'ex_expected'] <-
  lifetables$simulation[,,,,,'ex','projected']

arriaga_cntfc$simulation[,,,,,'ex_actual_minus_expected'] <-
  arriaga_cntfc$simulation[,,,,,'ex_actual'] - arriaga_cntfc$simulation[,,,,,'ex_expected']

arriaga_cntfc$simulation[,,,,,'e0_cntrb_d'] <-
  (
    arriaga_cntfc$simulation[,,,,,'nLx_actual']/
      arriaga_cntfc$simulation[,,,,,'lx_actual'] -
      arriaga_cntfc$simulation[,,,,,'nLx_expected']/
      arriaga_cntfc$simulation[,,,,,'lx_expected']
  ) * arriaga_cntfc$simulation[,,,,,'lx_actual']

arriaga_cntfc$simulation[,,,,,'e0_cntrb_i'] <-
  (
    arriaga_cntfc$simulation[,,,,,'lx_expected']/
      arriaga_cntfc$simulation[,,,,,'lx_actual']-
      apply(arriaga_cntfc$simulation[,,,,,'lx_expected'],
            # apply function to vector of data by age
            2:5, function (x) c(x[-1], 0))/
      apply(arriaga_cntfc$simulation[,,,,,'lx_actual'],
            # apply function to vector of data by age
            2:5, function (x) c(x[-1], 0))
  ) * apply(arriaga_cntfc$simulation[,,,,,'Tx_actual'],
            # apply function to vector of data by age
            2:5, function (x) c(x[-1], 0))
arriaga_cntfc$simulation[86,,,,,'e0_cntrb_i'] <- 0

arriaga_cntfc$simulation[,,,,,'e0_cntrb_t'] <-
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_d'] +
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_i']

# Calculates CI over simulations ----------------------------------

# ci's for life tables
lifetables$ci <-
  apply(
    lifetables$simulation[
      ,,,,-1,
      c('nmx', 'npx', 'nqx', 'lx', 'ex', 'ex_diff', 'ex_diff_2_year', 'bbi',
        'e0_cntrb_t', 'e0_cntrb_t_covid', 'e0_cntrb_t_noncovid'),
    ],
    -5,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:7,1))
V <- dimnames(lifetables$ci)
names(attr(lifetables$ci, 'dim'))[7] <- 'quantile'
dimnames(lifetables$ci) <- V

# ci's for sex-differences of life table statistics
sexdiff$ci <-
  apply(
    sexdiff$simulation[
      ,,,-1,
      c('nmx', 'npx', 'nqx', 'lx', 'ex', 'ex_diff', 'bbi', 'ex_diff_sign',
        'ex_diff_change_from_2019',
        'ex_diff_drop_from_2019_flag',
        'ex_diff_rise_from_2019_flag'),
    ],
    -4,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:6,1))
V <- dimnames(sexdiff$ci)
names(attr(sexdiff$ci, 'dim'))[6] <- 'quantile'
dimnames(sexdiff$ci) <- V

# ci's for average e0 annual change
e0avgdiff$ci <-
  apply(
    e0avgdiff$simulation[,,,-1,],
    -4,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:5,1))
V <- dimnames(e0avgdiff$ci)
names(attr(e0avgdiff$ci, 'dim'))[5] <- 'quantile'
dimnames(e0avgdiff$ci) <- V

# ci's for cause contributions to annual life expectancy change
codecomp$ci <-
  apply(
    codecomp$simulation[,,,-1,,],
    -4,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:6,1))
V <- dimnames(codecomp$ci)
names(attr(codecomp$ci, 'dim'))[6] <- 'quantile'
dimnames(codecomp$ci) <- V

# ci's for age specific contributions to e0 deviations from expectation
arriaga_cntfc$ci <-
  apply(
    arriaga_cntfc$simulation[,,,,-1,],
    -5,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:6,1))
V <- dimnames(arriaga_cntfc$ci)
names(attr(arriaga_cntfc$ci, 'dim'))[6] <- 'quantile'
dimnames(arriaga_cntfc$ci) <- V

# Transform to data frame -----------------------------------------

lifetables$ci_df <-
  as.data.frame.table(lifetables$ci, stringsAsFactors = FALSE)
names(lifetables$ci_df) <-
  c(names(attr(lifetables$ci, 'dim')), 'value')
lifetables$ci_df <-
  lifetables$ci_df %>%
  as_tibble() %>%
  pivot_wider(id_cols = c(age, year, region_iso, sex, projected),
              names_from = c(var_id, quantile),
              values_from = value) %>%
  mutate(across(c(age, year), ~as.integer(.x)))

sexdiff$ci_df <-
  as.data.frame.table(sexdiff$ci, stringsAsFactors = FALSE)
names(sexdiff$ci_df) <-
  c(names(attr(sexdiff$ci, 'dim')), 'value')
sexdiff$ci_df <-
  sexdiff$ci_df %>%
  as_tibble() %>%
  pivot_wider(id_cols = c(age, year, region_iso, projected),
              names_from = c(var_id, quantile),
              values_from = value) %>%
  mutate(across(c(age, year), ~as.integer(.x)))

e0avgdiff$ci_df <-
  as.data.frame.table(e0avgdiff$ci, stringsAsFactors = FALSE)
names(e0avgdiff$ci_df) <-
  c(names(attr(e0avgdiff$ci, 'dim')), 'value')
e0avgdiff$ci_df <-
  e0avgdiff$ci_df %>%
  as_tibble() %>%
  pivot_wider(id_cols = c(age, sex, region_iso, projected),
              names_from = c(quantile),
              values_from = value) %>%
  mutate(across(c(age), ~as.integer(.x)))

codecomp$ci_df <-
  as.data.frame.table(codecomp$ci, stringsAsFactors = FALSE)
names(codecomp$ci_df) <-
  c(names(attr(codecomp$ci, 'dim')), 'value')
codecomp$ci_df <-
  codecomp$ci_df %>%
  as_tibble() %>%
  pivot_wider(id_cols = c(sex, region_iso, year, projected),
              names_from = c(var_id, quantile),
              values_from = value)

arriaga_cntfc$ci_df <-
  as.data.frame.table(arriaga_cntfc$ci, stringsAsFactors = FALSE)
names(arriaga_cntfc$ci_df) <-
  c(names(attr(arriaga_cntfc$ci, 'dim')), 'value')
arriaga_cntfc$ci_df <-
  arriaga_cntfc$ci_df %>%
  filter(year >= 2020) %>% 
  as_tibble() %>%
  pivot_wider(id_cols = c(sex, region_iso, year, age),
              names_from = c(var_id, quantile),
              values_from = value)

# Test ------------------------------------------------------------

lifetables$simulation['0','2021',,'T',1,'population_py','projected'] ==
  lt_input$openage_85 %>%
  filter(age_start == 0, year == 2021, sex == 'Total') %>%
  select(population_py)

# check: cause specific contribution to life expectancy change by age
# sum to age specific contributions
lifetables$ci_df %>%
  mutate(
    diff = abs(e0_cntrb_t_mean - (e0_cntrb_t_covid_mean + e0_cntrb_t_noncovid_mean)),
    reldiff = diff/e0_cntrb_t_mean
  ) %>%
  filter(year %in% 2020:2021) %>%
  arrange(-reldiff) %>%
  pull(reldiff) %>% summary()
  
# check: cause specific contribution to life expectancy change by age
# sum to overall life expectancy change
lifetables$ci_df %>%
  group_by(region_iso, sex, year, projected) %>%
  summarise(
    e0_diff_1 = ex_diff_mean[1],
    e0_diff_2 = sum(e0_cntrb_t_covid_mean + e0_cntrb_t_noncovid_mean),
    reldiff = abs(e0_diff_2-e0_diff_1)/abs(e0_diff_1)
  ) %>%
  ungroup() %>%
  filter(year %in% 2020:2021) %>%
  arrange(-reldiff) %>%
  pull(reldiff) %>% summary()

# Export ----------------------------------------------------------

saveRDS(lifetables$ci_df, paths$output$lifetables)
saveRDS(sexdiff$ci_df, paths$output$sexdiff)
saveRDS(e0avgdiff$ci_df, paths$output$e0avgdiff)
saveRDS(codecomp$ci_df, paths$output$codecomp)
saveRDS(arriaga_cntfc$ci_df, paths$output$arriaga_cntfc)
