# Estimate life tables and associated statistics with CIs
#
# Simulation based life table inference based on Poisson samples of
# observed death counts. Implemented in a huge array.

# Init ------------------------------------------------------------

library(glue); library(yaml)
library(dplyr); library(tidyr); library(readr)

# Constants -------------------------------------------------------

# randomness in Poisson sampling
set.seed(1987)

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  lt_input = './out/30-lt_input.rds',
  figspecs = './cfg/figure_specification.R'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  lifetables = './out/40-lifetables.rds',
  lifetables_csv = './tmp/40-lifetables.csv',
  sexdiff = './out/40-sexdiff.rds',
  sexdiff_csv = './tmp/40-sexdiff.csv',
  e0avgdiff = './out/40-e0avgdiff.rds',
  e0avgdiff_csv = './tmp/40-e0avgdiff.csv',
  codecomp = './out/40-codecomp.rds',
  codecomp_csv = './tmp/40-codecomp.csv',
  arriaga_cntfc = './out/40-arriaga_cntfc.rds',
  arriaga_cntfc_csv = './tmp/40-arriaga_cntfc.csv'
)

# global configuration
config <- read_yaml(paths$input$config)

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  regions_for_analysis = config$regions_for_all_cause_analysis
  # number of Poisson life-table replicates
  n_sim = 100
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
  x <- x[!(is.na(x)|is.nan(x)|is.infinite(x))]
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

# lt_input$openage_85 <-
#   lt_input$openage_100 %>%
#   # for each life-table input stratum create age group 85+
#   group_by(region_iso, sex, year) %>%
#   group_modify(~{
#     
#     input_sorted <- arrange(.x, age_start)
#     lt85 <- filter(input_sorted, age_start <= 85)
#     lt85p <- filter(input_sorted, age_start >= 85)
#     nage <- nrow(lt85)
#     
#     lt85[nage, 'age_width'] <- Inf
#     lt85[nage, 'death_total'] <- sum(lt85p$death_total)
#     lt85[nage, 'population_midyear'] <- sum(lt85p$population_midyear)
#     lt85[nage, 'population_py'] <- sum(lt85p$population_py)
#     lt85[nage, 'death_covid'] <- sum(lt85p$death_covid)
#     # counterfactual death rate age 85+
#     lx <- head(cumprod(c(1, exp(-input_sorted$nmx_cntfc))), -1)
#     ndx <- c(-diff(lx), tail(lx, 1))
#     nLx <- lx-0.5*ndx # fix later
#     nLx[length(nLx)] <- 1/input_sorted$nmx_cntfc[length(input_sorted$nmx_cntfc)]*lx[length(lx)]
#     Tx <- rev(cumsum(rev(nLx)))
#     lt85[nage, 'nmx_cntfc'] <- lx[86]/Tx[86]
# 
#     return(lt85)
#   }) %>%
#   ungroup() %>%
#   # harmonize order of variables
#   select(all_of(names(lt_input$openage_100))) %>%
#   arrange(sex, region_iso, year, age_start)

lt_input$openage_85 <- lt_input$openage_100

# Create Poisson replicates of counts -----------------------------

# create life table replicates by region, sex, and year
# based on repeatedly sampling death counts from a Poisson
# distribution with mean equal to estimated mean from PCLM

lifetables <- list()
lifetables$cnst <- list(
  nage = length(unique(lt_input$openage_85$age_start)),
  nyears = config$skeleton$year$max-config$skeleton$year$min+1,
  nregions = length(config$skeleton$region)
)

lifetables$input <-
  lt_input$openage_85 %>%
  select(id, age_start, age_width, population_py,
         death_total, death_covid, nmx_cntfc,
         death_total_prop_q1, death_total_prop_q2,
         death_total_prop_q3, death_total_prop_q4,
         death_expected_prop_q1, death_expected_prop_q2,
         death_expected_prop_q3, death_expected_prop_q4,
         nweeks_year)
# vector ordered by sex, region, year, age distributed into
# 8D array [age, year, region_id, sex, sim_id, var_id, projected, quarter]
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
# (28) <ex_deficit>         Life-expectancy difference between projected and observed
lifetables$simulation <-
  array(
    dim = c(
      age = lifetables$cnst$nage,
      year = lifetables$cnst$nyears,
      region_iso = lifetables$cnst$nregions,
      sex = 3,
      sim_id = cnst$n_sim+1,
      var_id = 28,
      projected = 2,
      quarter = 5
    ),
    dimnames = list(
      0:(lifetables$cnst$nage-1),
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
        'e0_cntrb_t_covid', 'e0_cntrb_t_noncovid', 'bbi', 'ex_diff_2_year',
        'ex_deficit'),
      c('actual', 'projected'),
      c('annual', 'Q1', 'Q2', 'Q3', 'Q4')
    )
  )

# actual observables
lifetables$simulation[,,,,,'population_py','actual','annual'] <-
  lifetables$input$population_py
lifetables$simulation[,,,,,'population_py','actual','Q1'] <-
  lifetables$input$population_py*
  (length(config$weeks_in_quarter$q1)/lifetables$input$nweeks_year)
lifetables$simulation[,,,,,'population_py','actual','Q2'] <-
  lifetables$input$population_py*
  (length(config$weeks_in_quarter$q2)/lifetables$input$nweeks_year)
lifetables$simulation[,,,,,'population_py','actual','Q3'] <-
  lifetables$input$population_py*
  (length(config$weeks_in_quarter$q3)/lifetables$input$nweeks_year)
lifetables$simulation[,,,,,'population_py','actual','Q4'] <-
  lifetables$input$population_py*(
    ifelse(lifetables$input$nweeks_year==53, 14, 13)/
      lifetables$input$nweeks_year
  )

lifetables$simulation[,,,,1,'death_total','actual','annual'] <-
  lifetables$input$death_total
lifetables$simulation[,,,,1,'death_total','actual','Q1'] <-
  lifetables$input$death_total*lifetables$input$death_total_prop_q1
lifetables$simulation[,,,,1,'death_total','actual','Q2'] <-
  lifetables$input$death_total*lifetables$input$death_total_prop_q2
lifetables$simulation[,,,,1,'death_total','actual','Q3'] <-
  lifetables$input$death_total*lifetables$input$death_total_prop_q3
lifetables$simulation[,,,,1,'death_total','actual','Q4'] <-
  lifetables$input$death_total*lifetables$input$death_total_prop_q4

# we don't analyze the within year covid death distribution
# so the quarter-by-quarter fields remain NA
lifetables$simulation[,,,,1,'death_covid','actual','annual'] <-
  lifetables$input$death_covid

# projected observables based on 5 year average nmx change
# here we store all-cause death counts how we would
# expect them under continuing pre-pandemic trends
lifetables$simulation[,,,,,'population_py','projected','annual'] <-
  lifetables$simulation[,,,,,'population_py','actual','annual']
lifetables$simulation[,,,,,'population_py','projected','Q1'] <-
  lifetables$simulation[,,,,,'population_py','actual','Q1']
lifetables$simulation[,,,,,'population_py','projected','Q2'] <-
  lifetables$simulation[,,,,,'population_py','actual','Q2']
lifetables$simulation[,,,,,'population_py','projected','Q3'] <-
  lifetables$simulation[,,,,,'population_py','actual','Q3']
lifetables$simulation[,,,,,'population_py','projected','Q4'] <-
  lifetables$simulation[,,,,,'population_py','actual','Q4']

lifetables$simulation[,,,,1,'death_total','projected','annual'] <-
  lifetables$input$nmx_cntfc*
  lifetables$simulation[,,,,1,'population_py','projected','annual']
lifetables$simulation[,,,,1,'death_total','projected','Q1'] <-
  # distribute annual expected deaths to quarters by multiplying
  # with expected proportion of deaths by quarter as calculated
  # in a prior step from pre-covid distribution of deaths within a year
  # by age, sex, and region
  c(lifetables$simulation[,,,,1,'death_total','projected','annual'])*
  lifetables$input$death_expected_prop_q1
lifetables$simulation[,,,,1,'death_total','projected','Q2'] <-
  c(lifetables$simulation[,,,,1,'death_total','projected','annual'])*
  lifetables$input$death_expected_prop_q2
lifetables$simulation[,,,,1,'death_total','projected','Q3'] <-
  c(lifetables$simulation[,,,,1,'death_total','projected','annual'])*
  lifetables$input$death_expected_prop_q3
lifetables$simulation[,,,,1,'death_total','projected','Q4'] <-
  c(lifetables$simulation[,,,,1,'death_total','projected','annual'])*
  lifetables$input$death_expected_prop_q4

lifetables$simulation[,,,,,'death_covid','projected','annual'] <- 0

# simulate total death counts
lifetables$simulation[,,,,-1,'death_total',,] <-
  apply(lifetables$simulation[,,,,1,'death_total',,],
        MARGIN = 1:6, function (lambda) rpois(n = cnst$n_sim, lambda),
        simplify = TRUE) %>%
  aperm(c(2,3,4,5,1,6,7))
# simulate covid death counts
lifetables$simulation[,,,,-1,'death_covid','actual','annual'] <-
  apply(lifetables$simulation[,,,,1,'death_covid','actual','annual'],
        MARGIN = 1:4, function (lambda) rpois(n = cnst$n_sim, lambda),
        simplify = TRUE) %>%
  aperm(c(2,3,4,5,1))

# Calculate lifetables over simulated counts ----------------------

# nmx
lifetables$simulation[,,,,,'nmx',,] <-
  lifetables$simulation[,,,,,'death_total',,] /
  lifetables$simulation[,,,,,'population_py',,]

# npx, using constant hazard assumption,
# i.e. npx = exp(-nmx)) for single year age groups
lifetables$simulation[,,,,,'npx',,] <-
  exp(-lifetables$simulation[,,,,,'nmx',,])
lifetables$simulation[lifetables$cnst$nage,,,,,'npx',,] <- 0

# no need for fancy nax adjustment, I checked, virtually
# same results as with PWE
# nax <- lifetables$simulation[,,,,,'nmx',]
# nax[1:101,,,,,] <- 0.5
# I <- lifetables$simulation[1,,,'F',,'nmx',]<0.107
# I[is.na(I)] <- FALSE
# nax[1,,,'F',,][I] <-
#   0.053+2.8*lifetables$simulation[1,,,'F',,'nmx',][I]
# nax[1,,,'F',,][!I] <- 0.350
# I <- lifetables$simulation[1,,,'M',,'nmx',]<0.107
# I[is.na(I)] <- FALSE
# nax[1,,,'M',,][I] <-
#   0.045+2.684*lifetables$simulation[1,,,'M',,'nmx',][I]
# nax[1,,,'M',,][!I] <- 0.330
# lifetables$simulation[,,,,,'nqx',] <- 
#   lifetables$simulation[,,,,,'nmx',] /
#   (1+(1-nax)*lifetables$simulation[,,,,,'nmx',])
# lifetables$simulation[lifetables$cnst$nage,,,,,'nqx',] <- 1
# lifetables$simulation[,,,,,'npx',] <-
#   1-lifetables$simulation[,,,,,'nqx',]

# nqx
lifetables$simulation[,,,,,'nqx',,] <-
  1-lifetables$simulation[,,,,,'npx',,]
# lx
lifetables$simulation[,,,,,'lx',,] <-
  apply(
    lifetables$simulation[,,,,,'npx',,],
    # apply function to vector of data by age
    MARGIN = 2:7, function (npx) head(cumprod(c(1, npx)), -1)
  )
# ndx
lifetables$simulation[,,,,,'ndx',,] <-
  apply(
    lifetables$simulation[,,,,,'lx',,],
    # apply function to vector of data by age
    MARGIN = 2:7, function (lx) c(-diff(lx), tail(lx, 1))
  )
# nLx = ifelse(mx==0, lx*nx, ndx/nmx)
lifetables$simulation[,,,,,'nLx',,] <-
  lifetables$simulation[,,,,,'ndx',,]/lifetables$simulation[,,,,,'nmx',,]
tmp$I <- compareNA(lifetables$simulation[,,,,,'nmx',,],0)
lifetables$simulation[,,,,,'nLx',,][tmp$I] <-
  lifetables$simulation[,,,,,'lx',,][tmp$I]
# Tx = rev(cumsum(rev(nLx)))
lifetables$simulation[,,,,,'Tx',,] <-
  apply(
    lifetables$simulation[,,,,,'nLx',,],
    # apply function to vector of data by age
    MARGIN = 2:7, function (nLx) rev(cumsum(rev(nLx)))
  )
# ex = Tx/lx
lifetables$simulation[,,,,,'ex',,] <-
  lifetables$simulation[,,,,,'Tx',,] /
  lifetables$simulation[,,,,,'lx',,]

# Calculate ex deficit --------------------------------------------

lifetables$simulation[,,,,,'ex_deficit','actual',] <-
  lifetables$simulation[,,,,,'ex','projected',] -
  lifetables$simulation[,,,,,'ex','actual',]
lifetables$simulation[,,,,,'ex_deficit','projected',] <-
  lifetables$simulation[,,,,,'ex','projected',] -
  lifetables$simulation[,,,,,'ex','actual',]

# Calculate annual ex change --------------------------------------

lifetables$simulation[,,,,,'ex_lag',,] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'ex',,]
lifetables$simulation[,,,,,'ex_diff',,] <-
  lifetables$simulation[,,,,,'ex',,]-lifetables$simulation[,,,,,'ex_lag',,]
lifetables$simulation[,,,,,'ex_diff_lag',,] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'ex_diff',,]

# Calculate Arriaga decomposition ---------------------------------

# decompose annual changes in e0 into age specific mortality changes

# see Arriaga (1984)
# Measuring and explaining the change in life expectancies
# DOI 10.2307/2061029

lifetables$simulation[,,,,,'lx_lag',,] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'lx',,]
lifetables$simulation[,,,,,'nLx_lag',,] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'nLx',,]
lifetables$simulation[,,,,,'Tx_lag',,] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'Tx',,]

lifetables$simulation[,,,,,'e0_cntrb_d',,] <-
  (
    lifetables$simulation[,,,,,'nLx',,]/lifetables$simulation[,,,,,'lx',,]-
      lifetables$simulation[,,,,,'nLx_lag',,]/lifetables$simulation[,,,,,'lx_lag',,]
  ) * lifetables$simulation[,,,,,'lx',,]

lifetables$simulation[,,,,,'e0_cntrb_i',,] <-
  (
    lifetables$simulation[,,,,,'lx_lag',,]/
      lifetables$simulation[,,,,,'lx',,]-
      apply(lifetables$simulation[,,,,,'lx_lag',,],
            # apply function to vector of data by age
            2:7, function (x) c(x[-1], 0))/
      apply(lifetables$simulation[,,,,,'lx',,],
            # apply function to vector of data by age
            2:7, function (x) c(x[-1], 0))
  ) * apply(lifetables$simulation[,,,,,'Tx',,],
            # apply function to vector of data by age
            2:7, function (x) c(x[-1], 0))
lifetables$simulation[lifetables$cnst$nage,,,,,'e0_cntrb_i',,] <- 0

lifetables$simulation[,,,,,'e0_cntrb_t',,] <-
  lifetables$simulation[,,,,,'e0_cntrb_d',,] +
  lifetables$simulation[,,,,,'e0_cntrb_i',,]

# Calculate cause of death decomposition --------------------------

lifetables$simulation[,,,,,'nmx_lag',,] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'nmx',,]

# R_covid = D_covid / D_total
lifetables$simulation[,,,,,'R_covid',,] <-
  lifetables$simulation[,,,,,'death_covid',,]/lifetables$simulation[,,,,,'death_total',,]
# R_noncovid = 1 - R_covid
lifetables$simulation[,,,,,'R_covid_lag',,] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'R_covid',,]

lifetables$simulation[,,,,,'e0_cntrb_t_covid',,] <-
  lifetables$simulation[,,,,,'e0_cntrb_t',,] * (
    (
      lifetables$simulation[,,,,,'R_covid',,] *
        lifetables$simulation[,,,,,'nmx',,] -
        lifetables$simulation[,,,,,'R_covid_lag',,] *
        lifetables$simulation[,,,,,'nmx_lag',,]
    ) /
      (lifetables$simulation[,,,,,'nmx',,] -
         lifetables$simulation[,,,,,'nmx_lag',,])
  )

lifetables$simulation[,,,,,'e0_cntrb_t_noncovid',,] <-
  lifetables$simulation[,,,,,'e0_cntrb_t',,] -
  lifetables$simulation[,,,,,'e0_cntrb_t_covid',,]

# Calculate additional indicators ---------------------------------

# bounce back indicator
lifetables$simulation[,,,,,'bbi',,] <-
  1 - (
    lifetables$simulation[,,,,,'ex_diff_lag',,] +
      lifetables$simulation[,,,,,'ex_diff',,]
  ) /
  lifetables$simulation[,,,,,'ex_diff_lag',,]
tmp$I <- compareNA(sign(lifetables$simulation[,,,,,'ex_diff_lag',,]), 1)
lifetables$simulation[,,,,,'bbi',,][tmp$I] <-
  -lifetables$simulation[,,,,,'bbi',,][tmp$I]

# 2 year life expectancy difference
lifetables$simulation[,,,,,'ex_diff_2_year',,] <-
  lifetables$simulation[,,,,,'ex_diff',,] +
  lifetables$simulation[,,,,,'ex_diff_lag',,]

# Calculate sex differences ---------------------------------------

sexdiff <- list()
sexdiff$simulation <-
  lifetables$simulation[,,,'F',,,,]-
  lifetables$simulation[,,,'M',,,,]

# add additional variables to array
# var_id:
# (29) <ex_diff_sign>                sign of sex difference in ex change
# (30) <ex_diff_change_from_2019>    change in sex difference from 2019
# (31) <ex_diff_drop_from_2019_flag> drop in sex difference from 2019
# (32) <ex_diff_rise_from_2019_flag> rise in sex difference from 2019
D <- dim(sexdiff$simulation)
D[5] <- D[5]+4
Dn <- dimnames(sexdiff$simulation)
Dn[[5]] <- c(
  Dn[[5]], 'ex_diff_sign', 'ex_diff_change_from_2019',
  'ex_diff_drop_from_2019_flag',
  'ex_diff_rise_from_2019_flag'
)
sexdiff$temp <- array(NA, D, Dn)
sexdiff$temp[,,,,-(29:32),,] <- sexdiff$simulation
sexdiff$simulation <- sexdiff$temp

# +2: ex increase for both sexes
# -2: ex decrease for both sexes
# 0: mixed increase-decrease by sex
sexdiff$simulation[,,,,'ex_diff_sign',,] <-
  sign(lifetables$simulation[,,,'F',,'ex_diff',,])+
  sign(lifetables$simulation[,,,'M',,'ex_diff',,])

sexdiff$simulation[,,,,'ex_diff_change_from_2019',,] <-
  sexdiff$simulation[,,,,'ex_diff',,] -
  sexdiff$simulation[,rep('2019', dim(sexdiff$simulation)[2]),,,'ex_diff',,]

sexdiff$simulation[,,,,'ex_diff_drop_from_2019_flag',,] <-
  sexdiff$simulation[,,,,'ex_diff_change_from_2019',,] < 0

sexdiff$simulation[,,,,'ex_diff_rise_from_2019_flag',,] <-
  sexdiff$simulation[,,,,'ex_diff_change_from_2019',,] >= 0

# Calculate average ex change 2016 to 2019 ------------------------

e0avgdiff <- list()

e0avgdiff$simulation <- apply(
  lifetables$simulation[,2:5,,,,'ex_diff',,],
  # apply function to vector of data by year
  MARGIN = c(1, 3:7), mean
) %>%
  aperm(c(1,3,2,4,5,6))

# Calculate cause contribution to annual e0 change ----------------

# simply sum the age specific cause contributions

codecomp <- list()

codecomp$simulation <-
  apply(
    lifetables$simulation[
      ,,,,,
      c('e0_cntrb_t_covid', 'e0_cntrb_t_noncovid', 'e0_cntrb_t'),,
    ],
    # apply function to vector of data by age
    MARGIN = c(2:8), FUN = sum
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
      var_id = 22,
      quarter = 5
    ),
    dimnames = list(
      0:(lifetables$cnst$nage-1),
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
        'e0_cntrb_t_covid', 'e0_cntrb_t_noncovid'),
      c('annual', 'Q1', 'Q2', 'Q3', 'Q4')
    )
  )

arriaga_cntfc$simulation[,,,,,'death_total_actual',] <-
  lifetables$simulation[,,,,,'death_total','actual',]
arriaga_cntfc$simulation[,,,,,'death_total_expected',] <-
  lifetables$simulation[,,,,,'death_total','projected',]
arriaga_cntfc$simulation[,,,,,'death_covid_actual',] <-
  lifetables$simulation[,,,,,'death_covid','actual',]
arriaga_cntfc$simulation[,,,,,'death_covid_expected',] <-
  lifetables$simulation[,,,,,'death_covid','projected',]

arriaga_cntfc$simulation[,,,,,'nmx_actual',] <-
  lifetables$simulation[,,,,,'nmx','actual',]
arriaga_cntfc$simulation[,,,,,'nmx_expected',] <-
  lifetables$simulation[,,,,,'nmx','projected',]

arriaga_cntfc$simulation[,,,,,'lx_actual',] <-
  lifetables$simulation[,,,,,'lx','actual',]
arriaga_cntfc$simulation[,,,,,'lx_expected',] <-
  lifetables$simulation[,,,,,'lx','projected',]

arriaga_cntfc$simulation[,,,,,'nLx_actual',] <-
  lifetables$simulation[,,,,,'nLx','actual',]
arriaga_cntfc$simulation[,,,,,'nLx_expected',] <-
  lifetables$simulation[,,,,,'nLx','projected',]

arriaga_cntfc$simulation[,,,,,'Tx_actual',] <-
  lifetables$simulation[,,,,,'Tx','actual',]
arriaga_cntfc$simulation[,,,,,'Tx_expected',] <-
  lifetables$simulation[,,,,,'Tx','projected',]

arriaga_cntfc$simulation[,,,,,'ex_actual',] <-
  lifetables$simulation[,,,,,'ex','actual',]
arriaga_cntfc$simulation[,,,,,'ex_expected',] <-
  lifetables$simulation[,,,,,'ex','projected',]

arriaga_cntfc$simulation[,,,,,'ex_actual_minus_expected',] <-
  arriaga_cntfc$simulation[,,,,,'ex_actual',] -
  arriaga_cntfc$simulation[,,,,,'ex_expected',]

arriaga_cntfc$simulation[,,,,,'e0_cntrb_d',] <-
  (
    arriaga_cntfc$simulation[,,,,,'nLx_actual',]/
      arriaga_cntfc$simulation[,,,,,'lx_actual',] -
      arriaga_cntfc$simulation[,,,,,'nLx_expected',]/
      arriaga_cntfc$simulation[,,,,,'lx_expected',]
  ) * arriaga_cntfc$simulation[,,,,,'lx_actual',]

arriaga_cntfc$simulation[,,,,,'e0_cntrb_i',] <-
  (
    arriaga_cntfc$simulation[,,,,,'lx_expected',]/
      arriaga_cntfc$simulation[,,,,,'lx_actual',]-
      apply(arriaga_cntfc$simulation[,,,,,'lx_expected',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))/
      apply(arriaga_cntfc$simulation[,,,,,'lx_actual',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))
  ) * apply(arriaga_cntfc$simulation[,,,,,'Tx_actual',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))
arriaga_cntfc$simulation[lifetables$cnst$nage,,,,,'e0_cntrb_i',] <- 0

arriaga_cntfc$simulation[,,,,,'e0_cntrb_t',] <-
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_d',] +
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_i',]

arriaga_cntfc$simulation[,,,,,'R_covid_actual',] <-
  lifetables$simulation[,,,,,'R_covid','actual',]
arriaga_cntfc$simulation[,,,,,'R_covid_expected',] <-
  lifetables$simulation[,,,,,'R_covid','projected',]

arriaga_cntfc$simulation[,,,,,'e0_cntrb_t_covid',] <-
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_t',] * (
    (
      arriaga_cntfc$simulation[,,,,,'R_covid_actual',] *
        arriaga_cntfc$simulation[,,,,,'nmx_actual',] -
        arriaga_cntfc$simulation[,,,,,'R_covid_expected',] *
        arriaga_cntfc$simulation[,,,,,'nmx_expected',]
    ) /
      (arriaga_cntfc$simulation[,,,,,'nmx_actual',] -
         arriaga_cntfc$simulation[,,,,,'nmx_expected',])
  )

arriaga_cntfc$simulation[,,,,,'e0_cntrb_t_noncovid',] <-
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_t',] -
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_t_covid',]

# Calculates CI over simulations ----------------------------------

# ci's for life tables
lifetables$ci <-
  apply(
    lifetables$simulation[
      ,,,,-1,
      c('nmx', 'npx', 'nqx', 'lx', 'ex', 'ex_diff', 'ex_diff_2_year', 'bbi',
        'e0_cntrb_t', 'e0_cntrb_t_covid', 'e0_cntrb_t_noncovid', 'ex_deficit'),,
    ],
    -5,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:8,1))
V <- dimnames(lifetables$ci)
names(attr(lifetables$ci, 'dim'))[8] <- 'quantile'
dimnames(lifetables$ci) <- V

# ci's for sex-differences of life table statistics
sexdiff$ci <-
  apply(
    sexdiff$simulation[
      ,,,-1,
      c('nmx', 'npx', 'nqx', 'lx', 'ex', 'ex_diff', 'bbi', 'ex_diff_sign',
        'ex_diff_change_from_2019',
        'ex_diff_drop_from_2019_flag',
        'ex_diff_rise_from_2019_flag'),,
    ],
    -4,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:7,1))
V <- dimnames(sexdiff$ci)
names(attr(sexdiff$ci, 'dim'))[7] <- 'quantile'
dimnames(sexdiff$ci) <- V

# ci's for average e0 annual change
e0avgdiff$ci <-
  apply(
    e0avgdiff$simulation[,,,-1,,],
    -4,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:6,1))
V <- dimnames(e0avgdiff$ci)
names(attr(e0avgdiff$ci, 'dim'))[6] <- 'quantile'
dimnames(e0avgdiff$ci) <- V

# ci's for cause contributions to annual life expectancy change
codecomp$ci <-
  apply(
    codecomp$simulation[,,,-1,,,],
    -4,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:7,1))
V <- dimnames(codecomp$ci)
names(attr(codecomp$ci, 'dim'))[7] <- 'quantile'
dimnames(codecomp$ci) <- V

# ci's for age specific contributions to e0 deviations from expectation
arriaga_cntfc$ci <-
  apply(
    arriaga_cntfc$simulation[,,,,-1,,],
    -5,
    QuantileWithMean,
    simplify = TRUE
  ) %>%
  aperm(c(2:7,1))
V <- dimnames(arriaga_cntfc$ci)
names(attr(arriaga_cntfc$ci, 'dim'))[7] <- 'quantile'
dimnames(arriaga_cntfc$ci) <- V

# Transform to data frame -----------------------------------------

lifetables$ci_df <-
  as.data.frame.table(lifetables$ci, stringsAsFactors = FALSE)
names(lifetables$ci_df) <-
  c(names(attr(lifetables$ci, 'dim')), 'value')
lifetables$ci_df <-
  lifetables$ci_df %>%
  as_tibble() %>%
  pivot_wider(id_cols = c(age, year, region_iso, sex, projected, quarter),
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
  pivot_wider(id_cols = c(age, year, region_iso, projected, quarter),
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
  pivot_wider(id_cols = c(age, sex, region_iso, projected, quarter),
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
  pivot_wider(id_cols = c(sex, region_iso, year, projected, quarter),
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
  pivot_wider(id_cols = c(sex, region_iso, year, age, quarter),
              names_from = c(var_id, quantile),
              values_from = value)

# Test ------------------------------------------------------------

lifetables$simulation['0','2021',,'T',1,'population_py','projected','annual'] ==
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
  group_by(region_iso, sex, year, projected, quarter) %>%
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
lifetables$ci_df %>%
  mutate(across(.cols = where(is.numeric), .fns = ~round(.x,6))) %>%
  write_csv(paths$output$lifetables_csv)
saveRDS(sexdiff$ci_df, paths$output$sexdiff)
sexdiff$ci_df %>%
  mutate(across(.cols = where(is.numeric), .fns = ~round(.x,6))) %>%
  write_csv(paths$output$sexdiff_csv)
saveRDS(e0avgdiff$ci_df, paths$output$e0avgdiff)
e0avgdiff$ci_df %>%
  mutate(across(.cols = where(is.numeric), .fns = ~round(.x,6))) %>%
  write_csv(paths$output$e0avgdiff_csv)
saveRDS(codecomp$ci_df, paths$output$codecomp)
codecomp$ci_df %>%
  mutate(across(.cols = where(is.numeric), .fns = ~round(.x,6))) %>%
  write_csv(paths$output$codecomp_csv)
saveRDS(arriaga_cntfc$ci_df, paths$output$arriaga_cntfc)
arriaga_cntfc$ci_df %>%
  mutate(across(.cols = where(is.numeric), .fns = ~round(.x,6))) %>%
  write_csv(paths$output$arriaga_cntfc_csv)
