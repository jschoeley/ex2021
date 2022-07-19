# Check specific statements made in paper and quantify uncertainty

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
  lifetables_sim = './tmp/40-lifetables_sim.rds',
  sexdiff_sim = './tmp/40-sexdiff_sim.rds',
  cntfc_lt_debug = './tmp/24-counterfactual_lt_debug.rds',
  skeleton = './tmp/10-harmonized_skeleton.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  statements = './out/92-statements.rds'
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
statements <- list()

# Functions -------------------------------------------------------

# figure specifications
source(paths$input$figspec)

# Calculate two sided P values around bootstrapped
# effect estimates based on shifted empirical bootstrap
# distribution
# x: vector of bootstrapped effect estimates, first element
# of vector is the mean effect estimate
BootH0ECDF <- function (x, tail = 'two', null_value = 0) {
  # assume the the distribution of the null hypothesis
  # is the shifted bootstrap distribution of the estimate
  test <- (x-x[1]+null_value)[-1]
  p_ <- if (!any(is.na(test))) {
    stats::ecdf(test)(x[1]) 
  } else NA
  if (identical(tail, 'two')) {
    p <- 2*min(p_, 1-p_)
  }
  if (identical(tail, 'right')) {
    p <- 1-p_
  }
  if (identical(tail, 'left')) {
    p <- p_
  }
  
  return(p)
}
# Calculate 95% CI around bootstrapped effect estimates
# based on empirical bootstrap distribution
# x: vector of bootstrapped effect estimates, first element
# of vector is the mean effect estimate
Boot95CI <- function (x) {
  quantile(x[-1], probs = c(0.025, 0.975), na.rm = TRUE)
}

# Given a table with columns
# 'region', 'effect', 'lo', 'hi', 'p', create strings of text
# describing the regional estimates
DescribeRegionalEstimates <- function (tab, digits = 1) {
  paste0(
    tab$region, ' (',
    formatC(tab$effect, format = 'f', digits = digits, flag = '+'),
    ', CI ',
    formatC(tab$lo, format = 'f', digits = digits, flag = '+'), ' to ',
    formatC(tab$hi, format = 'f', digits = digits, flag = '+'), ', $p',
    ifelse(
      tab$p < 0.001, '< 0.001',
      paste0('=', formatC(tab$p, format = 'f', digits = 3))
    ),
    '$)'
  )
}

# Load data -------------------------------------------------------

# life tables
dat$lifetables <-
  readRDS(paths$input$lifetables) %>%
  left_join(region_meta, by = c(region_iso = 'region_code_iso3166_2')) %>%
  filter(region_iso %in% cnst$regions_for_analysis, quarter == 'annual')
dat$cntfc_lt_debug <- readRDS(paths$input$cntfc_lt_debug)
dat$skeleton <- readRDS(paths$input$skeleton)
# bootstrap simulations of life tables
dat$lifetables_sim <- readRDS(paths$input$lifetables_sim)
dat$lifetables_sim <-
  dat$lifetables_sim[,,cnst$regions_for_analysis,,,,,]

# bootstrap simulations of sex differences
dat$sexdiff_sim <- readRDS(paths$input$sexdiff_sim)

# Check -----------------------------------------------------------

# life expectancy loss in 2020
statements$regions_with_le_loss_in_2020 <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', year == 2020, projected == 'actual',
         quarter == 'annual', ex_diff_mean < 0) %>%
  select(region_name, region_iso)

# life expectancy loss in 2021
statements$regions_with_le_loss_in_2021 <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', year == 2021, projected == 'actual',
         quarter == 'annual', ex_diff_mean < 0) %>%
  select(region_name, region_iso)

# 2 subsequent losses
statements$regions_with_two_subsequent_le_losses <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual', quarter == 'annual') %>%
  select(region_name, region_iso, year, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_diff_mean)) %>%
  filter(`2020` < 0, `2021` < 0) %>%
  select(region_name, region_iso)

# incomplete bounce-back
statements$regions_with_incomplete_le_bounce_back <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual', quarter == 'annual') %>%
  select(region_name, region_iso, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_diff_mean_2020 < 0, ex_diff_mean_2021 > 0, ex_mean_2019 > ex_mean_2021) %>%
  select(region_name, region_iso)

# exceeding pre-pandemic LE in 2021
statements$regions_exceeding_pre_pandemic_le_in_2021 <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual', quarter == 'annual') %>%
  select(region_name, region_iso, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_mean_2019 < ex_mean_2021) %>%
  select(region_name, region_iso)

# still loss in 2021 compared with 2019
statements$regions_not_exceeding_pre_pandemic_le_in_2021 <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual', quarter == 'annual') %>%
  select(region_name, region_iso, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_mean_2019 > ex_mean_2021) %>%
  select(region_name, region_iso)

# first drop below 2019 in 2021
statements$regions_with_first_le_drop_below_2019_in_2021 <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual', quarter == 'annual') %>%
  select(region_name, region_iso, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_diff_mean_2020 >= 0, ex_diff_mean_2021 < 0, ex_mean_2019 > ex_mean_2021) %>%
  select(region_name, region_iso)

# lower LE in 2021 than expected
statements$regions_with_lower_than_expected_le_in_2021 <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', quarter == 'annual') %>%
  select(region_name, region_iso, year, projected, ex_mean) %>%
  pivot_wider(names_from = c(projected, year), values_from = c(ex_mean)) %>%
  filter(projected_2021 > actual_2021) %>%
  select(region_name, region_iso)

# any bounce back
statements$regions_with_any_le_bounce_back_in_2021 <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', projected == 'actual', quarter == 'annual') %>%
  select(region_name, region_iso, year, ex_mean, ex_diff_mean) %>%
  pivot_wider(names_from = year, values_from = c(ex_mean, ex_diff_mean)) %>%
  filter(ex_diff_mean_2020 < 0, ex_diff_mean_2021 > 0) %>%
  select(region_name, region_iso)

# younger age group contributed more than older age group to
# e0 losses in 2021 given losses in 2021
statements$regions_where_younger_age_group_contributed_most_to_le_losses <-
  dat$lifetables %>%
  filter(sex == 'T', projected == 'actual', quarter == 'annual') %>%
  mutate(
    age_group = cut(age, c(0, 60, Inf), right = FALSE,
                    labels = c('<60', '60+'))
  ) %>%
  select(
    region_name, region_iso, age, age_group, year, e0_cntrb_t_mean
  ) %>%
  group_by(region_name, region_iso, year, age_group) %>%
  summarise(e0_cntrb_t_mean = sum(e0_cntrb_t_mean)) %>%
  ungroup() %>%
  filter(year > 2019) %>%
  pivot_wider(names_from = c(year, age_group), values_from = e0_cntrb_t_mean) %>%
  filter((`2021_<60` + `2021_60+` < 0), `2021_<60` < 0, `2021_<60` < `2020_<60`) %>%
  select(region_name, region_iso)

# Changes in life expectancy since 2019 ---------------------------

# lifetable_sim is an
# 8D array [age, year, region_id, sex, sim_id, var_id, projected, quarter]

statements$le_changes <- list()
statements$le_changes$since_2019 <-
  dat$lifetables %>%
  filter(age == 0, year == 2021, sex == 'T', projected == 'actual',
         quarter == 'annual') %>%
  select(region_name, year, starts_with('ex_diff'))

# Among the 29 countries analyzed, 8 countries saw significant LE
# bounce-backs from 2020 losses:
# "Belgium (+10.8, CI +9.7 to +11.9, $p< 0.001$)"
# "Switzerland (+7.7, CI +6.4 to +8.8, $p< 0.001$)"
# "Spain (+7.6, CI +7.1 to +8.1, $p< 0.001$)"
# "France (+5.0, CI +4.4 to +5.6, $p< 0.001$)"
# "England and Wales (+2.1, CI +1.6 to +2.7, $p< 0.001$)"
# "Italy (+5.1, CI +4.6 to +5.5, $p< 0.001$)"
# "Sweden (+7.5, CI +6.0 to +8.6, $p< 0.001$)"
# "Slovenia (+3.1, CI +0.4 to +5.7, $p=0.010$)"
statements$le_changes$bounce_backs <- list(table = NA, sim = NA)
statements$le_changes$bounce_backs$table <-
  statements$le_changes$since_2019 %>%
  filter(
    ex_diff_mean > 0,
    region_name %in% statements$regions_with_le_loss_in_2020$region_name
  ) %>%
  select(region = region_name, effect = ex_diff_mean,
         lo = ex_diff_q0.025, hi = ex_diff_q0.975)
statements$le_changes$bounce_backs$sim <-
  dat$lifetables_sim[as.character(0),as.character(c(2021)),
                     statements$regions_with_le_loss_in_2020$region_iso,'T',,
                  'ex_diff',
                  'actual','annual']
statements$le_changes$bounce_backs$sim <-
  statements$le_changes$bounce_backs$sim[
    statements$le_changes$bounce_backs$sim[,1]>0,]
statements$le_changes$bounce_backs$table$p <-
  apply(statements$le_changes$bounce_backs$sim, MARGIN = 1,
        BootH0ECDF, tail = 'right')
statements$le_changes$bounce_backs$table <-
  statements$le_changes$bounce_backs$table %>%
  mutate(across(c(effect, lo, hi), ~.*12))
statements$le_changes$bounce_backs$text <-
  statements$le_changes$bounce_backs$table %>%
  filter(p<=0.05) %>%
  DescribeRegionalEstimates()

# Significant life expectancy losses in two subsequent years were
# suffered by
# Bulgaria (-25.1, CI -26.6 to -23.4, $p< 0.001$)
# Chile (-8.0, CI -9.0 to -7.0, $p< 0.001$)
# Czech Republic (-10.4, CI -11.5 to -9.4, $p< 0.001$)
# Germany (-3.1, CI -3.5 to -2.7, $p< 0.001$)
# Estonia (-21.5, CI -25.1 to -17.6, $p< 0.001$)
# Greece (-12.4, CI -13.8 to -11.0, $p< 0.001$)
# Croatia (-11.6, CI -13.3 to -9.7, $p< 0.001$)
# Hungary (-16.4, CI -17.6 to -15.3, $p< 0.001$)
# Lithuania (-7.9, CI -10.5 to -5.4, $p< 0.001$)
# Poland (-12.1, CI -12.7 to -11.3, $p< 0.001$)
# Slovakia (-23.9, CI -25.7 to -22.3, $p< 0.001$)
# USA (-2.7, CI -3.1 to -2.2, $p< 0.001$)
statements$le_changes$compound_loss <- list(table = NA, sim = NA)
statements$le_changes$compound_loss$table <-
  statements$le_changes$since_2019 %>%
  filter(
    region_name %in% statements$regions_with_two_subsequent_le_losses$region_name
  ) %>%
  select(region = region_name, effect = ex_diff_mean,
         lo = ex_diff_q0.025, hi = ex_diff_q0.975)
statements$le_changes$compound_loss$sim <-
  dat$lifetables_sim[as.character(0),as.character(c(2021)),
                     statements$regions_with_two_subsequent_le_losses$region_iso,
                     'T',,'ex_diff',
                     'actual','annual']
statements$le_changes$compound_loss$table$p <-
  apply(statements$le_changes$compound_loss$sim, MARGIN = 1,
        BootH0ECDF, tail = 'left')
statements$le_changes$compound_loss$table <-
  statements$le_changes$compound_loss$table %>%
  mutate(across(c(effect, lo, hi), ~.*12))
statements$le_changes$compound_loss$text <-
  statements$le_changes$compound_loss$table %>%
  filter(p<=0.05) %>%
  DescribeRegionalEstimates()

statements$le_changes$since_2019_complete <-
  dat$lifetables %>%
  filter(age == 0, year %in% 2019:2021, sex == 'T', projected == 'actual',
         quarter == 'annual') %>%
  select(region_name, year, starts_with('ex')) %>%
  pivot_wider(names_from = year, values_from = starts_with('ex')) %>%
  mutate(across(starts_with('ex'), ~.*12))

# Age contributions to LE changes ---------------------------------

statements$age <- list()

# Mortality increases among ages below 60 contributed LE losses of
# -7.2 months (CI -7.0 to -7.4, $H_0: \mu \geq 0$, $p < 0.001$) in 2021
# compared to 2020.
statements$age$us_le_losses_grew_due_to_younger_age_group <-
  list(sim = NA, effect = NA, CI = NA, p = NA)
statements$age$us_le_losses_grew_due_to_younger_age_group$sim <-
  dat$lifetables_sim[as.character(0:60), as.character(2021),'US','T',,
                     'e0_cntrb_t',
                     'actual','annual'] %>%
  colSums() %>% `*`(12)
statements$age$us_le_losses_grew_due_to_younger_age_group$effect <-
  statements$age$us_le_losses_grew_due_to_younger_age_group$sim[1]
statements$age$us_le_losses_grew_due_to_younger_age_group$CI <- Boot95CI(
  statements$age$us_le_losses_grew_due_to_younger_age_group$sim
)
statements$age$us_le_losses_grew_due_to_younger_age_group$p <- BootH0ECDF(
  statements$age$us_le_losses_grew_due_to_younger_age_group$sim,
  tail = 'left'
)

# Excess mortality among under-60s explained more than half of the
# loss in [USA] LE in 2021 compared to 2019
# 58.9 CI: 57.9 to 59.8, p < 0.001
statements$age$us_le_losse_since19_were_driven_by_young <-
  list(sim = NA, effect = NA, CI = NA, p = NA)
statements$age$us_le_losse_since19_were_driven_by_young$sim_0_60 <-
  dat$lifetables_sim[as.character(0:60), as.character(c(2020, 2021)),
                     'US','T',,
                     'e0_cntrb_t',
                     'actual','annual'] %>%
  colSums() %>% colSums() %>% `*`(12)
statements$age$us_le_losse_since19_were_driven_by_young$sim_total <-
  dat$lifetables_sim[, as.character(c(2020, 2021)),
                     'US','T',,
                     'e0_cntrb_t',
                     'actual','annual'] %>%
  colSums() %>% colSums() %>% `*`(12)
statements$age$us_le_losse_since19_were_driven_by_young$sim_share <-
  statements$age$us_le_losse_since19_were_driven_by_young$sim_0_60/
  statements$age$us_le_losse_since19_were_driven_by_young$sim_total
statements$age$us_le_losse_since19_were_driven_by_young$effect <-
  statements$age$us_le_losse_since19_were_driven_by_young$sim_share[1]
statements$age$us_le_losse_since19_were_driven_by_young$CI <-
  Boot95CI(statements$age$us_le_losse_since19_were_driven_by_young$sim_share)
statements$age$us_le_losse_since19_were_driven_by_young$p <-
  BootH0ECDF(
    statements$age$us_le_losse_since19_were_driven_by_young$sim_share,
    tail = 'right', null_value = 0
  )

# In 11 out of 16 countries with LE losses in 2021, the under-60 age
# groups contributed significantly more to LE loss in 2021 than in 2020
# ($p<0.05$).
statements$age$increasing_young_contribution_to_le_loss <-
  list(p = NA, sim = NA)
statements$age$increasing_young_contribution_to_le_loss$sim <-
  dat$lifetables_sim[as.character(0:60), as.character(c(2020,2021)),
                     statements$regions_with_le_loss_in_2021$region_iso,
                     'T',,
                     'e0_cntrb_t',
                     'actual','annual'] %>%
  colSums() %>%
  apply(c(2,3), function (X) X['2021']-X['2020'])
statements$age$increasing_young_contribution_to_le_loss$p <-
  statements$age$increasing_young_contribution_to_le_loss$sim %>%
  apply(1, BootH0ECDF, tail = 'left')
  
# Among the 13 countries which partially or completely bounced-back from
# their LE losses in 2020, 10
# (Austria, Belgium, Switzerland, Spain, France, England \& Wales,
# Italy, Netherlands, Sweden, Slovenia)
# achieved the bounce-back primarily or solely due to normalizing
# mortality among the older population
# ($H_0: \mu_\text{60+} \leq \mu_\text{<60}$, $p<0.05$).
statements$age$higher_old_contribution_to_bounce_back <-
  list(p = NA, sim = NA)
statements$age$higher_old_contribution_to_bounce_back$sim_0_60 <-
  dat$lifetables_sim[
    as.character(0:60), as.character(c(2021)),
    statements$regions_with_any_le_bounce_back_in_2021$region_iso,
    'T',,
    'e0_cntrb_t',
    'actual','annual'
  ] %>%
  colSums()
statements$age$higher_old_contribution_to_bounce_back$sim_60plus <-
  dat$lifetables_sim[
    as.character(61:100), as.character(c(2021)),
    statements$regions_with_any_le_bounce_back_in_2021$region_iso,
    'T',,
    'e0_cntrb_t',
    'actual','annual'] %>%
  colSums()
statements$age$higher_old_contribution_to_bounce_back$sim <-
  statements$age$higher_old_contribution_to_bounce_back$sim_60plus -
  statements$age$higher_old_contribution_to_bounce_back$sim_0_60
statements$age$higher_old_contribution_to_bounce_back$p <-
  apply(
    statements$age$higher_old_contribution_to_bounce_back$sim, 1,
    BootH0ECDF, tail = 'right'
  )

# LE dropped in 28 out of the 29 countries analyzed from 2019 to 2021,
# with only Norway exceeding the 2019 levels.
# Excess mortality in ages 60+ was the main or sole contributor to these
# losses in 19 out of 28 countries ($H_0: \mu \leq 50%$, $p<0.05$),
# with the USA being the prominent exception.
statements$age$le_since19_driven_by_old <-
  list(sim = NA, effect = NA, CI = NA, p = NA)
statements$age$le_since19_driven_by_old$sim_60plus <-
  dat$lifetables_sim[
    as.character(61:100), as.character(c(2020, 2021)),
    statements$regions_not_exceeding_pre_pandemic_le_in_2021$region_iso,
    'T',,
    'e0_cntrb_t',
    'actual','annual'] %>%
  colSums() %>% colSums()
statements$age$le_since19_driven_by_old$sim_total <-
  dat$lifetables_sim[
    , as.character(c(2020, 2021)),
    statements$regions_not_exceeding_pre_pandemic_le_in_2021$region_iso,
    'T',,
    'e0_cntrb_t',
    'actual','annual'] %>%
  colSums() %>% colSums()
# values above 100 mean that the LE drop was solely the result of
# increasing mortality among old population, whereas negative
# values mean that the LE drop was solely the result of increasing
# mortality among the younger population
# values <50 mean that the LE drop was mostly the result of increasing
# mortality among the young
statements$age$le_since19_driven_by_old$sim_share <-
  statements$age$le_since19_driven_by_old$sim_60plus/
  statements$age$le_since19_driven_by_old$sim_total*100
statements$age$le_since19_driven_by_old$effect <-
  statements$age$le_since19_driven_by_old$sim_share[,1]
statements$age$le_since19_driven_by_old$p <-
  apply(
    statements$age$le_since19_driven_by_old$sim_share, 1,
    BootH0ECDF, tail = 'right', null = 50
  )
statements$age$le_since19_driven_by_old$CI <-
  apply(
    statements$age$le_since19_driven_by_old$sim_share, 1,
    Boot95CI
  )
statements$age$le_since19_driven_by_old$table <-
  tibble(
    region = statements$regions_not_exceeding_pre_pandemic_le_in_2021$region_name,
    effect = statements$age$le_since19_driven_by_old$effect,
    lo = statements$age$le_since19_driven_by_old$CI[1,],
    hi = statements$age$le_since19_driven_by_old$CI[2,],
    p = statements$age$le_since19_driven_by_old$p
  )
statements$age$le_since19_driven_by_old$table %>%
  filter(p <= 0.05) %>%
  DescribeRegionalEstimates()

# Sex differences -------------------------------------------------

statements$sex <- list()

statements$sex$change_in_mortality_sex_diff_since_2019 <-
  list(sim = NA, effect = NA, CI95 = NA, p = NA)
statements$sex$change_in_mortality_sex_diff_since_2019$sim <-
  dat$sexdiff_sim[as.character(0),as.character(c(2019, 2021)),
                  cnst$regions_for_analysis,,
                  'ex',
                  'actual','annual'] %>%
  # change in LE sex difference 2019 to 2021 (months)
  apply(c(2,3), function(x) (x[2]-x[1]))
statements$sex$change_in_mortality_sex_diff_since_2019$effect <-
  statements$sex$change_in_mortality_sex_diff_since_2019$sim[,1]
statements$sex$change_in_mortality_sex_diff_since_2019$p <-
  apply(statements$sex$change_in_mortality_sex_diff_since_2019$sim,
        MARGIN = 1, BootH0ECDF)
statements$sex$change_in_mortality_sex_diff_since_2019$CI95 <-
  apply(statements$sex$change_in_mortality_sex_diff_since_2019$sim,
        MARGIN = 1, Boot95CI)
statements$sex$change_in_mortality_sex_diff_since_2019$table <-
  data.frame(
    effect = statements$sex$change_in_mortality_sex_diff_since_2019$effect,
    lo = statements$sex$change_in_mortality_sex_diff_since_2019$CI95[1,],
    hi = statements$sex$change_in_mortality_sex_diff_since_2019$CI95[2,],
    p = statements$sex$change_in_mortality_sex_diff_since_2019$p
  ) %>%
  as_tibble(rownames = 'region')

# However, our results show that the female advantage in LE significantly
# (p<0.05) increased in 16 of the 29 countries during the pandemic,
# thereby widening the sex gap (Figure 3).
statements$sex$regions_with_increasing_female_le_advantage_since_2019 <- list()
statements$sex$regions_with_increasing_female_le_advantage_since_2019$table <-
  statements$sex$change_in_mortality_sex_diff_since_2019$table %>%
  filter(effect > 0, p <= 0.05) %>%
  arrange(effect)

# The narrowing sex gaps observed in 6 countries were not significant.
statements$sex$regions_with_narrowing_sex_gap_since_2019 <- list()
statements$sex$regions_with_narrowing_sex_gap_since_2019$table <-
  statements$sex$change_in_mortality_sex_diff_since_2019$table %>%
  filter(effect < 0) %>%
  arrange(effect)

# The magnitude of the gap in 2021 varied from 3.17 years in Norway
# (CI 2.95 to 3.37, p < 0.001) to more than 9.65 years in Lithuania
# (CI 9.20 to 9.90, p < 0.001).
statements$sex$sex_diff_in_2021 <- list()
statements$sex$sex_diff_in_2021$sim <-
  dat$sexdiff_sim[as.character(0),as.character(c(2021)),
                  cnst$regions_for_analysis,,
                  'ex',
                  'actual','annual']
statements$sex$sex_diff_in_2021$effect <-
  statements$sex$sex_diff_in_2021$sim[,1]
statements$sex$sex_diff_in_2021$p <-
  apply(statements$sex$sex_diff_in_2021$sim, MARGIN = 1, BootH0ECDF)
statements$sex$sex_diff_in_2021$CI95 <-
  apply(statements$sex$sex_diff_in_2021$sim, MARGIN = 1, Boot95CI)
statements$sex$sex_diff_in_2021$table <- data.frame(
  effect = statements$sex$sex_diff_in_2021$effect,
  lo = statements$sex$sex_diff_in_2021$CI95[1,],
  hi = statements$sex$sex_diff_in_2021$CI95[2,],
  p = statements$sex$sex_diff_in_2021$p
) %>%
  as_tibble(rownames = 'region')

# The biggest increase in the sex gap was observed in the USA, where the
# gap increased by almost a year from 5.72 to 6.69 years
# (+0.97 years, CI 0.90 to 1.04, p < 0.001).
statements$sex$largest_increases_in_sex_gap_since_2019 <- list()
statements$sex$largest_increases_in_sex_gap_since_2019$table <-
  statements$sex$change_in_mortality_sex_diff_since_2019$table %>%
  filter(effect > 0) %>%
  arrange(-effect)

# Export ----------------------------------------------------------

saveRDS(statements, paths$output$statements)
