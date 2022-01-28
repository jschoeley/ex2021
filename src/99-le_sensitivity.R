# This R script is for sensitivity analysis. We examined the sensitivity of our
# lifetable estimation to the denominators. We mainly use a self-defined function
# to read and run R codes from 40-lifetable_estimation.R and 50-figure_e0diff.R.
# The former is used for lifetable estimation whereas the latter the visualization
# of the changes in life expectancy.

library(here); library(stringr)
setwd(here())

# Step 1. Define a function for reading R codes from part of a file ------------
source_part <- function(file, from, to) {
  source = scan(
    file,
    what = character(),
    sep = '\n',
    skip = from - 1, n = to - from + 1,
    encoding = 'UTF-8',
    quiet = TRUE,
    blank.lines.skip = FALSE
  )
  withVisible(eval.parent(parse(text = source)))
}



# Step 2. Lifetable estimation -------------------------------------------------

# (1) Init, constants, and functions
le <- 'src/40-lifetable_estimation.R'
source_part(le, 1, 100)


# (2) Use NSO pop instead of the WPP data for lifetable estimation
## load NSO pop data 
pop_nso <- readRDS('tmp/pop_nso.rds')
lt_input <- readRDS('out/lt_input.rds')

## metadata on the completeness of the 2021 death data
completeness <-
  lt_input %>%
  filter(year == 2021, age_start == 0, sex == 'Male') %>%
  select(year, region_iso, death_total_nweeksmiss) %>%
  mutate(as_of = 52-death_total_nweeksmiss)

## adjust exposures for the length of the observation period within 2021
pop_alt <- pop_nso %>%
  rename(region_iso = region_code_iso3166_2) %>%
  left_join(completeness, by = c('region_iso', 'year')) %>%
  mutate(
    population_py = ifelse(year==2021, midpop*as_of/52, midpop)
  ) 

## calculate the exposures for both sexes combined
pop_total <- pop_alt %>% 
  group_by(country, year, age_start) %>%
  summarise(population_py = sum(population_py)) %>% 
  ungroup() %>%
  mutate(sex='Total') %>% 
  left_join(pop_alt %>% 
              mutate(id = paste0(str_sub(id,1,6), 'T', str_sub(id,8,14))) %>%
              select(-sex, -contains('pop')),
            by = c('country', 'year', 'age_start')
  ) %>% 
  unique() %>%
  rbind(pop_alt %>% select(-pop, -midpop) ) %>%
  arrange(country, year, sex)

## use the adjusted NSO exposures for lifetable estimation 
lifetables$input <- lifetables$input %>%
  select(-population_py) %>%
  left_join(pop_total %>% 
              filter(!(country %in% c('Bulgaria','Greece')),
                     !(country=='Switzerland' & year==2021)) %>%
              select('id',population_py), by='id')


# (3) Create Poisson replicates of counts
source_part(le, 125, 165)


# (4) Calculate lifetables over simulated counts
source_part(le, 169, 206)


# (5) Calculate annual ex change
source_part(le, 210, 228)


# (6) Calculate Arriaga decomposition
source_part(le, 232, 256)


# (7) Calculate sex differences
source_part(le, 260, 279)


# (8) Calculates CI over simulations
source_part(le, 283, 313)


# (9) Test
lifetables$simulation['0','2021',,'T',1,'population_py'] ==
  lifetables$input %>%
  filter(str_sub(id, 7, 14)=='T2021000') %>%
  select(population_py)


# (10) Transform output to data frame
source_part(le, 324, 348)


# (11) Export
saveRDS(lifetables$ci_df, 'out/lifetables_alt.rds')
saveRDS(sexdiff$ci_df, 'out/sexdiff_alt.rds')



# Step 3. Ex0diff visualization ------------------------------------------------

# (1) Init, constants, and functions
fig_e0diff <- 'src/50-figure_e0diff.R'
source_part(fig_e0diff, 1, 118)


# (2) Load and prepare data
dat$lifetables <- lifetables$ci_df
dat$lt_input <- lt_input
dat$completeness <- completeness
source_part(fig_e0diff, 134, 153)


# (3) Parameterize figure
source_part(fig_e0diff, 157, 183)


# (4) Plot
source_part(fig_e0diff, 185, 490)


# (5) Arriaga decomposition
source_part(fig_e0diff, 494, 542)


# (6) Export
## export the ex0diff figure
fig_spec$ExportFigure(
  fig$e0diff$plot, device = 'svg',
  filename = 'e0diff_alt',
  path = paths$output$out,
  width = 180, height = 200
)

## export the arriaga decomposition figure
fig_spec$ExportFigure(
  fig$arriaga$plot, device = 'svg',
  filename = 'arriaga_alt',
  path = paths$output$out,
  width = 180, height = 200
)

