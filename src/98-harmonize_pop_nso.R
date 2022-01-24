# This script is for appending and harmonizing population data 
# (estimates and/or projections) from national statistical offices (NSOs).

# Average pop counts are calculated for countries who report pop data as of Jan
# 1st or Dec 31st. This is to facilitate the subsequent data quality check (i.e.,
# comparing the NSO data with the WPP data).

# For Switzerland, the 2015-2020 pop estimates are reported in single years of
# age, while the 2021 pop projection is reported in 5-year age groups, so extra
# harmonization is performed for calculating the average pop counts.

# The prepared data contain 11 columns: 
#   id, which contains information on country, year, sex, and age
#   pop: original pop counts from NSOs
#   source: pe = population estimates, pp = population projection
#   region_code_iso3166_2, country, year, sex, age_start, age_width
#   date: reference date; start=Jan 1st, midyear=Jun 30th or Jul 1st, end=Dec 31st
#   midpop: calculated average pop counts 

library(here); library(glue); library(readxl); library(tidyverse); library(stringr)

# 1. Import data ---------------------------------------------------------------
wd <- here()
setwd(glue('{wd}/dat/nso_raw'))

# population estimates
pp <- list.files(pattern='*_pp.xlsx')
# population projections
pe <- list.files(pattern='*.xlsx')[!(list.files(pattern='*.xlsx') %in% pp)]
# metadata on reference date
pop_info <- read_csv('pop_info.csv') 
# metadata on regions
region_meta <- read_csv(glue('{wd}/cfg/region_metadata.csv'), na = '.') %>%
  select(region_code_iso3166_2, country=region_name)

# import raw data into a list
raw <- lapply(1:2, function(x) lapply(list(pe,pp)[[x]], read_excel)) %>%
  set_names('pe', 'pp')
raw$pe <- raw$pe %>% set_names(pe)
raw$pp <- raw$pp %>% set_names(pp)


# 2. Harmonize data ------------------------------------------------------------
cleaned <- raw %>% 
  map(~{.x <- .x %>% 
    map(~{ .x <- .x %>%
      transmute(
        year = as.numeric(Year),
        sex = ifelse(Sex=='M',1,2) %>% factor(labels = c('Male','Female')),
        # open age group at 85+
        age_start = ifelse(as.numeric(Age)>=85, 85, as.numeric(Age)),
        pop = as.numeric(Pop)
      ) %>%
      group_by(year, sex, age_start) %>%
      summarise(pop=sum(pop))
    }) %>%
    bind_rows(.id='country')
  }) %>%
  bind_rows(.id='source') %>%
  group_by(country,year,sex) %>%
  mutate(
    country = case_when(
      source=='pe' ~ substr(country,1,nchar(country)-5),
      source=='pp' ~ substr(country,1,nchar(country)-8)
    ),
    age_width = c(diff(age_start), Inf),
    source = factor(source)
  ) %>% 
  left_join(pop_info, by='country') %>%
  group_by(country, sex, age_start) %>%
  mutate(
    country = str_to_title(country),
    country = case_when(
      country=='England And Wales' ~ 'England and Wales',
      country=='Usa' ~ 'USA',
      TRUE ~ country
    ) ,
    date = factor(date),
    # create mid-year/average pop counts
    midpop = case_when(
      date=='start' ~ 1/2*(pop+lead(pop)),
      date=='end' ~ 1/2*(pop+lag(pop)), 
      date=='midyear' ~ pop)
  ) %>%
  left_join(region_meta, by='country') %>%
  mutate(
    id = GenerateRowID(region_code_iso3166_2, sex, year, age_start)
  ) %>%
  relocate(id, pop, .before = source) %>%
  relocate(region_code_iso3166_2, .before = country) %>%
  drop_na() %>%
  ungroup()

# collapse 2021 pop estimates for Switzerland to 5-year age groups
# then calculate the average pop counts for 2021
cleaned[cleaned$country=='Switzerland' & cleaned$year==2021,] <- cleaned %>%
  filter(country=='Switzerland', year==2020) %>%
  mutate(
    age_start = rep(
      c(rep(seq(0,80,5), each=5), 85),
      2),
    id = case_when(
      age_start<10 ~ paste0(str_sub(id,1,12),0,age_start),
      age_start>=10 ~ paste0(str_sub(id,1,12),age_start)
    )
  ) %>% 
  group_by(sex,age_start) %>%
  summarise(pop2=sum(pop)) %>%
  left_join(cleaned %>% filter(country=='Switzerland', year==2021),
            by = c('sex', 'age_start')) %>%
  mutate(midpop=(pop+pop2)/2) %>%
  select(-pop2) %>%
  relocate(sex, age_start, .before = age_width) %>%
  ungroup()


# 3. Export data ---------------------------------------------------------------
saveRDS(cleaned, glue('{wd}/tmp/pop_nso.rds'))