# This script is for comparing the pop data from WPP to those from the national 
# statistical offices (NSOs). The NSO of all but two countries report pop data  
# in single years of age. For Bulgaria, the NSO pop data are reported in 0, 1-4,
# and 5-year age groups. For Greece, 5-year age groups are used. For Switzerland, 
# 2015-2020 data are reported in single years of age, while 2021 pop projection
# is reported in 5-year age groups.

# Comparison is based on mid-year/average population

library(here); library(glue); library(tidyverse)
library(ggplot2); library(gridExtra)

# 1. Import data ---------------------------------------------------------------
wd <- here()
setwd(glue('{wd}/tmp'))

pop <- readRDS('harmonized_population.rds')
pop_nso <- readRDS('pop_nso.rds') 


# 2. Prepare WPP data ----------------------------------------------------------
pop_wpp <- pop %>%
  drop_na() %>%
  mutate(
    country = str_sub(id,1,6),
    year = as.numeric(str_sub(id,8,11)),
    sex = str_sub(id,7,7),
    sex = ifelse(sex=='M',1,2) %>% factor(labels = c('Male', 'Female')),
    age_start = as.numeric(str_sub(id,12,14)),
    age_start = ifelse(age_start>=85,85,age_start),
  )

# collapse data for Bulgaria to 0, 1-4, and 5-year age groups
pop_wpp[pop_wpp$country=='BG----',] <- pop_wpp[pop_wpp$country=='BG----',] %>%
  mutate(
    age_start = rep(
      c(0, rep(1,4), rep(seq(5,80,5), each=5), rep(85,16)),
      14)
  )

# collapse data for Greece to 5-year age groups
pop_wpp[pop_wpp$country=='GR----',] <- pop_wpp[pop_wpp$country=='GR----',] %>%
  mutate(
    age_start = rep(
      c(rep(seq(0,80,5), each=5), rep(85,16)),
      14)
  )

# collapse 2021 data for Switzerland to 5-year age groups
pop_wpp[pop_wpp$country=='CH----' & pop_wpp$year==2021,] <- 
  pop_wpp[pop_wpp$country=='CH----' & pop_wpp$year==2021,] %>%
  mutate(
    age_start = rep(
      c(rep(seq(0,80,5), each=5), rep(85,16)),
      2)
  ) 

pop_wpp <- pop_wpp %>%
  mutate(
    id = case_when(
      age_start>=10 ~ paste0(str_sub(id,1,11),0,age_start),
      age_start<10 ~ paste0(str_sub(id,1,11),'00',age_start)
      )
    ) %>%
  group_by(id, country, year, sex, age_start, population_source) %>%
  summarise(population_midyear=sum(population_midyear)) %>%
  ungroup() %>%
  select(id, population_midyear, population_source)


# 3. Compare pop data ----------------------------------------------------------
# join data
compare <- left_join(pop_wpp, pop_nso, by='id') %>%
  drop_na() %>%
  group_by(country, year, population_source, source, date, age_start, age_width) %>%
  summarise(
    pop1 = sum(population_midyear),
    pop2 = sum(midpop)
  ) %>%
  ungroup()

# calculate correlations
cor <- lapply(compare %>% 
                group_by(country, year) %>%
                group_split(), 
              function(x) cor(x$pop1, x$pop2)) %>%
  unlist() 

cor_dat <- data.frame(cor=cor) %>%
  cbind(compare %>% select(country, year, population_source, source, date) %>% unique()) %>%
  relocate(cor, .after = last_col())

# plot the correlations between WPP and NSO data across years by country
cor_plot <- cor_dat %>%
  ggplot(aes(year, cor)) +
  facet_wrap(facets = vars(country), nrow = 6, ncol = 5) +
  geom_point(color='red') +
  geom_line(color='grey') +
  ggtitle('WPP-NSO correlations for all age groups over 2015-2021') +
  labs( x = 'Year', y = 'Correlation') +
  scale_x_continuous(limits = c(2016,2021)) +
  scale_y_continuous(limits = c(0.94,1),
                     breaks = seq(0.94,1.00,0.02)) +
  theme(axis.text.x=element_text(angle=45,hjust = 1))

# plot the WPP data against the NSO data by year for each country
compare_list <- compare %>% 
  group_by(country) %>%
  group_split() %>%
  set_names(unique(compare$country)) 

plot_list <- list()
for (i in seq(compare_list)) {
  plot_list[[i]] <- compare_list[[i]] %>%
    ggplot(aes(pop1, pop2)) +
    facet_wrap(facets = vars(year)) +
    geom_point(color='red') +
    geom_line(color='grey') +
    geom_abline(slope = 1) +
    labs( x = 'WPP pop', y = 'NSO pop') +
    ggtitle(paste0('WPP-NSO pop data comparison: ', names(compare_list[i])))
}


# 3. Explot
pdf(glue('{wd}/out/pop_diagnostics.pdf'))
pdf.options(width = 10, height = 8)
print(cor_plot)
for (i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()

