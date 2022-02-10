library(tidyverse)

lt <- readRDS('./out/lifetables.rds')

lt %>%
  select(year, region_iso, sex, age, projected, ex_mean) %>%
  filter(age == 0, sex == 'T') %>%
  #pivot_wider(names_from = projected, values_from = ex_mean) %>%
  ggplot(aes(x = year)) +
  geom_point(aes(y = ex_mean, color = projected)) +
  facet_wrap(~region_iso)
