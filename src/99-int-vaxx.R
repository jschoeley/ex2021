#===============================================================================
# 2022-01-07 -- ex 2021
# try integral vaccination
# Ilya Kashnitsky, ilya.kashnitsky@gmail.com, @ikashnitsky
#===============================================================================

library(tidyverse)
library(magrittr)
library(lubridate)
library(ggrepel)


c19 <- covid19()

owd <- read.csv("dat/owd/owid-covid-data.csv.gz")
ids <- read.csv("cfg/region_metadata.csv")

# quick own function to take the diff of any vector of an arbitrary length
diff_first_last <- function(vec) {
    clean <- vec[!is.na(vec)]
    d <- last(clean) %>% subtract(first(clean))
    return(d)
}

# subset 2021 data and calculate integral vaxx measures
eu21 <- owd %>%
    filter(year(date) == "2021", iso_code %in% ids$region_code_iso3166_1_alpha3) %>%
    group_by(iso_code, location) %>%
    summarise(
        vac_int = people_fully_vaccinated_per_hundred %>% divide_by(100) %>% mean(na.rm = T),
        vac_int_w = new_cases_per_million %>% prop.table %>%
            multiply_by(people_fully_vaccinated_per_hundred %>% divide_by(100)) %>% sum(na.rm = T),
        ex_d = excess_mortality_cumulative_per_million %>% diff_first_last
    )

eu21 %>%
    ggplot(aes(vac_int, ex_d))+
    geom_text_repel(aes(label = iso_code), size = 3, color = "#bababa")+
    geom_point()+
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme_minimal()+
    theme(plot.background = element_rect(color = NA, fill = "#fafafa"))+
    labs(
        x = "Integral vaccination measure [0, 1]",
        y = "Excess deaths per million in 2021"
    )

ggsave("out/90-int-vac.png", width = 5, height = 4, type = "cairo-png")

eu21 %>%
    ggplot(aes(vac_int_w, ex_d))+
    geom_point()+
    geom_text_repel(aes(label = iso_code), size = 3, color = "#bababa")+
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme_minimal()+
    theme(plot.background = element_rect(color = NA, fill = "#fafafa"))+
    labs(
        x = "Integral vaccination measure, weighted by new cases [0, 1]",
        y = "Excess deaths per million in 2021"
    )

ggsave("out/90-int-vac-w.png", width = 5, height = 4, type = "cairo-png")

