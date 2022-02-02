# Generate table of qualitative Arriaga decomposition results

# Init ------------------------------------------------------------

library(yaml); library(tidyverse); library(gt)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  region_metadata = './cfg/region_metadata.csv',
  figspecs = './cfg/figure_specification.R',
  lifetables = './out/lifetables.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  tab_arriaga = './out/tab_arriaga.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# meta data on regions
region_meta <- read_csv(paths$input$region_meta, na = '.')

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  regions_for_analysis = config$regions_for_all_cause_analysis
})

# list containers for analysis artifacts
dat <- list()
tab <- list()
fig <- list()

# Functions -------------------------------------------------------

# figure specifications
source(paths$input$figspec)

FormatTable <- function (x, scaler = 1, digits = 1) {
  lab <- formatC(x*scaler, digits = digits, format = 'f', flag = '+')
  ifelse(lab == 'NA', '.', lab)
}

FormatCI <- function (lo, hi, scaler = 1, digits = 1) {
  low <- formatC(lo*scaler, digits = digits, format = 'f')
  high <- formatC(hi*scaler, digits = digits, format = 'f')
  lab <- paste0('(', low, ',', high, ')')
  ifelse(lo == 'NA', '.', lab)
}

# Load ex differences ---------------------------------------------

dat$lifetables <- readRDS(paths$input$lifetables)

# Arriaga table ---------------------------------------------------

dat$subset <-
  dat$lifetables %>%
  filter(sex == 'T', region_iso %in% cnst$regions_for_analysis,
         projected == 'actual',
         year >= 2020)

dat$table_input <-
  dat$subset %>%
  mutate(
    age_group = cut(age, c(0, 60, Inf), right = FALSE,
                    labels = c('<60', '60+'))
  ) %>%
  select(
    region_iso, age, age_group, year, e0_cntrb_t_mean
  ) %>%
  group_by(region_iso, year, age_group) %>%
  summarise(e0_cntrb_t_mean = sum(e0_cntrb_t_mean)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = c(year),
    values_from = c(e0_cntrb_t_mean),
    names_prefix = 'e0_diff_cntrb_'
  ) %>%
  mutate(
    e0_diff_cntrb_since_2019 =
      e0_diff_cntrb_2020 + e0_diff_cntrb_2021
  ) %>%
  pivot_wider(names_from = age_group, values_from = c(
    e0_diff_cntrb_2020, e0_diff_cntrb_2021, e0_diff_cntrb_since_2019
  )) %>%
  left_join(
    dat$subset %>% filter(age == 0) %>% select(
      region_iso, year,
      e0_diff = ex_diff_mean,
      e0_diff_lo = ex_diff_q0.025,
      e0_diff_hi = ex_diff_q0.975,
      bbi = bbi_mean,
      bbi_lo = bbi_q0.025,
      bbi_hi = bbi_q0.975
    ) %>%
      pivot_wider(id = region_iso, names_from = year, values_from = c(
        e0_diff, e0_diff_lo, e0_diff_hi, bbi, bbi_lo, bbi_hi
      ))
  ) %>%
  mutate(
    e0_diff_since_2019 = e0_diff_2020 + e0_diff_2021
  ) %>%
  left_join(
    region_meta, by = c('region_iso' = 'region_code_iso3166_2')
  )

tab$arriaga <- list()

tab$arriaga$cnst <-
  list(
    primary_cutoff = 0.6,
    glyph_overall_old = '$\\blacktriangleright$',
    glyph_overall_young = '$\\vartriangleleft$',
    glyph_solely_old = '$\\blacktriangleright\\blacktriangleright$',
    glyph_solely_young = '$\\vartriangleleft\\vartriangleleft$'
    #glyph_overall_old = '\u25B7',
    #glyph_overall_young = '\u25C1',
    #glyph_solely_old = '\u25B6',
    #glyph_solely_young = '\u25C0'
  )

tab$arriaga$data <- list()

tab$arriaga$data$hypothesis <-
  dat$table_input %>%
  mutate(
    # 1 countries which lost e0 2019 through 2021
    h1 = e0_diff_since_2019 < 0,
    # 1.1 countries which lost e0 2019 through 2021
    # due to overall elevated mortality
    h1.1 = e0_diff_since_2019 < 0 &
      `e0_diff_cntrb_since_2019_<60` < 0 &
      `e0_diff_cntrb_since_2019_60+` < 0,
    # 1.1a countries for which e0 losses 2019 through 2021 were
    # primarily due elevated mortality in age groups 60+
    h1.1a = e0_diff_since_2019 < 0 &
      `e0_diff_cntrb_since_2019_<60` < 0 &
      `e0_diff_cntrb_since_2019_60+` < 0 &
      `e0_diff_cntrb_since_2019_<60` > `e0_diff_cntrb_since_2019_60+`,
    # 1.1b countries for which e0 losses 2019 through 2021
    # were primarily due to elevated mortality in age groups <60
    h1.1b = e0_diff_since_2019 < 0 &
      `e0_diff_cntrb_since_2019_<60` < 0 &
      `e0_diff_cntrb_since_2019_60+` < 0 &
      `e0_diff_cntrb_since_2019_<60` < `e0_diff_cntrb_since_2019_60+`,
    # 1.2 countries which lost e0 2019 through 2021 solely
    # due to elevated mortality ages 60+
    h1.2 = e0_diff_since_2019 < 0 &
      `e0_diff_cntrb_since_2019_<60` >= 0 &
      `e0_diff_cntrb_since_2019_60+` < 0,    
    # 1.3 countries which lost e0 2019 through 2021 solely
    # due to elevated mortality ages <60
    h1.3 = e0_diff_since_2019 < 0 &
      `e0_diff_cntrb_since_2019_<60` < 0 &
      `e0_diff_cntrb_since_2019_60+` >= 0,
    
    # 2 countries which lost e0 2020
    h2 = e0_diff_2020 < 0,
    # 2.1 countries which lost e0 2020
    # due to overall elevated mortality
    h2.1 = e0_diff_2020 < 0 &
      `e0_diff_cntrb_2020_60+` < 0 &
      `e0_diff_cntrb_2020_<60` < 0,
    # 2.1a countries which lost e0 2020
    # primarily due to elevated mortality ages 60+
    h2.1a = e0_diff_2020 < 0 &
      `e0_diff_cntrb_2020_<60` < 0 &
      `e0_diff_cntrb_2020_60+` < 0 &
      `e0_diff_cntrb_2020_<60` > `e0_diff_cntrb_2020_60+`,
    # 2.1b countries which lost e0 2020
    # primarily due to elevated mortality ages <60
    h2.1b = e0_diff_2020 < 0 &
      `e0_diff_cntrb_2020_<60` < 0 &
      `e0_diff_cntrb_2020_60+` < 0 &
      `e0_diff_cntrb_2020_<60` < `e0_diff_cntrb_2020_60+`,
    # 2.2 countries which lost e0 2020
    # solely due to elevated mortality ages 60+
    h2.2 = e0_diff_2020 < 0 &
      `e0_diff_cntrb_2020_<60` >= 0 &
      `e0_diff_cntrb_2020_60+` < 0,
    # 2.3 countries which lost e0 2019 through 2020
    # solely due to elevated mortality ages <60
    h2.3 = e0_diff_2020 < 0 &
      `e0_diff_cntrb_2020_<60` < 0 &
      `e0_diff_cntrb_2020_60+` >= 0,
    
    # 3 countries which lost e0 2021
    h3 = e0_diff_2021 < 0,
    # 3.1 countries which lost e0 2021
    # due to overall elevated mortality
    h3.1 = e0_diff_2021 < 0 &
      `e0_diff_cntrb_2021_60+` < 0 &
      `e0_diff_cntrb_2021_<60` < 0,
    # 3.1a countries which lost e0 2021
    # primarily due to elevated mortality ages 60+
    h3.1a = e0_diff_2021 < 0 &
      `e0_diff_cntrb_2021_<60` < 0 &
      `e0_diff_cntrb_2021_60+` < 0 &
      `e0_diff_cntrb_2021_<60` > `e0_diff_cntrb_2021_60+`,
    # 3.1b countries which lost e0 2021
    # primarily due to elevated mortality ages <60
    h3.1b = e0_diff_2021 < 0 &
      `e0_diff_cntrb_2021_<60` < 0 &
      `e0_diff_cntrb_2021_60+` < 0 &
      `e0_diff_cntrb_2021_<60` < `e0_diff_cntrb_2021_60+`,
    # 3.2 countries which lost e0 2021
    # solely due to elevated mortality ages 60+
    h3.2 = e0_diff_2021 < 0 &
      `e0_diff_cntrb_2021_<60` >= 0 &
      `e0_diff_cntrb_2021_60+` < 0,
    # 3.3 countries which lost e0 2021
    # solely due to elevated mortality ages <60
    h3.3 = e0_diff_2021 < 0 &
      `e0_diff_cntrb_2021_<60` < 0 &
      `e0_diff_cntrb_2021_60+` >= 0,
    
    # 4 countries with compound e0 losses since 2019
    h4 = e0_diff_2020 < 0 & e0_diff_2021 < 0,
    # 4.1 countries with compound e0 losses since 2019
    # due to overall elevated mortality
    h4.1 = e0_diff_2020 < 0 & e0_diff_2021 < 0 &
      `e0_diff_cntrb_since_2019_60+` < 0 &
      `e0_diff_cntrb_since_2019_<60` < 0,
    # 4.1a countries with compound e0 losses since 2019
    # primarily due to elevated mortality ages 60+
    h4.1a = e0_diff_2020 < 0 & e0_diff_2021 < 0 &
      `e0_diff_cntrb_since_2019_60+` < 0 &
      `e0_diff_cntrb_since_2019_<60` < 0 &
      `e0_diff_cntrb_since_2019_<60` > `e0_diff_cntrb_since_2019_60+`,
    # 4.1b countries with compound e0 losses since 2019
    # primarily due to elevated mortality ages <60
    h4.1b = e0_diff_2020 < 0 & e0_diff_2021 < 0 &
      `e0_diff_cntrb_since_2019_60+` < 0 &
      `e0_diff_cntrb_since_2019_<60` < 0 &
      `e0_diff_cntrb_since_2019_<60` < `e0_diff_cntrb_since_2019_60+`,
    # 4.2 countries with compound e0 losses since 2019
    # solely due to elevated mortality ages 60+
    h4.2 = e0_diff_2020 < 0 & e0_diff_2021 < 0 &
      `e0_diff_cntrb_since_2019_<60` >= 0 &
      `e0_diff_cntrb_since_2019_60+` < 0,
    # 4.3 countries with compound e0 losses since 2019
    # solely due to elevated mortality ages <60
    h4.3 = e0_diff_2020 < 0 & e0_diff_2021 < 0 &
      `e0_diff_cntrb_since_2019_<60` < 0 &
      `e0_diff_cntrb_since_2019_60+` >= 0,
    
    # 5 countries which bounced-back in 2021
    h5 = e0_diff_2020 < 0 & e0_diff_2021 > 0,
    # 5.1 countries which bounced-back in 2021
    # due to overall normalizing mortality
    h5.1 = e0_diff_2020 < 0 & e0_diff_2021 > 0 &
      `e0_diff_cntrb_2021_60+` > 0 &
      `e0_diff_cntrb_2021_<60` > 0,
    # 5.1a countries which bounced-back in 2021
    # primarily due to normalizing mortality ages 60+
    h5.1a = e0_diff_2020 < 0 & e0_diff_2021 > 0 &
      `e0_diff_cntrb_2021_60+` > 0 &
      `e0_diff_cntrb_2021_<60` > 0 &
      `e0_diff_cntrb_2021_<60` < `e0_diff_cntrb_2021_60+`,
    # 5.1b countries which bounced-back in 2021
    # primarily due to normalizing mortality ages <60
    h5.1b = e0_diff_2020 < 0 & e0_diff_2021 > 0 &
      `e0_diff_cntrb_2021_60+` > 0 &
      `e0_diff_cntrb_2021_<60` > 0 &
      `e0_diff_cntrb_2021_<60` > `e0_diff_cntrb_2021_60+`,
    # 5.2 countries which bounced-back in 2021
    # solely due to normalizing mortality ages 60+
    h5.2 = e0_diff_2020 < 0 & e0_diff_2021 > 0 &
      `e0_diff_cntrb_since_2019_<60` <= 0 &
      `e0_diff_cntrb_2021_60+` > 0,
    # 5.3 countries which bounced-back in 2021
    # solely due to normalizing mortality ages <60
    h5.3 = e0_diff_2020 < 0 & e0_diff_2021 > 0 &
      `e0_diff_cntrb_2021_<60` > 0 &
      `e0_diff_cntrb_2021_60+` <= 0,
  )

tab$arriaga$data$table <-
  tab$arriaga$data$hypothesis %>%
  mutate(
    h1_contribution = case_when(
      h1 & h1.1 ~ 'overall',
      h1 & (h1.2 | h1.3) ~ 'solely',
      TRUE ~ 'none'
    ),
    h1_age = case_when(
      h1_contribution == 'overall' & h1.1a ~ 'old',
      h1_contribution == 'overall' & h1.1b ~ 'young',
      h1_contribution == 'solely' & h1.2 ~ 'old',
      h1_contribution == 'solely' & h1.3 ~ 'young',
      TRUE ~ 'none'
    ),
    h2_contribution = case_when(
      h2 & h2.1 ~ 'overall',
      h2 & (h2.2 | h2.3) ~ 'solely',
      TRUE ~ 'none'
    ),
    h2_age = case_when(
      h2_contribution == 'overall' & h2.1a ~ 'old',
      h2_contribution == 'overall' & h2.1b ~ 'young',
      h2_contribution == 'solely' & h2.2 ~ 'old',
      h2_contribution == 'solely' & h2.3 ~ 'young',
      TRUE ~ 'none'
    ),
    h3_contribution = case_when(
      h3 & h3.1 ~ 'overall',
      h3 & (h3.2 | h3.3) ~ 'solely',
      TRUE ~ 'none'
    ),
    h3_age = case_when(
      h3_contribution == 'overall' & h3.1a ~ 'old',
      h3_contribution == 'overall' & h3.1b ~ 'young',
      h3_contribution == 'solely' & h3.2 ~ 'old',
      h3_contribution == 'solely' & h3.3 ~ 'young',
      TRUE ~ 'none'
    ),
    h4_contribution = case_when(
      h4 & h4.1 ~ 'overall',
      h4 & (h4.2 | h4.3) ~ 'solely',
      TRUE ~ 'none'
    ),
    h4_age = case_when(
      h4_contribution == 'overall' & h4.1a ~ 'old',
      h4_contribution == 'overall' & h4.1b ~ 'young',
      h4_contribution == 'solely' & h4.2 ~ 'old',
      h4_contribution == 'solely' & h4.3 ~ 'young',
      TRUE ~ 'none'
    ),
    h5_contribution = case_when(
      h5 & h5.1 ~ 'overall',
      h5 & (h5.2 | h5.3) ~ 'solely',
      TRUE ~ 'none'
    ),
    h5_age = case_when(
      h5_contribution == 'overall' & h5.1a ~ 'old',
      h5_contribution == 'overall' & h5.1b ~ 'young',
      h5_contribution == 'solely' & h5.2 ~ 'old',
      h5_contribution == 'solely' & h5.3 ~ 'young',
      TRUE ~ 'none'
    )
  ) %>%
  mutate(
    h1_glyph = case_when(
      h1_contribution == 'none' ~ '',
      h1_contribution == 'overall' & h1_age == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h1_contribution == 'overall' & h1_age == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h1_contribution == 'solely' & h1_age == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h1_contribution == 'solely' & h1_age == 'young' ~ tab$arriaga$cnst$glyph_overall_young
    ),
    h1_estimate = FormatTable(ifelse(h1, e0_diff_since_2019, NA), 12),
    h1_ci = '.',
    h2_glyph = case_when(
      h2_contribution == 'none' ~ '',
      h2_contribution == 'overall' & h2_age == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h2_contribution == 'overall' & h2_age == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h2_contribution == 'solely' & h2_age == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h2_contribution == 'solely' & h2_age == 'young' ~ tab$arriaga$cnst$glyph_solely_young
    ),
    h2_estimate = FormatTable(ifelse(h2, e0_diff_2020, NA), 12),
    h2_ci = ifelse(h2, FormatCI(e0_diff_lo_2020, e0_diff_hi_2020, 12), '.'),
    h3_glyph = case_when(
      h3_contribution == 'none' ~ '',
      h3_contribution == 'overall' & h3_age == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h3_contribution == 'overall' & h3_age == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h3_contribution == 'solely' & h3_age == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h3_contribution == 'solely' & h3_age == 'young' ~ tab$arriaga$cnst$glyph_solely_young
    ),
    h3_estimate = FormatTable(ifelse(h3, e0_diff_2021, NA), 12),
    h3_ci = FormatCI(e0_diff_lo_2021, e0_diff_hi_2021, 12),
    h5_glyph = case_when(
      h5_contribution == 'none' ~ '',
      h5_contribution == 'overall' & h5_age == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h5_contribution == 'overall' & h5_age == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h5_contribution == 'solely' & h5_age == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h5_contribution == 'solely' & h5_age == 'young' ~ tab$arriaga$cnst$glyph_solely_young
    ),
    h5_estimate = FormatTable(ifelse(h5, bbi_2021, NA), scaler = 100, digits = 0),
    h5_ci = ifelse(h5, FormatCI(bbi_lo_2021, bbi_hi_2021, scaler = 100, digits = 0), '.')
  ) %>%
  select(
    region_name_short,
    h1_glyph, h1_estimate, h1_ci,
    h2_glyph, h2_estimate, h2_ci,
    h3_glyph, h3_estimate, h3_ci,
    h5_glyph, h5_estimate, h5_ci
  )

tab$arriaga$table <-
  tab$arriaga$data$table %>%
  gt(rowname_col = 'region_name_short') %>%
  cols_label(
    h1_glyph = 'AT', h1_estimate = 'ES', h1_ci = 'CI',
    h2_glyph = 'AT', h2_estimate = 'ES', h2_ci = 'CI',
    h3_glyph = 'AT', h3_estimate = 'ES', h3_ci = 'CI',
    h5_glyph = 'AT', h5_estimate = 'ES', h5_ci = 'CI'
  ) %>%
  tab_spanner(columns = 2:4, 'Net LE losses 2019 to 2021', id = 'h1') %>%
  tab_spanner(columns = 5:7, 'LE losses 2020', id = 'h2') %>%
  tab_spanner(columns = 8:10, 'LE losses 2021', id = 'h3') %>%
  tab_spanner(columns = 11:13, 'LE bounce-back 2021', id = 'h5') %>%
  tab_footnote(
    locations = cells_column_labels('h1_glyph'),
    footnote = paste0(
      'Attribution, due to mortality changes among',
      ' primarily 60+', tab$arriaga$cnst$glyph_overall_old,
      ', solely 60+', tab$arriaga$cnst$glyph_solely_old,
      ', primarily <60', tab$arriaga$cnst$glyph_overall_young,
      ', solely <60', tab$arriaga$cnst$glyph_solely_young
    )
  ) %>%
  tab_footnote(
    locations = cells_column_labels('h1_estimate'),
    footnote = 'Central estimate in months'
  ) %>%
  tab_footnote(
    locations = cells_column_labels('h1_ci'),
    footnote = '95% confidence interval'
  )

tab$arriaga$table
as_latex(tab$arriaga$table) %>% as.character() %>% clipr::write_clip()

# tab$arriaga <-
#   tab$arriaga %>%
#   select(region_name = region_name_short, h1:h4.3) %>%
#   pivot_longer(cols = c(h1:h4.3)) %>%
#   mutate(value = ifelse(value, '\u2B24', '')) %>%
#   pivot_wider(
#     id_cols = c(region_name, name),
#     names_from = region_name,
#     values_from = value
#   ) %>%
#   mutate(name = str_sub(name, 2)) %>%
#   mutate(
#     description = case_when(
#       name == '1' ~ 'Net LE loss 2019 through 2021',
#       name == '1.1' ~ 'due to overall elevated mortality',
#       name == '1.1a' ~ 'primarily due to elevated mortality ages 60+',
#       name == '1.1b' ~ 'primarily due to elevated mortality ages <60',
#       name == '1.2' ~ 'solely due to elevated mortality ages 60+',
#       name == '1.3' ~ 'solely due to elevated mortality ages <60',
#       name == '2' ~ '2020 LE loss',
#       name == '2.1' ~ 'due to overall elevated mortality',
#       name == '2.1a' ~ 'primarily due to elevated mortality ages 60+',
#       name == '2.1b' ~ 'primarily due to elevated mortality ages <60',
#       name == '2.2' ~ 'solely due to elevated mortality ages 60+',
#       name == '2.3' ~ 'solely due to elevated mortality ages <60',
#       name == '3' ~ 'Compound 2021 LE loss',
#       name == '3.1' ~ 'due to overall elevated mortality',
#       name == '3.1a' ~ 'primarily due to elevated mortality ages 60+',
#       name == '3.1b' ~ 'primarily due to elevated mortality ages <60',
#       name == '3.2' ~ 'solely due to elevated mortality ages 60+',
#       name == '3.3' ~ 'solely due to elevated mortality ages <60',
#       name == '4' ~ '2021 LE bounce-back',
#       name == '4.1' ~ 'due to overall normalizing mortality',
#       name == '4.1a' ~ 'primarily due to normalizing mortality ages 60+',
#       name == '4.1b' ~ 'primarily due to normalizing mortality ages <60',
#       name == '4.2' ~ 'solely due to normalizing mortality ages 60+',
#       name == '4.3' ~ 'solely due to normalizing mortality ages <60',
#       TRUE ~ ''
#     )
#   ) %>%
#   relocate(description, .after = name) %>%
#   gt(rowname_col = 'name') %>% cols_label(description = '') %>%
#   cols_align('center', -(1:2)) %>%
#   tab_style(
#     cell_text(size = 'medium', weight = 'bold', align = 'left'),
#     cells_body(rows = str_length(name) == 1)
#   ) %>%
#   tab_style(
#     cell_text(size = 'small', align = 'center'),
#     cells_body(rows = str_length(name) == 3)
#   ) %>%
#   tab_style(
#     cell_text(size = 'x-small', style = 'italic', align = 'right'),
#     cells_body(rows = str_length(name) == 4)
#   ) %>%
#   tab_style(
#     cell_borders(sides = c('left', 'right'), color = 'grey90'),
#     cells_body(columns = -1)
#   ) %>%
#   tab_source_note(
#     source_note =
#       paste0('(', region_meta$region_name_short, ') ', region_meta$region_name,
#              collapse = ', ')
#   )

#as_latex(tab$arriaga) %>% as.character()  %>% clipr::write_clip()

# Export ----------------------------------------------------------

saveRDS(tab$arriaga, file = paths$output$tab_arriaga)
