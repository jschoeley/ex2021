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
  lab <- paste0('(', low, ';', high, ')')
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
      e0_diff_2_year = ex_diff_2_year_mean,
      e0_diff_2_year_lo = ex_diff_2_year_q0.025,
      e0_diff_2_year_hi = ex_diff_2_year_q0.975,
      bbi = bbi_mean,
      bbi_lo = bbi_q0.025,
      bbi_hi = bbi_q0.975
    ) %>%
      pivot_wider(id = region_iso, names_from = year, values_from = c(
        e0_diff, e0_diff_lo, e0_diff_hi,
        e0_diff_2_year, e0_diff_2_year_lo, e0_diff_2_year_hi,
        bbi, bbi_lo, bbi_hi
      ))
  ) %>%
  mutate(
    e0_diff_since_2019 = e0_diff_2_year_2021,
    e0_diff_since_2019_lo = e0_diff_2_year_lo_2021,
    e0_diff_since_2019_hi = e0_diff_2_year_hi_2021
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

tab$arriaga$data$hypotheses <-
  dat$table_input %>%
  transmute(
    region_name_short = region_name_short,
    # h1. LE change 2019 through 2021
    h1_effect =
      e0_diff_since_2019,
    h1_sign =
      sign(h1_effect),
    h1_lo =
      e0_diff_since_2019_lo,
    h1_hi =
      e0_diff_since_2019_hi,
    h1_significance =
      ifelse(sign(h1_lo) == sign(h1_hi), TRUE, FALSE),
    h1_attribution_1 =
      ifelse(
        sign(`e0_diff_cntrb_since_2019_<60`) == sign(`e0_diff_cntrb_since_2019_60+`),
        'primarily', 'solely'
      ),
    h1_attribution_2 =
      case_when(
        h1_sign > 0 &
          (`e0_diff_cntrb_since_2019_60+` > `e0_diff_cntrb_since_2019_<60`) ~ 'old',
        h1_sign > 0 &
          (`e0_diff_cntrb_since_2019_60+` < `e0_diff_cntrb_since_2019_<60`) ~ 'young',
        h1_sign < 0 &
          (`e0_diff_cntrb_since_2019_60+` < `e0_diff_cntrb_since_2019_<60`) ~ 'old',
        h1_sign < 0 &
          (`e0_diff_cntrb_since_2019_60+` > `e0_diff_cntrb_since_2019_<60`) ~ 'young'
      ),
    # h2. LE change 2019 to 2020
    h2_effect =
      e0_diff_2020,
    h2_sign =
      sign(h2_effect),
    h2_lo =
      e0_diff_lo_2020,
    h2_hi =
      e0_diff_hi_2020,
    h2_significance =
      ifelse(sign(h2_lo) == sign(h2_hi), TRUE, FALSE),
    h2_attribution_1 =
      ifelse(
        sign(`e0_diff_cntrb_2020_<60`) == sign(`e0_diff_cntrb_2020_60+`),
        'primarily', 'solely'
      ),
    h2_attribution_2 =
      case_when(
        h2_sign > 0 &
          (`e0_diff_cntrb_2020_60+` > `e0_diff_cntrb_2020_<60`) ~ 'old',
        h2_sign > 0 &
          (`e0_diff_cntrb_2020_60+` < `e0_diff_cntrb_2020_<60`) ~ 'young',
        h2_sign < 0 &
          (`e0_diff_cntrb_2020_60+` < `e0_diff_cntrb_2020_<60`) ~ 'old',
        h2_sign < 0 &
          (`e0_diff_cntrb_2020_60+` > `e0_diff_cntrb_2020_<60`) ~ 'young'
      ),
    # h3. LE change 2020 to 2021
    h3_effect =
      e0_diff_2021,
    h3_sign =
      sign(h3_effect),
    h3_lo =
      e0_diff_lo_2021,
    h3_hi =
      e0_diff_hi_2021,
    h3_significance =
      ifelse(sign(h3_lo) == sign(h3_hi), TRUE, FALSE),
    h3_attribution_1 =
      ifelse(
        sign(`e0_diff_cntrb_2021_<60`) == sign(`e0_diff_cntrb_2021_60+`),
        'primarily', 'solely'
      ),
    h3_attribution_2 =
      case_when(
        h3_sign > 0 &
          (`e0_diff_cntrb_2021_60+` > `e0_diff_cntrb_2021_<60`) ~ 'old',
        h3_sign > 0 &
          (`e0_diff_cntrb_2021_60+` < `e0_diff_cntrb_2021_<60`) ~ 'young',
        h3_sign < 0 &
          (`e0_diff_cntrb_2021_60+` < `e0_diff_cntrb_2021_<60`) ~ 'old',
        h3_sign < 0 &
          (`e0_diff_cntrb_2021_60+` > `e0_diff_cntrb_2021_<60`) ~ 'young'
      )
  )

tab$arriaga$data$table <-
  tab$arriaga$data$hypotheses %>%
  transmute(
    region_name_short = region_name_short,
    # h1
    h1_glyph = case_when(
      h1_attribution_1 == 'primarily' & h1_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h1_attribution_1 == 'primarily' & h1_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h1_attribution_1 == 'solely' & h1_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h1_attribution_1 == 'solely' & h1_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young
    ),
    h1_glyph = case_when(
      h1_sign < 0 & h1_significance ~ paste0('\\color{negativesig}', h1_glyph),
      h1_sign < 0 & !h1_significance ~ paste0('\\color{negativenonsig}', h1_glyph),
      h1_sign > 0 & h1_significance ~ paste0('\\color{positivesig}', h1_glyph),
      h1_sign > 0 & !h1_significance ~ paste0('\\color{positivenonsig}', h1_glyph)
    ),
    h1_estimate = FormatTable(h1_effect, 12),
    h1_estimate = case_when(
      h1_sign < 0 & h1_significance ~ paste0('\\color{negativesig}', h1_estimate),
      h1_sign < 0 & !h1_significance ~ paste0('\\color{negativenonsig}', h1_estimate),
      h1_sign > 0 & h1_significance ~ paste0('\\color{positivesig}', h1_estimate),
      h1_sign > 0 & !h1_significance ~ paste0('\\color{positivenonsig}', h1_estimate)
    ),
    h1_ci_lo = FormatTable(h1_lo, 12),
    h1_ci_lo = case_when(
      h1_sign < 0 & h1_significance ~ paste0('\\color{negativesig}', h1_ci_lo, '{to}'),
      h1_sign < 0 & !h1_significance ~ paste0('\\color{negativenonsig}', h1_ci_lo, '{to}'),
      h1_sign > 0 & h1_significance ~ paste0('\\color{positivesig}', h1_ci_lo, '{to}'),
      h1_sign > 0 & !h1_significance ~ paste0('\\color{positivenonsig}', h1_ci_lo, '{to}')
    ),
    h1_ci_hi = FormatTable(h1_hi, 12),
    h1_ci_hi = case_when(
      h1_sign < 0 & h1_significance ~ paste0('\\color{negativesig}', h1_ci_hi),
      h1_sign < 0 & !h1_significance ~ paste0('\\color{negativenonsig}', h1_ci_hi),
      h1_sign > 0 & h1_significance ~ paste0('\\color{positivesig}', h1_ci_hi),
      h1_sign > 0 & !h1_significance ~ paste0('\\color{positivenonsig}', h1_ci_hi)
    ),
    # h2
    h2_glyph = case_when(
      h2_attribution_1 == 'primarily' & h2_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h2_attribution_1 == 'primarily' & h2_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h2_attribution_1 == 'solely' & h2_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h2_attribution_1 == 'solely' & h2_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young
    ),
    h2_glyph = case_when(
      h2_sign < 0 & h2_significance ~ paste0('\\color{negativesig}', h2_glyph),
      h2_sign < 0 & !h2_significance ~ paste0('\\color{negativenonsig}', h2_glyph),
      h2_sign > 0 & h2_significance ~ paste0('\\color{positivesig}', h2_glyph),
      h2_sign > 0 & !h2_significance ~ paste0('\\color{positivenonsig}', h2_glyph)
    ),
    h2_estimate = FormatTable(h2_effect, 12),
    h2_estimate = case_when(
      h2_sign < 0 & h2_significance ~ paste0('\\color{negativesig}', h2_estimate),
      h2_sign < 0 & !h2_significance ~ paste0('\\color{negativenonsig}', h2_estimate),
      h2_sign > 0 & h2_significance ~ paste0('\\color{positivesig}', h2_estimate),
      h2_sign > 0 & !h2_significance ~ paste0('\\color{positivenonsig}', h2_estimate)
    ),
    h2_ci_lo = FormatTable(h2_lo, 12),
    h2_ci_lo = case_when(
      h2_sign < 0 & h2_significance ~ paste0('\\color{negativesig}', h2_ci_lo, '{to}'),
      h2_sign < 0 & !h2_significance ~ paste0('\\color{negativenonsig}', h2_ci_lo, '{to}'),
      h2_sign > 0 & h2_significance ~ paste0('\\color{positivesig}', h2_ci_lo, '{to}'),
      h2_sign > 0 & !h2_significance ~ paste0('\\color{positivenonsig}', h2_ci_lo,'{to}')
    ),
    h2_ci_hi = FormatTable(h1_hi, 12),
    h2_ci_hi = case_when(
      h2_sign < 0 & h2_significance ~ paste0('\\color{negativesig}', h2_ci_hi),
      h2_sign < 0 & !h2_significance ~ paste0('\\color{negativenonsig}', h2_ci_hi),
      h2_sign > 0 & h2_significance ~ paste0('\\color{positivesig}', h2_ci_hi),
      h2_sign > 0 & !h2_significance ~ paste0('\\color{positivenonsig}', h2_ci_hi)
    ),
    # h3
    h3_glyph = case_when(
      h3_attribution_1 == 'primarily' & h3_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h3_attribution_1 == 'primarily' & h3_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h3_attribution_1 == 'solely' & h3_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h3_attribution_1 == 'solely' & h3_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young
    ),
    h3_glyph = case_when(
      h3_sign < 0 & h3_significance ~ paste0('\\color{negativesig}', h3_glyph),
      h3_sign < 0 & !h3_significance ~ paste0('\\color{negativenonsig}', h3_glyph),
      h3_sign > 0 & h3_significance ~ paste0('\\color{positivesig}', h3_glyph),
      h3_sign > 0 & !h3_significance ~ paste0('\\color{positivenonsig}', h3_glyph)
    ),
    h3_estimate = FormatTable(h3_effect, 12),
    h3_estimate = case_when(
      h3_sign < 0 & h3_significance ~ paste0('\\color{negativesig}', h3_estimate),
      h3_sign < 0 & !h3_significance ~ paste0('\\color{negativenonsig}', h3_estimate),
      h3_sign > 0 & h3_significance ~ paste0('\\color{positivesig}', h3_estimate),
      h3_sign > 0 & !h3_significance ~ paste0('\\color{positivenonsig}', h3_estimate)
    ),
    h3_ci_lo = FormatTable(h1_lo, 12),
    h3_ci_lo = case_when(
      h3_sign < 0 & h3_significance ~ paste0('\\color{negativesig}', h3_ci_lo, '{to}'),
      h3_sign < 0 & !h3_significance ~ paste0('\\color{negativenonsig}', h3_ci_lo, '{to}'),
      h3_sign > 0 & h3_significance ~ paste0('\\color{positivesig}', h3_ci_lo, '{to}'),
      h3_sign > 0 & !h3_significance ~ paste0('\\color{positivenonsig}', h3_ci_lo, '{to}')
    ),
    h3_ci_hi = FormatTable(h1_hi, 12),
    h3_ci_hi = case_when(
      h3_sign < 0 & h3_significance ~ paste0('\\color{negativesig}', h3_ci_hi),
      h3_sign < 0 & !h3_significance ~ paste0('\\color{negativenonsig}', h3_ci_hi),
      h3_sign > 0 & h3_significance ~ paste0('\\color{positivesig}', h3_ci_hi),
      h3_sign > 0 & !h3_significance ~ paste0('\\color{positivenonsig}', h3_ci_hi)
    ),
  ) %>%
  select(
    region_name_short,
    h1_glyph, h1_estimate, h1_ci_lo, h1_ci_hi,
    h2_glyph, h2_estimate, h2_ci_lo, h2_ci_hi,
    h3_glyph, h3_estimate, h3_ci_lo, h3_ci_hi
  )

tab$arriaga$table <-
  tab$arriaga$data$table %>%
  gt(rowname_col = 'region_name_short') %>%
  cols_label(
    h1_glyph = 'AT', h1_estimate = 'ES', h1_ci_lo = 'CI', h1_ci_hi = '',
    h2_glyph = 'AT', h2_estimate = 'ES', h2_ci_lo = 'CI', h2_ci_hi = '',
    h3_glyph = 'AT', h3_estimate = 'ES', h3_ci_lo = 'CI', h3_ci_hi = ''
  ) %>%
  tab_spanner(columns = 2:5, 'Net LE losses 2019 to 2021', id = 'h1') %>%
  tab_spanner(columns = 6:9, 'LE losses 2020', id = 'h2') %>%
  tab_spanner(columns = 10:13, 'LE losses 2021', id = 'h3') %>%
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
    locations = cells_column_labels('h1_ci_lo'),
    footnote = '95% confidence interval'
  )

tab$arriaga$table
as_latex(tab$arriaga$table) %>% as.character() %>%
  gsub( x = ., pattern = '\\$', replacement = '$', fixed = TRUE) %>%
  gsub( x = ., pattern = '\\textbackslash ', replacement = '\\', fixed = TRUE) %>%
  gsub( x = ., pattern = '\\textbackslash{}', replacement = '\\', fixed = TRUE) %>%
  gsub( x = ., pattern = '\\{', replacement = '{', fixed = TRUE) %>%
  gsub( x = ., pattern = '\\}', replacement = '}', fixed = TRUE) %>%
  clipr::write_clip()

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
