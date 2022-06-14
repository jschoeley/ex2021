# Generate table of qualitative Arriaga decomposition results

# Init ------------------------------------------------------------

library(yaml); library(tidyverse); library(gt); library(purrr)

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
  arriaga_cntfc = './out/40-arriaga_cntfc.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  tab_arriaga = './out/54-tab_arriaga.rds',
  csv_arriaga = './tmp/54-tab_arriaga.csv'
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
  ifelse(is.na(x), '', lab)
}

FormatCI <- function (lo, hi, scaler = 1, digits = 1) {
  low <- formatC(lo*scaler, digits = digits, format = 'f')
  high <- formatC(hi*scaler, digits = digits, format = 'f')
  lab <- paste0('(', low, ';', high, ')')
  ifelse(is.na(lo), '', lab)
}

# Load ex differences ---------------------------------------------

dat$lifetables <- readRDS(paths$input$lifetables)
dat$arriaga_cntfc <- readRDS(paths$input$arriaga_cntfc)

# Arriaga table ---------------------------------------------------

dat$subset_lt <-
  dat$lifetables %>%
  filter(region_iso %in% cnst$regions_for_analysis, year >= 2020,
         quarter == 'annual')

dat$subset_arriaga_cntfc <-
  dat$arriaga_cntfc %>%
  filter(region_iso %in% cnst$regions_for_analysis, year >= 2020,
         quarter == 'annual')

dat$table_contributions <-
  dat$subset_lt %>%
  mutate(
    age_group = cut(age, c(0, 60, Inf), right = FALSE,
                    labels = c('<60', '60+'))
  ) %>%
  select(
    projected, region_iso, sex, age, age_group, year, e0_cntrb_t_mean
  ) %>%
  group_by(projected, region_iso, sex, year, age_group) %>%
  summarise(e0_cntrb_t_mean = sum(e0_cntrb_t_mean)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = c(year, projected),
    values_from = c(e0_cntrb_t_mean),
    names_prefix = 'e0_diff_cntrb_'
  ) %>%
  mutate(
    e0_diff_cntrb_since_2019_actual =
      e0_diff_cntrb_2020_actual + e0_diff_cntrb_2021_actual,
    e0_diff_cntrb_since_2019_projected =
      e0_diff_cntrb_2020_projected + e0_diff_cntrb_2021_projected
  ) %>%
  pivot_wider(names_from = age_group, values_from = c(
    e0_diff_cntrb_2020_actual, e0_diff_cntrb_2021_actual, e0_diff_cntrb_since_2019_actual,
    e0_diff_cntrb_2020_projected, e0_diff_cntrb_2021_projected, e0_diff_cntrb_since_2019_projected
  ))

dat$table_le_differences <-
  dat$subset_lt %>%
  filter(age == 0) %>%
  select(
    projected, region_iso, sex, year,
    e0_diff = ex_diff_mean,
    e0_diff_lo = ex_diff_q0.025,
    e0_diff_hi = ex_diff_q0.975,
    e0_diff_2_year = ex_diff_2_year_mean,
    e0_diff_2_year_lo = ex_diff_2_year_q0.025,
    e0_diff_2_year_hi = ex_diff_2_year_q0.975,
    bbi = bbi_mean,
    bbi_lo = bbi_q0.025,
    bbi_hi = bbi_q0.975,
    ex_deficit = ex_deficit_mean,
    ex_deficit_lo = ex_deficit_q0.025,
    ex_deficit_hi = ex_deficit_q0.975
  ) %>%
  pivot_wider(
    id_cols = c(region_iso, sex),
    names_from = c(projected, year), values_from = c(
      e0_diff, e0_diff_lo, e0_diff_hi,
      e0_diff_2_year, e0_diff_2_year_lo, e0_diff_2_year_hi,
      bbi, bbi_lo, bbi_hi, ex_deficit, ex_deficit_lo, ex_deficit_hi
    ))

dat$table_arriaga_cntfc_contributions <-
  dat$subset_arriaga_cntfc %>%
  mutate(
    age_group = cut(as.integer(age), c(0, 60, Inf), right = FALSE,
                    labels = c('<60', '60+'))
  ) %>%
  select(
    region_iso, sex, age, age_group, year, e0_cntrb_t_mean
  ) %>%
  group_by(region_iso, sex, year, age_group) %>%
  summarise(e0_deficit_cntrb = sum(e0_cntrb_t_mean)) %>%
  ungroup() %>%
  pivot_wider(names_from = c(year, age_group), values_from = c(
    e0_deficit_cntrb
  ), names_prefix = 'e0_deficit_cntrb_')

dat$table_input <-
  dat$table_contributions %>%
  left_join(dat$table_le_differences) %>%
  left_join(dat$table_arriaga_cntfc_contributions) %>%
  mutate(
    e0_diff_since_actual_2019 = e0_diff_2_year_actual_2021,
    e0_diff_since_2019_actual_lo = e0_diff_2_year_lo_actual_2021,
    e0_diff_since_2019_actual_hi = e0_diff_2_year_hi_actual_2021,
    e0_diff_since_projected_2019 = e0_diff_2_year_projected_2021,
    e0_diff_since_2019_projected_lo = e0_diff_2_year_lo_projected_2021,
    e0_diff_since_2019_projected_hi = e0_diff_2_year_hi_projected_2021
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
    glyph_solely_young = '$\\vartriangleleft\\vartriangleleft$',
    color_negative_significant_glyph = '\\color{negativesig}',
    color_negative_nonsignificant_glyph = '\\color{negativenonsig}',
    color_positive_significant_glyph = '\\color{positivesig}',
    color_positive_nonsignificant_glyph = '\\color{positivenonsig}',
    color_negative_significant_ci = '',
    color_negative_nonsignificant_ci = '',
    color_positive_significant_ci = '',
    color_positive_nonsignificant_ci = '',
    ci_start = '{[}',
    ci_end = '{]}'
  )

tab$arriaga$data <- list()

tab$arriaga$data$hypotheses <-
  dat$table_input %>%
  transmute(
    region_name_short = region_name_short, sex = sex,
    # h1. LE change 2019 through 2021
    h1_effect =
      e0_diff_since_actual_2019,
    h1_sign =
      sign(h1_effect),
    h1_lo =
      e0_diff_since_2019_actual_lo,
    h1_hi =
      e0_diff_since_2019_actual_hi,
    h1_significance =
      ifelse(sign(h1_lo) == sign(h1_hi), TRUE, FALSE),
    h1_attribution_1 =
      ifelse(
        sign(`e0_diff_cntrb_since_2019_actual_<60`) == sign(`e0_diff_cntrb_since_2019_actual_60+`),
        'primarily', 'solely'
      ),
    h1_attribution_2 =
      case_when(
        h1_sign > 0 &
          (`e0_diff_cntrb_since_2019_actual_60+` > `e0_diff_cntrb_since_2019_actual_<60`) ~ 'old',
        h1_sign > 0 &
          (`e0_diff_cntrb_since_2019_actual_60+` < `e0_diff_cntrb_since_2019_actual_<60`) ~ 'young',
        h1_sign < 0 &
          (`e0_diff_cntrb_since_2019_actual_60+` < `e0_diff_cntrb_since_2019_actual_<60`) ~ 'old',
        h1_sign < 0 &
          (`e0_diff_cntrb_since_2019_actual_60+` > `e0_diff_cntrb_since_2019_actual_<60`) ~ 'young'
      ),
    # h2. LE change 2019 to 2020
    h2_effect =
      e0_diff_actual_2020,
    h2_sign =
      sign(h2_effect),
    h2_lo =
      e0_diff_lo_actual_2020,
    h2_hi =
      e0_diff_hi_actual_2020,
    h2_significance =
      ifelse(sign(h2_lo) == sign(h2_hi), TRUE, FALSE),
    h2_attribution_1 =
      ifelse(
        sign(`e0_diff_cntrb_2020_actual_<60`) == sign(`e0_diff_cntrb_2020_actual_60+`),
        'primarily', 'solely'
      ),
    h2_attribution_2 =
      case_when(
        h2_sign > 0 &
          (`e0_diff_cntrb_2020_actual_60+` > `e0_diff_cntrb_2020_actual_<60`) ~ 'old',
        h2_sign > 0 &
          (`e0_diff_cntrb_2020_actual_60+` < `e0_diff_cntrb_2020_actual_<60`) ~ 'young',
        h2_sign < 0 &
          (`e0_diff_cntrb_2020_actual_60+` < `e0_diff_cntrb_2020_actual_<60`) ~ 'old',
        h2_sign < 0 &
          (`e0_diff_cntrb_2020_actual_60+` > `e0_diff_cntrb_2020_actual_<60`) ~ 'young'
      ),
    # h3. LE change 2020 to 2021
    h3_effect =
      e0_diff_actual_2021,
    h3_sign =
      sign(h3_effect),
    h3_lo =
      e0_diff_lo_actual_2021,
    h3_hi =
      e0_diff_hi_actual_2021,
    h3_significance =
      ifelse(sign(h3_lo) == sign(h3_hi), TRUE, FALSE),
    h3_attribution_1 =
      ifelse(
        sign(`e0_diff_cntrb_2021_actual_<60`) == sign(`e0_diff_cntrb_2021_actual_60+`),
        'primarily', 'solely'
      ),
    h3_attribution_2 =
      case_when(
        h3_sign > 0 &
          (`e0_diff_cntrb_2021_actual_60+` > `e0_diff_cntrb_2021_actual_<60`) ~ 'old',
        h3_sign > 0 &
          (`e0_diff_cntrb_2021_actual_60+` < `e0_diff_cntrb_2021_actual_<60`) ~ 'young',
        h3_sign < 0 &
          (`e0_diff_cntrb_2021_actual_60+` < `e0_diff_cntrb_2021_actual_<60`) ~ 'old',
        h3_sign < 0 &
          (`e0_diff_cntrb_2021_actual_60+` > `e0_diff_cntrb_2021_actual_<60`) ~ 'young'
      ),
    # h4. LE bounce back 2021
    h4_effect =
      ifelse(bbi_actual_2021>0,bbi_actual_2021,NA),
    h4_sign =
      sign(h4_effect),
    h4_lo =
      ifelse(bbi_lo_actual_2021>0,bbi_lo_actual_2021,NA),
    h4_hi =
      ifelse(bbi_hi_actual_2021>0,bbi_hi_actual_2021,NA),
    h4_significance =
      ifelse(sign(h4_lo) == sign(h4_hi), TRUE, FALSE),
    # h5. LE deficit 2021
    h5_effect =
      -ex_deficit_actual_2021,
    h5_sign =
      sign(h5_effect),
    h5_lo =
      -ex_deficit_lo_actual_2021,
    h5_hi =
      -ex_deficit_hi_actual_2021,
    h5_significance =
      ifelse(sign(h5_lo) == sign(h5_hi), TRUE, FALSE),
    h5_attribution_1 =
      ifelse(
        sign(`e0_deficit_cntrb_2021_<60`) == sign(`e0_deficit_cntrb_2021_60+`),
        'primarily', 'solely'
      ),
    h5_attribution_2 =
      case_when(
        h5_sign > 0 &
          (`e0_deficit_cntrb_2021_60+` > `e0_deficit_cntrb_2021_<60`) ~ 'old',
        h5_sign > 0 &
          (`e0_deficit_cntrb_2021_60+` < `e0_deficit_cntrb_2021_<60`) ~ 'young',
        h5_sign < 0 &
          (`e0_deficit_cntrb_2021_60+` < `e0_deficit_cntrb_2021_<60`) ~ 'old',
        h5_sign < 0 &
          (`e0_deficit_cntrb_2021_60+` > `e0_deficit_cntrb_2021_<60`) ~ 'young'
      )
  )

tab$arriaga$data$table <-
  tab$arriaga$data$hypotheses %>%
  transmute(
    region_name_short = region_name_short, sex = sex,
    # h1
    h1_glyph = case_when(
      h1_attribution_1 == 'primarily' & h1_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h1_attribution_1 == 'primarily' & h1_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h1_attribution_1 == 'solely' & h1_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h1_attribution_1 == 'solely' & h1_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young
    ),
    h1_glyph = case_when(
      h1_sign < 0 & h1_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_glyph, h1_glyph),
      h1_sign < 0 & !h1_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_glyph, h1_glyph),
      h1_sign > 0 & h1_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_glyph, h1_glyph),
      h1_sign > 0 & !h1_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_glyph, h1_glyph)
    ),
    h1_estimate = FormatTable(h1_effect, 12),
    h1_estimate = case_when(
      h1_sign < 0 & h1_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h1_estimate),
      h1_sign < 0 & !h1_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h1_estimate),
      h1_sign > 0 & h1_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h1_estimate),
      h1_sign > 0 & !h1_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h1_estimate)
    ),
    h1_ci_lo = FormatTable(h1_lo, 12),
    h1_ci_lo = case_when(
      h1_sign < 0 & h1_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_significant_ci, h1_ci_lo, '{;}'),
      h1_sign < 0 & !h1_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_nonsignificant_ci, h1_ci_lo, '{;}'),
      h1_sign > 0 & h1_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_significant_ci, h1_ci_lo, '{;}'),
      h1_sign > 0 & !h1_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_nonsignificant_ci, h1_ci_lo, '{;}')
    ),
    h1_ci_hi = FormatTable(h1_hi, 12),
    h1_ci_hi = case_when(
      h1_sign < 0 & h1_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h1_ci_hi, tab$arriaga$cnst$ci_end),
      h1_sign < 0 & !h1_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h1_ci_hi, tab$arriaga$cnst$ci_end),
      h1_sign > 0 & h1_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h1_ci_hi, tab$arriaga$cnst$ci_end),
      h1_sign > 0 & !h1_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h1_ci_hi, tab$arriaga$cnst$ci_end)
    ),
    # h2
    h2_glyph = case_when(
      h2_attribution_1 == 'primarily' & h2_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h2_attribution_1 == 'primarily' & h2_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h2_attribution_1 == 'solely' & h2_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h2_attribution_1 == 'solely' & h2_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young
    ),
    h2_glyph = case_when(
      h2_sign < 0 & h2_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_glyph, h2_glyph),
      h2_sign < 0 & !h2_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_glyph, h2_glyph),
      h2_sign > 0 & h2_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_glyph, h2_glyph),
      h2_sign > 0 & !h2_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_glyph, h2_glyph)
    ),
    h2_estimate = FormatTable(h2_effect, 12),
    h2_estimate = case_when(
      h2_sign < 0 & h2_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h2_estimate),
      h2_sign < 0 & !h2_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h2_estimate),
      h2_sign > 0 & h2_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h2_estimate),
      h2_sign > 0 & !h2_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h2_estimate)
    ),
    h2_ci_lo = FormatTable(h2_lo, 12),
    h2_ci_lo = case_when(
      h2_sign < 0 & h2_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_significant_ci, h2_ci_lo, '{;}'),
      h2_sign < 0 & !h2_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_nonsignificant_ci, h2_ci_lo, '{;}'),
      h2_sign > 0 & h2_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_significant_ci, h2_ci_lo, '{;}'),
      h2_sign > 0 & !h2_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_nonsignificant_ci, h2_ci_lo,'{;}')
    ),
    h2_ci_hi = FormatTable(h2_hi, 12),
    h2_ci_hi = case_when(
      h2_sign < 0 & h2_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h2_ci_hi, tab$arriaga$cnst$ci_end),
      h2_sign < 0 & !h2_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h2_ci_hi, tab$arriaga$cnst$ci_end),
      h2_sign > 0 & h2_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h2_ci_hi, tab$arriaga$cnst$ci_end),
      h2_sign > 0 & !h2_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h2_ci_hi, tab$arriaga$cnst$ci_end)
    ),
    # h3
    h3_glyph = case_when(
      h3_attribution_1 == 'primarily' & h3_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h3_attribution_1 == 'primarily' & h3_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h3_attribution_1 == 'solely' & h3_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h3_attribution_1 == 'solely' & h3_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young
    ),
    h3_glyph = case_when(
      h3_sign < 0 & h3_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_glyph, h3_glyph),
      h3_sign < 0 & !h3_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_glyph, h3_glyph),
      h3_sign > 0 & h3_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_glyph, h3_glyph),
      h3_sign > 0 & !h3_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_glyph, h3_glyph)
    ),
    h3_estimate = FormatTable(h3_effect, 12),
    h3_estimate = case_when(
      h3_sign < 0 & h3_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h3_estimate),
      h3_sign < 0 & !h3_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h3_estimate),
      h3_sign > 0 & h3_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h3_estimate),
      h3_sign > 0 & !h3_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h3_estimate)
    ),
    h3_ci_lo = FormatTable(h3_lo, 12),
    h3_ci_lo = case_when(
      h3_sign < 0 & h3_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_significant_ci, h3_ci_lo, '{;}'),
      h3_sign < 0 & !h3_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_nonsignificant_ci, h3_ci_lo, '{;}'),
      h3_sign > 0 & h3_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_significant_ci, h3_ci_lo, '{;}'),
      h3_sign > 0 & !h3_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_nonsignificant_ci, h3_ci_lo, '{;}')
    ),
    h3_ci_hi = FormatTable(h3_hi, 12),
    h3_ci_hi = case_when(
      h3_sign < 0 & h3_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h3_ci_hi, tab$arriaga$cnst$ci_end),
      h3_sign < 0 & !h3_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h3_ci_hi, tab$arriaga$cnst$ci_end),
      h3_sign > 0 & h3_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h3_ci_hi, tab$arriaga$cnst$ci_end),
      h3_sign > 0 & !h3_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h3_ci_hi, tab$arriaga$cnst$ci_end)
    ),
    # h4
    h4_estimate = FormatTable(h4_effect, 100),
    h4_estimate = case_when(
      h4_sign < 0 & h4_significance ~ paste0(tab$arriaga$cnst$color_negative_significant, h4_estimate),
      h4_sign < 0 & !h4_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant, h4_estimate),
      h4_sign > 0 & h4_significance ~ paste0(tab$arriaga$cnst$color_positive_significant, h4_estimate),
      h4_sign > 0 & !h4_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant, h4_estimate)
    ),
    h4_ci_lo = FormatTable(h4_lo, 100),
    h4_ci_lo = case_when(
      h4_sign < 0 & h4_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_significant_ci, h4_ci_lo, '{;}'),
      h4_sign < 0 & !h4_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_nonsignificant_ci, h4_ci_lo, '{;}'),
      h4_sign > 0 & h4_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_significant_ci, h4_ci_lo, '{;}'),
      h4_sign > 0 & !h4_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_nonsignificant_ci, h4_ci_lo, '{;}')
    ),
    h4_ci_hi = FormatTable(h4_hi, 100),
    h4_ci_hi = case_when(
      h4_sign < 0 & h4_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h4_ci_hi, tab$arriaga$cnst$ci_end),
      h4_sign < 0 & !h4_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h4_ci_hi, tab$arriaga$cnst$ci_end),
      h4_sign > 0 & h4_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h4_ci_hi, tab$arriaga$cnst$ci_end),
      h4_sign > 0 & !h4_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h4_ci_hi, tab$arriaga$cnst$ci_end)
    ),
    # h5
    h5_glyph = case_when(
      h5_attribution_1 == 'primarily' & h5_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_overall_old,
      h5_attribution_1 == 'primarily' & h5_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young,
      h5_attribution_1 == 'solely' & h5_attribution_2 == 'old' ~ tab$arriaga$cnst$glyph_solely_old,
      h5_attribution_1 == 'solely' & h5_attribution_2 == 'young' ~ tab$arriaga$cnst$glyph_overall_young
    ),
    h5_glyph = case_when(
      h5_sign < 0 & h5_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_glyph, h5_glyph),
      h5_sign < 0 & !h5_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_glyph, h5_glyph),
      h5_sign > 0 & h5_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_glyph, h5_glyph),
      h5_sign > 0 & !h5_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_glyph, h5_glyph)
    ),
    h5_estimate = FormatTable(h5_effect, 12),
    h5_estimate = case_when(
      h5_sign < 0 & h5_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h5_estimate),
      h5_sign < 0 & !h5_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h5_estimate),
      h5_sign > 0 & h5_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h5_estimate),
      h5_sign > 0 & !h5_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h5_estimate)
    ),
    h5_ci_lo = FormatTable(h5_lo, 12),
    h5_ci_lo = case_when(
      h5_sign < 0 & h5_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_significant_ci, h5_ci_lo, '{;}'),
      h5_sign < 0 & !h5_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_negative_nonsignificant_ci, h5_ci_lo, '{;}'),
      h5_sign > 0 & h5_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_significant_ci, h5_ci_lo, '{;}'),
      h5_sign > 0 & !h5_significance ~ paste0(tab$arriaga$cnst$ci_start, tab$arriaga$cnst$color_positive_nonsignificant_ci, h5_ci_lo, '{;}')
    ),
    h5_ci_hi = FormatTable(h5_hi, 12),
    h5_ci_hi = case_when(
      h5_sign < 0 & h5_significance ~ paste0(tab$arriaga$cnst$color_negative_significant_ci, h5_ci_hi, tab$arriaga$cnst$ci_end),
      h5_sign < 0 & !h5_significance ~ paste0(tab$arriaga$cnst$color_negative_nonsignificant_ci, h5_ci_hi, tab$arriaga$cnst$ci_end),
      h5_sign > 0 & h5_significance ~ paste0(tab$arriaga$cnst$color_positive_significant_ci, h5_ci_hi, tab$arriaga$cnst$ci_end),
      h5_sign > 0 & !h5_significance ~ paste0(tab$arriaga$cnst$color_positive_nonsignificant_ci, h5_ci_hi, tab$arriaga$cnst$ci_end)
    )
  ) %>%
  select(
    region_name_short, sex,
    h1_glyph, h1_estimate, h1_ci_lo, h1_ci_hi,
    h2_glyph, h2_estimate, h2_ci_lo, h2_ci_hi,
    h3_glyph, h3_estimate, h3_ci_lo, h3_ci_hi,
    #          h4_estimate, h4_ci_lo, h4_ci_hi,
    h5_glyph, h5_estimate, h5_ci_lo, h5_ci_hi
  )

name <- 'table'
strata <- c('T', 'F', 'M')
tab$arriaga$table <-
  map(strata, ~{
    table <-
      tab$arriaga$data$table %>%
      filter(sex == .x) %>%
      select(-sex) %>%
      gt(rowname_col = 'region_name_short') %>%
      cols_label(
        h1_glyph = 'AT', h1_estimate = 'ES', h1_ci_lo = 'CI', h1_ci_hi = '',
        h2_glyph = 'AT', h2_estimate = 'ES', h2_ci_lo = 'CI', h2_ci_hi = '',
        h3_glyph = 'AT', h3_estimate = 'ES', h3_ci_lo = 'CI', h3_ci_hi = '',
        #                 h4_estimate = 'ES', h4_ci_lo = 'CI', h4_ci_hi = ''
        h5_glyph = 'AT', h5_estimate = 'ES', h5_ci_lo = 'CI', h5_ci_hi = ''
      ) %>%
      tab_spanner(columns = 2:5, 'Net LE losses 2019 to 2021', id = 'h1') %>%
      tab_spanner(columns = 6:9, 'LE losses 2020', id = 'h2') %>%
      tab_spanner(columns = 10:13, 'LE losses 2021', id = 'h3') %>%
      tab_spanner(columns = 14:17, 'LE deficit 2021', id = 'h5') %>%
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
    table
  })
names(tab$arriaga$table) <- paste0(name, strata)

tab$arriaga$table$tableT

as_latex(tab$arriaga$table$tableT) %>% as.character() %>%
  gsub( x = ., pattern = '\\$', replacement = '$', fixed = TRUE) %>%
  gsub( x = ., pattern = '\\textbackslash ', replacement = '\\', fixed = TRUE) %>%
  gsub( x = ., pattern = '\\textbackslash{}', replacement = '\\', fixed = TRUE) %>%
  gsub( x = ., pattern = '\\{', replacement = '{', fixed = TRUE) %>%
  gsub( x = ., pattern = '\\}', replacement = '}', fixed = TRUE) %>%
  clipr::write_clip()

# Export ----------------------------------------------------------

saveRDS(tab$arriaga, file = paths$output$tab_arriaga)
tab$arriaga$data$hypotheses %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_arriaga)
