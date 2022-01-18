# Generate figure of sex differences in life expectancy changes

# Init ------------------------------------------------------------

library(yaml); library(tidyverse); library(prismatic)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  region_metadata = './cfg/region_metadata.csv',
  figspecs = './cfg/figure_specification.R',
  sexdiff = './out/sexdiff.rds',
  lifetables = './out/lifetables.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  out = './out'
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
fig <- list()

# Functions -------------------------------------------------------

# figure specifications
source(paths$input$figspec)

# Load ex differences ---------------------------------------------

dat$sexdiff <- readRDS(paths$input$sexdiff)
dat$e0diff <- readRDS(paths$input$lifetables)

# Sex differences in life expectancy ------------------------------

fig$sexdiff <- list()
fig$sexdiff$cnst <- list(); fig$sexdiff$cnst <- within(fig$sexdiff$cnst, {
  font_size_annotation = 2
  font_size_annotation_numeric = font_size_annotation*0.8
  lighten_uncertainty = 0.4
  width_uncertainty = 1
  font_family_label = 'robotocondensed'
  font_family_theme = 'robotocondensed'
  font_size_theme = 6
  font_color_country = 'grey50'
  color_positive = '#005784'
  color_negative = '#B70D0D'
  color_positive_light =
    clr_mix(color_positive, mix_in = '#FFFFFF', ratio = lighten_uncertainty)
  color_negative_light =
    clr_mix(color_negative, mix_in = '#FFFFFF', ratio = lighten_uncertainty)
  color_estimate = 'black'
  color_estimate_light =
    clr_mix(color_estimate, mix_in = '#FFFFFF', ratio = lighten_uncertainty)
  nudge_value = 0.25
  nudge_country = 0.6
})
fig$sexdiff$func$FormatLabel <- function (x) {
  lab <- formatC(x, digits = 2, format = 'f')
  str_replace_all(lab, '(^[+-])(0)(.+)', '\\1\\3')
}

fig$sexdiff$dat <- list()
fig$sexdiff$dat <-
  dat$sexdiff %>%
  filter(age == 0, year %in% c(2019, 2021),
         region_iso %in% cnst$regions_for_analysis) %>%
  mutate(
    # if the sex difference change from 2019 is
    # positive|negative, is it so with a probability of 95% or higher
    ex_sexdiff_change_sign_sig = ifelse(
      ex_diff_change_from_2019_mean < 0,
      ex_diff_drop_from_2019_flag_mean >= 0.95,
      ex_diff_rise_from_2019_flag_mean >= 0.95
    )
  ) %>%
  select(region_iso, year, age,
         ex_sexdiff_q500 = ex_q0.5,
         ex_sexdiff_q025 = ex_q0.025,
         ex_sexdiff_q975 = ex_q0.975,
         ex_sexdiff_change_sign_sig) %>%
  mutate(
    region_pos =
      as.integer(fct_reorder(region_iso, ex_sexdiff_q500))
  ) %>%
  pivot_wider(
    names_from = year,
    values_from = c(ex_sexdiff_q500, ex_sexdiff_q025, ex_sexdiff_q975, ex_sexdiff_change_sign_sig)
  ) %>%
  mutate(
    ex_sexdiff_change_sign = sign(ex_sexdiff_q500_2021 - ex_sexdiff_q500_2019)
  ) %>%
  left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2'))
fig$sexdiff$dat <-
  fig$sexdiff$dat %>%
  mutate(region_name = case_when(
    region_name == 'England and Wales' ~ 'England\nand Wales',
    region_name == 'Northern Ireland' ~ 'Northern\nIreland',
    TRUE ~ region_name
  ))

fig$sexdiff$plot <- list()
fig$sexdiff$plot <-
  ggplot(fig$sexdiff$dat) +
  # grid
  geom_hline(yintercept = 0, color = 'grey', size = 1) +
  # central estimates
  geom_segment(
    aes(
      x = region_pos, xend = region_pos,
      y = ex_sexdiff_q500_2019, yend = ex_sexdiff_q500_2021,
      color = as.character(ex_sexdiff_change_sign),
      alpha = ex_sexdiff_change_sign_sig_2021
    ),
    arrow = arrow(length = unit(1, 'mm'))
  ) +
  # numeric labels
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q500_2021,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q500_2021)
    ),
    color = fig$sexdiff$cnst$color_negative_light,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = -fig$sexdiff$cnst$nudge_value),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(ex_sexdiff_change_sign < 0)
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q500_2021,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q500_2021)
    ),
    color = fig$sexdiff$cnst$color_positive_light,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = +fig$sexdiff$cnst$nudge_value),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(ex_sexdiff_change_sign >= 0)
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q500_2019,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q500_2019)
    ),
    color = fig$sexdiff$cnst$color_negative_light,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = +fig$sexdiff$cnst$nudge_value),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(ex_sexdiff_change_sign < 0)
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q500_2019,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q500_2019)
    ),
    color = fig$sexdiff$cnst$color_positive_light,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = -fig$sexdiff$cnst$nudge_value),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(ex_sexdiff_change_sign >= 0)
  ) +
  # country labels
  geom_text(
    aes(
      x = region_pos,
      y = pmin(ex_sexdiff_q500_2019, ex_sexdiff_q500_2021),
      label = region_name,
    ),
    color = fig$sexdiff$cnst$font_color_country,
    family = fig$sexdiff$cnst$font_family_label,
    lineheight = 1,
    position = position_nudge(y = -fig$sexdiff$cnst$nudge_country),
    data =
      . %>% filter(region_pos%%2==1),
    size = fig$sexdiff$cnst$font_size_annotation
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = pmax(ex_sexdiff_q500_2019, ex_sexdiff_q500_2021),
      label = region_name,
    ),
    lineheight = 1,
    color = fig$sexdiff$cnst$font_color_country,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = fig$sexdiff$cnst$nudge_country),
    data =
      . %>% filter(region_pos%%2==0),
    size = fig$sexdiff$cnst$font_size_annotation
  ) +
  # meta plot
  scale_color_manual(
    values = c('-1' = fig$sexdiff$cnst$color_negative,
               '1' = fig$sexdiff$cnst$color_positive)
  ) +
  scale_x_continuous(breaks = 1:30, labels = NULL) +
  scale_y_continuous(breaks = 0:10) +
  scale_alpha_manual(values = c(0.3, 1)) +
  fig_spec$MyGGplotTheme(
    show_legend = FALSE, grid = 'y', axis = '',
    family = fig$sexdiff$cnst$font_family_theme,
    size = fig$sexdiff$cnst$font_size_theme
  ) +
  theme(panel.grid.major.y =
          element_line(color = 'grey90', size = 0.3, linetype = 1)) +
  coord_cartesian(ylim = c(0, 11), xlim = c(0,30), expand = c(0,0)) +
  labs(x = NULL, y = 'Female life expectancy advantage in years')

fig$sexdiff$plot

# Export ----------------------------------------------------------

fig_spec$ExportFigure(
  fig$sexdiff$plot, device = 'pdf',
  filename = 'sexdiff',
  path = paths$output$out,
  width = fig_spec$width, height = 100
)
