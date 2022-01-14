# Generate figure of sex differences in life expectancy changes

# Init ------------------------------------------------------------

library(yaml); library(tidyverse); library(prismatic)
#extrafont::font_import(prompt = FALSE)

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
  lighten_uncertainty = 0.8
  width_uncertainty = 1
  font_family_label = 'roboto_condensed'
  font_family_theme = 'roboto_condensed'
  font_size_theme = 6
  font_color_country = 'grey50'
  color_estimate = 'black'
  color_estimate_light =
    clr_mix(color_estimate, mix_in = '#FFFFFF', ratio = lighten_uncertainty)
  color_positive_light2 =
    clr_mix(color_estimate, mix_in = '#FFFFFF', ratio = lighten_uncertainty*0.5)
})
fig$sexdiff$func$FormatLabel <- function (x) {
  lab <- formatC(x, digits = 2, format = 'f', flag = '+')
  str_replace_all(lab, '(^[+-])(0)(.+)', '\\1\\3')
}

fig$sexdiff$dat <- list()
fig$sexdiff$dat$e0 <-
  dat$sexdiff %>%
  filter(age == 0, year %in% c(2019, 2021),
         region_iso %in% cnst$regions_for_analysis) %>%
  select(region_iso, year, age,
         ex_sexdiff_q500 = ex_q0.5,
         ex_sexdiff_q025 = ex_q0.025,
         ex_sexdiff_q975 = ex_q0.975) %>%
  mutate(
    region_pos =
      as.integer(fct_reorder(region_iso, ex_sexdiff_q500))
  ) %>%
  left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2'))

fig$sexdiff$plot <- list()
fig$sexdiff$plot$e0 <-
  ggplot(fig$sexdiff$dat$e0) +
  # grid
  geom_hline(yintercept = 0, color = 'grey', size = 1) +
  # central estimates & uncertainty
  geom_linerange(
    aes(
      x = region_pos,
      ymin = ex_sexdiff_q025,
      ymax = ex_sexdiff_q975,
      color = as.character(sign(ex_sexdiff_q500)+3)
    ),
    size = fig$sexdiff$cnst$width_uncertainty,
    data = . %>% filter(year == 2021)
  ) +
  geom_point(
    aes(
      x = region_pos,
      y = ex_sexdiff_q500,
      color = as.character(sign(ex_sexdiff_q500))
    ),
    data = . %>% filter(year == 2021)
  ) +
  geom_point(
    aes(
      x = region_pos,
      y = ex_sexdiff_q500,
      color = as.character(sign(ex_sexdiff_q500))
    ),
    shape = 21,
    data = . %>% filter(year == 2019)
  ) +
  # numeric labels
  # white offset
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q025,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q025),
    ),
    family = fig$sexdiff$cnst$font_family_label, color = 'white',
    position = position_nudge(y = -0.055, x = -0.005),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(year == 2021)
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q975,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q975)
    ),
    color = 'white',
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = 0.055, x = -0.005),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(year == 2021)
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q500,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q500)
    ),
    color = 'white',
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(x = 0.42, y = -0.005),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(year == 2021)
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q025,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q025),
      color = as.character(sign(ex_sexdiff_q500)+6)
    ),
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = -0.05),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(year == 2021)
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q975,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q975),
      color = as.character(sign(ex_sexdiff_q500)+6)
    ),
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = 0.05),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(year == 2021)
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q500,
      label = fig$sexdiff$func$FormatLabel(ex_sexdiff_q500),
      color = as.character(sign(ex_sexdiff_q500)+6)
    ),
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(x = 0.45),
    size = fig$sexdiff$cnst$font_size_annotation_numeric,
    data = . %>% filter(year == 2021)
  ) +
  # country labels
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q025,
      label = region_name,
    ),
    color = fig$sexdiff$cnst$font_color_country,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = -0.15),
    data =
      . %>% filter(region_pos%%2==1, year == 2021),
    size = fig$sexdiff$cnst$font_size_annotation
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_sexdiff_q975,
      label = region_name,
    ),
    color = fig$sexdiff$cnst$font_color_country,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = 0.15),
    data =
      . %>% filter(region_pos%%2==0, year == 2021),
    size = fig$sexdiff$cnst$font_size_annotation
  ) +
  # meta plot
  scale_color_manual(
    values = c('1' = fig$sexdiff$cnst$color_estimate,
               '4' = fig$sexdiff$cnst$color_estimate_light,
               '7' = fig$sexdiff$cnst$color_estimate_light2)
  ) +
  scale_x_continuous(breaks = 1:30, labels = NULL) +
  fig_spec$MyGGplotTheme(
    show_legend = FALSE, grid = 'y', axis = '',
    family = fig$sexdiff$cnst$font_family_theme,
    size = fig$sexdiff$cnst$font_size_theme
  ) +
  theme(panel.grid.major.y =
          element_line(color = 'grey70', size = 0.1, linetype = 1)) +
  coord_cartesian(ylim = c(0, 11), xlim = c(0,30), expand = c(0,0)) +
  labs(x = NULL, y = NULL)

# Prepare data ----------------------------------------------------

fig$sexdiff <- list()
fig$sexdiff$data <-
  dat$sexdiff %>%
  filter(age == 0, year %in% 2021,
         # +2: ex increase for both sexes
         # -2: ex decrease for both sexes
         # 0: mixed increase-decrease by sex
         ex_diff_sign_q0.5 == -2,
         region_iso %in% cnst$regions_for_analysis) %>%
  select(region_iso, year, age,
         ex_diff_sexdiff_q500 = ex_diff_q0.5,
         ex_diff_sexdiff_q025 = ex_diff_q0.025,
         ex_diff_sexdiff_q975 = ex_diff_q0.975) %>%
  mutate(
    region_pos =
      as.integer(fct_reorder(region_iso, ex_diff_sexdiff_q500))
  ) %>%
  left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2'))

# Parameterize figure ---------------------------------------------

fig$sexdiff$cnst <- list(); fig$sexdiff$cnst <- within(fig$sexdiff$cnst, {
  
  font_size_annotation = 2
  font_size_annotation_numeric = font_size_annotation*0.8
  lighten_uncertainty = 0.8
  width_uncertainty = 2.5
  font_family_label = 'roboto'
  font_family_theme = 'roboto_condensed'
  font_size_theme = 6
  font_color_country = 'grey50'
  color_positive = '#005784'
  color_positive_light =
    clr_mix(color_positive, mix_in = '#FFFFFF', ratio = lighten_uncertainty)
  color_positive_light2 =
    clr_mix(color_positive, mix_in = '#FFFFFF', ratio = lighten_uncertainty*0.5)
  color_negative = '#B70D0D'
  color_negative_light =
    clr_mix(color_negative, mix_in = '#FFFFFF', ratio = lighten_uncertainty)
  color_negative_light2 =
    clr_mix(color_negative, mix_in = '#FFFFFF', ratio = lighten_uncertainty*0.5)
  
})

fig$sexdiff$func$FormatLabel <- function (x) {
  lab <- formatC(x, digits = 2, format = 'f', flag = '+')
  str_replace_all(lab, '(^[+-])(0)(.+)', '\\1\\3')
}

fig$sexdiff$plot <-
  fig$sexdiff$data %>%
  ggplot() +
  # grid
  geom_hline(yintercept = 0, color = 'grey', size = 1) +
  # central estimates & uncertainty
  geom_linerange(
    aes(
      x = region_pos,
      ymin = ex_diff_sexdiff_q025,
      ymax = ex_diff_sexdiff_q975,
      color = as.character(sign(ex_diff_sexdiff_q500)+3)
    ),
    size = fig$sexdiff$cnst$width_uncertainty
  ) +
  geom_point(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q500,
      color = as.character(sign(ex_diff_sexdiff_q500))
    )
  ) +
  # numeric labels
  # white offset
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q025,
      label = fig$sexdiff$func$FormatLabel(ex_diff_sexdiff_q025),
    ),
    family = fig$sexdiff$cnst$font_family_label, color = 'white',
    position = position_nudge(y = -0.055, x = -0.005),
    size = fig$sexdiff$cnst$font_size_annotation_numeric
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q975,
      label = fig$sexdiff$func$FormatLabel(ex_diff_sexdiff_q975)
    ),
    color = 'white',
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = 0.055, x = -0.005),
    size = fig$sexdiff$cnst$font_size_annotation_numeric
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q500,
      label = fig$sexdiff$func$FormatLabel(ex_diff_sexdiff_q500)
    ),
    color = 'white',
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(x = 0.42, y = -0.005),
    size = fig$sexdiff$cnst$font_size_annotation_numeric
  ) +
  # colored text
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q025,
      label = fig$sexdiff$func$FormatLabel(ex_diff_sexdiff_q025),
      color = as.character(sign(ex_diff_sexdiff_q500)+6)
    ),
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = -0.05),
    size = fig$sexdiff$cnst$font_size_annotation_numeric
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q975,
      label = fig$sexdiff$func$FormatLabel(ex_diff_sexdiff_q975),
      color = as.character(sign(ex_diff_sexdiff_q500)+6)
    ),
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = 0.05),
    size = fig$sexdiff$cnst$font_size_annotation_numeric
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q500,
      label = fig$sexdiff$func$FormatLabel(ex_diff_sexdiff_q500),
      color = as.character(sign(ex_diff_sexdiff_q500)+6)
    ),
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(x = 0.45),
    size = fig$sexdiff$cnst$font_size_annotation_numeric
  ) +
  # country labels
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q025,
      label = region_name,
    ),
    color = fig$sexdiff$cnst$font_color_country,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = -0.15),
    data =
      fig$sexdiff$data %>% filter(region_pos%%2==1, region_name != 'Russia'),
    size = fig$sexdiff$cnst$font_size_annotation
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q975,
      label = region_name,
    ),
    color = fig$sexdiff$cnst$font_color_country,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = 0.15),
    data =
      fig$sexdiff$data %>% filter(region_pos%%2==0, region_name != 'Russia'),
    size = fig$sexdiff$cnst$font_size_annotation
  ) +
  geom_text(
    aes(
      x = region_pos,
      y = ex_diff_sexdiff_q025,
      label = region_name,
    ),
    color = fig$sexdiff$cnst$font_color_country,
    family = fig$sexdiff$cnst$font_family_label,
    position = position_nudge(y = -0.15),
    data =
      fig$sexdiff$data %>% filter(region_name == 'Russia'),
    size = fig$sexdiff$cnst$font_size_annotation*1.2,
    fontface = 'bold'
  ) +
  # annotation
  annotate(
    'text', x = max(fig$sexdiff$data$region_pos),
    y = 1.15, label = 'Larger life expectancy drop for men',
    color = fig$sexdiff$cnst$color_positive,
    hjust = 1, vjust = 0, alpha = 0.5,
    family = fig$sexdiff$cnst$font_family_label,
    size = fig$sexdiff$cnst$font_size_annotation
  ) +
  annotate(
    'text', x = 0, y = -0.9, label = 'Larger life expectancy drop for women',
    color = fig$sexdiff$cnst$color_negative,
    hjust = 0, vjust = 1, alpha = 0.5,
    family = fig$sexdiff$cnst$font_family_label,
    size = fig$sexdiff$cnst$font_size_annotation
  ) +
  annotate(
    'text', x = 0, y = 1,
    label = 'Sex difference in years of life expectancy drop',
    hjust = 0, vjust = -1, fontface = 'bold',
    size = fig$sexdiff$cnst$font_size_annotation*1.3
  ) +
  # meta plot
  scale_color_manual(
    values = c('1' = fig$sexdiff$cnst$color_positive,
               '-1' = fig$sexdiff$cnst$color_negative,
               '2' = fig$sexdiff$cnst$color_negative_light,
               '4' = fig$sexdiff$cnst$color_positive_light,
               '5' = fig$sexdiff$cnst$color_negative_light2,
               '7' = fig$sexdiff$cnst$color_positive_light2)
  ) +
  scale_x_continuous(breaks = 1:30, labels = NULL) +
  fig_spec$MyGGplotTheme(
    show_legend = FALSE, grid = 'y', axis = '',
    family = fig$sexdiff$cnst$font_family_theme,
    size = fig$sexdiff$cnst$font_size_theme
  ) +
  theme(panel.grid.major.y =
          element_line(color = 'grey70', size = 0.1, linetype = 1)) +
  coord_cartesian(ylim = c(-1.4, 1.4), expand = c(0,0)) +
  labs(x = NULL, y = NULL)
fig$sexdiff$plot


# -----------------------------------------------------------------

dat$e0diff %>%
  filter(age == 0, year %in% 2019:2021) %>%
  select(region_iso, year, sex, ex_q0.5) %>%
  pivot_wider(names_from = sex, values_from = ex_q0.5) %>%
  ggplot(aes(x = F, y = M)) +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_abline(color = 'grey') +
  geom_path(aes(group = region_iso)) +
  coord_equal(xlim = c(65, 85), ylim = c(65, 85)) +
  fig_spec$MyGGplotTheme(family = 'roboto_condensed', axis = '', grid = 'xy')

# Export ----------------------------------------------------------

extrafont::loadfonts()
fig_spec$ExportFigure(
  fig$e0diff$plot, device = 'svg',
  filename = 'e0diff',
  path = paths$output$out,
  #width = 164, height = 195
  width = 1.3*164, height = 150
)
