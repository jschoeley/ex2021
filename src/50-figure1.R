# Generate Figure of annual e0 changes

# Init ------------------------------------------------------------

library(yaml); library(tidyverse)
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
  lifetables = './out/lifetables.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  data = './dat/output_data.rds',
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

# copied from ggplot2 `geom_curve`
geom_curve2 <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        curvature = 0.5,
                        angle = 90,
                        ncp = 5,
                        arrow = NULL,
                        lineend = "butt",
                        inflect = FALSE,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomCurve2, # call `GeomCurve2` instead of `GeomCurve`
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      arrow = arrow,
      curvature = curvature,
      angle = angle,
      ncp = ncp,
      lineend = lineend,
      inflect = inflect,
      na.rm = na.rm,
      ...
    )
  )
}

# copied from ggplot2 `GeomCurve`
GeomCurve2 <-
  ggproto(
    "GeomCurve2", GeomSegment,
    # the following `default_aes =` statement is missing in ggplot2 `GeomCurve`
    default_aes = aes(colour = "black", fill = "black", size = 0.5, linetype = 1, alpha = NA),
    draw_panel = function(data, panel_params, coord, curvature = 0.5, angle = 90,
                          ncp = 5, arrow = NULL, lineend = "butt", inflect = FALSE, na.rm = FALSE) {
      if (!coord$is_linear()) {
        warning("geom_curve is not implemented for non-linear coordinates",
                call. = FALSE)
      }
      trans <- coord$transform(data, panel_params)
      
      grid::curveGrob(
        trans$x, trans$y, trans$xend, trans$yend,
        default.units = "native",
        curvature = curvature, angle = angle, ncp = ncp,
        square = FALSE, squareShape = 1, inflect = inflect, open = TRUE,
        gp = grid::gpar(
          col = alpha(trans$colour, trans$alpha),
          # the following `fill = ` statement is missing in ggplot2 `GeomCurve`
          fill = alpha(trans$fill, trans$alpha),
          lwd = trans$size * .pt,
          lty = trans$linetype,
          lineend = lineend),
        arrow = arrow
      )
    }
  )

BounceBack <- function (delta1, delta2) {
  bb <- (1 - (delta1 + delta2) / delta1)*100
  ifelse(delta1 > 0 & delta2 > 0, NA, bb)
}

FormatTable <- function (x) {
  lab <- formatC(x, digits = 1, format = 'f', flag = '+')
  ifelse(lab == 'NA', '\u00b7', lab)
}

# Load ex differences ---------------------------------------------

dat$lifetables <- readRDS(paths$input$lifetables)

# Prepare data ----------------------------------------------------

fig$e0diff <- list()
fig$e0diff$data <-
  dat$lifetables %>%
  filter(age == 0, sex == 'T', year %in% 2020:2021,
         region_iso %in% cnst$regions_for_analysis) %>%
  select(region_iso, sex, year, age,
         ex_diff_q0.5, ex_diff_q0.025, ex_diff_q0.975,
         bbi_q0.5, bbi_q0.025, bbi_q0.975) %>%
  pivot_wider(
    id_cols = c(region_iso, sex, age),
    names_from = year,
    values_from = c(ex_diff_q0.5, ex_diff_q0.025, ex_diff_q0.975,
                    bbi_q0.5, bbi_q0.025, bbi_q0.975)
  ) %>%
  mutate(
    region_position =
      as.integer(fct_reorder(region_iso, -(ex_diff_q0.5_2020+ex_diff_q0.5_2021)))
  ) %>%
  left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2'))

# Parameterize figure ---------------------------------------------

fig$e0diff$cnst <- list()
fig$e0diff$cnst <- within(fig$e0diff$cnst, {
  n_countries = length(unique(fig$e0diff$data$region_iso))
  vertical_gap = 0.3
  curvature = 0.7
  size_vline = 3
  ribbons_size = 7*0.7
  segment_size = 1.2
  segment_nudge_y = -0.12
  color_vline = '#FFFFFF'
  color_line = '#FF007C'
  color_positive = '#005784'
  color_positive_light = '#92CFEC'
  color_negative = '#B70D0D'
  color_negative_light = '#FFB3AB'
  color_ribbon = '#E5E5E5'
  text_x_position_min = 1.0
  text_x_position_max = 1.7
  text_x_position1 = text_x_position_min
  text_x_position2 = text_x_position_min+(text_x_position_max-text_x_position_min)/2
  text_x_position3 = text_x_position_max
  text_x_position4 = 2.1
  font_table = 'Roboto Condensed'
  font_countries = 'Roboto Condensed'
  font_xaxis = 'Roboto Condensed'
  fontsize_table = 3.3
})

# Create figure ---------------------------------------------------

fig$e0diff$plot <-
  fig$e0diff$data %>%
  ggplot() +
  
# Grid lines ------------------------------------------------------

  geom_hline(
    yintercept =
      seq(1, fig$e0diff$cnst$n_countries, 2),
    color = fig$e0diff$cnst$color_ribbon,
    size = fig$e0diff$cnst$ribbons_size
  ) +
  geom_hline(
    yintercept = 0,
    color = fig$e0diff$cnst$color_ribbon, size = 1
  ) +
  geom_vline(
    xintercept = seq(-3, 1, 0.5),
    color = '#FFFFFF', size = 0.5
  ) +

# e0 diff 19/20 ---------------------------------------------------
# 
#   geom_segment(
#     aes(
#       y = region_position-
#         fig$e0diff$cnst$segment_nudge_y,
#       yend = region_position-
#         fig$e0diff$cnst$segment_nudge_y,
#       x = ex_diff_q0.025_2020, xend = ex_diff_q0.975_2020,
#       color = as.character(sign(ex_diff_q0.5_2020)+3)
#     ), size = 1.8
#   ) +
  geom_segment(
    aes(
      y = region_position-fig$e0diff$cnst$segment_nudge_y,
      yend = region_position-fig$e0diff$cnst$segment_nudge_y,
      x = 0, xend = ex_diff_q0.5_2020
    ),
    size = fig$e0diff$cnst$segment_size,
    color = fig$e0diff$cnst$color_positive,
    data =
      filter(fig$e0diff$data, ex_diff_q0.5_2020 > 0)
  ) +
  geom_segment(
    aes(
      y = region_position-fig$e0diff$cnst$segment_nudge_y,
      yend = region_position-fig$e0diff$cnst$segment_nudge_y,
      x = 0, xend = ex_diff_q0.5_2020
    ),
    size = fig$e0diff$cnst$segment_size,
    color = fig$e0diff$cnst$color_negative,
    data =
      filter(fig$e0diff$data, ex_diff_q0.5_2020 <= 0)
  ) +
  
  # Curve connector -------------------------------------------------

  geom_curve2(
    aes(
      y = region_position-fig$e0diff$cnst$segment_nudge_y,
      yend = region_position - fig$e0diff$cnst$segment_nudge_y -
        fig$e0diff$cnst$vertical_gap,
      x = ex_diff_q0.5_2020, xend = ex_diff_q0.5_2020
    ),
    inflect = FALSE, curvature = 1*fig$e0diff$cnst$curvature,
    size = fig$e0diff$cnst$segment_size,
    color = fig$e0diff$cnst$color_positive,
    data =
      filter(fig$e0diff$data, ex_diff_q0.5_2020 <= 0, ex_diff_q0.5_2021 > 0)
  ) +
  geom_curve2(
    aes(
      y = region_position - fig$e0diff$cnst$segment_nudge_y,
      yend = region_position - fig$e0diff$cnst$vertical_gap -
        fig$e0diff$cnst$segment_nudge_y,
      x = ex_diff_q0.5_2020, xend = ex_diff_q0.5_2020
    ),
    inflect = TRUE, curvature = 1*fig$e0diff$cnst$curvature,
    size = fig$e0diff$cnst$segment_size,
    color = fig$e0diff$cnst$color_negative,
    data =
      filter(fig$e0diff$data, ex_diff_q0.5_2020 <= 0, ex_diff_q0.5_2021 <= 0)
  ) +
  geom_curve2(
    aes(
      y = region_position - fig$e0diff$cnst$segment_nudge_y,
      yend = region_position - fig$e0diff$cnst$vertical_gap -
        fig$e0diff$cnst$segment_nudge_y,
      x = ex_diff_q0.5_2020, xend = ex_diff_q0.5_2020
    ),
    inflect = FALSE, curvature = -1*fig$e0diff$cnst$curvature,
    size = fig$e0diff$cnst$segment_size,
    color = fig$e0diff$cnst$color_negative,
    data =
      filter(fig$e0diff$data, ex_diff_q0.5_2020 > 0, ex_diff_q0.5_2021 <= 0)
  ) +
  geom_curve2(
    aes(
      y = region_position - fig$e0diff$cnst$segment_nudge_y,
      yend = region_position - fig$e0diff$cnst$vertical_gap -
        fig$e0diff$cnst$segment_nudge_y,
      x = ex_diff_q0.5_2020, xend = ex_diff_q0.5_2020
    ),
    inflect = TRUE, curvature = -1*fig$e0diff$cnst$curvature,
    size = fig$e0diff$cnst$segment_size,
    color = fig$e0diff$cnst$color_positive,
    data =
      filter(fig$e0diff$data, ex_diff_q0.5_2020 > 0, ex_diff_q0.5_2021 > 0)
  ) +
  geom_segment(
    x = 0, xend = 0, y = -0.1, yend = fig$e0diff$cnst$n_countries+1,
    color = fig$e0diff$cnst$color_vline,
    size = fig$e0diff$cnst$size_vline
  ) +
  geom_vline(
    xintercept = 0,
    color = fig$e0diff$cnst$color_ribbon,
    size = 1, lty = 3
  ) +
  
# e0 diff 20/21 ---------------------------------------------------

  geom_segment(
    aes(
      y = region_position - fig$e0diff$cnst$vertical_gap -
        fig$e0diff$cnst$segment_nudge_y,
      yend = region_position - fig$e0diff$cnst$vertical_gap -
        fig$e0diff$cnst$segment_nudge_y,
      x = ex_diff_q0.5_2020, xend = ex_diff_q0.5_2020 + ex_diff_q0.5_2021
    ),
    lineend = 'butt',
    arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
    size = fig$e0diff$cnst$segment_size,
    color = fig$e0diff$cnst$color_positive,
    data =
      filter(fig$e0diff$data, ex_diff_q0.5_2021 > 0)
  ) +
  geom_segment(
    aes(
      y = region_position - fig$e0diff$cnst$vertical_gap -
        fig$e0diff$cnst$segment_nudge_y,
      yend = region_position - fig$e0diff$cnst$vertical_gap -
        fig$e0diff$cnst$segment_nudge_y,
      x = ex_diff_q0.5_2020, xend = ex_diff_q0.5_2020 + ex_diff_q0.5_2021
    ),
    lineend = 'butt',
    arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
    size = fig$e0diff$cnst$segment_size,
    color = fig$e0diff$cnst$color_negative,
    data =
      filter(fig$e0diff$data, ex_diff_q0.5_2021 <= 0)
  ) +
  
# Table -----------------------------------------------------------

  geom_text(
    aes(
      y = region_position,
      label = FormatTable(ex_diff_q0.5_2020*12),
      color = as.character(sign(ex_diff_q0.5_2020))
    ),
    x = fig$e0diff$cnst$text_x_position2,
    family = fig$e0diff$cnst$font_table, hjust = 1,
    size = fig$e0diff$cnst$fontsize_table,
    size = fig$e0diff$cnst$fontsize_table
  ) +
  geom_text(
    aes(
      y = region_position,
      label = FormatTable(ex_diff_q0.5_2021*12),
      color = as.character(sign(ex_diff_q0.5_2021))
    ),
    x = fig$e0diff$cnst$text_x_position3,
    family = fig$e0diff$cnst$font_table, hjust = 1,
    size = fig$e0diff$cnst$fontsize_table
  ) +
  geom_text(
    aes(
      y = region_position,
      label = FormatTable((ex_diff_q0.5_2020+ex_diff_q0.5_2021)*12),
      color = as.character(sign(ex_diff_q0.5_2020+ex_diff_q0.5_2021))
    ),
    x = fig$e0diff$cnst$text_x_position1,
    family = fig$e0diff$cnst$font_table,
    hjust = 1, fontface = 'bold',
    size = fig$e0diff$cnst$fontsize_table
  ) +
  geom_text(
    aes(
      y = region_position,
      label = FormatTable(bbi_q0.5_2021*100),
      color = as.character(sign(bbi_q0.5_2021))
    ),
    x = fig$e0diff$cnst$text_x_position4,
    family = fig$e0diff$cnst$font_table, hjust = 1,
    fontface = 'bold',
    size = fig$e0diff$cnst$fontsize_table
  ) +

# Table header ----------------------------------------------------

  annotate(
    'text',
    x = fig$e0diff$cnst$text_x_position2,
    y = fig$e0diff$cnst$n_countries+1.7,
    label = '\u0394~e[0]',
    family = fig$e0diff$cnst$font_table,
    color = 'black', parse = TRUE, hjust = 1,
    fontface = 'bold', size = fig$e0diff$cnst$fontsize_table
  ) +
  geom_segment(
    y = fig$e0diff$cnst$n_countries+1.4,
    yend = fig$e0diff$cnst$n_countries+1.4,
    x = fig$e0diff$cnst$text_x_position1-0.25,
    xend = fig$e0diff$cnst$text_x_position3,
    size = 0.2, size = fig$e0diff$cnst$fontsize_table
  ) +
  annotate(
    'text',
    x = fig$e0diff$cnst$text_x_position1,
    y = fig$e0diff$cnst$n_countries+1,
    label = '19/21',
    color = 'black', parse = FALSE,
    family = fig$e0diff$cnst$font_table, hjust = 1,
    fontface = 'bold', size = fig$e0diff$cnst$fontsize_table
  ) +
  annotate(
    'text',
    x = fig$e0diff$cnst$text_x_position2,
    y = fig$e0diff$cnst$n_countries+1,
    label = '19/20',
    color = 'black', parse = FALSE,
    family = fig$e0diff$cnst$font_table, hjust = 1,
    size = fig$e0diff$cnst$fontsize_table
  ) +
  annotate(
    'text',
    x = fig$e0diff$cnst$text_x_position3,
    y = fig$e0diff$cnst$n_countries+1,
    label = '20/21',
    color = 'black', parse = FALSE,
    family = fig$e0diff$cnst$font_table, hjust = 1,
    size = fig$e0diff$cnst$fontsize_table
  ) +
  annotate(
    'text',
    x = fig$e0diff$cnst$text_x_position4,
    y = fig$e0diff$cnst$n_countries+1,
    label = 'BBI%',
    color = 'black', parse = FALSE,
    family = fig$e0diff$cnst$font_table, hjust = 1,
    fontface = 'bold',
    size = fig$e0diff$cnst$fontsize_table
  ) +
  
# Scales and labels -----------------------------------------------

  scale_x_continuous(
    breaks = seq(-3, 1, 0.5),
    labels = c(
      '-36 months', '-30', '-24', '-18', '-12', '-6',
      'Life\nexpectancy\nin 2019',
      '+6', '+12 months'
    )
  ) +
  scale_y_continuous(
    breaks = unique(fig$e0diff$data$region_position),
    labels = unique(fig$e0diff$data$region_name),
    expand = c(0,0.3)
  ) +
  scale_color_manual(
    values = c(`-1` = fig$e0diff$cnst$color_negative,
               `1` = fig$e0diff$cnst$color_positive,
               `2` = fig$e0diff$cnst$color_negative_light,
               `3` = fig$e0diff$cnst$color_positive_light)
  ) +
  fig_spec$MyGGplotTheme(
    grid = '', family = fig$e0diff$cnst$font_xaxis, axis = '',
    size = 10, show_legend = FALSE
  ) +
  theme(plot.background = element_rect(colour = NA, fill = '#FFFFFF')) +
  coord_cartesian(
    xlim = c(-36/12, 23/12),
    ylim = c(0, fig$e0diff$cnst$n_countries+1.7)
  ) +
  labs(
    title = 'Bounce backs amid continued losses',
    subtitle = 'Life expectancy changes since COVID-19',
    caption = '@jschoeley',
    y = NULL, x = NULL
  )
fig$e0diff$plot

# Export ----------------------------------------------------------

extrafont::loadfonts()
fig_spec$ExportFigure(
  fig$e0diff$plot, device = 'svg',
  filename = 'e0diff',
  path = paths$output$out,
  #width = 164, height = 195
  width = 1.3*164, height = 150
)
