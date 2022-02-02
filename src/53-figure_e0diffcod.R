# Generate Figure of annual e0 changes

# Init ------------------------------------------------------------

library(yaml); library(tidyverse)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  region_metadata = './cfg/region_metadata.csv',
  figspecs = './cfg/figure_specification.R',
  lifetables = './out/codecomp.rds',
  e0avgdiff = './out/e0avgdiff.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  data = './dat/output_data.rds',
  fig_e0diff = './out',
  rds_e0diff = './out/rds_e0diffcod.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# meta data on regions
region_meta <- read_csv(paths$input$region_meta, na = '.')

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  regions_for_analysis = config$regions_for_covid_cause_analysis
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
                        lineend = "round",
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
                          ncp = 5, arrow = NULL, lineend = "round", inflect = FALSE, na.rm = FALSE) {
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
dat$e0avgdiff <- readRDS(paths$input$e0avgdiff)

# completeness of data
dat$lt_input <- readRDS('./out/lt_input.rds')
dat$completeness <-
  dat$lt_input %>%
  filter(year == 2021, age_start == 0, sex == 'Male') %>%
  select(year, region_iso, death_total_nweeksmiss) %>%
  mutate(as_of = 52-death_total_nweeksmiss)

# Prepare data ----------------------------------------------------

fig$e0diffcod <- list()
fig$e0diffcod$data <- list()
fig$e0diffcod$data$e0diff2021 <-
  dat$lifetables %>%
  filter(sex == 'T', year %in% 2020:2021,
         region_iso %in% cnst$regions_for_analysis) %>%
  select(region_iso, sex, year,
         e0_cntrb_t_covid_q0.5, e0_cntrb_t_covid_q0.025, e0_cntrb_t_covid_q0.975,
         e0_cntrb_t_noncovid_q0.5, e0_cntrb_t_noncovid_q0.025, e0_cntrb_t_noncovid_q0.975) %>%
  pivot_wider(
    id_cols = c(region_iso, sex),
    names_from = year,
    values_from = c(e0_cntrb_t_covid_q0.5, e0_cntrb_t_covid_q0.025, e0_cntrb_t_covid_q0.975,
                    e0_cntrb_t_noncovid_q0.5, e0_cntrb_t_noncovid_q0.025, e0_cntrb_t_noncovid_q0.975)
  ) %>%
  mutate(
    region_position =
      as.integer(fct_reorder(region_iso, -(e0_cntrb_t_covid_q0.5_2020+e0_cntrb_t_noncovid_q0.5_2020)))
  ) %>%
  left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2')) %>%
  left_join(dat$completeness)

fig$e0diffcod$data$e0diff2021 <-
  fig$e0diffcod$data$e0diff2021 %>%
  left_join(
    dat$e0avgdiff %>%
      filter(age == 0, sex == 'T') %>%
      select(region_iso, e0avgdiff1619_q0.5 = q0.5,
             e0avgdiff1619_q0.025 = q0.025,
             e0avgdiff1619_q0.975 = q0.975)
  )

# Parameterize figure ---------------------------------------------

fig$e0diffcod$cnst <- list()
fig$e0diffcod$cnst <- within(fig$e0diffcod$cnst, {
  n_countries = length(unique(fig$e0diffcod$data$e0diff2021$region_iso))
  vertical_gap = 0.3
  curvature = 0.6
  size_vline = 0.4
  ribbons_size = 7*1
  segment_size = 3
  segment_nudge_y = -0.12
  color_vline = '#4A4A4A'
  color_line = '#FF007C'
  color_covid = '#005e59'
  color_covid_light = '#C341BA'
  color_noncovid = '#1A1A1A'
  color_noncovid_light = '#5C5C5C'
  color_ribbon = '#E5E5E5'
  fill_avge0diff = 'grey40'
  color_avge0diff = 'grey40'
  text_x_position_min = 1.0
  text_x_position_max = 2.3
  text_x_position1 = text_x_position_min
  text_x_position2 = text_x_position_min+(text_x_position_max-text_x_position_min)/3
  text_x_position3 = text_x_position_min+2*(text_x_position_max-text_x_position_min)/3
  text_x_position4 = 2.9
  text_x_position5 = text_x_position_max
  font_table = 'robotocondensed'
  font_countries = 'robotocondensed'
  font_xaxis = 'robotocondensed'
  fontsize_table = 3.3
})

# Create figure ---------------------------------------------------

fig$e0diffcod$plot <-
  fig$e0diffcod$data$e0diff2021 %>%
  ggplot() +
  
  # grid lines

  geom_vline(
    xintercept = seq(-3.5, 0.5, 0.5),
    color = '#FFFFFF', size = 0.5
  ) +
  
  # e0 diff 19/20 covid

  geom_segment(
    aes(
      y = 4, yend = 4,
      x = 0, xend = e0_cntrb_t_covid_q0.5_2020
    ),
    lineend = 'butt', linejoin = 'mitre',
    size = fig$e0diffcod$cnst$segment_size,
    color = fig$e0diffcod$cnst$color_covid
  ) +
  
  # central line
  
  geom_segment(
    x = 0, xend = 0, y = -0.1, yend = 4.3,
    color = fig$e0diffcod$cnst$color_vline,
    size = fig$e0diffcod$cnst$size_vline
  ) +
  
  # e0 diff 19/20 covid-noncovid connector

  geom_segment(
    aes(
      y = 4, yend = 3,
      x = e0_cntrb_t_covid_q0.5_2020,
      xend = e0_cntrb_t_covid_q0.5_2020
    ),
    lineend = 'butt', linejoin = 'mitre',
    size = 0.1*fig$e0diffcod$cnst$segment_size,
    color = fig$e0diffcod$cnst$color_covid
  ) +
  
  # e0 diff 19/20 noncovid
  
  geom_segment(
    aes(
      y = 3, yend = 3,
      x = e0_cntrb_t_covid_q0.5_2020,
      xend = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020
    ),
    lineend = 'butt', linejoin = 'mitre',
    size = fig$e0diffcod$cnst$segment_size,
    color = fig$e0diffcod$cnst$color_noncovid
  ) +

  
  # e0 diff 19/20 noncovid 20/21 covid connector
  
  geom_segment(
    aes(
      y = 3, yend = 2,
      x = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020,
      xend = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020
    ),
    lineend = 'butt', linejoin = 'mitre',
    size = 0.1*fig$e0diffcod$cnst$segment_size,
    color = fig$e0diffcod$cnst$color_noncovid
  ) +
    
  # e0 diff 20/21 covid
  
  geom_segment(
    aes(
      y = 2, yend = 2,
      x = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020,
      xend = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020 +
        e0_cntrb_t_covid_q0.5_2021
    ),
    lineend = 'butt', linejoin = 'mitre',
    size = fig$e0diffcod$cnst$segment_size,
    color = fig$e0diffcod$cnst$color_covid
  ) +
  
  # e0 diff 20/21 covid - noncovid connector
  
  geom_segment(
    aes(
      y = 2, yend = 1,
      x = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020 + e0_cntrb_t_covid_q0.5_2021,
      xend = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020 + e0_cntrb_t_covid_q0.5_2021
    ),
    lineend = 'butt', linejoin = 'mitre',
    size = 0.1*fig$e0diffcod$cnst$segment_size,
    color = fig$e0diffcod$cnst$color_covid
  ) +
  
  # e0 diff 20/21 non-covid
  
  geom_segment(
    aes(
      y = 1, yend = 1,
      x = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020 + e0_cntrb_t_covid_q0.5_2021,
      xend = e0_cntrb_t_covid_q0.5_2020 + e0_cntrb_t_noncovid_q0.5_2020 +
        e0_cntrb_t_covid_q0.5_2021 + e0_cntrb_t_noncovid_q0.5_2021
    ),
    lineend = 'butt', linejoin = 'mitre',
    arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
    size = fig$e0diffcod$cnst$segment_size,
    color = fig$e0diffcod$cnst$color_noncovid
  ) +
  
  # scales and labels
  
  scale_x_continuous(
    breaks = seq(-3, 0.5, 0.5),
    labels = c(
      '-36m', '', '-24', '', '-12', '',
      'Life\nexpectancy\nin 2019',
      ''
    )
  ) +
  scale_y_continuous(
    breaks = NULL,
    labels = NULL,
    expand = c(0,0.3)
  ) +
  scale_color_manual(
    values = c(`-1` = fig$e0diffcod$cnst$color_negative,
               `1` = fig$e0diffcod$cnst$color_positive,
               `2` = fig$e0diffcod$cnst$color_negative_light,
               `3` = fig$e0diffcod$cnst$color_positive_light)
  ) +
  facet_wrap(~region_name, ncol = 4) +
  fig_spec$MyGGplotTheme(
    grid = '', family = fig$e0diffcod$cnst$font_xaxis, axis = '',
    size = 10, show_legend = FALSE
  ) +
  theme(
    panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    #strip.background = element_rect(fill = 'grey90', color = 'grey95')
  ) +
  coord_cartesian(
    xlim = c(-36/12, 32/12),
    ylim = c(1, 4)
  ) +
  labs(
    #title = 'Bounce backs amid continued losses',
    #subtitle = 'Life expectancy changes since COVID-19. Estimates for 2021 preliminary and adjusted for missing data.',
    #caption = '@jschoeley @jm_aburto @ridhikash07 @MaxiKniffka',
    y = NULL, x = NULL
  )
fig$e0diffcod$plot

# Export ----------------------------------------------------------

fig_spec$ExportFigure(
  fig$e0diffcod$plot, device = 'pdf',
  filename = 'e0diffcod',
  path = paths$output$fig_e0diff,
  width = fig_spec$width, height = fig_spec$width
)

saveRDS(fig$e0diffcod, file = paths$output$rds_e0diff)
