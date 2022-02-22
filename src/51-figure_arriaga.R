# Generate Figure of arriaga decompositions of e0 changes

# Init ------------------------------------------------------------

library(yaml); library(tidyverse); library(purrr)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  region_metadata = './cfg/region_metadata.csv',
  figspecs = './cfg/figure_specification.R',
  lifetables = './out/40-lifetables.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  fig_arriaga = './out',
  rds_arriagaT = './out/51-arriagaT.rds',
  csv_arriagaT = './tmp/51-arriagaT.csv',
  rds_arriagaF = './out/51-arriagaF.rds',
  csv_arriagaF = './tmp/51-arriagaF.csv',
  rds_arriagaM = './out/51-arriagaM.rds',
  csv_arriagaM = './tmp/51-arriagaM.csv'
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

# Load ex differences ---------------------------------------------

dat$lifetables <- readRDS(paths$input$lifetables)

# Plot ------------------------------------------------------------

name <- 'arriaga_'
strata <- c('T', 'F', 'M')

fig <- map(strata, ~{
  
  const <- list(); const <- within(const, {
    age_breaks = c(0, 20, 40, 60, 80, Inf)
    age_names = c('0-19', '20-39', '40-59', '60-79', '80+')
    segment_nudge_y = -0.12
    segment_size = 0.8
    color_positive = '#555555'
    color_negative = '#B70D0D'
    color_vline = 'white'
    vertical_gap = 0.3
    curvature = 0.7
    arrow_length = 0.5
    n_age = length(age_breaks)
    color_ribbon = 'white'
    size_ribbon = 6
    font_plot = 'roboto'
  })
  
  data <-
    dat$lifetables %>%
    filter(sex == .x, region_iso %in% cnst$regions_for_analysis,
           year >= 2020, projected == 'actual') %>%
    mutate(
      age_group = cut(as.integer(age),
                      const$age_breaks,
                      right = FALSE,
                      labels = const$age_names),
      age_position = as.integer(age_group)
    ) %>%
    select(region_iso, age, age_group, sex, year, age_position, e0_cntrb_t_q0.5) %>%
    group_by(region_iso, sex, year, age_group, age_position) %>%
    summarise(e0_cntrb_t_q0.5 = sum(e0_cntrb_t_q0.5)) %>%
    ungroup() %>%
    group_by(region_iso, sex, year) %>%
    mutate(e0_diff = sum(e0_cntrb_t_q0.5)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = year,
      values_from = c(e0_cntrb_t_q0.5, e0_diff)
    ) %>%
    ungroup() %>%
    left_join(
      region_meta, by = c('region_iso' = 'region_code_iso3166_2')
    )

  plot <-
    data %>%
    ggplot(aes(y = age_position, yend = age_position)) +
    geom_hline(
      yintercept =
        seq(1, const$n_age, 2),
      color = const$color_ribbon,
      size = const$size_ribbon
    ) +
    geom_vline(
      xintercept = seq(-2, 0.5, 0.5)*12,
      color = '#FFFFFF', size = 0.3
    ) +
    geom_segment(
      aes(
        y = age_position-const$segment_nudge_y,
        yend = age_position-const$segment_nudge_y,
        x = 0, xend = e0_cntrb_t_q0.5_2020*12
      ),
      size = const$segment_size,
      color = const$color_positive,
      data =
        . %>% filter(e0_cntrb_t_q0.5_2020 > 0)
    ) +
    geom_segment(
      aes(
        y = age_position-const$segment_nudge_y,
        yend = age_position-const$segment_nudge_y,
        x = 0, xend = e0_cntrb_t_q0.5_2020*12
      ),
      size = const$segment_size,
      color = const$color_negative,
      data =
        . %>% filter(e0_cntrb_t_q0.5_2020 <= 0)
    ) +
    # central vline
    geom_vline(
      xintercept = 0,
      color = const$color_vline,
      size = 0.8, lty = 1
    ) +
    geom_segment(
      aes(
        y = age_position - const$vertical_gap -
          const$segment_nudge_y,
        yend = age_position - const$vertical_gap -
          const$segment_nudge_y,
        x = e0_cntrb_t_q0.5_2020*12,
        xend = (e0_cntrb_t_q0.5_2020 + e0_cntrb_t_q0.5_2021)*12
      ),
      lineend = 'round',
      arrow = arrow(length = unit(const$arrow_length, 'mm'), angle = 30),
      size = const$segment_size,
      color = const$color_positive,
      data =
        . %>% filter(e0_cntrb_t_q0.5_2021 > 0)
    ) +
    geom_segment(
      aes(
        y = age_position - const$vertical_gap -
          const$segment_nudge_y,
        yend = age_position - const$vertical_gap -
          const$segment_nudge_y,
        x = e0_cntrb_t_q0.5_2020*12,
        xend = (e0_cntrb_t_q0.5_2020 + e0_cntrb_t_q0.5_2021)*12
      ),
      lineend = 'round',
      arrow = arrow(length = unit(const$arrow_length, 'mm'), angle = 30),
      size = const$segment_size,
      color = const$color_negative,
      data =
        . %>% filter(e0_cntrb_t_q0.5_2021 <= 0)
    ) +
    # bounce back
    geom_curve2(
      aes(
        y = age_position-const$segment_nudge_y,
        yend = age_position - const$segment_nudge_y -
          const$vertical_gap,
        x = e0_cntrb_t_q0.5_2020*12, xend = (e0_cntrb_t_q0.5_2020+0.01)*12
      ),
      inflect = FALSE, curvature = 1*const$curvature,
      size = const$segment_size,
      color = const$color_positive,
      data =
        . %>% filter(e0_cntrb_t_q0.5_2020 <= 0, e0_cntrb_t_q0.5_2021 > 0)
    ) +
    # compound losses
    geom_curve2(
      aes(
        y = age_position - const$segment_nudge_y,
        yend = age_position - const$vertical_gap -
          const$segment_nudge_y,
        x = e0_cntrb_t_q0.5_2020*12, xend = (e0_cntrb_t_q0.5_2020+0.01)*12
      ),
      inflect = TRUE, curvature = 1*const$curvature,
      size = const$segment_size,
      color = const$color_negative,
      data =
        . %>% filter(e0_cntrb_t_q0.5_2020 <= 0, e0_cntrb_t_q0.5_2021 <= 0)
    ) +
    # late losses
    geom_curve2(
      aes(
        y = age_position - const$segment_nudge_y,
        yend = age_position - const$vertical_gap -
          const$segment_nudge_y,
        x = e0_cntrb_t_q0.5_2020*12, xend = (e0_cntrb_t_q0.5_2020+0.01)*12
      ),
      inflect = FALSE, curvature = -1*const$curvature,
      size = const$segment_size,
      color = const$color_negative,
      data =
        . %>% filter(e0_cntrb_t_q0.5_2020 > 0, e0_cntrb_t_q0.5_2021 <= 0)
    ) +
    # compound gains
    geom_curve2(
      aes(
        y = age_position - const$segment_nudge_y,
        yend = age_position - const$vertical_gap -
          const$segment_nudge_y,
        x = e0_cntrb_t_q0.5_2020*12, xend = (e0_cntrb_t_q0.5_2020+0.01)*12
      ),
      inflect = TRUE, curvature = -1*const$curvature,
      size = const$segment_size,
      color = const$color_positive,
      data =
        . %>% filter(e0_cntrb_t_q0.5_2020 > 0, e0_cntrb_t_q0.5_2021 > 0)
    ) +
    scale_x_continuous(
      breaks = seq(-3, 0.5, 0.5)*12,
      labels = c(
        '-36', '-30' , '-24', '-18', '-12', '-6',
        'LE\nin 2019',
        '+6'
      )
    ) +
    scale_y_continuous(
      breaks = unique(data$age_position),
      labels = const$age_names,
      expand = c(0,0.3)
    ) +
    coord_cartesian(xlim = c(NA, 0.2*12)) +
    facet_wrap(~region_name, ncol = 4) +
    fig_spec$MyGGplotTheme(
      grid = '', family = const$font_plot, axis = '',
      size = 10, show_legend = FALSE
    ) +
    theme(
      panel.background = element_rect(fill = 'grey95', color = 'grey95'),
      strip.background = element_rect(fill = 'grey90', color = 'grey95')
    ) +
    labs(
      x = 'Agewise contributions to months of life expectancy change since 2019',
      y = 'Age group'
    )

  list(cnst = const, data = data, plot = plot)
      
})
names(fig) <- paste0(name, strata)

# Export ----------------------------------------------------------

fig_spec$ExportFigure(
  fig$arriaga_T$plot, device = 'pdf',
  filename = '51-arriagaT',
  path = paths$output$fig_arriaga,
  width = fig_spec$width, height = 200, scale = 1.2
)
saveRDS(fig$arriaga_T, file = paths$output$rds_arriagaT)
fig$arriaga_T$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_arriagaT)

fig_spec$ExportFigure(
  fig$arriaga_F$plot, device = 'pdf',
  filename = '51-arriagaF',
  path = paths$output$fig_arriaga,
  width = fig_spec$width, height = 200, scale = 1.2
)
saveRDS(fig$arriaga_F, file = paths$output$rds_arriagaF)
fig$arriaga_F$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_arriagaF)

fig_spec$ExportFigure(
  fig$arriaga_M$plot, device = 'pdf',
  filename = '51-arriagaM',
  path = paths$output$fig_arriaga,
  width = fig_spec$width, height = 200, scale = 1.2
)
saveRDS(fig$arriaga_M, file = paths$output$rds_arriagaM)
fig$arriaga_M$data %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_arriagaM)
