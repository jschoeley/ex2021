# Generate Figure of annual e0 changes

# Init ------------------------------------------------------------

library(yaml); library(tidyverse); library(readr)

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
  e0avgdiff = './out/40-e0avgdiff.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  fig_e0diff = './out',
  rds_e0diff = './out',
  csv_e0diffT = './tmp/50-e0diffT.csv',
  csv_e0diffF = './tmp/50-e0diffF.csv',
  csv_e0diffM = './tmp/50-e0diffM.csv'
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

FormatTable <- function (x) {
  lab <- formatC(x, digits = 1, format = 'f', flag = '+')
  ifelse(lab == 'NA', '\u00b7', lab)
}

# Load ex differences ---------------------------------------------

dat$lifetables <- readRDS(paths$input$lifetables)
dat$e0avgdiff <- readRDS(paths$input$e0avgdiff)

# completeness of data
dat$lt_input <- readRDS('./out/30-lt_input.rds')
dat$completeness <-
  dat$lt_input %>%
  filter(year == 2021, age_start == 0, sex == 'Male') %>%
  select(year, region_iso, death_total_nweeksmiss) %>%
  mutate(as_of = 52-death_total_nweeksmiss)

# Prepare data ----------------------------------------------------

name <- 'e0diff_'
strata <- c('T', 'F', 'M')
fig <- map(strata, ~{

  data <- list()
  data$e0diff2021 <-
    dat$lifetables %>%
    filter(age == 0, sex == .x, year %in% 2020:2021,
           region_iso %in% cnst$regions_for_analysis,
           projected == 'actual', quarter == 'annual') %>%
    select(region_iso, sex, year, age,
           ex_diff_mean, ex_diff_q0.025, ex_diff_q0.975,
           bbi_q0.5, bbi_q0.025, bbi_q0.975) %>%
    pivot_wider(
      id_cols = c(region_iso, sex, age),
      names_from = year,
      values_from = c(ex_diff_mean, ex_diff_q0.025, ex_diff_q0.975,
                      bbi_q0.5, bbi_q0.025, bbi_q0.975)
    ) %>%
    mutate(
      region_position =
        as.integer(fct_reorder(region_iso, -(ex_diff_mean_2020+ex_diff_mean_2021)))
    ) %>%
    left_join(region_meta, by = c('region_iso' = 'region_code_iso3166_2')) %>%
    left_join(dat$completeness)
  
  data$e0diff2021 <-
    data$e0diff2021 %>%
    left_join(
      dat$e0avgdiff %>%
        filter(age == 0, sex == .x, projected == 'actual',
               quarter == 'annual') %>%
        select(region_iso, e0avgdiff1619_q0.5 = q0.5,
               e0avgdiff1619_q0.025 = q0.025,
               e0avgdiff1619_q0.975 = q0.975)
    )
  
  # Parameterize figure ---------------------------------------------
  
  cnst <- list()
  cnst <- within(cnst, {
    n_countries = length(unique(data$e0diff2021$region_iso))
    vertical_gap = 0.3
    curvature = 0.7
    size_vline = 3
    ribbons_size = 7*1
    segment_size = 1.2
    segment_nudge_y = -0.12
    color_vline = '#FFFFFF'
    color_line = '#FF007C'
    color_positive = '#005784'
    color_positive_light = '#92CFEC'
    color_negative = '#B70D0D'
    color_negative_light = '#FFB3AB'
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
  
  plot <-
    data$e0diff2021 %>%
    ggplot() +
    
    # grid lines
    
    geom_hline(
      yintercept =
        seq(1, cnst$n_countries, 2),
      color = cnst$color_ribbon,
      size = cnst$ribbons_size
    ) +
    geom_hline(
      yintercept = 0,
      color = cnst$color_ribbon, size = 1
    ) +
    geom_vline(
      xintercept = seq(-3.5, 0.5, 0.5),
      color = '#FFFFFF', size = 0.5
    ) +
    
    # e0 diff 19/20
    
    geom_segment(
      aes(
        y = region_position-cnst$segment_nudge_y,
        yend = region_position-cnst$segment_nudge_y,
        x = 0, xend = ex_diff_mean_2020
      ),
      size = cnst$segment_size,
      color = cnst$color_positive,
      data =
        . %>% filter(ex_diff_mean_2020 > 0)
    ) +
    geom_segment(
      aes(
        y = region_position-cnst$segment_nudge_y,
        yend = region_position-cnst$segment_nudge_y,
        x = 0, xend = ex_diff_mean_2020
      ),
      size = cnst$segment_size,
      color = cnst$color_negative,
      data =
        . %>% filter(ex_diff_mean_2020 <= 0)
    ) +
    
    # curve connector
    
    geom_curve2(
      aes(
        y = region_position-cnst$segment_nudge_y,
        yend = region_position - cnst$segment_nudge_y -
          cnst$vertical_gap,
        x = ex_diff_mean_2020, xend = ex_diff_mean_2020
      ),
      inflect = FALSE, curvature = 1*cnst$curvature,
      size = cnst$segment_size,
      color = cnst$color_positive,
      data =
        . %>% filter(ex_diff_mean_2020 <= 0, ex_diff_mean_2021 > 0)
    ) +
    geom_curve2(
      aes(
        y = region_position - cnst$segment_nudge_y,
        yend = region_position - cnst$vertical_gap -
          cnst$segment_nudge_y,
        x = ex_diff_mean_2020, xend = ex_diff_mean_2020
      ),
      inflect = TRUE, curvature = 1*cnst$curvature,
      size = cnst$segment_size,
      color = cnst$color_negative,
      data =
        . %>% filter(ex_diff_mean_2020 <= 0, ex_diff_mean_2021 <= 0)
    ) +
    geom_curve2(
      aes(
        y = region_position - cnst$segment_nudge_y,
        yend = region_position - cnst$vertical_gap -
          cnst$segment_nudge_y,
        x = ex_diff_mean_2020, xend = ex_diff_mean_2020
      ),
      inflect = FALSE, curvature = -1*cnst$curvature,
      size = cnst$segment_size,
      color = cnst$color_negative,
      data =
        . %>% filter(ex_diff_mean_2020 > 0, ex_diff_mean_2021 <= 0)
    ) +
    geom_curve2(
      aes(
        y = region_position - cnst$segment_nudge_y,
        yend = region_position - cnst$vertical_gap -
          cnst$segment_nudge_y,
        x = ex_diff_mean_2020, xend = ex_diff_mean_2020
      ),
      inflect = TRUE, curvature = -1*cnst$curvature,
      size = cnst$segment_size,
      color = cnst$color_positive,
      data =
        . %>% filter(ex_diff_mean_2020 > 0, ex_diff_mean_2021 > 0)
    ) +
    geom_segment(
      x = 0, xend = 0, y = -0.1, yend = cnst$n_countries+1,
      color = cnst$color_vline,
      size = cnst$size_vline
    ) +
    geom_vline(
      xintercept = 0,
      color = cnst$color_ribbon,
      size = 0.5, lty = 1
    ) +
    
    # e0 diff 20/21
    
    geom_segment(
      aes(
        y = region_position - cnst$vertical_gap -
          cnst$segment_nudge_y,
        yend = region_position - cnst$vertical_gap -
          cnst$segment_nudge_y,
        x = ex_diff_mean_2020, xend = ex_diff_mean_2020 + ex_diff_mean_2021
      ),
      lineend = 'round',
      arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
      size = cnst$segment_size,
      color = cnst$color_positive,
      data =
        . %>% filter(ex_diff_mean_2021 > 0)
    ) +
    geom_segment(
      aes(
        y = region_position - cnst$vertical_gap -
          cnst$segment_nudge_y,
        yend = region_position - cnst$vertical_gap -
          cnst$segment_nudge_y,
        x = ex_diff_mean_2020, xend = ex_diff_mean_2020 + ex_diff_mean_2021
      ),
      lineend = 'round',
      arrow = arrow(length = unit(1.2, 'mm'), angle = 30),
      size = cnst$segment_size,
      color = cnst$color_negative,
      data =
        . %>% filter(ex_diff_mean_2021 <= 0)
    ) +
    
    # table
    
    geom_text(
      aes(
        y = region_position,
        label = FormatTable(ex_diff_mean_2020*12),
        color = as.character(sign(ex_diff_mean_2020))
      ),
      x = cnst$text_x_position2,
      family = cnst$font_table, hjust = 1,
      size = cnst$fontsize_table,
      size = cnst$fontsize_table
    ) +
    geom_text(
      aes(
        y = region_position,
        label = FormatTable(ex_diff_mean_2021*12),
        color = as.character(sign(ex_diff_mean_2021))
      ),
      x = cnst$text_x_position3,
      family = cnst$font_table, hjust = 1,
      size = cnst$fontsize_table
    ) +
    geom_text(
      aes(
        y = region_position,
        label = FormatTable((ex_diff_mean_2020+ex_diff_mean_2021)*12),
        color = as.character(sign(ex_diff_mean_2020+ex_diff_mean_2021))
      ),
      x = cnst$text_x_position1,
      family = cnst$font_table,
      hjust = 1, fontface = 'bold',
      size = cnst$fontsize_table
    ) +
    # geom_text(
    #   aes(
    #     y = region_position,
    #     label = as_of
    #   ),
    #   color = '#666666',
    #   x = cnst$text_x_position4,
    #   family = cnst$font_table, hjust = 1,
    #   #fontface = 'bold',
    #   size = cnst$fontsize_table
    # ) +
    
    # avg e0 diff 2016-19
    
    geom_segment(
      aes(
        x = e0avgdiff1619_q0.025, xend = e0avgdiff1619_q0.975,
        y = region_position, yend = region_position
      ),
      size = 0.5,
      color = cnst$fill_avge0diff
    ) +
    geom_point(
      aes(
        x = e0avgdiff1619_q0.5,
        y = region_position
      ),
      size = 1,
      color = cnst$color_avge0diff
    ) + 
    geom_text(
      aes(
        y = region_position,
        label = FormatTable(e0avgdiff1619_q0.5*12),
      ),
      color = cnst$color_avge0diff,
      x = cnst$text_x_position5,
      family = cnst$font_table,
      hjust = 1,
      size = cnst$fontsize_table
    ) +
    
    # table header
    
    annotate(
      'text',
      x = (cnst$text_x_position2+cnst$text_x_position3)/2,
      y = cnst$n_countries+1.7,
      label = '\u0394~e[0]',
      family = cnst$font_table,
      color = '#666666', parse = TRUE, hjust = 1,
      fontface = 'bold', size = cnst$fontsize_table
    ) +
    geom_segment(
      y = cnst$n_countries+1.4,
      yend = cnst$n_countries+1.4,
      x = cnst$text_x_position1-0.25,
      xend = cnst$text_x_position5,
      size = 0.2, size = cnst$fontsize_table
    ) +
    annotate(
      'text',
      x = cnst$text_x_position1,
      y = cnst$n_countries+1,
      label = '19/21',
      color = '#666666', parse = FALSE,
      family = cnst$font_table, hjust = 1,
      fontface = 'bold', size = cnst$fontsize_table
    ) +
    annotate(
      'text',
      x = cnst$text_x_position2,
      y = cnst$n_countries+1,
      label = '19/20',
      color = '#666666', parse = FALSE,
      family = cnst$font_table, hjust = 1,
      size = cnst$fontsize_table
    ) +
    annotate(
      'text',
      x = cnst$text_x_position3,
      y = cnst$n_countries+1,
      label = '20/21',
      color = '#666666', parse = FALSE,
      family = cnst$font_table, hjust = 1,
      size = cnst$fontsize_table
    ) +
    annotate(
      'text',
      x = cnst$text_x_position5,
      y = cnst$n_countries+1,
      label = 'Av[16-19]',
      color = '#666666', parse = TRUE,
      family = cnst$font_table, hjust = 1,
      #fontface = 'bold',
      size = cnst$fontsize_table
    ) +
    # annotate(
    #   'text',
    #   x = cnst$text_x_position4,
    #   y = cnst$n_countries+1,
    #   label = 'As of\nweek', vjust = 0,
    #   color = '#666666', parse = FALSE, lineheight = 0.7,
    #   family = cnst$font_table, hjust = 1,
    #   #fontface = 'bold',
    #   size = cnst$fontsize_table
    # ) +
    
    # scales and labels
    
    scale_x_continuous(
      breaks = seq(-3.5, 0.5, 0.5),
      labels = c(
        '-42 months', '-36', '-30', '-24', '-18', '-12', '-6',
        'Life\nexpectancy\nin 2019',
        '+6 months'
      )
    ) +
    scale_y_continuous(
      breaks = unique(data$e0diff2021$region_position),
      labels = unique(data$e0diff2021$region_name),
      expand = c(0,0.3)
    ) +
    scale_color_manual(
      values = c(`-1` = cnst$color_negative,
                 `1` = cnst$color_positive,
                 `2` = cnst$color_negative_light,
                 `3` = cnst$color_positive_light)
    ) +
    fig_spec$MyGGplotTheme(
      grid = '', family = cnst$font_xaxis, axis = '',
      size = 10, show_legend = FALSE
    ) +
    theme(plot.background = element_rect(colour = NA, fill = '#FFFFFF')) +
    coord_cartesian(
      xlim = c(-40/12, 26/12),
      ylim = c(0, cnst$n_countries+1.8)
    ) +
    labs(
      #title = 'Bounce backs amid continued losses',
      #subtitle = 'Life expectancy changes since COVID-19. Estimates for 2021 preliminary and adjusted for missing data.',
      #caption = '@jschoeley @jm_aburto @ridhikash07 @MaxiKniffka',
      y = NULL, x = NULL
    )
  #fig[[paste0(name, .x)]]$plot
  list(cnst = cnst, data = data, plot = plot)
})
names(fig) <- paste0(name, strata)

# Export ----------------------------------------------------------

fig_spec$ExportFigure(
  fig$e0diff_T$plot, device = 'pdf',
  filename = '50-e0diffT',
  path = paths$output$fig_e0diff,
  width = fig_spec$width, height = 200
)
saveRDS(fig$e0diff_T, file = paste0(paths$output$rds_e0diff, '/50-e0diffT.rds'))
fig$e0diff_T$data$e0diff2021 %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0diffT)

fig_spec$ExportFigure(
  fig$e0diff_F$plot, device = 'pdf',
  filename = '50-e0diffF',
  path = paths$output$fig_e0diff,
  width = fig_spec$width, height = 200
)
saveRDS(fig$e0diff_F, file = paste0(paths$output$rds_e0diff, '/50-e0diffF.rds'))
fig$e0diff_F$data$e0diff2021 %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0diffF)

fig_spec$ExportFigure(
  fig$e0diff_M$plot, device = 'pdf',
  filename = '50-e0diffM',
  path = paths$output$fig_e0diff,
  width = fig_spec$width, height = 200
)
saveRDS(fig$e0diff_M, file = paste0(paths$output$rds_e0diff, '/50-e0diffM.rds'))
fig$e0diff_M$data$e0diff2021 %>%
  mutate(across(where(is.numeric), ~round(.x, 6))) %>%
  write_csv(paths$output$csv_e0diffM)
