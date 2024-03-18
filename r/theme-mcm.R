library(ggthemes)

theme_mcm <- function() {
  font <- "Noto Sans"
  dgy <- "grey15"
  lgy <- "#cacaca"
  mrg <- 2
  theme_few() +
    theme(
      plot.title = element_text(
        family = font,
        size = 9,
        face = 'bold',
        margin = margin(b = mrg)
      ),
      
      plot.subtitle = element_text(family = font,
                                   size = 8),
      
      plot.caption = element_text(family = font,
                                  size = 8),
      
      panel.border = element_rect(
        color = dgy,
        size = 0.25,
        fill = NA
      ),
      
      
      axis.title = element_text(
        family = font,
        size = 9,
        color = dgy,
      ),
      axis.title.y.left = element_text(margin = margin(r = mrg)),
      axis.title.y.right = element_text(margin = margin(l = mrg)),
      axis.title.x.bottom =  element_text(margin = margin(t = mrg)),
      axis.title.x.top =  element_text(margin = margin(b = mrg)),
      
      axis.text = element_text(
        family = font,
        size = 8,
        color = dgy,
      ),
      
      axis.line = element_line(color = 'black',
                               size = 0),
      axis.ticks = element_line(color = 'black',
                                size = 0.5),
      
      text = element_text(size = 8, family = font),
      
      strip.text = element_text(
        family = font,
        size = 8,
        face = 'bold',
        hjust = 0,
        margin = margin(b = mrg)
      ),
      
      legend.text = element_text(size = 8, color = dgy),
      
      legend.title = element_text(size = 9),
      
      panel.grid.major.y = element_line(
        color = lgy,
        size = 0.15,
        linetype = 1
      ),
      
      panel.grid.major.x = element_line(
        color = lgy,
        size = 0.15,
        linetype = 1
      ),
      line = element_line(size = 0.25),
    )
}

# Set theme globally
theme_set(theme_mcm())