#these functions were copied from 
#https://github.com/clauswilke/gridtext/blob/master/examples/ggplot2-integration.R
#by Claus Wilke, in order to render titles in bold.

# define new theme element that inherits from `element_text()`
element_markdown <- function(family = NULL, face = NULL, colour = NULL, size = NULL,
                          hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
                          color = NULL, margin = NULL,
                          debug = FALSE, inherit.blank = FALSE) {
  if (!is.null(color))
    colour <- color
  structure(
    list(
      family = family, face = face, colour = colour,
      size = size, hjust = hjust, vjust = vjust, angle = angle,
      lineheight = lineheight, margin = margin, debug = debug,
      inherit.blank = inherit.blank),
    class = c("element_markdown", "element_text", "element")
  )
}

# rendering of the theme element is handled by `richtext_grob()`
element_grob.element_markdown <- function(element, label = "", x = NULL, y = NULL,
                                          family = NULL, face = NULL, colour = NULL, size = NULL,
                                          hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
                                          margin = NULL, margin_x = FALSE, margin_y = FALSE, ...) {
  if (is.null(label))
    return(ggplot2:::zeroGrob())

  vj <- vjust %||% element$vjust
  hj <- hjust %||% element$hjust
  margin <- margin %||% element$margin %||% ggplot2::margin(0, 0, 0, 0)
  angle <- angle %||% element$angle %||% 0

  x <- x %||% hj
  if (!is.unit(x))
    x <- unit(x, "npc")
  y <- y %||% vj
  if (!is.unit(y))
    y <- unit(y, "npc")

  # The gp settings can override element_gp
  gp <- gpar(
    fontsize = size %||% element$size,
    col = colour %||% element$colour,
    fontfamily = family %||% element$family,
    fontface = face %||% element$face,
    lineheight = lineheight %||% element$lineheight
  )

  richtext_grob(
    label, x = x, y = y, hjust = hj, vjust = vj, rot = angle,
    padding = margin, gp = gp
  )
}
