# ==============================================================================
# 00_setup.R
# libraries & defining parameters
# ==============================================================================




# libraries ====================================================================

library(tidyverse)
library(readxl)

library(scales)
library(ComplexHeatmap)
library(cowplot)
library(circlize)
library(ggbeeswarm)



# palettes =====================================================================

palette_vaccine <-
  c(RZV = "#969696", ZVL = "#CFCFCF")




# custom functions =============================================================

# format a continuous variable as mean (SD) ------------------------------------
# (looks excessive, but copied over from personal utilities package 'choomisc')

summaryMeanSD <-
  function(x, digits = 2, na.rm = F) {
    
    if ( !(typeof(x) %in% c("double", "integer")) ) {
      warning("input vector not numeric")
      x <- as.numeric(x)
    }
    tot <- length(x)
    n_NA <- sum(is.na(x))
    
    if (n_NA != 0) {
      warning(paste0(n_NA,
                     " values are missing/NA in input vector"))
    }
    
    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }
    
    if (na.rm == T) {
      x <- x[!is.na(x)]
    }
    
    output <- c(mean(x), sd(x))
    output <- formatC(output, format = "f", digits = digits)
    
    output <- paste0( output[1], " (", output[2], ")")
    
    return(output)
    
  }

# custom boxplot style ---------------------------------------------------------
# boxplot style by Dr. Kong, via https://github.com/rpkgs/Ipaper
# copied code here because library(Ipaper) failed installation

geom_boxplot2 <-
  function(mapping = NULL, data = NULL,
           stat = "boxplot", position = "dodge2",
           ...,
           outlier.colour = NULL,
           outlier.color = NULL,
           outlier.fill = NULL,
           outlier.shape = 19,
           outlier.size = 1.5,
           outlier.stroke = 0.5,
           outlier.alpha = NULL,
           show.errorbar = TRUE,
           width.errorbar = 0.7,
           notch = FALSE,
           notchwidth = 0.5,
           varwidth = FALSE,
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = TRUE) {
    if (is.character(position)) {
      if (varwidth == TRUE) position <- position_dodge2(preserve = "single")
    } else {
      if (identical(position$preserve, "total") & varwidth == TRUE) {
        warning("Can't preserve total widths when varwidth = TRUE.",
                call. = FALSE)
        position$preserve <- "single"
      }
    }

    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomBoxplot2,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        outlier.colour = outlier.color %||% outlier.colour,
        outlier.fill = outlier.fill,
        outlier.shape = outlier.shape,
        outlier.size = outlier.size,
        outlier.stroke = outlier.stroke,
        outlier.alpha = outlier.alpha,
        show.errorbar = show.errorbar,
        width.errorbar = width.errorbar,
        notch = notch,
        notchwidth = notchwidth,
        varwidth = varwidth,
        na.rm = na.rm,
        ...
      )
    )
  }

GeomBoxplot2 <-
  ggproto("GeomBoxplot2", Geom,
    extra_params = c("na.rm", "width"),
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)

      data$outliers <- NULL
      if (!is.null(data$outliers)) {
        suppressWarnings({
          out_min <- vapply(data$outliers, min, numeric(1))
          out_max <- vapply(data$outliers, max, numeric(1))
        })

        data$ymin_final <- pmin(out_min, data$ymin)
        data$ymax_final <- pmax(out_max, data$ymax)
      }

      if (is.null(params) || is.null(params$varwidth) ||
          !params$varwidth || is.null(data$relvarwidth)) {
        data$xmin <- data$x - data$width / 2
        data$xmax <- data$x + data$width / 2
      } else {
        data$relvarwidth <- data$relvarwidth / max(data$relvarwidth)
        data$xmin <- data$x - data$relvarwidth * data$width / 2
        data$xmax <- data$x + data$relvarwidth * data$width / 2
      }
      if (!is.null(data$relvarwidth)) data$relvarwidth <- NULL
      data
    },
    draw_group = function(data, panel_params, coord, fatten = 2,
                          outlier.colour = NULL, outlier.fill = NULL,
                          outlier.shape = 19,
                          outlier.size = 1.5, outlier.stroke = 0.5,
                          outlier.alpha = NULL,
                          show.errorbar = TRUE,
                          width.errorbar = 0.7,
                          notch = FALSE, notchwidth = 0.5, varwidth = FALSE) {
      common <- list(
        colour = data$colour,
        size = data$size,
        linetype = data$linetype,
        fill = alpha(data$fill, data$alpha),
        group = data$group
      )

      whiskers <- new_data_frame(c(
        list(
          x = c(data$x, data$x),
          xend = c(data$x, data$x),
          y = c(data$upper, data$lower),
          yend = c(data$ymax, data$ymin),
          alpha = c(NA_real_, NA_real_)
        ),
        common
      ), n = 2)

      box <- new_data_frame(c(
        list(
          xmin = data$xmin,
          xmax = data$xmax,
          ymin = data$lower,
          y = data$middle,
          ymax = data$upper,
          ynotchlower = ifelse(notch, data$notchlower, NA),
          ynotchupper = ifelse(notch, data$notchupper, NA),
          notchwidth = notchwidth,
          alpha = data$alpha
        ),
        common
      ))

      errorbar <- new_data_frame(c(
        list(
          xmin = data$x - width.errorbar / 2,
          xmax = data$x + width.errorbar / 2,
          x = data$x,
          ymin = data$ymin,
          ymax = data$ymax,
          alpha = data$alpha
        ),
        common
      ))

      grob_whiskers <- GeomSegment$draw_panel(whiskers, panel_params, coord)
      grob_errorbar <- NULL

      if (show.errorbar) {
        grob_errorbar <- GeomErrorbar$draw_panel(errorbar, panel_params, coord)
      }
      ggplot2:::ggname("geom_boxplot2", grobTree(
        grob_errorbar,
        GeomCrossbar$draw_panel(box, fatten = fatten, panel_params, coord)
      ))
    },
    draw_key = draw_key_boxplot,
    default_aes = aes(
      weight = 1, colour = "grey20", fill = "white", size = 0.5,
      alpha = NA, shape = 19, linetype = "solid"
    ),
    required_aes = c("x", "lower", "upper", "middle", "ymin", "ymax")
  )

`%||%` <- function(a, b) {
  if (!is.null(a)) {
    a
  } else {
    b
  }
}

new_data_frame <-
  function(x = list(), n = NULL) {
    if (length(x) != 0 && is.null(names(x)))
      stop("Elements must be named", call. = FALSE)
    lengths <- vapply(x, length, integer(1))
    if (is.null(n)) {
      n <- if (length(x) == 0) 0 else max(lengths)
    }
    for (i in seq_along(x)) {
      if (lengths[i] == n) next
      if (lengths[i] != 1)
        stop("Elements must equal the number of rows or 1", call. = FALSE)
      x[[i]] <- rep(x[[i]], n)
    }

    class(x) <- "data.frame"

    attr(x, "row.names") <- .set_row_names(n)
    x
  }



