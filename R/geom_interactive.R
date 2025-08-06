#' @export
geom_segment2_interactive <- function(...){
  rlang::check_installed('ggiraph', "for `geom_segment2_interactive()`.")
  layer_interactive(geom_segment2, interactive_geom = GeomInteractiveSegmentGGtree, ...)
}


#' @export
geom_text2_interactive <- function(...){
  layer_interactive(geom_text2, interactive_geom = GeomInteractiveTextGGtree, ...)
}

#' @export
geom_label2_interactive <- function(...){
  layer_interactive(geom_label2, interactive_geom = GeomInteractiveLabelGGtree, ...)
}

#' @export 
geom_point2_interactive <- function(...){
  layer_interactive(geom_point2, interactive_geom = GeomInteractivePointGGtree, ...)
}

#' @export
geom_hilight_rect2_interactive <- function(...){
  layer_interactive(geom_hilight_rect2, interactive_geom = GeomInteractiveHilightRect, ...)
}

#' @export
geom_hilight_encircle2_interactive <- function(...){
  layer_interactive(geom_hilight_encircle2, interactive_geom = GeomInteractiveHilightEncircle, ...)
}

#' @export
geom_curvelink_interactive <- function(...){
  layer_interactive(geom_curvelink, interactive_geom = GeomInteractiveCurvelink, ...)
}

#' @export
geom_shadowtext_interactive <- function(...){
  rlang::check_installed(c('ggiraph', 'shadowtext'), "for `geom_shadowtext_interactive()`.")
  layer_interactive(geom_shadowtext, interactive_geom = GeomInteractiveShadowtext,...)
}


#' @importFrom ggimage geom_image
#' @export
geom_image_interactive <- function(...){
  rlang::check_installed(c('ggiraph', 'ggimage'), "for `geom_image_interactive()`.")
  layer_interactive(geom_image, interactive_geom = GeomInteractiveImage,...)
}

#' @importFrom ggimage geom_phylopic
#' @export
geom_phylopic_interactive <- function(...){
  rlang::check_installed(c('ggiraph', 'ggimage'), "for `geom_phylopic_interactive()`.")
  layer_interactive(geom_phylopic, interactive_geom = GeomInteractiveImage,...)
}


# the internal functions of ggiraph
layer_interactive <- getFromNamespace("layer_interactive", "ggiraph")
add_default_interactive_aes <- getFromNamespace("add_default_interactive_aes", "ggiraph")
interactive_geom_parameters <- getFromNamespace("interactive_geom_parameters", "ggiraph")
interactive_geom_draw_key <- getFromNamespace("interactive_geom_draw_key", "ggiraph")
IPAR_NAMES <- getFromNamespace("IPAR_NAMES", "ggiraph")
add_interactive_attrs <- getFromNamespace("add_interactive_attrs", "ggiraph")

GeomSegmentGGtree <- getFromNamespace("GeomSegmentGGtree", "ggtree")
GeomTextGGtree <- getFromNamespace("GeomTextGGtree", "ggtree")
GeomPointGGtree <- getFromNamespace("GeomPointGGtree", "ggtree")
GeomLabelGGtree <- getFromNamespace("GeomLabelGGtree", "ggtree")
GeomCurvelink <- getFromNamespace("GeomCurvelink", "ggtree")
GeomHilightRect <- getFromNamespace("GeomHilightRect", "ggtree")
GeomHilightEncircle <- getFromNamespace("GeomHilightEncircle", "ggtree")

geom_curvelink <- getFromNamespace("geom_curvelink", "ggtree")
geom_hilight_rect2 <- getFromNamespace("geom_hilight_rect2", "ggtree")
geom_hilight_encircle2 <- getFromNamespace("geom_hilight_encircle2", "ggtree")


GeomImage <- getFromNamespace("GeomImage", "ggimage")

#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @importFrom ggiraph GeomInteractiveSegment
#' @importFrom ggplot2 ggproto
#' @export
GeomInteractiveSegmentGGtree <- ggproto(
  "GeomInteractiveSegmentGGtree",
  GeomSegmentGGtree,
  default_aes = add_default_interactive_aes(GeomSegmentGGtree),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(data, ..., nudge_x = 0, .ipar = IPAR_NAMES){
      data$x <- data$x + nudge_x
      GeomInteractiveSegment$draw_panel(data, ..., .ipar = .ipar)
  }
)


#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @importFrom ggiraph GeomInteractiveText
#' @importFrom ggplot2 ggproto
#' @export
GeomInteractiveTextGGtree <- ggproto(
  "GeomInteractiveTextGGtree",
  GeomTextGGtree,
  default_aes = add_default_interactive_aes(GeomTextGGtree),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(..., .ipar = IPAR_NAMES){
    GeomInteractiveText$draw_panel(..., .ipar = .ipar)
  }
)

#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @importFrom ggiraph GeomInteractivePoint
#' @importFrom ggplot2 ggproto
#' @export
GeomInteractivePointGGtree <- ggproto(
  "GeomInteractivePointGGtree",
  GeomPointGGtree,
  default_aes = add_default_interactive_aes(GeomPointGGtree),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(..., .ipar = IPAR_NAMES){
    GeomInteractivePoint$draw_panel(..., .ipar = IPAR_NAMES)
  }
)


#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @importFrom ggiraph GeomInteractiveLabel
#' @importFrom ggplot2 ggproto
#' @export
GeomInteractiveLabelGGtree <- ggproto(
  "GeomInteractiveLabelGGtree",
  GeomLabelGGtree,
  default_aes = add_default_interactive_aes(GeomLabelGGtree),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(..., .ipar = IPAR_NAMES){
    GeomInteractiveLabel$draw_panel(..., .ipar = .ipar)
  }
)

#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @importFrom shadowtext GeomShadowtext
#' @importFrom ggplot2 ggproto
#' @export
GeomInteractiveShadowtext <- ggproto(
  "GeomInteractiveShadowtext",
  GeomShadowtext,
  default_aes = add_default_interactive_aes(GeomShadowtext),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(data, panel_params, coord, ..., .ipar = IPAR_NAMES){
    gr <- GeomShadowtext$draw_panel(data, panel_params, coord, ...)
    coords <- coord$transform(data, panel_params)
    gr$children[[1]] <- add_interactive_attrs(gr$children[[1]], coords, ipar=.ipar)
    gr
  }
)


#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @importFrom ggplot2 ggproto
#' @export
GeomInteractiveImage <- ggproto(
  "GeomInteractiveImage",
  GeomImage,
  default_aes = add_default_interactive_aes(GeomImage),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(data, panel_params, coord, ..., .ipar = IPAR_NAMES){
    gr <- GeomImage$draw_panel(data, panel_params, coord, ...)
    coords <- coord$transform(data, panel_params)
    for (i in seq_along(gr$children)){
       gr$children[[i]] <- add_interactive_attrs(gr$children[[i]], coords[i,], ipar=.ipar)
    }
    gr
  }
)

#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @importFrom ggplot2 ggproto
#' @importFrom grid gTree gpar curveGrob
#' @export
GeomInteractiveCurvelink <- ggproto(
  "GeomInteractiveCurvelink",
  GeomCurvelink,
  default_aes = add_default_interactive_aes(GeomCurvelink),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(data, panel_params, coord, shape=0.5, outward=TRUE,
                        arrow = NULL, arrow.fill=NULL, lineend = "butt", 
                        na.rm = FALSE, .ipar = IPAR_NAMES) {
    if (!coord$is_linear()) {
        tmpgroup <- data$group
        starts <- subset(data, select = c(-xend, -yend))
        starts$group <- 1
        ends <- rename(subset(data, select = c(-x, -y)), c("x" = "xend", "y" = "yend"))
        ends$group <- 2
        pieces <- rbind(starts, ends)

        trans <- coord$transform(pieces, panel_params)
        starts <- trans[trans$group==1, ,drop=FALSE]
        ends <- trans[trans$group==2, ,drop=FALSE]
        if (outward){
            curvature <- unlist(mapply(generate_curvature2, starttheta=starts$theta,
                                       endtheta=ends$theta, hratio=starts$hratio, ncp=starts$ncp,
                                       SIMPLIFY=FALSE))
        }else{
            curvature <- unlist(mapply(generate_curvature, starttheta=starts$theta,
                                       endtheta=ends$theta, hratio=starts$hratio, ncp=starts$ncp,
                                       SIMPLIFY=FALSE))
        }
        ends <- rename(subset(ends, select=c(x, y)), c("xend"="x", "yend"="y"))
        trans <- cbind(starts, ends)
        trans$group <- tmpgroup
        trans$curvature <- curvature
    }else{
        trans <- coord$transform(data, panel_params)
        if (inherits(coord, 'CoordFlip')){
            trans$curvature <- -1 * trans$curvature
        }
    }
    arrow.fill <- arrow.fill %|||% trans$colour

    grobs <- lapply(seq_len(nrow(trans)), function(i){
                        subgrob <- curveGrob(
                              trans$x[i], trans$y[i], trans$xend[i], trans$yend[i],
                              default.units = "native",
                              curvature = trans$curvature[i], angle = trans$curveangle[i], ncp = trans$ncp[i],
                              square = trans$square[i], squareShape = 1, inflect = FALSE, open = TRUE,
                              gp = gpar(col = alpha(trans$colour[i], trans$alpha[i]),
                                        fill = alpha(arrow.fill[i], trans$alpha[i]),
                                        lwd = trans$linewidth[i] * ggplot2::.pt,
                                        lty = trans$linetype[i],
                                        lineend = lineend),
                              arrow = arrow,
                              shape = shape)
                        add_interactive_attrs(subgrob, trans[i,], ipar = .ipar)
                        })
    class(grobs) <- "gList"
    return(ggname("geom_curvelink_interactive", gTree(children=grobs))) 
  }
)

#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @importFrom ggplot2 ggproto
#' @importFrom grid gTree gpar curveGrob
#' @export
GeomInteractiveHilightRect <- ggproto(
  "GeomInteractiveHilightRect",
  GeomHilightRect,
  default_aes = add_default_interactive_aes(GeomHilightRect),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(data, panel_params, coord, ..., .ipar = IPAR_NAMES){
     gr <- GeomHilightRect$draw_panel(data, panel_params, coord, ...)
     coords <- coord$transform(data, panel_params)
     gr <- add_interactive_attrs(gr, coords, ipar=.ipar)
     gr
  } 
)


#' @title ggproto classes for ggiraph
#' @description
#' ggproto classes for ggiraph
#' @format NULL
#' @usage NULL
#' @export
GeomInteractiveHilightEncircle <- ggproto(
  "GeomInteractiveHilightEncircle",
  GeomHilightEncircle,
  default_aes = add_default_interactive_aes(GeomHilightEncircle),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(data, panel_params, coord, ..., .ipar = IPAR_NAMES){
      gr <- GeomHilightEncircle$draw_panel(data, panel_params, coord, ...)
      coords <- coord$transform(data, panel_params)
      index <- coords[,!names(coords) %in% c('x', 'y')]  |> duplicated() 
      coords <- coords[!index, , drop=FALSE]
      for (i in seq_along(gr$children)){
        gr$children[[i]] <- add_interactive_attrs(gr$children[[i]], coords[i,], ipar=.ipar)
      }
      gr
  }
)


