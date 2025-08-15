#' @title Create interactive line segments of ggtree  
#' @description
#' The geometry is based on `geom_segment2()`.
#' See the documentation for those functions for more details.
#' @param ... see also the parameters of `geom_segment2()` of `ggtree`
#' @export
geom_segment2_interactive <- function(...){
  rlang::check_installed('ggiraph', "for `geom_segment2_interactive()`.")
  layer_interactive(geom_segment2, interactive_geom = GeomInteractiveSegmentGGtree, ...)
}

#' @title Create interactive text of ggtree
#' @description
#' The geometry is based on `geom_text2()`.
#' See the documentation for those functions for more details.
#' @param ... see also the parameters of `geom_text2()` of `ggtree`
#' @export
geom_text2_interactive <- function(...){
  layer_interactive(geom_text2, interactive_geom = GeomInteractiveTextGGtree, ...)
}


#' @title Create interactive label of ggtree
#' @description
#' The geometry is based on `geom_label2()`.
#' See the documentation for those functions for more details.
#' @param ... see also the parameters of `geom_label2()` of `ggtree`
#' @export
geom_label2_interactive <- function(...){
  layer_interactive(geom_label2, interactive_geom = GeomInteractiveLabelGGtree, ...)
}

#' @title Create interactive point of ggtree
#' @description
#' The geometry is based on `geom_point2()`.
#' See the documentation for those functions for more details.
#' @param ... see also the parameters of `geom_point2()` of `ggtree`
#' @export 
geom_point2_interactive <- function(...){
  layer_interactive(geom_point2, interactive_geom = GeomInteractivePointGGtree, ...)
}


#' @keywords internal
geom_hilight_rect2_interactive <- function(...){
  layer_interactive(geom_hilight_rect2, interactive_geom = GeomInteractiveHilightRect, ...)
}

#' @keywords internal
geom_hilight_encircle2_interactive <- function(...){
  layer_interactive(geom_hilight_encircle2, interactive_geom = GeomInteractiveHilightEncircle, ...)
}



#' @keywords internal
geom_curvelink_interactive <- function(...){
  layer_interactive(geom_curvelink, interactive_geom = GeomInteractiveCurvelink, ...)
}

#' @title Create interactive shadow text
#' @description
#' The geometry is based on `geom_shadowtext()`.
#' See the documentation for those functions for more details.
#' @param ... see also the parameters of `geom_shadowtext()` of `shadowtext`
#' @export
geom_shadowtext_interactive <- function(...){
  rlang::check_installed(c('ggiraph', 'shadowtext'), "for `geom_shadowtext_interactive()`.")
  layer_interactive(shadowtext::geom_shadowtext, interactive_geom = GeomInteractiveShadowtext,...)
}


#' @title Create interactive image of ggimage 
#' @description
#' The geometry is based on `geom_image()`.
#' See the documentation for those functions for more details.
#' @param ... see also the parameters of `geom_image()` of `ggimage`
#' @export
geom_image_interactive <- function(...){
  rlang::check_installed(c('ggiraph', 'ggimage'), "for `geom_image_interactive()`.")
  layer_interactive(ggimage::geom_image, interactive_geom = GeomInteractiveImage,...)
}


#' @title Create interactive phylopic of ggimage
#' @description
#' The geometry is based on `geom_phylopic()`.
#' See the documentation for those functions for more details.
#' @param ... see also the parameters of `geom_phylopic()` of `ggimage`
#' @export
geom_phylopic_interactive <- function(...){
  rlang::check_installed(c('ggiraph', 'ggimage'), "for `geom_phylopic_interactive()`.")
  layer_interactive(ggimage::geom_phylopic, interactive_geom = GeomInteractiveImage,...)
}

#' @importFrom ggplot2 aes_
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


geom_segment2 <- getFromNamespace("geom_segment2", "ggtree")
geom_text2 <- getFromNamespace("geom_text2", "ggtree")
geom_label2 <- getFromNamespace("geom_label2", "ggtree")
geom_point2 <- getFromNamespace("geom_point2", "ggtree")
geom_curvelink <- getFromNamespace("geom_curvelink", "ggtree")
geom_hilight_rect2 <- getFromNamespace("geom_hilight_rect2", "ggtree")
geom_hilight_encircle2 <- getFromNamespace("geom_hilight_encircle2", "ggtree")


generate_curvature <- getFromNamespace("generate_curvature", "ggtree")
generate_curvature2 <- getFromNamespace("generate_curvature2", "ggtree")
build_align_data <- getFromNamespace("build_align_data", "ggtree")
rect_to_poly <- getFromNamespace("rect_to_poly", "ggtree")

GeomImage <- getFromNamespace("GeomImage", "ggimage")

GeomShadowtext <- getFromNamespace("GeomShadowtext", "shadowtext")

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
#' @importFrom grid gTree gpar curveGrob grobTree
#' @importFrom cli cli_alert_warning
#' @importFrom ggiraph GeomInteractivePolygon
#' @export
GeomInteractiveHilightRect <- ggproto(
  "GeomInteractiveHilightRect",
  GeomHilightRect,
  default_aes = add_default_interactive_aes(GeomHilightRect),
  parameters = interactive_geom_parameters,
  draw_key = interactive_geom_draw_key,
  draw_panel = function(data, panel_params, coord, align='none', 
                        gradient = FALSE, roundrect = FALSE, ..., 
                        .ipar = IPAR_NAMES){
     if (coord$is_linear()){
        gr <- GeomHilightRect$draw_panel(data, panel_params, coord, gradient = gradient, roundrect = roundrect, align = align, ...)
        coords <- coord$transform(data, panel_params)
        gr <- add_interactive_attrs(gr, coords, ipar=.ipar)
     }else{
        data$xmax <- data$xmax + data$extend
        if (!any(is.null(data$extendto)) && !any(is.na(data$extendto))){
            # check whether the x of tree is reversed.
            flag1 <- data$xmin < data$xmax
            # check whether extendto is more than xmax
            flag2 <- data$extendto < data$xmax
            flag <- flag1 == flag2
            if (all(flag1) && any(flag)){
                cli_alert_warning(c("{.code extendto} ", paste0(data$extendto[flag], collapse="; "),
                             ifelse(length(data$extendto[flag])>1, " are", " is")," too small for node: ",
                             paste0(data$clade_root_node[flag], collapse="; "),", keep the original xmax value(s): ",
                             paste0(data$xmax[flag], collapse="; "), "."), wrap = TRUE)
                data$xmax[!flag] <- data$extendto[!flag]
            }else if(!all(flag1) && any(flag)){
                cli_alert_warning(c("{.code extendto} ", paste0(data$extendto[flag], collapse="; "),
                             ifelse(length(data$extendto[flag])>1, " are", " is"), " too big for node: ",
                             paste0(data$clade_root_node[flag], collapse="; "), ", keep the original xmax value(s): ",
                             paste0(data$xmax[flag], collapse="; "), "."), wrap = TRUE)
                data$xmax[!flag] <- data$extendto[!flag]
            }else{
                data$xmax <- data$extendto
            }
        }
        data <- build_align_data(data=data, align=align)
        if (gradient){
            cli_alert_warning("The gradient color hight light layer only presents in
                              rectangular, ellipse, roundrect layouts.", wrap = TRUE)
        }
        if (roundrect){
            cli_alert_warning("The round rectangular hight light layer only presents in
                              rectangular, ellipse, roundrect layouts.", wrap =TRUE)
        }
        aesthetics <- setdiff(colnames(data), c("xmin", "xmax", "ymin", "ymax", "clade_root_node")) 
        gr <- lapply(split(data, seq_len(nrow(data))), function(row) {
                      poly <- rect_to_poly(row$xmin, row$xmax, row$ymin, row$ymax)
                      aes <- row[rep(1,5), aesthetics]
                      GeomInteractivePolygon$draw_panel(vctrs::vec_cbind(poly, aes), panel_params, coord, ..., .ipar=.ipar)
                      })
        gr <- ggname("geom_hilight_rect2_interactive", do.call("grobTree", gr)) 
     }
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


