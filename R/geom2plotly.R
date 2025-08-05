#' @importFrom plotly to_basic
#' @method to_basic GeomSegmentGGtree
#' @export
to_basic.GeomSegmentGGtree <- getFromNamespace("to_basic.GeomSegment", asNamespace("plotly"))


#' @method to_basic GeomPointGGtree
#' @export
to_basic.GeomPointGGtree <- function(data, prestats_data, layout, params, p, ...){
   prefix_class(data, "GeomPoint")
}

prefix_class <- function(x, y) {
  structure(x, class = unique(c(y, class(x))))
}

## #' @importFrom plotly geom2trace
## #' @method geom2trace GeomPointGGtree
## #' @export
## geom2trace.GeomPointGGtree <- getFromNamespace("geom2trace.GeomPoint", "plotly")


#' @method to_basic GeomTextGGtree
#' @export
to_basic.GeomTextGGtree <- function(data, prestats_data, layout, params, p, ...){
   prefix_class(data, "GeomText")
}

#' @method to_basic GeomHilightRect
#' @export
to_basic.GeomHilightRect <- getFromNamespace("to_basic.GeomRect", "plotly")


#' @method to_basic GeomHilightEncircle
#' @export
to_basic.GeomHilightEncircle <- getFromNamespace("to_basic.GeomRect", "plotly")

