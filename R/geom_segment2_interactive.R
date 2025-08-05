#' @export
geom_segment2_interactive <- function(...){
  rlang::check_installed('ggiraph', "for `geom_segment2_interactive()`.")
  layer_interactive(geom_segment2, ...)
}

# the internal functions of ggiraph
layer_interactive <- getFromNamespace("layer_interactive", "ggiraph")
add_default_interactive_aes <- getFromNamespace("add_default_interactive_aes", "ggiraph")
interactive_geom_parameters <- getFromNamespace("interactive_geom_parameters", "ggiraph")
interactive_geom_draw_key <- getFromNamespace("interactive_geom_draw_key", "ggiraph")
IPAR_NAMES <- getFromNamespace("IPAR_NAMES", "ggiraph")
add_interactive_attrs <- getFromNamespace("add_interactive_attrs", "ggiraph")

GeomSegmentGGtree <- getFromNamespace("GeomSegmentGGtree", "ggtree")

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
