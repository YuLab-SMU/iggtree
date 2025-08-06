#' @importFrom grid grobName
ggname <- function (prefix, grob){
  grob$name <- grobName(grob, prefix)
  grob
}

`%|||%` <- function(x, y){
  if (is.null(x)) {
      return(y)
  }
  if (is.null(y)) {
      return(x)
  }
  if (length(x) < length(y)) {
      return(y)
  }
  else {
      return(x)
  }
}
