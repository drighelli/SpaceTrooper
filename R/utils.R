#' .getActiveGeometryName
#'
#' @param sf
#'
#' @return
#' @export
#'
#' @examples
.getActiveGeometryName <- function(sf)
{
    stopifnot(is(sf, "sf"))
    cn = attr(sf, "sf_column");
    return(cn);
}
