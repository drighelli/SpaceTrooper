#' .getActiveGeometryName
#'
#' @param sf an sf object
#'
#' @return character with the name of the active geometry
#' @export
#'
#' @examples
.getActiveGeometryName <- function(sf)
{
    stopifnot(is(sf, "sf"))
    cn = attr(sf, "sf_column")
    return(cn)
}

#' .setActiveGeometry
#'
#' @param sf an sf object
#' @param name character for the geometry to activate
#'
#' @return an sf object
#' @export
#'
#' @examples
.setActiveGeometry <- function(sf, name)
{
    stopifnot(is(sf, "sf"), name %in% names(sf))
    sf::st_geometry(sf) <- name
    return(sf)
}

#' .renameGeometry
#'
#' @description renames the `from` to `to` geometry of the `sf` object.
#' If `activate` is `TRUE` it set as the active geometry the new geometry name.
#' Default behaviour is to check if the renamed geometry is already active and
#' leave it as active with the new name.
#'
#' @param sf an sf object with the `from` geometry
#' @param from character indicating the name of the geometry to change
#' @param to character indicating the new name of the geometry
#' @param activate logical indicating if the renamed geometry has to be activated
#'
#' @return an sf object
#' @export
#'
#' @examples
.renameGeometry <- function(sf, from, to, activate=FALSE)
{
    stopifnot(all(is(sf, "sf"), from %in% names(sf)))
    act <- .getActiveGeometryName(sf)

    names(sf)[which(names(sf)==from)] <- to
    if(from==act) { sf <- .setActiveGeometry(sf, to) }
    if(activate) { sf <- .setActiveGeometry(sf, to) }
    return(sf)
}
