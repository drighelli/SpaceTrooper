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

#' getFencesOutlier
#' Retrieve Threshold (Fence) Values from a SpatialExperiment Object
#'
#' This function extracts the threshold values, also known as fences,
#' from a specified column in the `colData` of a `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param fences_of A character string specifying the name of the column in
#' `colData(spe)` from which to extract the fence values. This column should
#' contain an `outlier.filter` object.
#' @param decimal.round An optional integer specifying the number of decimal
#' places to which the fence values should be rounded. If `NULL`, no rounding is
#'  applied. Default is `NULL`.
#'
#' @return A numeric vector containing the lower and upper threshold values
#' extracted from the specified column.
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' #TBD
getFencesOutlier <- function(spe, fences_of,
                    high_low=c("both", "lower", "higher"), decimal_round=NULL)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(fences_of %in% names(colData(spe)))
    stopifnot(is(colData(spe)[[fences_of]], "outlier.filter"))
    high_low <- match.arg(high_low)
    fences <- attr(colData(spe)[[fences_of]], "thresholds")
    if(!is.null(decimal_round)) {fences <- round(fences, decimal_round)}
    fences <- switch (high_low,
        both = {fences},
        lower = {fences[1]},
        higher = {fences[2]}
    )
    return(fences)
}
