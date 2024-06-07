#' readPolygonsCosMx
#'
#' @param polygonsfpattern
#' @param x
#' @param y
#' @param xloc
#' @param yloc
#' @param micronConvFact
#'
#' @return
#' @export
#' @importFrom data.table fread
#'
#' @examples
readPolygonsCosMx <- function(polygonsFile, x="x_global_px", y="y_global_px",
                            xloc="x_local_px", yloc="y_local_px",
                            micronConvFact=0.12)
{
    spat_obj <- fread(polygonsFile)
    spat_obj$cell_id <- paste0("f", spat_obj$fov, "_c",
                               spat_obj$cellID)
    spat_obj$cell_id <- as.factor(spat_obj$cell_id)

    polygons_glo <- .createPolygons(spat_obj, x=x, y=y, polygon_id="cell_id")[,-c(4,5)]
    polygons_loc <- .createPolygons(spat_obj, x=xloc, y=yloc, polygon_id="cell_id")[,-c(4,5)]

    polygons <- cbind(polygons_glo, polygons_loc$geometry)

    st_geometry(polygons) <- "geometry.1"
    st_geometry(polygons) <- "local"
    polygons <- .checkPolygonsValidity(polygons)
    st_geometry(polygons) <- "global"
    polygons <- .checkPolygonsValidity(polygons)
    if (sum(rownames(polygons) == colnames(spe)) != dim(spe)[2])
    {
        cd <- colData(spe)
        polygons <- merge(polygons, cd[,"cell_id"])
        cd <- merge(cd, polygons[,"cell_id"])
        spe <- spe[, rownames(cd)]
    }

    ### write polygons as parquet file

    ###############
    # um_area <- st_area(polygons)*(micronConvFact^2) # working on the global geometry
    # spe$um_area <- um_area
    # polygons$um_area <- um_area
    # aspect_ratio_list <- lapply(spe$polygons_whole_custom$geometry, function(x){
    #     (max(x[[1]][,2]) - min(x[[1]][,2]))/(max(x[[1]][,1]) - min(x[[1]][,1]))
    # })
    # spe$polygons_whole_custom$log2_AspectRatio <- log2(unlist(aspect_ratio_list))
    # spe$log2_AspectRatio <- log2(unlist(aspect_ratio_list))
    ################
}

#' .createPolygons
#'
#' @param spat_obj
#' @param x
#' @param y
#' @param polygon_id
#'
#' @return
#' @keywords internal
#' @importFrom sfheaders sf_polygon
#'
#' @examples
.createPolygons <- function(spat_obj, x, y, polygon_id)
{
    polygons <- sfheaders::sf_polygon(obj=spat_obj,
                                      x=x, y=y,
                                      polygon_id=polygon_id, keep=TRUE)
    polygons <- polygons[base::order(polygons$fov, polygons$cellID),]
    rownames(polygons) <- polygons$cell_id
    return(polygons)
}

#' Title
#' @description checks validity on the active geometry of `sf` parameter
#' @param sf
#'
#' @return
#' @keywords internal
#' @importFrom sf st_is_valid st_buffer
#'
#' @examples
.checkPolygonsValidity <- function(sf)
{
    stopifnot(is(sf, "sf"))
    sf_tf <- st_is_valid(sf) # to parallelize? how? split sf in multiple sf and parallelize on it?

    if(sum(sf_tf)!=dim(sf)[1]) sf <- st_buffer(sf, dist=0)
    gmn <- .getActiveGeometryName(sf)
    #    Subsetting to remove non-polygons
    cellids <- unlist(apply(sf, 1, function(geom)
    {
        if(attr(geom[[gmn]], "class")[2] == "MULTIPOLYGON") {
            return(geom$cell_id)
        }
    }))

    if(length(cellids)!=0) sf <- sf[sf$cell_id!=cellids,]
    return(sf)
}


