#' readPolygons
#'
#' @param polygonsFile
#' @param x
#' @param y
#' @param xloc
#' @param yloc
#' @param micronConvFact
#' @param keepMultiPol
#' @param verbose
#'
#' @return an sf object with the loaded and validated polygons
#'
#' @export
#' @importFrom data.table fread
#' @importFrom arrow read_parquet
#' @importFrom sf st_geometry
#'
#' @examples
# polygons <- readPolygonsCosMx("~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/Run5810_Case2-polygons.csv")
readPolygons <- function(polygonsFile, type=c("csv", "parquet"),
                            x=c("x_global_px", "vertex_x"),
                            y=c("y_global_px", "vertex_y"),
                            xloc="x_local_px", yloc="y_local_px",
                            #micronConvFact=0.12,
                            keepMultiPol=TRUE,
                            verbose=FALSE)
{

    stopifnot(file.exists(polygonsFile))
    type <- match.arg(type)
    # type <- grep("csv", polygonsFile)

    spat_obj <- switch(type, csv=fread(polygonsFile),
                            parquet=read_parquet(polygonsFile))
    if (! "cell_id" %in% colnames(spat_obj))
    {
        spat_obj$cell_id <- paste0("f", spat_obj$fov, "_c", spat_obj$cellID)
    }

    spat_obj$cell_id <- as.factor(spat_obj$cell_id)

    polygons <- .createPolygons(spat_obj, x=x, y=y,
                                    polygon_id="cell_id")
    polygons <- .renameGeometry(polygons, "geometry", "global")

    if(all(c(xloc, yloc) %in% colnames(spat_obj)))
    {
        idxs <- which(colnames(polygons) %in% c(xloc, yloc))
        polygons <- polygons[,-idxs]
        polygons_loc <- .createPolygons(spat_obj, x=xloc, y=yloc,
                                        polygon_id="cell_id")[,-c(4,5)]
        polygons <- cbind(polygons, polygons_loc$geometry)
        polygons <- .renameGeometry(polygons, "geometry.1", "local")
    }

    if(verbose) message("Polygons detected: ", dim(polygons)[1])#### otherwise
    #### will print number of columns next to the numer of rows

    polygons <- .checkPolygonsValidity(polygons, keepMultiPol=keepMultiPol,
                                        verbose=verbose)

    rownames(polygons) <- polygons$cell_id #### if not here, then the check
    #### in addPolygonsToSPE cannot be be done

    #### identical cannot work since coordinates are different, but the validity
    #### and the geometries as well should stay the same
    #### needs to be changed or ignored
    # if(!table(st_is_valid(polygons$global))==table(st_is_valid(polygons$global)))
    #     warning("Global and Local geometries are not identical")
    if(verbose) message("Polygons after validity: ", dim(polygons)[1])#### otherwise
    #### will print number of columns next to the numer of rows
    return(polygons)
    ### write polygons as parquet file

    ############### custom metrics computed on custom read polygons -> to implement in separate function(s)
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


#' .checkPolygonsValidity
#' @description checks validity on a geometry of `sf` object.
#' It removes multipolygons when `keepMultiPol` is `FALSE`
#'
#' @details
#' In case geometry is NULL validity is checked on the active geometry,
#' otherwise it is checked on the passed geometry without changing the active
#' geometry of the sf object.
#' In case of not valid polygons, these are removed.
#' If keeMultiPol is FALSE, possible detected multipolygons are removed.
#'
#'
#' @param sf an sf class object
#' @param geometry character for the geometry to check validity, if `NULL`
#' it checks the active geometry (default is `NULL`)
#' @param keepMultiPol logical for keeping/removing moltipolygons, if any
#' (default is `TRUE`, so keeping the multipolygons)
#' @param verbose logical to print verbose output (default is `FALSE`)
#'
#' @return
#' @keywords internal
#' @importFrom sf st_is_valid st_buffer
#'
#' @examples
.checkPolygonsValidity <- function(sf, geometry=NULL, keepMultiPol=TRUE,
                                    verbose=FALSE)
{
    stopifnot(is(sf, "sf"))
    if(is.null(geometry))
    {
        geometry <- .getActiveGeometryName(sf)
    } else {
        act <- .getActiveGeometryName(sf)
        sf <- .setActiveGeometry(sf, geometry)
    }
    sf_tf <- st_is_valid(sf) # to parallelize? how? split sf in multiple sf and parallelize on it?

    if(sum(sf_tf)!=dim(sf)[1]) sf <- st_buffer(sf, dist=0)

    # Subsetting to remove non-polygons
    cellids <- unlist(apply(sf, 1, function(geom)
    {
        # just in case of custom in cosmx but present in automatic in xenium
        if(attr(geom[[geometry]], "class")[2] == "MULTIPOLYGON")
        {
            return(geom$cell_id)
        }
    }))

    ## print an alert on the detection/removal of multipolygons
    if(length(cellids)!=0)
    {
        if(verbose) message("Detected ", length(cellids), " multipolygons.")
        sf$is_multi <- sf$cell_id %in% cellids
    }
    if(!keepMultiPol)
    {
        if(verbose) message("Removing ", sum(sf$is_multi), " multipolygons.")
        sf <- sf[!sf$is_multi,]
    }

    if(exists("act")) sf <- .setActiveGeometry(sf, act)

    return(sf)
}


#' addPolygonsToSPE
#'
#' @param spe
#' @param polygons
#'
#' @return
#' @importFrom dplyr left_join
#' @export
#'
#' @examples
addPolygonsToSPE <- function(spe, polygons)
{
    stopifnot(all(is(spe, "SpatialExperiment"), is(polygons, "sf")))

    if (sum(rownames(polygons) == colnames(spe)) != dim(spe)[2])
    {
        cd <- data.frame(colData(spe))
        polygons <- left_join(polygons, cd[, c("fov", "cellID")],
                              by=c("fov", "cellID"))
        cd <- left_join(cd, polygons[ , c("fov", "cellID")],
                        by=c("fov","cellID"))
        rownames(cd) <- cd$cell_id
        rownames(polygons) <- polygons$cell_id
        spe <- spe[, spe$cell_id %in% rownames(polygons)]
        polygons <- polygons[polygons$cell_id %in% rownames(cd),]
    }
    spe <- spe[, rownames(polygons)]
    colData(spe)$polygons <- polygons
    return(spe)
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
    # polygons <- polygons[order(polygons$cell_id),]
    return(polygons)
}

#' Title
#'
#' @param polygonsFile
#' @param type
#' @param x
#' @param y
#' @param keepMultiPol
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
readPolygonsXenium <- function(polygonsFile, type=c("parquet", "csv"),
                   x="vertex_x", y="vertex_y", keepMultiPol=TRUE,
                   verbose=FALSE)
{
    type <- match.arg(type)
    readPolygons(polygonsFile=polygonsFile, type=type, x=x, y=y,
        keepMultiPol=keepMultiPol, verbose=verbose)
}




