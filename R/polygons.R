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
readPolygons <- function(polygonsFile, type=c("csv", "parquet", "h5"),
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

    if(type=="h5")
    {
        # polfiles <- list.files(polygonsFolder, pattern=hdf5pattern,
        #                        full.names=TRUE)
        # dfsfl <- lapply(seq_along(polfiles), function(i)
        # {
        #     poll <- readh5polygons(pol_file=polfiles[i])
        #     df <- data.frame(cell_id=paste0("f", i-1, "_c", poll$ids),
        #                      cell_ID=poll$ids,
        #                      fov=i-1, geometry=sf::st_sfc(poll$g))
        #     dfsf <- sf::st_sf(df)
        # })
        # polygons <- do.call(rbind, dfsfl)
    } else {
        spat_obj <- switch(type, csv=fread(polygonsFile),
                                parquet=read_parquet(polygonsFile))
        if (! "cell_id" %in% colnames(spat_obj))
        {
            spat_obj$cell_id <- paste0("f", spat_obj$fov, "_c", spat_obj$cellID)
        }

        spat_obj$cell_id <- as.factor(spat_obj$cell_id)

        polygons <- .createPolygons(spat_obj, x=x, y=y,
                                        polygon_id="cell_id")
        polygons <- .renameGeometry(polygons, "geometry", "global") ## only for cosmx

        if(all(c(xloc, yloc) %in% colnames(spat_obj)))
        {
            idxs <- which(colnames(polygons) %in% c(xloc, yloc))
            polygons <- polygons[,-idxs]
            polygons_loc <- .createPolygons(spat_obj, x=xloc, y=yloc,
                                            polygon_id="cell_id")[,-c(4,5)]
            polygons <- cbind(polygons, polygons_loc$geometry)
            polygons <- .renameGeometry(polygons, "geometry", "local")
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
    }
    ### write polygons as parquet file
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
#' @importFrom sf st_is_valid st_buffer st_geometry_type
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

    ############ CHECKING MULTIPOLYGONS
    # Subsetting to remove non-polygons
    # cellids <- unlist(apply(sf, 1, function(geom)
    # {

        # just in case of custom in cosmx but present in automatic in xenium
        # if(attr(geom[[geometry]], "class")[2] == "MULTIPOLYGON")
        # {
        #     return(geom$cell_id)
        # }
    # }))
    is_multi <- sf::st_geometry_type(sf$global) == "MULTIPOLYGON"
    ## TRUE is merscope - FALSE is cosmx and xenium
    merscopeFl <- (sum(is_multi) == dim(sf)[1])
    funct <- ifelse(merscopeFl, "lengths", "length")
    ll <- lapply(sf[[geometry]], funct)
    sums <- lapply(ll, sum)
    idx <- which(sums!=1)
    sf$is_multi <- FALSE
    sf$multi_n <- 1
    if(length(idx)!=0)
    {
        sf$is_multi[idx] <- TRUE
        sf$multi_n[idx] <- unlist(sums[idx])

    }
    if( all(merscopeFl, length(idx)!=0) )
    {
        sf$global[-idx] <- st_cast(sf$global[-idx], "POLYGON")
    }

    if(verbose) message("Detected ", sum(sf$is_multi), " multipolygons.")

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
#' @importFrom sf st_as_sf
#'
#' @examples
.createPolygons <- function(spat_obj, x=NULL, y=NULL, polygon_id=NULL, geometry="Geometry")
{
    if(all(!is.null(x), !is.null(y)))
    {
        polygons <- sfheaders::sf_polygon(obj=spat_obj,
                                      x=x, y=y,
                                      polygon_id=polygon_id, keep=TRUE)
    } else {
        polygons <- sf::st_as_sf(as.data.frame(spat_obj))
        ## check structure of polygons
        # ll <- lapply(polygons$Geometry, lengths)
        # sums <- lapply(ll, sum)
        # idx <- which(sums==1)
        # polygons$is_multi <- FALSE
        # polygons$is_multi[!idx] <- TRUE
        # polygons$multi_n <- 1
        # polygons$multi_n[!idx] <- sums[!idx]
        # polygons[idx, ] <- sf::st_cast(polygons[idx, ], "POLYGON")
    }
    # polygons <- polygons[order(polygons$cell_id),]
    return(polygons)
}

#' readPolygonsCosmx
#' @description
#'
#' @param polygonsFile
#' @param type
#' @param x
#' @param y
#' @param xloc
#' @param yloc
#' @param keepMultiPol
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
readPolygonsCosmx <- function(polygonsFile, type=c("csv", "parquet"),
                              x="x_global_px",
                              y="y_global_px",
                              xloc="x_local_px",
                              yloc="y_local_px",
                              keepMultiPol=TRUE,
                              verbose=FALSE)
{
    type=match.arg(type)
    polygons <- readPolygons(polygonsFile, type=type, x=x, y=y, xloc=xloc,
                    yloc=yloc, keepMultiPol=keepMultiPol, verbose=verbose)
    polygons <- st_cast(polygons, "GEOMETRY")
    mandatory <- c("cell_id", "global", "is_multi", "multi_n")
    cnames <- colnames(polygons)[!colnames(polygons) %in% mandatory]
    polygons <- polygons[,c(mandatory, cnames)]
    return(polygons)

}

#' readPolygonsXenium
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
    polygons <- readPolygons(polygonsFile=polygonsFile, type=type, x=x, y=y,
        keepMultiPol=keepMultiPol, verbose=verbose)
    mandatory <- c("cell_id", "global", "is_multi", "multi_n")
    cnames <- colnames(polygons)[!colnames(polygons) %in% mandatory]
    polygons <- polygons[,c(mandatory, cnames)]
    return(polygons)
}

#' readPolygonsMerfish
#'
#' @param polygonsFolder
#' @param type
#' @param hdf5pattern
#'
#' @return
#' @export
#' @importFrom sf st_sf st_cast st_sfc
#' @importFrom arrow read_parquet
#'
#' @examples
readPolygonsMerfish <- function(polygonsFolder, type=c("HDF5", "parquet"),
                                keepMultiPol=TRUE, hdf5pattern="hdf5",
                                z_lev=3L, zcolumn="ZIndex",
                                geometry="Geometry",
                                verbose=FALSE)
{
    type <- match.arg(type)
    if (type=="HDF5")
    {
        polfiles <- list.files(polygonsFolder, pattern=hdf5pattern,
                                full.names=TRUE)
        dfsfl <- lapply(seq_along(polfiles), function(i)
        {
            poll <- readh5polygons(pol_file=polfiles[i])
            df <- data.frame(cell_id=paste0("f", i-1, "_c", poll$ids),
                             cell_ID=poll$ids,
                             fov=i-1, geometry=sf::st_sfc(poll$g)) ## geometry can be a simple column
            dfsf <- sf::st_sf(df)
        })
        polygons <- do.call(rbind, dfsfl)
        ## check if polygons geometry are polygons
        polygons <- .checkPolygonsValidity(polygons, keepMultiPol=keepMultiPol)
        return(polygons)
    } else { ## case parquet
        polfile <- list.files(polygonsFolder, pattern=type,
                              full.names=TRUE)
        polygons <- arrow::read_parquet(polfile, as_data_frame=TRUE)
        polygons <- polygons[polygons[[zcolumn]]==z_lev,]
        polygons$cell_id <- polygons$EntityID
        polygons <- .createPolygons(polygons)
        polygons <- .renameGeometry(polygons, geometry, "global")
        polygons <- .checkPolygonsValidity(polygons, keepMultiPol=keepMultiPol,
                                        verbose=verbose)

    }
    mandatory <- c("cell_id", "global", "is_multi", "multi_n")
    cnames <- colnames(polygons)[!colnames(polygons) %in% mandatory]
    polygons <- polygons[,c(mandatory, cnames)]
    return(polygons)
}

#' computeAreaFromPolygons
#'
#' @param polygons
#' @param coldata
#'
#' @return
#' @export
#'
#' @examples
computeAreaFromPolygons <- function(polygons, coldata)
{
    cd <- coldata
    # cd$Area <- NA
    area <- sf::st_area(polygons)
    # idx <- match(names(area), rownames(cd))
    cd$um_area <- unlist(area)
    return(cd)
}

#' computeAspectRatioFromPolygons
#'
#' @param polygons
#' @param coldata
#'
#' @return
#' @export
#'
#' @examples
computeAspectRatioFromPolygons <- function(polygons, coldata)
{
    cd <- coldata
    stopifnot("cell_id" %in% colnames(cd))
    # aspRatL <- list()
    aspRatL <- lapply(polygons$global[!polygons$is_multi], function(x) ## get active name of geometry instead of global
    {
        # xx <- polygons$global[!polygons$is_multi]
        # for(i in seq_along(xx))
        # {
            # print(i)
            # x <- xx[[i]]
            # aspRatL[[i]] <- (max(x[[1]][,2]) - min(x[[1]][,2]))/(max(x[[1]][,1]) - min(x[[1]][,1]))
        # }
        (max(x[[1]][,2]) - min(x[[1]][,2]))/(max(x[[1]][,1]) - min(x[[1]][,1]))
    }) ## to parallelize with bplapply
    names(aspRatL) <- polygons$cell_id[!polygons$is_multi]

    cd$AspectRatio <- NA
    posz <- match(names(aspRatL), cd$cell_id)
    cd$AspectRatio[posz] <- unlist(aspRatL)
    return(cd)
}


#' readh5polygons
#'
#' @param pol_file
#' @author Lambda Moses
#' @return
#' @export
#'
#' @examples
readh5polygons <- function(pol_file)
{
    l <- rhdf5::h5dump(pol_file)[[1]]
    cell_ids <- names(l)
    geometries <- lapply(l, function(m) {
        sf::st_polygon(list(t(m[["zIndex_0"]]$p_0$coordinates[,,1])))
    })
    return(list(g=geometries, ids=cell_ids))
}

