#' readPolygons
#' @description
#'  Read and Validate Polygons from a File
#'
#' This function reads polygon data from a specified file, validates the polygons,
#' and returns them as an `sf` object. It supports multiple file formats and can
#' handle both global and local coordinates.
#'
#' @param polygonsFile A character string specifying the path to the polygon
#' file.
#' @param type A character string specifying the file type. Supported types are
#' `"csv"`, `"parquet"`, and `"h5"`. Default is `"csv"`.
#' @param x A character vector specifying the column names for the x-coordinates
#' in the polygon data. Default is `c("x_global_px", "vertex_x")`.
#' @param y A character vector specifying the column names for the y-coordinates
#' in the polygon data. Default is `c("y_global_px", "vertex_y")`.
#' @param xloc A character string specifying the column name for the local
#' x-coordinates. Default is `"x_local_px"`.
#' @param yloc A character string specifying the column name for the local
#' y-coordinates. Default is `"y_local_px"`.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons
#' during validation. Default is `TRUE`.
#' @param verbose A logical value indicating whether to print additional
#' information during processing. Default is `FALSE`.
#'
#' @return An `sf` object with the loaded and validated polygons.
#'
#' @details The function reads polygon data from the specified file and formats.
#' It validates the polygons and handles both global and local coordinates if
#' provided. If the file type is `"h5"`, the function currently does not handle
#' the data, as this part of the code is not implemented.
#'
#' @importFrom data.table fread
#' @importFrom arrow read_parquet
#' @importFrom sf st_geometry
#' @export
#'
#' @examples
#' # Reading polygon data from a CSV file:
#' # polygons <- readPolygons("~/Downloads/CosMx_data/polygons.csv")
#'
#' # Reading polygon data from a Parquet file with verbose output:
#' # polygons <- readPolygons("~/Downloads/CosMx_data/polygons.parquet",
#' #                         type = "parquet", verbose = TRUE)
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
#' @param sf An `sf` class object containing the spatial data.
#' @param geometry character for the geometry to check validity, if `NULL`
#' it checks the active geometry (default is `NULL`)
#' @param keepMultiPol logical for keeping/removing moltipolygons, if any
#' (default is `TRUE`, so keeping the multipolygons)
#' @param verbose logical to print verbose output (default is `FALSE`)
#'
#' @return An `sf` object with valid geometries, possibly with multipolygons
#' removed.
#' @keywords internal
#' @importFrom sf st_is_valid st_buffer st_geometry_type
#'
#' @examples
#' # TBD
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

readAndAddPolygonsToSPE(spe, keepMultiPol=TRUE,
                    boundaries_type=c("HDF5", "parquet")
                    )
{
    boundaries_type<-match.arg(boundaries_type)
    stopifnot("technology" %in% names(metadata(spe)))
    tech <- metadata(spe)$technology
    if(is.null(polygons))
    {
        switch(tech,
               "Nanostring_CosMx"=
                   {
                       polygons <- readPolygonsCosmx(metadata(spe)$polygons)
                   },
               "Vizgen_MERFISH"=
                   {
                       # ifelse(boundaries_type=="HDF5", merpol=)
                       polygons <- readPolygonsMerfish(polygonsFolder,
                                        keepMultiPol=TRUE, type=boundaries_type)

                   },
               "10X_Xenium"=
                   {
                       polygons <- readPolygonsXenium(pol_file, keepMultiPol=TRUE)
                   },
               stop("Unrecognized technology, please use an SPE from one of ",
                    "Nanostring_CosMx, Vizgen_MERFISH and 10x_Xenium")
        )
    }
    spe <- addPolygonsToSPE(spe, polygons)
}



#' addPolygonsToSPE
#'
#' @description This function adds polygon data to a `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object to which polygons will be added.
#' @param polygons An `sf` object containing the polygon data.
#'
#' @return The `SpatialExperiment` object with polygons added to the `colData`.
#'
#' @importFrom dplyr left_join
#' @export
#'
#' @examples
#' # Assuming `spe` is a SpatialExperiment object and `polygons` is an sf object:
#' # spe <- addPolygonsToSPE(spe, polygons)
addPolygonsToSPE <- function(spe, polygons=NULL)
{
    stopifnot(all(is(spe, "SpatialExperiment"), is(polygons, "sf")))

    if (sum(rownames(polygons) == colnames(spe)) != dim(spe)[2])
    {
        cd <- data.frame(colData(spe))
        polygons <- left_join(polygons, cd[, c("fov", "cellID")],
                              by=c("fov", "cellID"))
        cd <- left_join(cd, polygons[ , c("fov", "cellID")],
                        by=c("fov","cellID"))
        polygons$cell_id <- paste0("f", polygons$fov, "_c", polygons$cellID)
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
#' @description This internal function creates polygons from a spatial data object.
#'
#' @param spat_obj A data frame or similar object containing spatial data.
#' @param x A character vector specifying the x-coordinates.
#' @param y A character vector specifying the y-coordinates.
#' @param polygon_id A character string specifying the polygon ID.
#'
#' @return An `sf` object containing the created polygons.
#' @keywords internal
#' @importFrom sfheaders sf_polygon
#' @importFrom sf st_as_sf
#'
#' @examples
#' #TBD
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
#'
#' @description This function reads polygon data specific to CosMx technology.
#'
#' @param polygonsFile A character string specifying the file path to the
#' polygon data.
#' @param type A character string specifying the file type ("csv" or "parquet").
#' @param x A character string specifying the x-coordinate column.
#' @param y A character string specifying the y-coordinate column.
#' @param xloc A character string specifying the local x-coordinate column.
#' @param yloc A character string specifying the local y-coordinate column.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons.
#' @param verbose A logical value indicating whether to print additional
#' information.
#'
#' @return An `sf` object containing the CosMx polygon data.
#' @export
#'
#' @examples
#' # Read CosMx polygon data from a CSV file:
#' # polygons <- readPolygonsCosmx("path/to/polygons.csv", type="csv")
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
#' @description This function reads polygon data specific to Xenium technology.
#'
#' @param polygonsFile A character string specifying the file path to the
#' polygon data.
#' @param type A character string specifying the file type ("parquet" or "csv").
#' @param x A character string specifying the x-coordinate column.
#' @param y A character string specifying the y-coordinate column.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons.
#' @param verbose A logical value indicating whether to print additional
#' information.
#'
#' @return An `sf` object containing the Xenium polygon data.
#' @export
#'
#' @examples
#' # Read Xenium polygon data from a Parquet file:
#' # polygons <- readPolygonsXenium("path/to/polygons.parquet", type="parquet")
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
#' @description This function reads polygon data specific to MERFISH technology.
#'
#' @param polygonsFolder A character string specifying the folder containing the
#' polygon data files.
#' @param type A character string specifying the file type ("HDF5" or "parquet").
#' @param hdf5pattern A character string specifying the pattern to match HDF5
#' files.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons.
#' @param z_lev An integer specifying the Z level to filter the data. Default is
#' `3L`.
#' @param zcolumn A character string specifying the column name for the Z index.
#' @param geometry A character string specifying the geometry column name.
#' @param verbose A logical value indicating whether to print additional
#' information.
#'
#' @return An `sf` object containing the MERFISH polygon data.
#' @export
#' @importFrom sf st_sf st_cast st_sfc
#' @importFrom arrow read_parquet
#'
#' @examples
#' # Read MERFISH polygon data from a Parquet file:
#' # polygons <- readPolygonsMerfish("path/to/polygon_folder", type="parquet")
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
#' @description This function computes the area from polygon data and adds it to the `colData`.
#'
#' @param polygons An `sf` object containing polygon data.
#' @param coldata A `DataFrame` containing the `colData` to which area information will be added.
#'
#' @return A `DataFrame` with the added area information.
#' @export
#'
#' @examples
#' # Assuming `polygons` is an sf object and `coldata` is a DataFrame:
#' # coldata <- computeAreaFromPolygons(polygons, coldata)
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
#' @description This function computes the aspect ratio from polygon data and adds it to the `colData`.
#'
#' @param polygons An `sf` object containing polygon data.
#' @param coldata A `DataFrame` containing the `colData` to which aspect ratio information will be added.
#'
#' @return A `DataFrame` with the added aspect ratio information.
#' @export
#'
#' @examples
#' # Assuming `polygons` is an sf object and `coldata` is a DataFrame:
#' # coldata <- computeAspectRatioFromPolygons(polygons, coldata)
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
#' @description This function reads polygon data from an HDF5 file.
#'
#' @param pol_file A character string specifying the file path to the HDF5 polygon data.
#'
#' @return A list containing the polygon geometries and their associated cell IDs.
#' @author Lambda Moses
#' @export
#'
#' @examples
#' # Read polygons from an HDF5 file:
#' # polygons <- readh5polygons("path/to/polygons.h5")
readh5polygons <- function(pol_file)
{
    l <- rhdf5::h5dump(pol_file)[[1]]
    cell_ids <- names(l)
    geometries <- lapply(l, function(m) {
        sf::st_polygon(list(t(m[["zIndex_0"]]$p_0$coordinates[,,1])))
    })
    return(list(g=geometries, ids=cell_ids))
}

