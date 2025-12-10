#' readPolygons
#' @name readPolygons
#' @rdname readPolygons
#' @description
#' Read and Validate Polygons from a File
#'
#' This function reads polygon data from a specified file, validates the
#' polygons, and returns them as an `sf` object. It supports multiple file
#' formats and can handle both global and local coordinates.
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
#' @importFrom dplyr group_by mutate ungroup filter n
#' @export
#' @examples
#' example(readCosmxSPE)
#' polygons <- readPolygons(metadata(spe)$polygons)
#' polygons
readPolygons <- function(polygonsFile, type=c("csv", "parquet", "h5"),
    x=c("x_global_px", "vertex_x"), y=c("y_global_px", "vertex_y"),
    xloc="x_local_px", yloc="y_local_px", #micronConvFact=0.12,
    keepMultiPol=TRUE, verbose=FALSE){
    stopifnot(file.exists(polygonsFile))
    type <- match.arg(type)
    x <- match.arg(x)
    y <- match.arg(y)

    if(type=="h5") {
        stop("h5 support is momentarily disabled, please report an issue on GH")
    } else {
        spat_obj <- switch(type, csv=fread(polygonsFile),
                                parquet=read_parquet(polygonsFile))
        if (! "cell_id" %in% colnames(spat_obj)) {
            spat_obj$cell_id <- paste0("f", spat_obj$fov, "_c", spat_obj$cellID)
        }
        spat_obj$cell_id <- as.factor(spat_obj$cell_id)
        # removing polygons with less than 4 points
        if(any(as.vector(table(spat_obj$cell_id)<4)))
        {
            warning(paste0("Removing ",table(table(spat_obj$cell_id)<4)["TRUE"],
                            " polygons with less than 4 points."))
            spat_obj <- spat_obj |> group_by(cell_id) |> mutate(n_points = n()) |>
                ungroup() |>  filter(n_points >= 4)
        }
        polygons <- .createPolygons(spat_obj, x=x, y=y,
                                        polygon_id="cell_id")
        polygons <- .renameGeometry(polygons, "geometry", "global") ##for cosmx
        if(all(c(xloc, yloc) %in% colnames(spat_obj))) {
            idxs <- which(colnames(polygons) %in% c(xloc, yloc))
            polygons <- polygons[,-idxs]
            polygons_loc <- .createPolygons(spat_obj, x=xloc, y=yloc,
                                            polygon_id="cell_id")[,-c(4,5)]
            polygons <- cbind(polygons, polygons_loc$geometry)
            polygons <- .renameGeometry(polygons, "geometry", "local")
        }
        if(verbose) message("Polygons detected: ", dim(polygons)[1])
        polygons <- .checkPolygonsValidity(polygons, keepMultiPol=keepMultiPol,
                                            verbose=verbose)
        rownames(polygons) <- polygons$cell_id #needed x check addPolygonsToSPE
        if(verbose) message("Polygons after validity: ", dim(polygons)[1])
        return(polygons)
    }
}


#' .checkPolygonsValidity
#' @name .checkPolygonsValidity
#' @rdname dot-checkPolygonsValidity
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
.checkPolygonsValidity <- function(sf, geometry=NULL, keepMultiPol=TRUE,
    verbose=FALSE) {
    stopifnot(is(sf, "sf"))
    if(is.null(geometry)) { geometry <- .getActiveGeometryName(sf)
    } else {
        act <- .getActiveGeometryName(sf)
        sf <- .setActiveGeometry(sf, geometry)
    }
    sf_tf <- st_is_valid(sf) # TODO: parallelize splitting sf ?
    if(sum(sf_tf)!=dim(sf)[1]) sf <- st_buffer(sf, dist=0)
    is_multi <- sf::st_geometry_type(sf[[geometry]]) == "MULTIPOLYGON" ## <<<<
    ## TRUE is merscope - FALSE is cosmx and xenium
    merscopeFl <- (sum(is_multi) == dim(sf)[1])
    funct <- ifelse(merscopeFl, "lengths", "length")
    ll <- lapply(sf[[geometry]], funct)
    sums <- lapply(ll, sum)
    idx <- which(sums!=1)
    sf$is_multi <- FALSE
    sf$multi_n <- 1
    if(length(idx)!=0) {
        sf$is_multi[idx] <- TRUE
        sf$multi_n[idx] <- unlist(sums[idx])
    }
    if( merscopeFl ) {
        if( length(idx)!=0 ) {
            sf$global[-idx] <- st_cast(sf$global[-idx], "POLYGON")
        } else {
            sf$global <- st_cast(sf$global, "POLYGON")
        }
    }

    if(verbose) message("Detected ", sum(sf$is_multi), " multipolygons.")
    if(!keepMultiPol) {
        if(verbose) message("Removing ", sum(sf$is_multi), " multipolygons.")
        sf <- sf[!sf$is_multi,]
    }
    if(exists("act")) sf <- .setActiveGeometry(sf, act)
    return(sf)
}

#' readAndAddPolygonsToSPE
#' @name readAndAddPolygonsToSPE
#' @rdname readAndAddPolygonsToSPE
#' @description Read and Add Polygons to a SpatialExperiment Object
#'
#' This function reads polygon boundary data based on the technology associated
#' with the provided SpatialExperiment (SPE) object and adds the polygons to the
#' SPE.
#'
#' @param spe A \code{SpatialExperiment} object. The object should contain
#' metadata with the field \code{"technology"}, specifying the technology used
#' (e.g., "Nanostring_CosMx", "Vizgen_MERFISH", or "10X_Xenium").
#' @param keepMultiPol Logical. If \code{TRUE}, multi-polygon features will be
#' kept when reading the boundary data. Defaults to \code{TRUE}.
#' @param boundariesType Character. Specifies the type of boundary file format
#' to read. Options are \code{"HDF5"} or \code{"parquet"}. Defaults to
#' \code{"HDF5"}.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#'
#' @return A \code{SpatialExperiment} object with the added polygon data.
#'
#' @export
#'
#' @examples
#' example(readCosmxSPE)
#' spe <- readAndAddPolygonsToSPE(spe)
#' colData(spe)
readAndAddPolygonsToSPE <- function(spe, polygonsCol="polygons",
    keepMultiPol=TRUE, boundariesType=c("csv", "HDF5", "parquet"))
{
    boundariesType <- match.arg(boundariesType)
    stopifnot("technology" %in% names(metadata(spe)))
    tech <- metadata(spe)$technology
    switch(tech,
            "Nanostring_CosMx_Protein"=,
            "Nanostring_CosMx"={
                polygons <- readPolygonsCosmx(metadata(spe)$polygons)
            },
            "Vizgen_MERFISH"={
                ##### NEED TO HANDLE THE DIFFERENCES BETWEEN HDF5 FILES
                ##### AND PARQUET, TO PROPAGATE TO READING FUNCTION
                # ifelse(boundariesType=="HDF5", merpol=)
                polygons <- readPolygonsMerfish(metadata(spe)$polygons,
                                keepMultiPol=TRUE, type=boundariesType)
            },
            "10X_Xenium"={
                polygons <- readPolygonsXenium(metadata(spe)$polygons,
                                            keepMultiPol=TRUE)
            },
            stop("Unrecognized technology, please use an SPE from one of ",
                "Nanostring_CosMx, Vizgen_MERFISH and 10x_Xenium")
    )
    spe <- addPolygonsToSPE(spe, polygons, polygonsCol=polygonsCol)
    return(spe)
}



#' .addPolygonsToCD
#' @name .addPolygonsToCD
#' @rdname dot-addPolygonsToCD
#' @description
#' This function enriches a DataFrame (e.g., from colData) with matching polygon
#' geometries.
#'
#' @param cd A DataFrame containing at least `fov` and `cellID` columns.
#' @param polygons An sf object with matching `fov` and `cellID` columns.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#' @return A DataFrame identical to `cd`, but row‑subset to cells present in
#' `polygons` and with a new `polygons` list‑column of sf geometries.
#' @keywords internal
.addPolygonsToCD <- function(cd, polygons, polygonsCol="polygons") {
    stopifnot(inherits(cd, "DataFrame"), inherits(polygons, "sf"))

    if("fov" %in% colnames(polygons)) ## case of CosMx
    {
        # Build consistent cell IDs
        cd$cell_id <- paste0("f", cd$fov, "_c", cd$cellID)
        polygons$cell_id <- paste0("f", polygons$fov, "_c", polygons$cellID)
        rownames(cd) <- cd$cell_id
        rownames(polygons) <- polygons$cell_id
    } else { ## valid for both xenium and merfish
        rownames(cd) <- cd$cell_id
        rownames(polygons) <- polygons$cell_id
    }

    # Identify matching cells
    common_ids <- intersect(rownames(cd), rownames(polygons))
    if (length(common_ids) == 0) {
        stop("No matching cell IDs between 'cd' and 'polygons'.")
    }
    if (length(common_ids) < nrow(cd)) {
        warning(sprintf(
            "%d/%d cells have polygon data; subsetting to matches.",
            length(common_ids), nrow(cd)
        ))
    }

    # Subset and attach
    cd <- cd[common_ids, , drop = FALSE]
    polygons <- polygons[common_ids, , drop = FALSE]
    cd[[polygonsCol]] <- polygons

    return(cd)
}


#' addPolygonsToSPE
#' @name addPolygonsToSPE
#' @rdname addPolygonsToSPE
#' @description This function adds polygon data to a `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object to which polygons will be added.
#' @param polygons An `sf` object containing the polygon data.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#' @return The `SpatialExperiment` object with polygons added to the `colData`.
#'
#' @importFrom dplyr left_join
#' @export
#'
#' @examples
#' example(readCosmxSPE)
#' polygons <- readPolygonsCosmx(metadata(spe)$polygons)
#' spe <- addPolygonsToSPE(spe, polygons)
#' spe$polygons
addPolygonsToSPE <- function(spe, polygons, polygonsCol="polygons") {
    stopifnot(inherits(spe, "SpatialExperiment"), inherits(polygons, "sf"))

    # Extract and enrich colData via DataFrame
    cd <- S4Vectors::DataFrame(colData(spe))
    cd <- .addPolygonsToCD(cd, polygons, polygonsCol)

    # Subset and reorder the SpatialExperiment
    spe <- spe[, which(colnames(spe) %in% cd$cell_id), drop = FALSE]

    # Replace colData with enriched DataFrame
    colData(spe) <- cd

    return(spe)
}


#' .createPolygons
#' @name .createPolygons
#' @rdname dot-createPolygons
#' @description This internal function creates polygons from a data.frame or
#' similar object.
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
.createPolygons <- function(spat_obj, x=NULL, y=NULL, polygon_id=NULL,
                            geometry="Geometry")
{
    if(all(!is.null(x), !is.null(y)))
    {
        polygons <- sfheaders::sf_polygon(obj=spat_obj, x=x, y=y,
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
#' @name readPolygonsCosmx
#' @rdname readPolygonsCosmx
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
#' example(readCosmxSPE)
#' polygons <- readPolygonsCosmx(metadata(spe)$polygons)
#' polygons
readPolygonsCosmx <- function(polygonsFile, type=c("csv", "parquet"),
                            x="x_global_px",
                            y="y_global_px",
                            xloc="x_local_px",
                            yloc="y_local_px",
                            keepMultiPol=TRUE,
                            verbose=FALSE)
{
    type <- match.arg(type)
    polygons <- readPolygons(polygonsFile, type=type, x=x, y=y, xloc=xloc,
                    yloc=yloc, keepMultiPol=keepMultiPol, verbose=verbose)
    polygons <- sf::st_cast(polygons, "GEOMETRY")
    mandatory <- c("cell_id", "global", "is_multi", "multi_n")
    cnames <- colnames(polygons)[!colnames(polygons) %in% mandatory]
    polygons <- polygons[,c(mandatory, cnames)]
    polygons <- .setActiveGeometry(polygons, "local")
    polygons <- sf::st_cast(polygons, "GEOMETRY")
    polygons <- .setActiveGeometry(polygons, "global")
    return(polygons)
}

#' readPolygonsXenium
#' @name readPolygonsXenium
#' @rdname readPolygonsXenium
#' @description This function reads polygon data specific to Xenium technology.
#'
#' @param polygonsFile A character string specifying the file path to the
#' polygon data.
#' @param type A character string specifying the file type ("parquet" or "csv").
#' Default is parquet.
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
#' example(readXeniumSPE)
#' polygons <- readPolygonsXenium(metadata(spe)$polygons, type="parquet")
#' polygons
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
#' @name readPolygonsMerfish
#' @rdname readPolygonsMerfish
#' @description This function reads polygon data specific to MERFISH technology.
#'
#' @param polygons A character string specifying the folder containing the
#' polygon data files in case of HDF5, or a path to a parquet file (see `type`).
#' @param type A character string specifying the file type("HDF5" or "parquet").
#' Default is parquet.
#' @param hdf5pattern A character string specifying the pattern to match HDF5
#' files.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons.
#' @param zLev An integer specifying the Z level to filter the data. Default is
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
#' example(readMerfishSPE)
#' polygons <- readPolygonsMerfish(metadata(spe)$polygons, type="parquet")
#' polygons
readPolygonsMerfish <- function(polygons, type=c("parquet", "HDF5"),
                                keepMultiPol=TRUE, hdf5pattern="hdf5",
                                zLev=3L, zcolumn="ZIndex",
                                geometry="Geometry",
                                verbose=FALSE)
{
    type <- match.arg(type)
    if (type=="HDF5")
    {
        polfiles <- list.files(polygons, pattern=hdf5pattern,
                                full.names=TRUE)
        dfsfl <- lapply(seq_along(polfiles), function(i)
        {
            poll <- readh5polygons(polFile=polfiles[i])
            df <- data.frame(cell_id=paste0("f", i-1, "_c", poll$ids),
                cell_ID=poll$ids, fov=i-1,
                geometry=sf::st_sfc(poll$g)) ## geometry can be a simple column
            dfsf <- sf::st_sf(df)
        })
        polygons <- do.call(rbind, dfsfl)
        ## check if polygons geometry are polygons
        polygons <- .checkPolygonsValidity(polygons, keepMultiPol=keepMultiPol)
        return(polygons)
    } else { ## case parquet
        # polfile <- list.files(polygonsFolder, pattern=type,
                            # full.names=TRUE)
        polfile <- polygons
        polygons <- arrow::read_parquet(polfile, as_data_frame=TRUE)
        polygons <- polygons[polygons[[zcolumn]]==zLev,]
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
#' computeCenterFromPolygons
#' @name computeCenterFromPolygons
#' @rdname computeCenterFromPolygons
#' @description This function computes the center coordinates on x and y axis
#' from polygon data and adds it to the `colData`. It is necessary only for
#' Merfish.
#'
#' @param polygons An `sf` object containing polygon data.
#' @param coldata A `DataFrame` containing the `colData` to which center
#' coordinates information will be added.
#'
#' @return A `DataFrame` with the added center information.
#' @export
#'
#' @examples
#' example(readPolygonsMerfish)
#' coldata <- computeCenterFromPolygons(polygons, colData(spe))
#' colData(spe) <- coldata
computeCenterFromPolygons <- function(polygons, coldata)
{
    # this needs to be changed to be independent from coldata
    cd <- coldata
    # cd$Area <- NA
    centroid <- sf::st_centroid(polygons)
    ## get active name of geometry instead of global
    center_x <- lapply(centroid$geometry, function(x)
    {
        x[1]
    })
    ## get active name of geometry instead of global
    center_y <- lapply(centroid$geometry, function(x)
    {
        x[2]
    })
    cd$center_x <- unlist(center_x)
    cd$center_y <- unlist(center_y)
    return(cd)
}


#' computeAreaFromPolygons
#' @name computeAreaFromPolygons
#' @rdname computeAreaFromPolygons
#' @description This function computes the area from polygon data.
#'
#' @param polygons An `sf` object containing polygon data.
#'
#' @return A `numeric` vector with the area information.
#' @export
#'
#' @examples
#' example(readPolygonsMerfish)
#' area <- computeAreaFromPolygons(polygons)
#' area
computeAreaFromPolygons <- function(polygons)
{
    stopifnot(is(polygons, "sf"))
    area <- sf::st_area(polygons)
    Area_um <- unlist(area)
    return(Area_um)
}

#' computeAspectRatioFromPolygons
#' @name computeAspectRatioFromPolygons
#' @rdname computeAspectRatioFromPolygons
#' @description This function computes the aspect ratio (width / height)
#' from polygon data.
#'
#' @param polygons An `sf` object containing polygon data.
#'
#' @return A `numeric` vector with the aspect ratio information.
#' @export
#' @importFrom sf st_geometry st_is_empty st_bbox
#' @examples
#' example(readPolygonsMerfish)
#' ar <- computeAspectRatioFromPolygons(polygons)
#' ar
computeAspectRatioFromPolygons <- function(polygons)
{
    stopifnot(inherits(polygons, "sf"))
    if (!"cell_id" %in% names(polygons)) {
        stop("'polygons' must contain a 'cell_id' column.")
    }

    g <- sf::st_geometry(polygons)
    n <- length(g)
    # counting multi polygons if present
    if ("is_multi" %in% names(polygons)) {
        is_multi <- polygons$is_multi %in% TRUE
        if (any(is_multi)) {
            warning(
                "Found ", sum(is_multi),
                " multi-polygons: returning NA aspect ratio for them."
            )
        }
    } else {
        is_multi <- rep(FALSE, n)
    }
    # helper fun: aspect ratio via bounding box (CosMx-style: width / height)
    .bbox_aspect_ratio <- function(geom) {
        # geom: sfg (POLYGON, MULTIPOLYGON, ecc.)
        bb <- sf::st_bbox(geom) # bounding box computes xmin, ymin, xmax, ymax
        width  <- unname(bb["xmax"] - bb["xmin"])
        height <- unname(bb["ymax"] - bb["ymin"])
        if (height == 0) {
            return(NA_real_)
        }
        width / height
    }

    ar <- rep(NA_real_, n)
    idx <- which(!is_multi & !sf::st_is_empty(g))
    if (length(idx) > 0) {
        ar[idx] <- vapply(g[idx], .bbox_aspect_ratio, numeric(1))
    }
    names(ar) <- polygons$cell_id
    return(ar)
}

#' readh5polygons
#' @name readh5polygons
#' @rdname readh5polygons
#' @description This function reads polygon data from an HDF5 file.
#'
#' @param polFile A character string specifying the file path to the HDF5
#' polygon data.
#'
#' @return A list containing the polygon geometries and their associated cell_id
#' @author Lambda Moses
#' @importFrom rhdf5 h5dump
#' @importFrom sf st_polygon
#' @keywords internal
readh5polygons <- function(polFile)
{
    l <- rhdf5::h5dump(polFile)[[1]]
    cell_ids <- names(l)
    geometries <- lapply(l, function(m) {
        sf::st_polygon(list(t(m[["zIndex_0"]]$p_0$coordinates[,,1])))
    })
    return(list(g=geometries, ids=cell_ids))
}

