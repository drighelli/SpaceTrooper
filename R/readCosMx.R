#' readCosmxSPE
#' @name readCosmxSPE
#' @rdname readCosmxSPE
#' @aliases readCosmxSPE readCosmxProteinSPE
#' @description
#' Read and Construct a SpatialExperiment Object from CosMx Data
#'
#' This function reads in data from Nanostring CosMx files and constructs a
#' `SpatialExperiment` object, optionally including polygons data.
#'
#' @param dirName A character string specifying the directory containing the
#' CosMx data files.
#' @param sampleName A character string specifying the sample name. Default is
#' `"sample01"`.
#' @param coordNames A character vector specifying the names of the spatial
#' coordinate columns in the data. Default is `c("CenterX_global_px",
#' "CenterY_global_px")`.
#' @param countMatFPattern A character string specifying the pattern to match
#' the count matrix file. Default is `"exprMat_file.csv"`.
#' @param metadataFPattern A character string specifying the pattern to match
#' the metadata file. Default is `"metadata_file.csv"`.
#' @param polygonsFPattern A character string specifying the pattern to match
#' the polygons file. Default is `"polygons.csv"`.
#' @param fovPosFPattern A character string specifying the pattern to match the
#' FOV positions file. Default is `"fov_positions_file.csv"`.
#' @param fovdims A named numeric vector specifying the dimensions of the FOV
#' in pixels. Default is `c(xdim=4256, ydim=4256)`.
#' @param keepPolygons Logical indicating if the polygons need to be loaded into
#' memory or not (Default is `FALSE`).
#' @param polygonsCol Character name of the column in colData to store polygons.
#'
#' @return A `SpatialExperiment` object containing the read CosMx data,
#' including count matrices, metadata, and optionally polygons.
#'
#' @details The function firstly relies on
#' \link[SpatialExperimentIO]{readCosmxSXE} to read in the specified files
#' for count matrices, metadata, and FOV positions constructing a
#' `SpatialExperiment` object.
#' Then it harmonizes the object to have the same metadata as for the other
#' technologies, setting the colData names as required in further QC analysis.
#'
#' Optionally, polygons data can be read and added to the object by seeting
#' the `keepPolygons` argument to `TRUE`, otherwise it only stores the
#' polygons file path into the object metadata.
#'
#'
#' readCosmxProteinSPE is a wrapper of readCosmxSPE, it changes the
#' technology metadata in Nanostring_CosMx_Protein.
#'
#' @author Dario Righelli, Benedetta Banzi
#'
#' @importFrom data.table fread merge.data.table
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr mutate
#' @importFrom SpatialExperimentIO readCosmxSXE
#' @export
#'
#' @examples
#' cospath <- system.file(file.path("extdata", "CosMx_DBKero_Tiny"),
#'    package="SpaceTrooper")
#' spe <- readCosmxSPE(cospath, sampleName="DBKero_Tiny")
readCosmxSPE <- function(dirName, sampleName="sample01",
    coordNames=c("CenterX_global_px", "CenterY_global_px"),
    countMatFPattern="exprMat_file.csv", metadataFPattern="metadata_file.csv",
    polygonsFPattern="polygons.csv", fovPosFPattern="fov_positions_file.csv",
    fovdims=c(xdim=4256, ydim=4256), keepPolygons=FALSE, polygonsCol="polygons") {

    stopifnot(all(names(fovdims) == c("xdim", "ydim"), file.exists(dirName)))

    spe <- SpatialExperimentIO::readCosmxSXE(dirName=dirName, returnType="SPE",
        countMatPattern=countMatFPattern, metaDataPattern=metadataFPattern,
        coordNames=coordNames, addFovPos=TRUE, fovPosPattern=fovPosFPattern,
        altExps=NULL, addParquetPaths=FALSE)

    spe <- .setupCosmxSPE(spe, dirName, sampleName, polygonsFPattern, fovdims,
        keepPolygons, polygonsCol)
    return(spe)
}


#' @export
readCosmxProteinSPE <- function(dirName, sampleName="sample01",
    coordNames=c("CenterX_global_px", "CenterY_global_px"),
    countMatFPattern="exprMat_file.csv", metadataFPattern="metadata_file.csv",
    polygonsFPattern="polygons.csv", fovPosFPattern="fov_positions_file.csv",
    fovdims=c(xdim=4256, ydim=4256), keepPolygons=FALSE, polygonsCol="polygons")
{
    spe <- readCosmxSPE(dirName, sampleName, coordNames, countMatFPattern,
        metadataFPattern, polygonsFPattern, fovPosFPattern, fovdims,
        keepPolygons, polygonsCol)

    metadata(spe)$technology <- "Nanostring_CosMx_Protein"
    return(spe)
}

#' .checkFovPositionVersion
#' @rdname dot-checkFovPositionVersion
#' @description
#'
#' Check and Standardize FOV Position Column Names
#'
#' This internal utility function standardizes column names of a data frame
#' containing Field of View (FOV) positional information.
#' It modifies column names to ensure compatibility with expected naming
#' conventions, including support for older formats.
#'
#' Specifically, it:
#' - Renames any column containing "FOV" to "fov"
#' - Converts columns with coordinates matching "X", "Y", or "Z" to lowercase
#' - Replaces suffix "_px" with "_global_px" for coordinate pixel columns
#' - If the input contains `x_mm` and `y_mm` columns, the function computes
#' corresponding `x_global_px` and `y_global_px` values by converting from
#' millimeters to pixels using a fixed resolution factor (0.12028 mm/pixel).
#' @param spe A `SpatialExperiment` containing FOV position information
#' in the metadata to be standardized.
#'
#' @return A `SpatialExperiment` with updated and standardized column names
#' for the metadata `fov_position` `data.frame`.
#'
#' @keywords internal
.checkFovPositionVersion <- function(spe)
{
    fovpos <- metadata(spe)$fov_positions
    fovcidx <- grep("FOV", colnames(fovpos)) # works also with older vers
    if(length(fovcidx)!=0) colnames(fovpos)[fovcidx] <- "fov"
    fovcrdx <- grep("[X|Y|Z]", colnames(fovpos))
    if(length(fovcrdx)!=0) colnames(fovpos)[fovcrdx] <-
        tolower(colnames(fovpos)[fovcrdx])
    fovccdx <- grep("[x|y]_px", colnames(fovpos))
    if(length(fovccdx)!=0) colnames(fovpos)[fovccdx] <-
        gsub("_px", "_global_px", colnames(fovpos)[fovccdx])

    if(length(grep("x_mm", colnames(fovpos))!=0)) {
        fovpos <- fovpos |>
            dplyr::mutate(x_global_px = x_mm/0.12028*10^3,
                          y_global_px = (y_mm/0.12028*10^3))
    }

    idx <- fovpos$fov %in% unique(spe$fov)
    fovpos <- fovpos[idx, ]

    fovpos <- fovpos[order(fovpos$fov), ]
    metadata(spe)$fov_positions <- fovpos
    return(spe)
}


#' updateCosmxSPE
#'
#' @description
#' Update a SpatialExperiment object derived from CosMx data by adding
#' polygons, FOV dimensions, standardized column names, and metadata.
#'
#' @param spe SpatialExperiment object.
#' @param dirName Directory containing CosMx output files (e.g., polygon CSVs).
#' @param sampleName Character scalar, sample identifier stored in
#'   \code{colData(spe)$sample_id}. Default \code{"sample01"}.
#' @param polygonsFPattern Character, pattern used by \code{list.files()} to
#'   locate polygon files. Default \code{"polygons.csv"}.
#' @param fovdims Named numeric vector with entries \code{xdim} and \code{ydim}
#'   representing the FOV dimensions in pixels.
#' @param keepPolygons Logical indicating if the polygons need to be loaded into
#' memory or not (Default is `FALSE`).
#' @param polygonsCol Character name of the column in colData to store polygons.
#' @details
#' The function standardizes CosMx SPE structure by:
#' - creating unique cell names of the form \code{f<fov>_c<cell_ID>};
#' - ensuring consistent cell identifiers and sample metadata;
#' - recording FOV dimensions, polygon paths, and technology type in
#'   \code{metadata(spe)}.
#'
#' @return A SpatialExperiment object with updated metadata and column names.
#'
#' @seealso readCosmxSPE, readCosmxProteinSPE
#'
#' @export
#' @examples
#' cospath <- system.file(file.path("extdata", "CosMx_DBKero_Tiny"),
#'     package="SpaceTrooper")
#' spe <- SpatialExperimentIO::readCosmxSXE(dirName=cospath,
#'     addParquetPaths=FALSE)
#' spe <- updateCosmxSPE(spe, dirName=cospath, sampleName="DBKero_Tiny")
updateCosmxSPE <- function(spe, dirName, sampleName="sample01",
                        polygonsFPattern="polygons.csv",
                        fovdims=c(xdim=4256, ydim=4256), keepPolygons=FALSE,
                        polygonsCol="polygons") {
    stopifnot("spe is not a SpatialExperiment"=is(spe, "SpatialExperiment"))
    stopifnot("fovdims not x|y dim"=all(names(fovdims) == c("xdim", "ydim")))
    stopifnot("dirName not exists"=file.exists(dirName))
    spe <- .setupCosmxSPE(spe, dirName, sampleName, polygonsFPattern, fovdims,
        keepPolygons, polygonsCol)
    return(spe)
}

#' .setupCosmxSPE
#'
#' @description
#' Internal helper to configure a SpatialExperiment from CosMx data:
#' standardizes column names, sets cell IDs, and updates metadata.
#'
#' @param spe SpatialExperiment object.
#' @param dirName Directory containing CosMx output files.
#' @param sampleName Character sample identifier. Default \code{"sample01"}.
#' @param polygonsFPattern Character, pattern for polygon files.
#' Default \code{"polygons.csv"}.
#' @param fovdims Named numeric vector \code{c(xdim=, ydim=)} giving
#' field-of-view size.
#' @param keepPolygons Logical indicating if the polygons need to be loaded into
#' memory or not (Default is `FALSE`).
#' @param polygonsCol Character name of the column in colData to store polygons.
#'
#' @return Updated SpatialExperiment object.
#'
#' @keywords internal
#' @noRd
.setupCosmxSPE <- function(spe, dirName, sampleName="sample01",
                            polygonsFPattern="polygons.csv",
                            fovdims=c(xdim=4256, ydim=4256),
                            keepPolygons=FALSE, polygonsCol="polygons")
{
    pol_file <- list.files(dirName, polygonsFPattern, full.names=TRUE)
    cn <- paste0("f", spe$fov, "_c", spe$cell_ID)
    colnames(spe) <- cn
    rownames(colData(spe)) <- cn
    if(length(grep("cell_id", colnames(colData)))!=0)
        warning("Overwriting existing cell_id column in colData")
    spe$cell_id <- cn
    spe <- .checkFovPositionVersion(spe)
    metadata(spe) <- list(fov_positions=metadata(spe)$fov_positions,
            fov_dim=fovdims, polygons=pol_file, technology="Nanostring_CosMx")

    names(colData(spe))[names(colData(spe)) == "cell_ID"] <- "cellID"
    spe$sample_id <- sampleName
    if (keepPolygons) spe <- readAndAddPolygonsToSPE(spe, polygonsCol=polygonsCol)
    return(spe)
}

#' updateCosmxProteinSPE
#'
#' @description
#' Update a SpatialExperiment object corresponding to Nanostring CosMx
#' Protein data by adding metadata identifying the technology and
#' optionally passing through file-location parameters.
#'
#' @param spe SpatialExperiment object.
#' @param dirName Directory containing CosMx Protein data files.
#' @param sampleName Character sample ID. Default \code{"sample01"}.
#' @param coordNames Character vector of length two indicating coordinate
#'   column names in the per-cell metadata. Default
#'   \code{c("CenterX_global_px","CenterY_global_px")}.
#' @param countMatFPattern Character pattern for the counts matrix file.
#' @param metadataFPattern Character pattern for the single-cell metadata file.
#' @param polygonsFPattern Character pattern for the polygon file(s).
#' @param fovPosFPattern Character pattern for the FOV positions file.
#' @param fovdims Named numeric vector with FOV size in pixels.
#' @param keepPolygons Logical indicating if the polygons need to be loaded into
#' memory or not (Default is `FALSE`).
#' @param polygonsCol Character name of the column in colData to store polygons.
#' @details
#' This function sets \code{metadata(spe)$technology <- "Nanostring_CosMx_Protein"}.
#' It does not modify other assay or metadata components.
#'
#' @return SpatialExperiment object with updated technology metadata.
#'
#' @seealso readCosmxProteinSPE, readCosmxSPE
#'
#' @export
#' @examples
#' protfolder <- system.file( "extdata", "S01_prot", package="SpaceTrooper")
#' spe <- SpatialExperimentIO::readCosmxSXE(dirName=protfolder,
#'     addParquetPaths=FALSE)
#' spe <- updateCosmxProteinSPE(spe, protfolder, sampleName="cosmx_prots")
updateCosmxProteinSPE <- function(spe, dirName, sampleName="sample01",
    coordNames=c("CenterX_global_px", "CenterY_global_px"),
    countMatFPattern="exprMat_file.csv", metadataFPattern="metadata_file.csv",
    polygonsFPattern="polygons.csv", fovPosFPattern="fov_positions_file.csv",
    fovdims=c(xdim=4256, ydim=4256), keepPolygons=FALSE, polygonsCol="polygons")
{

    stopifnot(is(spe, "SpatialExperiment"))
    spe <- updateCosmxSPE(spe, dirName, sampleName, polygonsFPattern, fovdims,
        keepPolygons, polygonsCol)
    metadata(spe)$technology <- "Nanostring_CosMx_Protein"
    return(spe)
}
