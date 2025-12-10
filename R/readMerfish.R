#' readMerfishSPE
#' @name readMerfishSPE
#' @rdname readMerfishSPE
#' @description
#' `readMerfishSPE()` imports MERFISH/Merscore outputs (counts, metadata,
#' and optionally cell boundary polygons) from a directory and builds a
#' SpatialExperiment object.
#'
#' @param dirName `character(1)`
#'   Path to a folder containing MERFISH output files.
#' @param sampleName `character(1)`
#'   Identifier to assign to the `sample_id` field in the returned object.
#'   Default: `"sample01"`.
#' @param computeMissingMetrics `logical(1)`
#'   If `TRUE`, compute area and aspect‐ratio metrics from the cell boundary
#'   polygons. Default: `TRUE`.
#' @param keepPolygons `logical(1)`
#'   If `TRUE`, attach raw polygon geometries as extra columns in `colData`.
#'   Default: `FALSE`.
#' @param boundariesType `character(1)`
#'   One of `"HDF5"` or `"parquet"`. If `"HDF5"`, uses a folder of HDF5
#'   polygon files; if `"parquet"`, reads a single Parquet file of boundaries.
#' @param countmatFPattern `character(1)`
#'   Regex passed to `list.files()` to find the count matrix CSV. Default:
#'   `"cell_by_gene.csv"`.
#' @param metadataFPattern `character(1)`
#'   Pattern to find the cell metadata CSV. Default: `"cell_metadata.csv"`.
#' @param polygonsFPattern `character(1)`
#'   Pattern to find the cell boundaries file. Default:
#'   `"cell_boundaries.parquet"`.
#' @param coordNames `character(2)`
#'   Names of the columns in `colData` that store X/Y spatial coordinates.
#'   Default: `c("center_x", "center_y")`.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#'
#' @return A `SpatialExperiment` object with:
#'   - `assays$counts`: gene × cell count matrix
#'   - `colData`: per‐cell metadata (including computed metrics)
#'   - spatial coordinates named by `coordNames`
#'   - `metadata$polygons`: path to the boundaries file
#'   - `metadata$technology`: `"Vizgen_MERFISH"`.
#' @author Dario Righelli, Benedetta Banzi
#' @export
#' @importFrom data.table fread
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr left_join
#' @importFrom SpatialExperiment SpatialExperiment
#' @examples
#' path <- system.file("extdata", "Merfish_Tiny",
#'                     package = "SpaceTrooper")
#' spe <- readMerfishSPE(
#'   dirName = path,
#'   sampleName = "Patient2",
#'   keepPolygons = TRUE,
#'   boundariesType = "parquet"
#' )
#' spe
readMerfishSPE <- function(dirName, sampleName="sample01",
    computeMissingMetrics=TRUE, keepPolygons=FALSE,
    boundariesType=c("HDF5", "parquet"), countmatFPattern="cell_by_gene.csv",
    metadataFPattern="cell_metadata.csv",
    polygonsFPattern="cell_boundaries.parquet",
    coordNames=c("center_x", "center_y"), polygonsCol="polygons") {

    countmat_file <- list.files(dirName, countmatFPattern, full.names=TRUE)
    metadata_file <- list.files(dirName, metadataFPattern, full.names=TRUE)
    countmat <- data.table::fread(countmat_file)
    # spe <- readMerscopeSXE(dirName=dirName, returnType="SPE",
    #     countMatPattern=countmatFPattern, metaDataPattern=metadataFPattern,
    #     coordNames=coordNames)
    pol_file <- list.files(dirName, polygonsFPattern, full.names=TRUE) #parquet?
    names(countmat)[names(countmat) %in% c("V1", "cell")] <- "cell_id"
    metadata <- data.table::fread(metadata_file) # cell metadata
    names(metadata)[names(metadata) %in% c("EntityID", "V1")] <- "cell_id"
    countmat <- left_join(countmat, metadata[, "cell_id"], by = "cell_id")
    cn <- countmat$cell_id
    counts <- subset(countmat, select = -cell_id)
    features <- colnames(counts)
    counts <- t(as.matrix(counts))
    rownames(counts) <- features
    colnames(counts) <- as.character(cn)
    # TODO: rowData (does not exist) tx file use readSparseCSV from sparseArray
    colData <- metadata[match(countmat$cell_id, metadata$cell_id), ]
    rownames(colData) <- colData$cell_id
    cd <- subset(colData, select = c(2,1,3:dim(colData)[2]))
    if (computeMissingMetrics) {
        message("Computing missing metrics, this could take a while...")
        cd <- computeMissingMetricsMerfish(pol_file, colData, boundariesType,
                                            keepPolygons, polygonsCol)
    }
    spe <- SpatialExperiment::SpatialExperiment(sample_id=sampleName,
        assays=list(counts=counts), colData=cd,
        spatialCoordsNames=coordNames,
        metadata=list(polygons=pol_file, technology="Vizgen_MERFISH"))
    colnames(spe) <- spe$cell_id
    return(spe)
}

#' computeMissingMetricsMerfish
#' @name computeMissingMetricsMerfish
#' @rdname computeMissingMetricsMerfish
#' @description
#' `computeMissingMetricsMerfish()` takes cell metadata and boundary
#' polygons, calculates per‐cell area and aspect‐ratio, and optionally
#' appends the raw polygon geometries.
#'
#' @param coldata `DataFrame` or `data.frame`
#' Cell metadata with at least a `cell_id` column.
#' @param boundariesType `character(1)`
#'   One of `"HDF5"` or `"parquet"`—passed on to
#'   `readPolygonsMerfish()`.
#' @param keepPolygons `logical(1)`
#'   If `TRUE`, cbinds the raw polygon `sf` columns onto `coldata`.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#' @param polFile path to the polygon file
#'
#' @return A `DataFrame` (or `data.frame`) with:
#'   - all columns of `coldata`
#'   - `um_area`: area of each cell’s polygon
#'   - `AspectRatio`: width/height aspect ratio
#'   - (optionally) the polygon geometries
#'
#' @export
#' @importFrom S4Vectors DataFrame
#' @examples
#' example(readMerfishSPE)
#' cd <- computeMissingMetricsMerfish(metadata(spe)$polygons, colData(spe),
#'     boundariesType="parquet")
#' colData(spe) <- cd
#' cd
computeMissingMetricsMerfish <- function(polFile, coldata,
    boundariesType=c("parquet","HDF5"), keepPolygons=FALSE,
    polygonsCol="polygons") {
    boundariesType <- match.arg(boundariesType)
    polygons <- readPolygonsMerfish(polFile, type=boundariesType)
    cd <- DataFrame(coldata)
    cd$Area_um <- computeAreaFromPolygons(polygons)
    cd$AspectRatio <- computeAspectRatioFromPolygons(polygons)
    if (keepPolygons) cd <- .addPolygonsToCD(cd, polygons, polygonsCol)
    return(cd)
}

