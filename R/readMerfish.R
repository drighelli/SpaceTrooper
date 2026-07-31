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
#'   In particular for area, if a "volume" column is present in the colData, it
#'   will be used as area value, otherwise area will be computed from polygons.
#'   This is relevant for MERFISH as the volume is computed on the entire 3D
#'   cell available data, while polygons are 2D sections.
#' @param keepPolygons `logical(1)`
#'   If `TRUE`, attach raw polygon geometries as extra columns in `colData`.
#'   Default: `FALSE`.
#' @param boundariesType `character(1)`
#' One of `"parquet"` or `"HDF5"`. If `"parquet"`, reads a single Parquet
#' file of boundaries; if `"HDF5"`, uses a folder of HDF5 polygon files.
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
#' @param useVolume `logical(1)` If `TRUE`, prefer a "volume" column in the
#' metadata for area (when present). Default: `TRUE`.
#' @return A `SpatialExperiment` object with:
#'   - `assays$counts`: gene × cell count matrix
#'   - `colData`: per‐cell metadata (including computed metrics)
#'   - spatial coordinates named by `coordNames`
#'   - `metadata$polygons`: path to the boundaries file
#'   - `metadata$technology`: `"Vizgen_MERFISH"`.
#'
#' @details
#' The function searches `dirName` for three resources using the provided
#' filename patterns: the count matrix (`countmatFPattern`), the cell metadata
#' (`metadataFPattern`), and the polygon resource (`polygonsFPattern`). The
#' count matrix and metadata are expected to contain a `cell_id` column
#' (aliases `V1`, `cell`, or `EntityID` are renamed). If `computeMissingMetrics`
#' is TRUE, `computeMissingMetricsMerfish()` is invoked and receives `useVolume`.
#'
#' The returned `SpatialExperiment` contains:
#' - `assays$counts`: gene × cell count matrix (rows = genes, cols = `cell_id`);
#' - `colData`: per-cell metadata (may be reordered by the function);
#' - spatial coordinates named by `coordNames`;
#' - `metadata$polygons`: path(s) matched by `polygonsFPattern`;
#' - `metadata$technology`: `"Vizgen_MERFISH"`.
#'
#' Side effects: the function may reorder `colData` columns, attach polygon
#' file paths to `metadata(spe)$polygons`, and (when requested) append
#' `Area_um` and `AspectRatio` to `colData`.
#'
#'
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
    boundariesType=c("parquet", "HDF5"), countmatFPattern="cell_by_gene.csv",
    metadataFPattern="cell_metadata.csv",
    polygonsFPattern="cell_boundaries.parquet",
    coordNames=c("center_x", "center_y"), polygonsCol="polygons",
    useVolume=TRUE) {

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
                                            keepPolygons, polygonsCol, useVolume)
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
#' @param polFile character for the path to the polygon file. Tipically a
#' parquet file or a folder of HDF5 files in the `metadata(spe)$polygons`.
#' @param coldata `DataFrame` or `data.frame`
#' Cell metadata with at least a `cell_id` column.
#' @param boundariesType `character(1)`
#' One of `"HDF5"` or `"parquet"`—passed on to `readPolygonsMerfish()`.
#' @param keepPolygons `logical(1)`
#' If `TRUE`, cbinds the raw polygon `sf` columns onto `coldata`.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#' @param useVolume `logical(1)` it assigns the area from the "volume" column.
#' If the column is not present it computes the area from the polygons.
#' Default: `TRUE`.
#'
#' @return A `DataFrame` (or `data.frame`) with:
#'   - all columns of `coldata`
#'   - `Area_um`: area/volume of each cell’s polygon
#'   - `AspectRatio`: width/height aspect ratio
#'   - (optionally) the polygon geometries
#'
#' @details
#' `computeMissingMetricsMerfish()` reads polygon geometries via
#' `readPolygonsMerfish(polFile, type=boundariesType)` where `polFile`
#' is a Parquet file (parquet) or a folder of HDF5 files (HDF5). It expects
#' `coldata` to contain at least `cell_id`.
#'
#' Behavior:
#' - If `useVolume=TRUE` and `coldata` contains a `volume` column, the
#'   returned `Area_um` is set to `volume`. Otherwise `Area_um` is computed
#'   from the polygons using `computeAreaFromPolygons()`.
#' - `AspectRatio` is computed from polygons using `computeAspectRatioFromPolygons()`.
#' - If `keepPolygons=TRUE`, geometries are attached to the returned
#'   `DataFrame` under the name given by `polygonsCol`.
#'
#' Important: polygons must align (in order or by matching logic) with the
#' rows of `coldata`; otherwise area/aspect values may be assigned incorrectly.
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
    boundariesType=c("parquet", "HDF5"), keepPolygons=FALSE,
    polygonsCol="polygons", useVolume=TRUE) {

    boundariesType <- match.arg(boundariesType)
    polygons <- readPolygonsMerfish(polFile, type=boundariesType)
    cd <- DataFrame(coldata)

    if (useVolume)
    {
        if(!"volume" %in% colnames(cd))
        {
            warning("Volume column not found in colData.\nComputing area from polygons instead.")
            area <- computeAreaFromPolygons(polygons)
        } else {
            warning("Volume is used to compute Quality Score for MERFISH technology.
                For simplicity, it is renamed as Area_um.")
            area <- cd$volume
        }
    }else {
        area <- computeAreaFromPolygons(polygons)
    }
    cd$Area_um <- area
    cd$AspectRatio <- computeAspectRatioFromPolygons(polygons)
    if (keepPolygons) cd <- .addPolygonsToCD(cd, polygons, polygonsCol)
    return(cd)
}

