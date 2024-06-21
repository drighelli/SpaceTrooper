#' readMerfishSPE
#'
#' @description
#'
#'
#' @param dirname
#' @param sample_name
#' @param countmatfpattern
#' @param metadatafpattern
#' @param coord_names
#' @param polygonsfpattern
#' @param fov_dims
#'
#' @return A SpatialExperiment object
#' @export
#'
#' @importFrom data.table fread merge
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr left_join
#' @importFrom SpatialExperiment SpatialExperiment
#' @examples
## for old fovs consider dimensions 5472 x 3648 pixels.
readMerfishSPE <- function(dirname,
                           countmatfpattern = "cell_by_gene.csv",
                           metadatafpattern = "cell_metadata.csv",
                           polygonsfpattern = "cell_boundaries.parquet",
                           coord_names = c("center_x", "center_y"))
{
    countmat_file <- list.files(dirname, countmatfpattern, full.names=TRUE)
    metadata_file <- list.files(dirname, metadatafpattern, full.names=TRUE)
    pol_file <- list.files(dirname, polygonsfpattern, full.names=TRUE) #check if parquet

    # stopifnot(all(file.exists(countmat_file), file.exists(metadata_file),
    #               file.exists(fovpos_file), file.exists(pol_file)))

    # Read in
    countmat <- data.table::fread(countmat_file)
    names(countmat)[names(countmat) == "V1"] <- "cell_id"
    metadata <- data.table::fread(metadata_file) # cell metadata
    names(metadata)[names(metadata) == "V1"] <- "cell_id"

    # Count matrix
    countmat <- left_join(countmat, metadata[, "cell_id"], by = "cell_id")
    cn <- countmat$cell_id
    counts <- subset(countmat, select = -cell_id)

    features <- colnames(counts)
    counts <- t(as.matrix(counts))

    rownames(counts) <- features
    colnames(counts) <- cn

    # rowData (does not exist)
    # To be associated to the tx file
    # use readSparseCSV sparseArray from harve pege

    # colData
    colData <- left_join(metadata, countmat[, "cell_id"], by = "cell_id")
    rownames(metadata) <- metadata$cell_id
    colData <- subset(colData, select = c(2,1,3:dim(colData)[2]))


    spe <- SpatialExperiment::SpatialExperiment(
        sample_id=sample_name,
        assays = list(counts = counts),
        # rowData = rowData,
        colData = colData,
        spatialCoordsNames = coord_names,
        metadata=list(polygons=pol_file, technology="Vizgen_MERFISH")
    )

    return(spe)

}
