#' readCosMx
#'
#' @description
#'
#' @param sample character indicating the name of the sample
#' @param path character indicating the path pointing to the sample
#' @param rm.unmpd.tx logical (default TRUE) indicating if the unmapped
#' molecules have to be discarded
#' @param fov optional numeric indicating a list of FoVs to consider,
#' if NULL (default) all FoVs are considered
#' @param discard.tx character indicating a list of molecules to discard,
#' if NULL (default) all genes are kept
#' @param scalefactor
#' @param force logical, if TRUE the function uses the `path` as passed from
#' the user
#' @param verbose
#'
#' @importFrom data.table fread merge
#' @importFrom S4Vectors DataFrame
#' @importFrom SparseArray readSparseCSV
#' @examples
## for old fovs consider dimensions 5472 x 3648 pixels.
readCosmxSPE <- function(dirname=dirname,
                        countmatfpattern="exprMat_file.csv",
                        metadatafpattern="metadata_file.csv",
                        coord_names=c("CenterX_global_px", "CenterY_global_px"),
                        polygonsfpattern="polygons.csv",
                        fovposfpattern="fov_positions_file.csv",
                        fov_dims=c(xdim=4256, ydim=4256))
{
    stopifnot(all(names(fov_dims) == c("xdim", "ydim")))
    countmat_file <- list.files(dirname, countmatfpattern, full.names=TRUE)
    metadata_file <- list.files(dirname, metadatafpattern, full.names=TRUE)
    fovpos_file <- list.files(dirname, fovposfpattern, full.names=TRUE)
    pol_file <- list.files(dirname, polygonsfpattern, full.names=TRUE)

    # Read in
    countmat <- data.table::fread(countmat_file) # cell count matrix
    metadata <- data.table::fread(metadata_file) # cell metadata

    # Count matrix
    counts <- merge(countmat, metadata[, c("fov", "cell_ID")])
    cn <- paste0("f", counts$fov, "_c", counts$cell_ID)
    counts <- subset(counts, select = -c(fov, cell_ID))
    # cell_codes <- rownames(counts)
    features <- colnames(counts)
    counts <- t(counts)

    rownames(counts) <- features
    colnames(counts) <- cn

    # rowData (does not exist)
    # To be associated to the tx file
    # use readSparseCSV from harve pege

    # colData
    colData <- DataFrame(merge(metadata, countmat[, c("fov", "cell_ID")]))
    rn <- paste0("f", colData$fov, "_c", colData$cell_ID)
    rownames(colData) <- rn
    ## multiply spatial coordinates in micron multiplying by 0.18 and store
    ## two additional columns depends by technology version

    fov_positions <- as.data.frame(data.table::fread(fovpos_file, header = T))

    ## patch for let this work also with older versions of CosMx fov position output file
    fovcidx <- grep("FOV", colnames(fov_positions))
    if(length(fovcidx)!=0) colnames(fov_positions)[fovcidx] <- "fov"
    fovccdx <- grep("[X|Y]_px", colnames(fov_positions))
    if(length(fovccdx)!=0) colnames(fov_positions)[fovccdx] <- gsub("_px", "_global_px", colnames(fov_positions)[fovccdx])

    ## tracking if one of more fov is not present in the metadata file
    idx <- fov_positions$fov %in% unique(metadata$fov)
    fov_positions <- fov_positions[idx,]
    fov_positions <- fov_positions[order(fov_positions$fov),]
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        # rowData = rowData,
        colData = colData,
        spatialCoordsNames = coord_names,
        metadata=list(fov_positions=fov_positions, fov_dim=fov_dims,
                        polygons=pol_file)
    )
    return(spe)
}
