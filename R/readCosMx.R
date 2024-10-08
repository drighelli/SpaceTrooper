#' readCosmxSPE
#' @description
#' Read and Construct a SpatialExperiment Object from CosMx Data
#'
#' This function reads in data from Nanostring CosMx files and constructs a
#' `SpatialExperiment` object, optionally including polygon data.
#'
#' @param dirname A character string specifying the directory containing the
#' CosMx data files.
#' @param sample_name A character string specifying the sample name. Default is
#' `"sample01"`.
#' @param coord_names A character vector specifying the names of the spatial
#' coordinate columns in the data. Default is `c("CenterX_global_px",
#' "CenterY_global_px")`.
#' @param countmatfpattern A character string specifying the pattern to match
#' the count matrix file. Default is `"exprMat_file.csv"`.
#' @param metadatafpattern A character string specifying the pattern to match
#' the metadata file. Default is `"metadata_file.csv"`.
#' @param polygonsfpattern A character string specifying the pattern to match
#' the polygons file. Default is `"polygons.csv"`.
#' @param keep_polygons A logical value indicating whether to include polygon
#' data in the resulting `SpatialExperiment` object. Default is `FALSE`.
#' @param fovposfpattern A character string specifying the pattern to match the
#' FOV positions file. Default is `"fov_positions_file.csv"`.
#' @param fov_dims A named numeric vector specifying the dimensions of the FOV
#' in pixels. Default is `c(xdim=4256, ydim=4256)`.
#'
#' @return A `SpatialExperiment` object containing the read CosMx data,
#' including count matrices, metadata, and optionally polygons.
#'
#' @details The function reads in the specified files for count matrices,
#' metadata, and FOV positions, and constructs a `SpatialExperiment` object.
#' Optionally, polygon data can be read and added to the object.
#'
#' @importFrom data.table fread merge.data.table
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' # Assuming the data files are located in "path/to/cosmx_data":
#' #spe <- readCosmxSPE(dirname = "path/to/cosmx_data")
## for old fovs consider dimensions 5472 x 3648 pixels.
readCosmxSPE <- function(dirname,
                        sample_name="sample01",
                        coord_names=c("CenterX_global_px", "CenterY_global_px"),
                        countmatfpattern="exprMat_file.csv",
                        metadatafpattern="metadata_file.csv",
                        polygonsfpattern="polygons.csv",
                        keep_polygons=FALSE,
                        fovposfpattern="fov_positions_file.csv",
                        fov_dims=c(xdim=4256, ydim=4256))
{
    stopifnot(all(names(fov_dims) == c("xdim", "ydim"), file.exists(dirname)))
    countmat_file <- list.files(dirname, countmatfpattern, full.names=TRUE)
    metadata_file <- list.files(dirname, metadatafpattern, full.names=TRUE)
    fovpos_file <- list.files(dirname, fovposfpattern, full.names=TRUE)
    pol_file <- list.files(dirname, polygonsfpattern, full.names=TRUE) #check if parquet

    # stopifnot(all(file.exists(countmat_file), file.exists(metadata_file),
    #               file.exists(fovpos_file), file.exists(pol_file)))

    # Read in
    countmat <- data.table::fread(countmat_file, showProgress=FALSE) # cell count matrix
    metadata <- data.table::fread(metadata_file, showProgress=FALSE) # cell metadata

    # Count matrix
    counts <- merge.data.table(countmat, metadata[, c("fov", "cell_ID")])
    cn <- paste0("f", counts$fov, "_c", counts$cell_ID)
    counts <- subset(counts, select = -c(fov, cell_ID))

    # cell_codes <- rownames(counts)
    features <- colnames(counts)
    counts <- t(as.matrix(counts)) #### faster when it comes to big numbers

    rownames(counts) <- features
    colnames(counts) <- cn

    # rowData (does not exist)
    # To be associated to the tx file
    # use readSparseCSV sparseArray from harve pege

    # colData
    colData <- DataFrame(merge.data.table(metadata,
                                        countmat[, c("fov", "cell_ID")]))
    rn <- paste0("f", colData$fov, "_c", colData$cell_ID)
    rownames(colData) <- rn

    if(length(grep("cell_id", colnames(colData)))!=0)
        message("Warning: overwriting existing cell_id column in colData")

    colData$cell_id <- rn
    colData <- colData[,c(1,2,dim(colData)[2], 3:(dim(colData)[2]-1))]
    ## multiply spatial coordinates in micron multiplying by 0.18 and store
    ## two additional columns depends by technology version

    fov_positions <- as.data.frame(data.table::fread(fovpos_file, header=TRUE))

    ## patch for let this work also with older versions of CosMx fov position
    ## output file
    fovcidx <- grep("FOV", colnames(fov_positions))
    if(length(fovcidx)!=0) colnames(fov_positions)[fovcidx] <- "fov"
    fovccdx <- grep("[X|Y]_px", colnames(fov_positions))
    if(length(fovccdx)!=0)
    {
        colnames(fov_positions)[fovccdx] <- gsub("_px", "_global_px",
                                            colnames(fov_positions)[fovccdx])
    }

    ## tracking if one of more fov is not present in the metadata file ##
    idx <- fov_positions$fov %in% unique(metadata$fov)
    fov_positions <- fov_positions[idx,]
    fov_positions <- fov_positions[order(fov_positions$fov),]
    ####
    spe <- SpatialExperiment::SpatialExperiment(
        sample_id=sample_name,
        assays = list(counts = counts),
        # rowData = rowData,
        colData = colData,
        spatialCoordsNames = coord_names,
        metadata=list(fov_positions=fov_positions, fov_dim=fov_dims,
                        polygons=pol_file, technology="Nanostring_CosMx")
        ## keep atomx versioning in metadata, if possible
    )

    #### CHANGE SPE constructor WITH COORDINATES IN COLDATA #########
    colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))

    # Polygons file has cellID instead of cell_ID and it distinguish better
    # when compared to our cell_id
    names(colData(spe))[names(colData(spe))=="cell_ID"] <- "cellID"
    if(keep_polygons)
    {
        polygons <- readPolygonsCosmx(metadata(spe)$polygons)
        spe <- addPolygonsToSPE(spe, polygons)
    }
    return(spe)
}

