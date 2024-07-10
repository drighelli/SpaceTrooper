#' @rdname readXeniumSPE
#'
#' @title Load data from a 10x Geonomics Xenium experiment
#'
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded unzipped Xenium
#' Output Bundle directory for 10x Genomics Xenium spatial gene expression data.
#'
#' @param dirname a directory path to Xenium Output Bundle download that contains
#' files of interest.
#' @param sample_name
#' @param type
#' @param coord_names
#' @param boundaries_type
#' @param compute_missing_metrics
#' @param countsfilepattern
#' @param metadatafpattern a filename pattern of the zipped .csv file that
#' contains cell metadata and spatial coords. Default value is \code{"cells.csv.gz"}, and there is no
#' need to change.
#' @param polygonsfpattern a vector of two strings specify the spatial coord names.
#' Default value is \code{c("x_centroid", "y_centroid")}, and there is no need to change.
#'
#' @details # WHAT ABOUT THE OTHER PARAMETERS/FILES? WHAT ABOUT THE OUTS FOLDER(added)?
#' The constructor assumes the downloaded unzipped Xenium Output Bundle has the
#' following structure, with mandatory file of cells.csv.gz and either folder
#' /cell_feature_matrix or .h5 file cell_feature_matrix.h5:
#' Xenium_unzipped \cr
#'    | - outs
#'        | — cell_feature_matrix.h5 \cr
#'        | — cell_feature_matrix \cr
#'            | - barcodes.tsv.gz \cr
#'            | - features.tsv.gz \cr
#'            | - matrix.mtx.gz \cr
#'        | — cells.csv.gz \cr
#'
#' @return a \code{\link{SpatialExperiment}} object
#'
#' @author Estella Yixing Dong
#'
#' @examples
#' \dontrun{
#' # Data download is from:
#' # https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_
#' # Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
#'
#' xepath <- system.file(
#'   file.path("extdata", "10xXenium"),
#'   package = "SpatialExperiment")
#'
#' list.files(xepath)
#'
#' # read the count matrix .h5 file
#' xe_spe <- readXeniumSPE(dirname = xepath,
#'                         countfname = "cell_feature_matrix.h5",
#'                         coordfpattern = "cells.csv.gz",
#'                         coord_names = c("x_centroid", "y_centroid"))
#'
#' # or read the count matrix folder
#' xe_spe <- readXeniumSPE(dirname = xepath,
#'                         countfname = "cell_feature_matrix",
#'                         coordfpattern = "cells.csv.gz",
#'                         coord_names = c("x_centroid", "y_centroid"))
#'
#' # Subset to no control genes
#' xe_spe <- xe_spe[rowData(xe_spe)$Type == "Gene Expression"]
#' }
#' @importFrom DropletUtils read10xCounts
#' @importFrom data.table fread
#' @importFrom SpatialExperiment SpatialExperiment
readXeniumSPE <- function(dirname,
                          sample_name="sample01",
                          type=c("HDF5", "sparse"),
                          coord_names=c("x_centroid", "y_centroid"),
                          boundaries_type=c("parquet", "csv"),
                          compute_missing_metrics=TRUE,
                          countsfilepattern="cell_feature_matrix",
                          metadatafpattern="cells",
                          polygonsfpattern="cell_boundaries")
{
    stopifnot(file.exists(dirname))
    type <- match.arg(type)
    boundaries_type <- match.arg(boundaries_type)

    # add "outs/" directory if not already included
    if(basename(dirname) != "outs")
    {
        dirbkup <- dirname
        dirname <- file.path(dirname, "outs")
        if (!file.exists(dirname))
        {
            dirname <- dirbkup
        } else {
            warning("automatically detected/added outs dir in the 10x filepath")
        }
    }

    cfm <- paste0(countsfilepattern, switch(type, HDF5=".h5", ""))
    counts <- file.path(dirname, cfm)

    metadata_file <- file.path(dirname, paste0(metadatafpattern, ".csv.gz"))
    pex <- paste0(polygonsfpattern, switch(boundaries_type,
                                            parquet=".parquet",
                                            csv=".csv.gz"))
    pol_file <- list.files(dirname, pex, full.names=TRUE)
    stopifnot(all(file.exists(c(metadata_file, pol_file))))

    # Count matrix + rowData
    sce <- DropletUtils::read10xCounts(counts, col.names=TRUE)

    # Spatial and colData
    cd <- DataFrame(fread(metadata_file, header=TRUE))
    rownames(cd) <- cd$cell_id

    if ( dim(sce)[2] != dim(cd)[1] )
    {
        sce <- sce[, colnames(sce) %in% rownames(cd)]
        cd <- cd[rownames(cd) %in% colnames(sce), ]
    }
    if (compute_missing_metrics)
    {
        message("Computing missing metrics, this could take some time...")
        cd <- computeMissingMetricsXenium(pol_file, cd)
    }
    # construct 'SpatialExperiment'
    spe <- SpatialExperiment::SpatialExperiment(
        sample_id=sample_name,
        assays = assays(sce),
        rowData = rowData(sce),
        colData = cd,
        spatialCoordsNames = coord_names,
        metadata=list(polygons=pol_file, technology="10X_Xenium")
    )
    return(spe)
}


computeMissingMetricsXenium <- function(pol_file, coldata)
{
    polygons <- readPolygonsXenium(pol_file, keepMultiPol=TRUE)
    cd <- computeAspectRatioFromPolygons(polygons, coldata)
    return(cd)
}

