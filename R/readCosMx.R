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
#' @importFrom data.table fread
#' @importFrom SparseArray readSparseCSV
#' @examples
# readCosMx <- function(sample, path,
#     rm.unmpd.tx=TRUE,
#     fov=NULL, discard.tx=NULL,
#     scalefactor=0.18, force=FALSE, verbose=FALSE)
# {
#     ###########
#     path="~/Downloads/CosMx_Data"
#     sample="Lung5_Rep1"
#     force=FALSE
#     verbose=FALSE
#     ##########
#     datadir <- paste0(sample, "-Flat_files_and_images")
#     if(!force)
#     {
#         if (length(grep(sample, path)) == 0)
#         {
#             if (length(grep(datadir, sample)) != 0)
#                 stop("Please provide a valid path or use the force argument")
#             path <- file.path(path, sample, datadir)
#         } else if (length(grep(datadir, path)) == 0) {
#             path <- file.path(path, datadir)
#         }
#     }
#
#     nanorules <- c("[_a-zA-Z0-9]*_exprMat_file.csv",
#                    "[_a-zA-Z0-9]*_metadata_file.csv",
#                    "[_a-zA-Z0-9]*_tx_file.csv")#, "[_a-zA-Z0-9]*-polygons.csv")
#     af <- list.files(path=path, recursive=TRUE, full.names=TRUE)
#     af <- af[unlist(lapply(nanorules, grep, af))]
#     names(af) <- c("counts", "molecules", "metadata")
#     if (length(af) > 3)
#         stop("Found more files than expected!\n",
#             "Check if there is more than one experiment in the ", path)
#     message("Reading input files from ", path)
#     #### code for sparse array is faster, but at the moment is under dedvelopment
#     # system.time({ml <- fread(af["counts"])})
#     # mls <- readSparseCSV(af["counts"])
#     # afo <- lapply(af[-which(names(af)=="counts")], fread, verbose=verbose)
#     # names(afo) <- gsub(".csv", "", basename(af[-which(names(af)=="counts")]))
#     afo <- lapply(af, fread, verbose=verbose)
#     names(afo) <- gsub(".csv", "", basename(af))
#     ########## All these filters need to be done on a complete object
#     ########## this should allow to do the subsetting only once on all the tables
#     txidx <- grep("tx_file", names(afo))
#     if (rm.unmpd.tx)
#     {
#         ## remove transcripts with cellID=0, they are not assigned to any cell
#         afo[[txidx]] <- afo[[txidx]][!afo[[txidx]]$cell_ID == 0,]
#     }
#     if (!is.null(fov)) ## subsetting FoVs based on the input fov parameter
#     {
#         afo[[txidx]] <- afo[[txidx]][afo[[txidx]]$fov == fov,]
#     }
#     ## the tx tibble needs to be mapped in a BumpyMatrix
#     if (!is.null(discard.tx))
#     {
#         afo[[txidx]] <- afo[[txidx]][afo[[txidx]]$target == discard.tx,]
#     }
#     ##########
#     mols <- rep_len(x = files[["molecules.file"]],
#                     length.out = length(x = mol.type))
#     names(x = mols) <- mol.type
#     files <- c(files, mols)
#     files <- files[setdiff(x = names(x = files), y = "molecules.file")]
#
#
#     # library(pryr)
#     # object_size(txdata)
#
#     # length(table(txdata$target))
#     # table(txdata$fov, txdata$cell_ID)
#
#
#
#     # mol <- BumpyMatrix::splitAsBumpyMatrix(
#     #     df[, c("x", "y")],
#     #     row = df$gene, column = df$cell)
# }
## for old fovs consider dimensions 5472 x 3648 pixels.
readCosmxSPE <- function(dirname=dirname,
                        countmatfpattern="exprMat_file.csv",
                        metadatafpattern="metadata_file.csv",
                        coord_names=c("CenterX_global_px", "CenterY_global_px"),
                        fovposfpattern="fov_positions_file.csv",
                        fov_dims=c(xdim=4256, ydim=4256)){
    stopifnot(all(names(fov_dims) == c("xdim", "ydim")))
    countmat_file <- file.path(dirname, list.files(dirname, countmatfpattern))
    metadata_file <- file.path(dirname, list.files(dirname, metadatafpattern))
    fovpos_file <- file.path(dirname, list.files(dirname, fovposfpattern))

    # Read in
    countmat <- data.table::fread(countmat_file)
    metadata <- data.table::fread(metadata_file)

    # Count matrix
    counts <- merge(countmat, metadata[, c("fov", "cell_ID")])
    counts <- subset(counts, select = -c(fov, cell_ID))
    cell_codes <- rownames(counts)
    features <- colnames(counts)
    counts <- t(counts)

    rownames(counts) <- features
    colnames(counts) <- cell_codes

    # rowData (does not exist)
    # To be associated to the tx file

    # colData
    colData <- merge(metadata, countmat[, c("fov", "cell_ID")])


    ############# NOTES #################
    # to add FoV positions file as metadata(spe)$fov_positions
    # to add fov dimesions (they changes based on the technology) as
    # metadata(spe)$fov_dim
    ## FOV positions file to load and check the presence of all FOVs in the
    ## colData

    fov_positions <- as.data.frame(data.table::fread(fovpos_file, header = T))
    idx <- fov_positions$fov %in% unique(metadata$fov)
    ## track if one of more fov is not present in the metadata file
    # metadata(spe)$fov_positions <- fov_positions[idx,]
    # metadata(spe)$fov_dim <- fov_dims
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        # rowData = rowData,
        colData = colData,
        spatialCoordsNames = coord_names,
        metadata=list(fov_positions=fov_positions[idx,], fov_dim=fov_dims)
    )
    return(spe)
}
