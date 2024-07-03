#' spatialPerCellQC
#'
#' @param spe
#' @param negProbList
#' @param micronConvFact
#'
#' @return
#' @export
#' @importFrom scater addPerCellQC
#' @importFrom S4Vectors cbind.DataFrame
#'
#' @examples
#' spe <- readCosmxSPE("~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/", sample_name="DBKero_BC")
#' spe <- spatialPerCellQC(spe)
spatialPerCellQC <- function(spe, micronConvFact=0.12,
    negProbList=c("NegPrb", "Negative", "SystemControl", # CosMx
        "NegControlProbe", "NegControlCodeWord", "UnassignedCodeWord", # Xenium
        "Blank" # MERFISH
        ))
{
    ## CHECK EXISTENCE OF POLYGONS/AREAS ETC -> create function for metrics creation
    stopifnot(is(object=spe, "SpatialExperiment"))
    idxlist <- lapply(negProbList, function(ng){
        grep(paste0("^", ng), rownames(spe))
    })
    names(idxlist) <- negProbList
    idxlist <- idxlist[which(lengths(idxlist)!=0)] ## lengths instead of length

    spe <- addPerCellQC(spe, subsets=idxlist)
    idx <- grep("^subsets_.*_sum$", colnames(colData(spe)))

    if ( length(idx) !=0 )
    {
        ## TO TEST -> BENEDETTA
        npc <- rowSums(as.matrix(colData(spe)[ , idx, drop=FALSE])) # meglio dataframe, perché rowsums non funziona
        ## getting detected probes as the column suddenly after the sum column
        npd <- rowSums(as.matrix(colData(spe)[ , idx+1, drop=FALSE])) # not robust at all!
    }

    spe$control_sum <- npc
    spe$control_detected <- npd
    spe$target_sum <- spe$sum - npc
    spe$target_detected <- spe$detected - npd

    #### CHANGE SPE constructor WITH COORDINATES IN COLDATA #########
    colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))

    if(metadata(spe)$technology == "Nanostring_CosMx")
    {
        spnc <- spatialCoords(spe) * micronConvFact
        colnames(spnc) <- gsub("px", "um", spatialCoordsNames(spe))
        colData(spe) <- cbind.DataFrame(colData(spe), spnc)
        spe$Area_um <- spe$Area * (micronConvFact^2)
    }

    #### compute AspectRatio for other technologies ####
    if ("AspectRatio" %in% colnames(colData(spe)))
    {
        warning(paste0("Missing aspect ratio in colData...\n",
                "NB: This could lead to additional warnings or errors.\n",
                "Missing AspectRatio can be computed by loading polygons."))
        spe$log2AspectRatio <- log2(spe$AspectRatio)
    }

    spe$ctrl_total_ratio <- spe$control_sum/spe$total
    spe$ctrl_total_ratio[which(is.na(spe$ctrl_total_ratio))] <- 0

    return(spe)
}
