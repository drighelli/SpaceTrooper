#' spatialPerCellQC
#'
#' @param spe
#'
#' @return
#' @export
#' @importFrom scater addPerCellQC
#' @importFrom S4Vectors cbind.DataFrame
#'
#' @examples
#' spe <- readCosmxSPE("~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/", sample_name="DBKero_BC")
spatialPerCellQC <- function(spe,
    negProbList=c("NegPrb", "Negative", "SystemControl", # CosMx
        "NegControlProbe", "NegControlCodeWord", "UnassignedCodeWord", # Xenium
        "Blank" # MERFISH
        ), micronConvFact=0.12)
{
    stopifnot(is(object=spe, "SpatialExperiment"))
    idxlist <- lapply(negProbList, function(ng){
        grep(paste0("^", ng), rownames(spe))
    })
    names(idxlist) <- negProbList
    idxlist <- idxlist[which(length(idxlist)!=0)]

    spe <- addPerCellQC(spe, subsets=idxlist)
    idx <- grep("^subsets_.*_sum$", colnames(colData(spe)))
    if(length(idx)>1)
    {   ## TO TEST -> BENEDETTA
        npc <- rowSums(matrix(colData(spe)[,idx]))
        ## getting detected probes as the column suddenly after the sum column
        npd <- rowSums(matrix(colData(spe)[,idx+1])) # not robust at all!
    } else{
        npc <- as.matrix(colData(spe)[, idx, drop=FALSE])
        ## getting detected probes as the column suddenly after the sum column
        npd <- as.matrix(colData(spe)[, idx+1, drop=FALSE]) # not robust at all!
    }

    spe$control_sum <- npc
    spe$control_detected <- npd
    spe$target_sum <- spe$sum - npc
    spe$target_detected <- spe$detected - npd

    #### CHANGE SPE WITH COORDINATES IN COLDATA
    colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))

    if(metadata(spe)$technology == "Nanostring_CosMx")
    {
        spnc <- spatialCoords(spe) * micronConvFact
        colnames(spnc) <- gsub("px", "um", spatialCoordsNames(spe))
        colData(spe) <- cbind.DataFrame(colData(spe), spnc)
    }






    return(spe)
}
