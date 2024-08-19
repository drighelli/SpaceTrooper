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
        spe$log2AspectRatio <- log2(spe$AspectRatio)
    } else {
        warning(paste0("Missing aspect ratio in colData...\n",
                   "NB: This could lead to additional warnings or errors.\n",
                   "Missing AspectRatio can be computed by loading polygons."))
    }

    spe$ctrl_total_ratio <- spe$control_sum/spe$total
    spe$ctrl_total_ratio[which(is.na(spe$ctrl_total_ratio))] <- 0

    return(spe)
}


#' computeQCScore
#' @description
#' Compute for each cell a flag score based on the target_counts, area in micron
#' and log2 of the aspect ratio.
#'
#' @param spe a SpatialExperiment object with target_counts, area in micron
#' and log2 of the aspect ratio in the `colData`, tipically computed with
#' \link{spatialPerCellQC}
#'
#' @return
#' @export
#'
#' @examples
#' #changed default weights because I did tests and found that these works pretty good on DBKero
computeQCScore <- function(spe, a=0.3, b=1)#, threshold=0.15)
{
    stopifnot(is(spe, "SpatialExperiment"))
    cd <- colData(spe)
    if(!all(c("target_sum", "Area_um", "log2AspectRatio") %in% colnames(cd)))
    {
        stop("One of target_sum, Area_um, log2AspectRatio is missing, did you run spatialPerCellQC?")
    }
    spe$flag_score <- 1/(1 + exp(-a*spe$target_sum/spe$Area_um + b*abs(spe$log2AspectRatio)))

    return(spe)
}


#' computeSpatialOutlier
#' @description
#' Identifies outliers for cell area and mean DAPI signal through a statistical
#' test on the basis of the variable skewness
#'
#' @param spe a SpatialExperiment object with area in micron
#' and mean DAPI signal in the `colData`, tipically computed with
#' \link{spatialPerCellQC}
#'
#' @return
#' @export
#' @importFrom scuttle isOutlier
#' @importFrom robustbase mc adjbox
#' @examples
#' #TBD
computeSpatialOutlier <- function(spe)
{
    # control on the presence of cell area
    message("Outlier detection on cell area")
    stopifnot(all(c("Area_um") %in% colnames(colData(spe))))
    if(!skewness(spe$Area_um)[1]>0)
    {
        # skewness value must be > 0 to use Medcouple, redirection to MAD
        warning("Skewness value is not > 0. Outlier detection will be ",
                "performed using MAD")
        spe$is_area_outlier <- scuttle::isOutlier(spe$Area_um, type="both")
    }
    # Medcouple must be between -0.6 and 0.6
    if(!mc(spe$Area_um)>-0.6 & mc(spe$Area_um)<0.6)
    {
        warning("Medcouple value is not suitable for boxplot adjustment. ",
                "Outlier detection will be performed using MAD")
        spe$is_area_outlier <- scuttle::isOutlier(spe$Area_um, type="both")
    }
    names(spe$Area_um) <- colnames(spe)
    out <- robustbase::adjbox(spe$Area_um, plot = FALSE)
    message(paste0(length(out$out)," outliers have been found"))
    #saving area thresholds in spe metadata for later
    metadata(spe)$area_fence <- out$fence

    area_outliers <- vector(length=dim(spe)[2])
    names(area_outliers) <- colnames(spe)
    area_outliers[which(names(area_outliers) %in% names(out$out))] <- TRUE
    spe$is_area_outlier <- area_outliers

    # repeating the entire code for DAPI variable, again not very optimized
    if("Mean.DAPI" %in% colnames(colData(spe)))
    { #### NON MI è CHIARO! si può fare questa detection solo sul dapi?
    # stop is inserted because we have DAPI only for CosMx,
    # the function should not continue at this point
        warning("Mean DAPI signal is not in colData. Outlier detection cannot be ",
            "performed on this variable.")
    } else {
        if("Mean.DAPI" %in% colnames(colData(spe)))
            message("Outlier detection on mean DAPI signal")
        if(!skewness(spe$Mean.DAPI)[1]>0)
        {
            warning("Skewness value is not > 0. Outlier detection will be ",
                    "performed using MAD")
            spe$is_area_outlier <- scuttle::isOutlier(spe$Mean.DAPI, type = "both")

        }
        if(!mc(spe$Mean.DAPI)>-0.6 & mc(spe$Mean.DAPI)<0.6)
        {
            warning("Medcouple value is not suitable for boxplot adjustment. ",
                    "Outlier detection will be performed using MAD")
            spe$is_area_outlier <- scuttle::isOutlier(spe$Mean.DAPI, type = "both")
        }
        names(spe$Mean.DAPI) <- colnames(spe)
        out <- adjbox(spe$Mean.DAPI, plot = FALSE)
        message(paste0(length(out$out)," outliers have been found"))
        metadata(spe)$dapi_fence <- out$fence

        dapi_outliers <- vector(length=dim(spe)[2])
        names(dapi_outliers) <- colnames(spe)
        dapi_outliers[which(names(dapi_outliers) %in% names(out$out))] <- TRUE
        spe$is_dapi_outlier <- dapi_outliers
    }
    return(spe)
}

#' computeFilterFlags
#' @description
#' defines flagged cells for variables with fixed thresholds
#'
#' @param spe a SpatialExperiment object with total probe counts, control probe
#' counts on total counts ratio and flag score in the `colData`,
#' tipically computed with \link{spatialPerCellQC}
#'
#' @return
#' @export
#' @examples
#' #TBD
computeFilterFlags <- function(spe, fs_threshold=0.6)
{
    stopifnot(all(c("total", "ctrl_total_ratio", "flag_score") %in%
                colnames(colData(spe))))
    #flagging cells with zero total counts
    spe$is_0counts <- ifelse(spe$total == 0, TRUE, FALSE)
    #flagging cells with probe counts on total counts ratio > 0.1
    spe$is_ctrl_tot_outlier <- ifelse(spe$ctrl_total_ratio > 0.1, TRUE, FALSE)

    spe$is_fscore_outlier <- ifelse(spe$flag_score < fs_threshold, TRUE, FALSE)
    ## add additional flags
    return(spe)
}

