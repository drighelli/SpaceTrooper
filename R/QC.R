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
#' #TBD
computeQCScore <- function(spe, a=0.3, b=0.8)
{
    stopifnot(is(spe, "SpatialExperiment"))
    cd <- colData(spe)
    if(!all(c("target_sum", "Area_um", "log2AspectRatio") %in% colnames(cd)))
    {
        stop("One of target_sum, Area_um, log2AspectRatio is missing, ",
        "did you run spatialPerCellQC?")
    }
    spe$flag_score <- 1/(1 + exp(-a*spe$target_sum/spe$Area_um + b*abs(spe$log2AspectRatio)))

    return(spe)
}


#' computeSpatialOutlier
#' @description
#' Computes outliers based on the Area (in micron) of the experiment.
#' It gives the possibility to choose between the medcouple (mc method argument)
#' and the MADs (scuttle method argument).
#'
#' @details
#' The medcouple method is a measure for the skeweness of univariate distribution
#' as described in Hubert M. et al. (2008).
#' In particular, the computed medcouple value must be in a range between -0.6
#' and 0.6 to computed adjusted boxplots and perform the outlier detection.
#' For median absolute deviations (MADs) method we just wrap the isOutlier
#' function in the scuttle package. Please see McCarthy DJ et al (2017)
#' for further details.
#'
#' @param spe a SpatialExperiment object with target_counts, area in micron
#' and log2 of the aspect ratio in the `colData`
#' @param compute_by
#' @param method one of `mc`, `scuttle`, `both`.
#' Use `mc` for medcouple, `scuttle` for median absolute deviations as computed
#' in `scuttle`, `both` for computing both of them.
#' @param mcDoScale logical indicating if the values to compute the medcouple
#' for the outlier detection should be scaled (default is FALSE, as suggested
#' by the original Medcouple authors.). See \link[robustbase]{mc} for further
#' readings.
#'
#' @return a SpatialExperiment object with additional column(s) (named as
#' the column name indicated in `column_by` followed by the outlier_sc/mc
#' nomenclature) with the outlier detection as `outlier.filter` logical class
#' object. This allows to store the thresholds as attributes of the column.
#' use attr(,"thresholds") to retrieve them.
#' @export
#' @importFrom robustbase mc adjbox
#' @importFrom e1071 skewness
#' @importFrom scuttle isOutlier outlier.filter
#'
#' @examples
#' TBD
computeSpatialOutlier <- function(spe, compute_by=NULL,
                                method=c("mc", "scuttle", "both"),
                                mcDoScale=FALSE,
                                scuttleType=c("both", "lower", "higher"))
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(!is.null(compute_by))
    stopifnot(compute_by %in% names(colData(spe)))

    options(mc_doScale_quiet=TRUE)

    method <- match.arg(method)
    scuttleType <- match.arg(scuttleType)
    cd <- colData(spe)
    cdcol <- cd[[compute_by]]
    mcfl=scuttlefl=FALSE
    switch(method,
           both={ mcfl=scuttlefl=TRUE },
           mc={ mcfl=TRUE },
           scuttle={ scuttlefl=TRUE },
           {stop("Method is not one of allowed methods")}
    )

    if (mcfl)
    {
        skw <- e1071::skewness(cdcol)
        if (skw<=0) warning("The distribution is symmetric. ",
                    "The medcouple is not suited for symmetric distributions. ",
                    "In your case we suggest to use the scater method instead.")

        mcval <- robustbase::mc(cdcol, doScale=mcDoScale)
        if ( any( (mcval <= -0.6), (mcval >= 0.6) ) )
            stop("Obtained medcouple value is: ", round(mcval, digits=4),
                 "\nIt doesn't meet the needed requirements for outlier",
                 " identification with this method.")
        names(cdcol) <- colnames(spe)
        outl <- robustbase::adjbox(cdcol, plot=FALSE)
        outsmc <- rep(FALSE, dim(cd)[1])
        outsmc[rownames(cd) %in% names(outl$out)] <- TRUE
        ## using scuttle outlier.filter class to store fences for each filtering
        ## in the attributes of the class
        outlier_mc <- scuttle::outlier.filter(outsmc)
        thrs <- as.numeric(outl$fence)
        names(thrs) <- c("lower", "higher")
        attr(outlier_mc, "thresholds") <- thrs
        cd$outlier_mc <- outlier_mc
        names(cd)[names(cd)=="outlier_mc"] <- paste0(compute_by, "_outlier_mc")
        # metadata(spe)$outlier_fences[[compute_by]] <- outl$fence ## DEPRECATED

        ## TODO: compute distributions in the adjusted boxplots to store them
        ## in the coldata
    }

    if (scuttlefl)
    {
        outssc <- scuttle::isOutlier(cdcol, type=scuttleType)
        cd$outlier_sc <- scuttle::outlier.filter(outssc)
        names(cd)[names(cd)=="outlier_sc"] <- paste0(compute_by, "_outlier_sc")
    }

    colData(spe) <- cd
    return(spe)
}

#' computeFilterFlags
#' @description
#' defines flagged cells for variables based on the thresholds passed as input
#'
#' @param spe a SpatialExperiment object with total probe counts, control probe
#' counts on total counts ratio and flag score in the `colData`,
#' tipically computed with \link{spatialPerCellQC}
#' @param fs_threshold
#' @param use_fs_quantiles
#' @param total_threshold
#' @param ctrl_tot_ratio_threshold
#'
#' @return
#' @export
#' @examples
computeFilterFlags <- function(spe, fs_threshold=0.5,
                        use_fs_quantiles=FALSE,
                        total_threshold=0, ctrl_tot_ratio_threshold=0.1)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("flag_score" %in% names(colData(spe)))
    stopifnot("total" %in% names(colData(spe)))
    stopifnot("ctrl_total_ratio" %in% names(colData(spe)))

    spe$is_zero_counts <- ifelse(spe$total == total_threshold, TRUE, FALSE)
    #flagging cells with probe counts on total counts ratio > 0.1
    spe$is_ctrl_tot_outlier <- ifelse(spe$ctrl_total_ratio >
                                ctrl_tot_ratio_threshold, TRUE, FALSE)
    if(!use_fs_quantiles)
    {
        spe$is_fscore_outlier <- ifelse(spe$flag_score <
                                quantile(spe$flag_score, probs=fs_threshold),
                                TRUE, FALSE)
    } else {
        spe$is_fscore_outlier <- ifelse(spe$flag_score < fs_threshold,
                                TRUE, FALSE)
    }

    spe$filter_out <- (spe$is_fscore_outlier & spe$is_ctrl_tot_outlier &
                           spe$is_zero_counts)

    return(spe)
}
