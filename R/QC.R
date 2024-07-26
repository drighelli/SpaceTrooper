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
computeQCScore <- function(spe, a=0.3, b=0.8, threshold=0.15)
{
    stopifnot(is(spe, "SpatialExperiment"))
    cd <- colData(spe)
    if(!all(c("target_sum", "Area_um", "log2AspectRatio") %in% colnames(cd)))
    {
        stop("One of target_sum, Area_um, log2AspectRatio is missing, did you run spatialPerCellQC?")
    }
    spe$flag_score <- 1/(1 + exp(-a*spe$target_sum/spe$Area_um + b*abs(spe$log2AspectRatio)))
    spe$is_fscore_outlier <- ifelse(spe$flag_score < quantile(spe$flag_score, probs=threshold), TRUE, FALSE)

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
#' @param method
#' @param mcDoScale logical indicating if the values to compute the medcouple
#' for the outlier detection should be scaled (default is FALSE, as suggested
#' by the original Medcouple authors.). See \link[robustbase]{mc} for further
#' readings.
#'
#' @return a SpatialExperiment object with "column name" with the outlier detection...
#' @export
#' @importFrom robustbase mc
#' @importFrom e1071 skewness adjbox
#' @importFrom scuttle isOutlier
#'
#' @examples
#' TBD
computeSpatialOutlier <- function(spe, method=c("mc", "scuttle", "both"),
                                mcDoScale=FALSE,
                                scuttleType=c("both", "lower", "higher"))
{
    stopifnot(is(spe, "SpatialExperiment"))
    method <- match.arg(method)
    scuttleType <- match.arg(scuttleType)
    cd <- colData(spe)

    mcfl=scuttlefl=FALSE
    switch(method,
           both={ mcfl=scuttlefl=TRUE },
           mc={ mcfl=TRUE },
           scuttle={ scuttlefl=TRUE },
           {stop("Method is not one of allowed methods")}
    )

    if (mcfl)
    {
        skw <- e1071::skewness(cd$Area_um)
        if (skw<=0) warning("The distribution is symmetric. ",
                    "The medcouple is not suited for symmetric distributions. ",
                    "In your case we suggest to use the scater method instead.")

        mcval <- robustbase::mc(cd$Area_um, doScale=mcDoScale)
        if ( any( (mcval <= -0.6), (mcval >= 0.6) ) )
            stop("Obtained medcouple value is: ", round(mcval, digits=4),
                 "\nIt doesn't meet the needed requirements for outlier",
                 " identification with this method.")
        names(cd$Area_um) <- colnames(spe)
        outl <- robustbase::adjbox(cd$Area_um, plot=FALSE)$out
        outsmc <- rep(FALSE, dim(cd)[1])
        outsmc[rownames(cd) %in% names(outl)] <- TRUE
        ## compute distributions in the adjusted boxplots to store them in the
        ## coldata
    }

    if (scuttlefl) { outssc <- scuttle::isOutlier(cd$Area_um, type=scuttleType) }

    if (method=="both")
    {
        cd$umarea_outlier_mc <- outsmc
        cd$umarea_outlier_sc <- outssc
    } else {
        cd$umarea_outlier <- ifelse(mcfl, outsmc, outssc)
    }

    colData(spe) <- cd
    return(spe)


    # area_fence <- out$fence
    # outlier_color <- colnames(spe)
    # for(i in 1:out$n){
    #     outlier_color[i] <- ifelse(outlier_color[i] %in% names(out$out), "red", "grey80")
    # }
    # spe$umarea_outlier <- ifelse(outlier_color == "red", TRUE, FALSE)
    # spe$polygons$umarea_outlier <- outlier_color

    ## scater

}

#
# computeFilterFlags <- function(spe, fs_threshold=0.15)
# {
#     stopifnot(is(spe, "SpatialExperiment"))
#     stopifnot(all(c("flag_score") %in% colnames(colData(spe))))
#     spe$is_fscore_outlier <- ifelse(spe$flag_score < quantile(spe$flag_score, probs=fs_threshold), TRUE, FALSE)
#     ## add additional flags
#     return(spe)
# }
