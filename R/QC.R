#' spatialPerCellQC
#' @name spatialPerCellQC
#' @rdname spatialPerCellQC
#' @description
#' Computes quality‐control metrics for each cell and adds them to `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial data.
#' @param micronConvFact Numeric factor to convert pixels to microns. Default
#'   `0.12`.
#' @param rmZeros logical for removing zero counts cells (default is TRUE).
#' @param negProbList Character vector of patterns to identify negative probes.
#'   Defaults include:
#'   - Nanostring CosMx: `"NegPrb"`, `"Negative"`, `"SystemControl"`
#'   - Xenium: `"NegControlProbe"`, `"NegControlCodeword"`,
#'     `"UnassignedCodeword"`
#'   - MERFISH: `"Blank"`
#' @param use_altexps logical for `use_altexps` in `scuttle` package.
#' If TRUE uses the altexps for computing some metrics on it.
#' Useful for interoperability with `SpatialExperimentIO`.
#' (See \link[scuttle]{addPerCellQC} for additional details).
#'
#' @return A `SpatialExperiment` object with added QC metrics in `colData`.
#'
#' @details
#' Calculates sums and detected counts for control and target probes,
#' computes ratio and count‐area metrics, converts coords to microns for
#' CosMx, and drops zero‐count cells.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom scater addPerCellQC
#' @importFrom S4Vectors cbind.DataFrame
#' @export
#' @examples
#' example(readCosmxSPE)
#' spe <- spatialPerCellQC(spe)
spatialPerCellQC <- function(spe, micronConvFact=0.12, rmZeros=TRUE,
    negProbList=c("NegPrb", "Negative", "SystemControl", "Ms IgG1", "Rb IgG",
    "BLANK_", "NegControlProbe", "NegControlCodeword", "UnassignedCodeword",
    "Blank"),
    use_altexps=NULL) {
    stopifnot(is(object=spe, "SpatialExperiment"))
    idxlist <- lapply(negProbList, function(ng) {
        grep(paste0("^", ng), rownames(spe))
    })
    names(idxlist) <- negProbList
    idxlist <- idxlist[which(lengths(idxlist)!=0)]
    spe <- addPerCellQC(spe, subsets=idxlist, use_altexps=use_altexps)
    idx <- grep("^subsets_.*_sum$", colnames(colData(spe)))
    npc <- npd <- 0
    if ( length(idx) !=0 ) {
        npc <- rowSums(as.matrix(colData(spe)[ , idx, drop=FALSE]))
        # TODO: not robust at all! the +1 is not a really good choice
        npd <- rowSums(as.matrix(colData(spe)[ , idx+1, drop=FALSE]))
    }
    spe$control_sum <- npc
    spe$control_detected <- npd
    spe$target_sum <- spe$sum - npc
    spe$target_detected <- spe$detected - npd
    if(!all(spatialCoordsNames(spe) %in% names(colData(spe)))) {
        # TODO: CHANGE SPE constructor WITH COORDINATES IN COLDATA
        colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))
    }

    spe$ctrl_total_ratio <- spe$control_sum/spe$total
    spe$ctrl_total_ratio[which(is.na(spe$ctrl_total_ratio))] <- 0
    spe$log2Ctrl_total_ratio <- log2(spe$ctrl_total_ratio+0.0001)
    if(metadata(spe)$technology == "Nanostring_CosMx_Protein") {
        idx <- which(names(colData(spe)) == "Area.um2")
        if(length(idx)!=0) { names(colData(spe))[idx] <- "Area_um" }
    }

    if(any(metadata(spe)$technology %in%
            c("Nanostring_CosMx", "Nanostring_CosMx_Protein"))) {
        spnc <- spatialCoords(spe) * micronConvFact
        colnames(spnc) <- gsub("px", "um", spatialCoordsNames(spe))
        colData(spe) <- cbind.DataFrame(colData(spe), spnc)
        spe$Area_um <- spe$Area * (micronConvFact^2)
        spe <- .computeBorderDistanceCosMx(spe)
    }

    if (metadata(spe)$technology == "10X_Xenium") {
        spe$Area_um <- spe$cell_area # standardized across other techs
    }
    if ("AspectRatio" %in% colnames(colData(spe))) {
        spe$log2AspectRatio <- log2(spe$AspectRatio) # not cosmx
    } else { warning("Missing aspect ratio in colData") }

    spe$CountArea <- spe$sum/spe$Area_um
    spe$log2CountArea <- log2(spe$CountArea)
    if (rmZeros) {
        if (sum(spe$sum==0) > 0) {
            message("Removing ", dim(spe[,spe$sum==0])[2],
                    " cells with 0 counts!")
            spe <- spe[,!spe$sum==0]
        }
    }
    return(spe)
}

#' .computeBorderDistanceCosMx
#' @name .computeBorderDistanceCosMx
#' @rdname dot-computeBorderDistanceCosMx
#' @description
#' Calculates the minimum distance of each cell to the field‐of‐view border
#' and adds it to `colData`.
#'
#' @param spe A `SpatialExperiment` object with CosMx data.
#' @param xwindim Width of FOV in x (default from `metadata(spe)$fov_dim`).
#' @param ywindim Height of FOV in y (default from `metadata(spe)$fov_dim`).
#'
#' @return A `SpatialExperiment` object with `dist_border` columns in
#' `colData`.
#'
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @keywords internal
.computeBorderDistanceCosMx <- function(spe,
                                    xwindim=metadata(spe)$fov_dim[["xdim"]],
                                    ywindim=metadata(spe)$fov_dim[["ydim"]]) {
    stopifnot(is(spe, "SpatialExperiment"))
    cd <- colData(spe)
    cdf <- left_join(as.data.frame(cd),metadata(spe)$fov_positions,by="fov")
    spcn <- spatialCoordsNames(spe)
    fovpn <- colnames(metadata(spe)$fov_positions)[colnames(
        metadata(spe)$fov_positions) %in% c("x_global_px", "y_global_px")]
    cd$dist_border_vert <- pmin(cdf[,spcn[1]] - cdf[,fovpn[1]],
                            (cdf[,fovpn[1]] + xwindim) - cdf[,spcn[1]])
    cd$dist_border_hor <- pmin(cdf[,spcn[2]] - cdf[,fovpn[2]],
                            (cdf[,fovpn[2]] + ywindim) - cdf[,spcn[2]])
    cd$dist_border <- pmin(cd$dist_border_vert, cd$dist_border_hor)
    colData(spe) <- cd
    return(spe)
}

#' computeSpatialOutlier
#' @name computeSpatialOutlier
#' @rdname computeSpatialOutlier
#' @description
#' Computes outliers based on the Area (in micron) of the experiment.
#' It gives the possibility to choose between the medcouple (mc method argument)
#' and the MADs (scuttle method argument).
#'
#' @details
#' The medcouple method is a measure for the skeweness of univariate
#' distribution as described in Hubert M. et al. (2008).
#' In particular, the computed medcouple value must be in a range between -0.6
#' and 0.6 to computed adjusted boxplots and perform the outlier detection.
#' For median absolute deviations (MADs) method we just wrap the isOutlier
#' function in the scuttle package. Please see McCarthy DJ et al (2017)
#' for further details.
#'
#' @param spe a SpatialExperiment object with target_counts, area in micron
#' and log2 of the aspect ratio in the `colData`.
#' @param computeBy character indicating a `colData` column name on which
#' compute the outlier.
#' @param method one of `mc`, `scuttle`, `both`.
#' Use `mc` for medcouple, `scuttle` for median absolute deviations as computed
#' in `scuttle`, `both` for computing both of them.
#' @param mcDoScale logical indicating if the values to compute the medcouple
#' for the outlier detection should be scaled (default is FALSE, as suggested
#' by the original Medcouple authors.). See \link[robustbase]{mc} for further
#' readings.
#' @param scuttleType One of `"both"`, `"lower"`, `"higher"` for scuttle method.
#'
#' @return a SpatialExperiment object with additional column(s) (named as
#' the column name indicated in `column_by` followed by the outlier_sc/mc
#' nomenclature) with the outlier detection as `outlier.filter` logical class
#' object. This allows to store the thresholds as attributes of the column.
#' use attr(,"thresholds") to retrieve them.
#'
#' @export
#' @importFrom robustbase mc adjbox
#' @importFrom e1071 skewness
#' @importFrom scuttle isOutlier outlier.filter
#'
#' @examples
#' example(spatialPerCellQC)
#' spe <- computeSpatialOutlier(spe, computeBy="log2CountArea", method="both")
#' table(spe$log2CountArea_outlier_mc)
#' table(spe$log2CountArea_outlier_sc)
computeSpatialOutlier <- function(spe, computeBy=NULL,
    method=c("mc", "scuttle", "both"), mcDoScale=FALSE,
    scuttleType=c("both", "lower", "higher")) {
    stopifnot(all(is(spe, "SpatialExperiment"), !is.null(computeBy)))
    stopifnot(computeBy %in% names(colData(spe)))
    options(mc_doScale_quiet=TRUE)
    method <- match.arg(method)
    scuttleType <- match.arg(scuttleType)
    cd <- colData(spe)
    cdcol <- cd[[computeBy]]
    mcfl <- scuttlefl <- FALSE
    switch(method, both={ mcfl <- scuttlefl <- TRUE },
            mc={ mcfl <- TRUE }, scuttle={ scuttlefl <- TRUE },
            {stop("Method is not one of allowed methods")} )
    if (mcfl) {
        skw <- e1071::skewness(cdcol, na.rm = TRUE) # NAs arise problems
        if (skw >- 1 & skw < 1) warning("Distribution is symmetric: ",
                "mc is for asymmetric distributions. Use scuttle instead.")
        mcval <- robustbase::mc(cdcol, doScale=mcDoScale, na.rm=TRUE)
        if ( any( (mcval <= -0.6), (mcval >= 0.6) ) )
            stop("mc is: ",round(mcval, digits=4),"outliers reqs not satisfied")
        names(cdcol) <- colnames(spe)
        outl <- robustbase::adjbox(cdcol, plot=FALSE)
        outsmc <- rep("NO", dim(cd)[1])
        outsmc[rownames(cd) %in% names(outl$out)] <-
            ifelse(outl$out <= outl$fence[1], "LOW", "HIGH")
        outlier_mc <- scuttle::outlier.filter(outsmc) #using scuttle class
        thrs <- as.numeric(outl$fence)
        names(thrs) <- c("lower", "higher")
        attr(outlier_mc, "thresholds") <- thrs
        cd$outlier_mc <- outlier_mc
        names(cd)[names(cd) =="outlier_mc"] <- paste0(computeBy, "_outlier_mc")
        # TODO: compute distributions in the adjusted boxplots to store in cd
    }
    if (scuttlefl) {
        outssc <- scuttle::isOutlier(cdcol, type=scuttleType)
        sctri <- rep("NO", dim(cd)[1])
        sctri <- ifelse(outssc == TRUE & cdcol <= attr(outssc, "thresholds")[1],
                        "LOW", sctri)
        outlier_sc <- ifelse(outssc == TRUE &
                            cdcol >= attr(outssc, "thresholds")[2],
                            "HIGH", sctri)
        outlier_sc <- scuttle::outlier.filter(outlier_sc)
        attr(outlier_sc, "thresholds") <- attr(outssc, "thresholds")
        cd$outlier_sc <- outlier_sc
        names(cd)[names(cd)=="outlier_sc"] <- paste0(computeBy, "_outlier_sc")
    }
    colData(spe) <- cd
    return(spe)
}

#' computeThresholdFlags
#' @name computeThresholdFlags
#' @rdname computeThresholdFlags
#' @description
#' Compute Flagged cells using fixed thresholds for SpatialExperiment.
#'
#' This function calculates flagged cells only for total counts and control on
#' total probe counts ratio using fixed thresholds for a `SpatialExperiment`
#' object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param totalThreshold A numeric value for the threshold of total counts to
#' identify cells with low counts. Default is `0`.
#' @param ctrlTotRatioThreshold A numeric value for the threshold of
#' control-to-total ratio to flag cells over a certain threshold. Default is
#' `0.1`.
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#' @details The function flags cells basing on zero counts and control-to-total
#' ratio to identify junk cells.
#' It also combines these flags into a single filter flag.
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' example(readCosmxSPE)
#' spe <- spatialPerCellQC(spe)
#' spe <- computeThresholdFlags(spe)
#' table(spe$threshold_flags)
computeThresholdFlags <- function(spe, totalThreshold=0,
                            ctrlTotRatioThreshold=0.1)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("total" %in% names(colData(spe)))
    stopifnot("ctrl_total_ratio" %in% names(colData(spe)))

    spe$is_zero_counts <- ifelse(spe$total == totalThreshold, TRUE, FALSE)
    #flagging cells with probe counts on total counts ratio > 0.1
    spe$is_ctrl_tot_outlier <- ifelse(spe$ctrl_total_ratio >
                                        ctrlTotRatioThreshold, TRUE, FALSE)

    spe$threshold_flags <- (spe$is_ctrl_tot_outlier &
                                spe$is_zero_counts)
    return(spe)
}

#' computeOutliersQCScore
#' @name computeOutliersQCScore
#' @rdname computeOutliersQCScore
#' @description
#' Compute outlier cells for each metric that can be used in QC score formula
#' for SpatialExperiment.
#'
#' This function calculates outlier cells for each variable specified in
#' `metric_list` for a `SpatialExperiment`. Log2CountArea must be present in
#' the `colData` of the `SpatialExperiment` object as a minimum requirement.
#' The user can choose which metrics to include among the following: Area_um,
#' log2Ctrl_total_ratio, log2AspectRatio. For Xenium and Merfish datasets,
#' log2AspectRatio is automatically removed from the formula.
#'
#' @param spe A `SpatialExperiment` object with spatial omics data.
#' @param metric_list A character vector specifying the metrics to include in
#' the QC score formula. Default is `c("log2CountArea", "Area_um",
#' "log2AspectRatio", "log2Ctrl_total_ratio")`.
#'
#' @return The `SpatialExperiment` object with added outlier variables in
#' `colData` and the temporary QCScore metric variables that in the
#' `metadata`.
#'
#' @details The function computes outliers for each specified metric after
#' automatically choosing the appropriate method according to the skewness of
#' the distribution.
#' Internally the function:
#' \enumerate{
#'    \item Calls \code{.checkSkw()} to choose the proper outlier detection
#'     method according to the variable skewness,
#'    \item Calls \code{computeSpatialOutlier()} on each included metric to get
#'    fences,
#'    \item Labels cells as “LOW”/“HIGH” outliers or “NO”
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr case_when
#' @importFrom scuttle outlier.filter
#' @importFrom stats quantile
#' @export
#' @examples
#' example(readCosmxSPE)
#' spe <- spatialPerCellQC(spe)
#' spe <- computeOutliersQCScore(spe)
#' table(spe$log2CountArea_outlier_train)

computeOutliersQCScore <- function(spe, metric_list = c("log2CountArea", "Area_um",
                                                        "log2AspectRatio", "log2Ctrl_total_ratio")){

    method <- .checkSkw(spe, metric_list)
    # log2CountArea
    if("log2CountArea" %in% names(method)){
        spe_temp <- computeSpatialOutlier(spe[,spe$total>0],
                                          computeBy="log2CountArea", method=method["log2CountArea"])

        out_var <- colnames(colData(spe_temp))[grep(colnames(colData(spe_temp)), pattern=paste0("log2CountArea_outlier_", method["log2CountArea"]))]

        if(getFencesOutlier(spe_temp, out_var, "lower") <
           min(spe_temp$log2CountArea)) {
            low_thr <- quantile(spe$log2CountArea, probs = 0.01)
        } else {
            low_thr <- getFencesOutlier(spe_temp, out_var,
                                        "lower")
        }

        high_thr <- getFencesOutlier(spe_temp, out_var, "higher")
        spe$log2CountArea_outlier_train <- case_when(spe$total==0 ~ "NO",
                                                     spe$log2CountArea<low_thr ~ "LOW", spe$log2CountArea>high_thr ~ "HIGH",
                                                     TRUE ~ "NO")
        spe$log2CountArea_outlier_train <- scuttle::outlier.filter(spe$log2CountArea_outlier_train)

        attr(spe$log2CountArea_outlier_train, "thresholds") <-
            getFencesOutlier(spe_temp, out_var)
        attr(spe$log2CountArea_outlier_train, "thresholds")[1] <- low_thr
    }

    # log2Ctrl_total_ratio
    if("log2Ctrl_total_ratio" %in% names(method)){
        spe_temp <- computeSpatialOutlier(spe[,spe$ctrl_total_ratio!=0], computeBy="log2Ctrl_total_ratio", method=method["log2Ctrl_total_ratio"])

        out_var <- colnames(colData(spe_temp))[grep(colnames(colData(spe_temp)), pattern=paste0("log2Ctrl_total_ratio_outlier_", method["log2Ctrl_total_ratio"]))]

        spe$log2Ctrl_total_ratio_outlier_train <- case_when(spe$ctrl_total_ratio==0 ~ "NO",
                                                            spe$log2Ctrl_total_ratio<getFencesOutlier(spe_temp, out_var, "lower") ~ "LOW",
                                                            spe$log2Ctrl_total_ratio>getFencesOutlier(spe_temp, out_var, "higher") ~ "HIGH",
                                                            TRUE ~ "NO")

        spe$log2Ctrl_total_ratio_outlier_train <- scuttle::outlier.filter(spe$log2Ctrl_total_ratio_outlier_train)

        attr(spe$log2Ctrl_total_ratio_outlier_train, "thresholds") <- getFencesOutlier(spe_temp, out_var)

    }

    submethod <- method[!names(method)%in%c("log2CountArea", "log2Ctrl_total_ratio")]

    for(j in names(submethod)){
        spe <- computeSpatialOutlier(spe, computeBy=j,method=submethod[j])
    }

    out_var <- paste0(names(method), "_outlier_", method)
    names(out_var) <- names(method)
    # gives warning if one of the variables is missing, but still works!
    out_var[names(out_var)%in%c("log2CountArea", "log2Ctrl_total_ratio")] <- c("log2CountArea_outlier_train", "log2Ctrl_total_ratio_outlier_train")

    metadata(spe)$formula_variables <- out_var

    return(spe)
}

#' .checkSkw
#' @name .checkSkw
#' @rdname dot-checkSkw
#' @description
#' Check skewness of metrics to choose outlier detection method for
#' `SpatialExperiment`.
#'
#' @param spe A `SpatialExperiment` object with spatial omics data.
#' @param metric_list A character vector specifying the metrics to include in
#' the QC score formula. Default is `c("log2CountArea", "Area_um",
#' "log2AspectRatio", "log2Ctrl_total_ratio")`.
#'
#' @return
#' A vector containing the list of chosen outlier detection method for each
#' metric.
#'
#' @examples
#' example(readCosmxSPE)
#' spe <- spatialPerCellQC(spe)
#' .checkSkw(spe, metric_list = c("log2CountArea", "Area_um", "log2AspectRatio",
#' "log2Ctrl_total_ratio")`.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom e1071 skewness
#' @export


.checkSkw <- function(spe = spe, metric_list = metric_list){
    cd <- colData(spe)
    method <- c()
    for(i in metric_list){
        cdcol <- cd[[i]]
        skw <- e1071::skewness(cdcol, na.rm = TRUE)
        method[i] <- ifelse((skw>-1 & skw<1), "sc", "mc")
    }
    if("log2CountArea" %in% names(method)){
        logca_skw <- e1071::skewness(spe[,spe$total>0]$log2CountArea, na.rm = TRUE)
        logca_method <- ifelse((logca_skw>-1 & logca_skw<1), "sc", "mc")
        if(method[names(method)== "log2CountArea"]!= logca_method |
            is.na(method[names(method)== "log2CountArea"])){
            method[names(method)== "log2CountArea"] <- logca_method
        }
    }
    if("log2Ctrl_total_ratio" %in% names(method)){
        logctr_skw <- e1071::skewness(spe[,spe$ctrl_total_ratio!=0]$log2Ctrl_total_ratio, na.rm = TRUE)
        logctr_method <- ifelse((logctr_skw>-1 & logctr_skw<1), "sc", "mc")
        if(method[names(method)== "log2Ctrl_total_ratio"]!= logctr_method |
           is.na(method[names(method)== "log2Ctrl_total_ratio"])){
            method[names(method)== "log2Ctrl_total_ratio"] <- logctr_method
        }
    }
    return(method)
}

#' checkOutliers
#' @name checkOutliers
#' @rdname checkOutliers
#' @description
#' Checks if computed outliers meet the minimum numerical requirement, being
#' at least 0.1% of total cells for each metric to be used in QC score formula.
#' If the requirement is not met, the variable is removed from the formula.
#'
#' @param spe A `SpatialExperiment` object with spatial omics data.
#' @param verbose Logical. If `TRUE`, prints how many outliers were found for
#' each metric.
#'
#' @return The `SpatialExperiment` object with added QCScore metric variables
#'  in the `metadata`.
#'
#' @details The function checks if computed outliers for each metric meet
#' the minimum number to get the metric included in the QC score formula.
#' If verbose is TRUE, it also prints how many outliers were found for each
#' metric.
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' example(computeOutliersQCScore)
#' spe <- checkOutliers(spe, verbose=TRUE)
#' metadata(spe)$formula_variables

checkOutliers <- function(spe, verbose = FALSE){
    out_var <- metadata(spe)$formula_variables
    cd <- colData(spe)
    if(verbose){
        for(i in names(out_var)){
            cd <- colData(spe)
            print(paste0("How many outliers were found for ", i))
            print(table(cd[[out_var[i]]]))
        }
    }
    stopifnot("log2CountArea is not included in the QC score formula. QC score cannot be computed"="log2CountArea"%in%names(out_var))

    if ("log2CountArea"%in%names(out_var)){
        if(table(cd[[out_var[grep(out_var, pattern="log2CountArea_outlier")]]])["LOW"]<dim(spe)[2]*0.001){
            stop("Not enough outlier cells were found for log2CountArea.
      QC score computation cannot be performed")

            out_var <- out_var[-grep(out_var, pattern = "log2CountArea_outlier_train")]
        }
    }
    if ("Area_um"%in%names(out_var)){
        if(table(cd[[out_var[grep(out_var, pattern="Area_um_outlier")]]])["HIGH"]<dim(spe)[2]*0.001){
            warning("Not enough outlier cells were found for Area_um.
      This variable will not be used in the final formula")

            out_var <- out_var[-grep(out_var, pattern = "Area_um_outlier")]
        }
    }
    if ("log2Ctrl_total_ratio"%in%names(out_var)){
        if(table(cd[[out_var[grep(out_var, pattern="log2Ctrl_total_ratio_outlier")]]])["HIGH"]<dim(spe)[2]*0.001){
            warning("Not enough outlier cells were found for log2Ctrl_total_ratio.
      This variable will not be used in the final formula")

            out_var <- out_var[-grep(out_var, pattern = "log2Ctrl_total_ratio_outlier_train")]
        }
    }
    if (metadata(spe)$technology %in% c("Nanostring_CosMx","Nanostring_CosMx_Protein")){
        if("log2AspectRatio"%in%names(out_var)){
            if(table(cd[[out_var[grep(out_var, pattern="log2AspectRatio_outlier")]]])["LOW"]<dim(spe)[2]*0.001 &
               table(cd[[out_var[grep(out_var, pattern="log2AspectRatio_outlier")]]])["HIGH"]<dim(spe)[2]*0.001){
                warning("Not enough outlier cells were found for log2AspectRatio.
      This variable will not be used in the final formula")
                out_var <- out_var[-grep(out_var, pattern = "log2AspectRatio_outlier")]
            }
        }
    }else{
        out_var <- out_var[-grep(out_var, pattern = "log2AspectRatio_outlier")]
    }
    metadata(spe)$formula_variables <- out_var

    return(spe)
}

#' getModelFormula
#' @name getModelFormula
#' @rdname getModelFormula
#' @description
#' Returns the right‐hand side of a model formula string based on formula variables
#' found in the `metadata` of a `SpatialExperiment` object.
#' @param spe A `SpatialExperiment` object with spatial omics data.
#' @param verbose Logical. If `TRUE`, prints the final formula used for QC score
#' @return \[character\]
#'   A one‐sided formula as a string (e.g. "~ log2CountArea + ...").
#' @export
#' @examples
#' example(checkOutliers)
#' getModelFormula(spe, verbose=TRUE)

getModelFormula <- function(spe, verbose = verbose)
{
    out_var <- metadata(spe)$formula_variables
    if("log2AspectRatio"%in%names(out_var)){
        names(out_var)[grep(out_var, pattern = "log2AspectRatio_outlier")] <- "I(abs(log2AspectRatio) * as.numeric(dist_border<50))"
    }
    model_formula <- paste0("~(", paste(names(out_var), collapse = " + "), ")^2", sep = "")

    if(verbose){
        message("Final formula used for QC score computation:")
        print(model_formula)
    }

    return(model_formula)
}

#' computeTrainDF
#' @name computeTrainDF
#' @rdname computeTrainDF
#' @description
#' Build a Balanced Training Data Frame from a SpatialExperiment
#'
#' \code{computeTrainDF} takes a \linkS4class{SpatialExperiment} object
#' and assembles a balanced training set of “good” vs “bad” cells for
#' subsequent model fitting.
#'
#' @param spe \linkS4class{SpatialExperiment}
#'   A SpatialExperiment containing at least:
#'   \itemize{
#'     \item assay(s) with nonzero \code{total} counts,
#'     \item \code{colData(spe)} columns including \code{log2CountArea},
#'     \code{Area_um}, \code{log2Ctrl_total_ratio}, etc.
#'   }
#'
#' @param verbose \[logical(1)\] (default \code{FALSE})
#'   If \code{TRUE}, prints the number of “bad” and “good” cells selected.
#'
#' @return
#' A \code{data.frame} with one row per cell, including:
#' \itemize{
#'   \item \code{qcscore_train} (0/1) indicating “bad” vs “good”,
#'   \item relevant \code{colData} columns used for modeling.
#'   \item Deduplicates and down-samples “good” cells to match the number of
#'    “bad” cells.
#' }
#'
#' @details The function builds a training set using the variables specified
#' in the `metadata` of the `SpatialExperiment` object.
#'
#' @examples
#' example(spatialPerCellQC)
#' df_train <- computeTrainDF(spe, verbose = TRUE)
#' table(df_train$qcscore_train)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter mutate distinct pull
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats as.formula model.matrix quantile predict
#'
#' @export

computeTrainDF <- function(spe, verbose = TRUE){
    out_var <- metadata(spe)$formula_variables

    train_bad <- data.frame(colData(spe))
    train_good <- data.frame(colData(spe))

    train_bad_var <- c()
    train_good_var <- c()

    stopifnot("log2CountArea is not included in the QC score formula. QC score cannot be computed"="log2CountArea"%in%names(out_var))

    if ("log2CountArea"%in%names(out_var)){
        train_bad_temp <- train_bad |> filter(log2CountArea_outlier_train == "LOW") |> dplyr::pull(cell_id)

        train_good_temp <- train_good |> filter((log2CountArea > quantile(log2CountArea, probs = 0.90) &
                                                     log2CountArea < quantile(log2CountArea, probs = 0.99))) |> dplyr::pull(cell_id)

        train_bad_var <- unique(c(train_bad_var, train_bad_temp))
        train_good_var <- unique(c(train_good_var, train_good_temp))
    }
    if ("Area_um"%in%names(out_var)){
        train_bad_temp <- train_bad |> filter(.data[[out_var[names(out_var)=="Area_um"]]] == "HIGH") |> dplyr::pull(cell_id)

        train_good_temp <- train_good |> filter(Area_um > quantile(Area_um, probs = 0.25) &
                                                    Area_um < quantile(Area_um, probs = 0.75)) |> dplyr::pull(cell_id)

        train_bad_var <- unique(c(train_bad_var, train_bad_temp))
        train_good_var <- unique(c(train_good_var, train_good_temp))
    }

    if ("log2Ctrl_total_ratio"%in%names(out_var)){
        train_bad_temp <- train_bad |> filter(.data[[out_var[names(out_var)=="log2Ctrl_total_ratio"]]] == "HIGH") |> dplyr::pull(cell_id)

        train_bad_var <- unique(c(train_bad_var, train_bad_temp))
    }

    if (metadata(spe)$technology %in% c("Nanostring_CosMx", "Nanostring_CosMx_Protein") &
        "log2AspectRatio"%in%names(out_var)){
        train_bad_temp <- train_bad |> filter((.data[[out_var[names(out_var)=="log2AspectRatio"]]] ==
                                                   "HIGH" & dist_border < 50) |
                                                  (.data[[out_var[names(out_var)=="log2AspectRatio"]]] ==
                                                       "LOW" & dist_border < 50)) |> dplyr::pull(cell_id)

        train_good_temp <- train_good |> filter(log2AspectRatio > quantile(log2AspectRatio, probs = 0.25) &
                                                    log2AspectRatio < quantile(log2AspectRatio, probs = 0.75) & dist_border > 50) |>
            dplyr::pull(cell_id)

        names(out_var)[grep(out_var, pattern = "log2AspectRatio_outlier")] <- "I(abs(log2AspectRatio) * as.numeric(dist_border<50))"

        train_bad_var <- unique(c(train_bad_var, train_bad_temp))
        train_good_var <- unique(c(train_good_var, train_good_temp))
    }

    train_bad <- train_bad |> filter(cell_id%in%train_bad_var) |> mutate(qcscore_train = 0)
    train_good <- train_good |> filter(cell_id%in%train_good_var) |> mutate(qcscore_train=1, is_a_bad_boy=cell_id%in%train_bad$cell_id)

    train_bad <- train_bad |> distinct(cell_id, .keep_all = TRUE)

    message(paste0("Chosen low quality examples: ", dim(train_bad)[1]))

    # good example duplicates removal without any warning to the user

    train_good <- train_good |> distinct(cell_id, .keep_all = TRUE)

    train_good <- train_good[!train_good$is_a_bad_boy,]
    train_good <- train_good[sample(rownames(train_good), dim(train_bad)[1], replace = FALSE),]

    if(verbose){
        message(paste0("Chosen good quality examples, (should be the same number
              of bad quality examples): ", dim(train_good)[1]))
    }

    train_good$is_a_bad_boy <- NULL
    train_df <- rbind(train_bad, train_good)

    train_df <- train_df |> distinct(cell_id, .keep_all = TRUE)

    return(train_df)

}

#' trainModel
#' @name trainModel
#' @rdname trainModel
#' @description
#' Fit a Ridge Logistic Regression Model
#'
#' \code{trainModel} fits an L2-regularized (ridge) logistic regression
#' using \pkg{glmnet}, given a design matrix and a training data frame.
#'
#' @param model_matrix \[matrix\]
#'   The design matrix of predictors (e.g. from \code{model.matrix()}).
#'
#' @param train_df \[data.frame\]
#'   A data frame containing at least the response column
#'   \code{qcscore_train}, coded as 0/1.
#'
#' @return
#' A \code{\link[glmnet]{glmnet}} model object fitted with
#' \code{family="binomial"}, \code{alpha=0} (ridge), and a sequence of
#' \eqn{\lambda} values.
#'
#' @examples
#' example(computeTrainDF)
#' model_formula <- getModelFormula(metadata(spe)$technology)
#' model_matrix <- model.matrix(as.formula(model_formula), data=df_train)
#' fit <- trainModel(model_matrix, df_train)
#' coef(fit, s = 0.01)
#'
#' @export

trainModel <- function(model_matrix, train_df){
    model <- glmnet(x=model_matrix, y=train_df$qcscore_train,
                    family="binomial", lambda=NULL, alpha=0)
    return(model)
}

#' computeLambda
#' @description
#' Compute Optimal Ridge Regularization Parameter \eqn{\lambda} via
#' Cross-Validation
#'
#' \code{computeLambda} performs ridge (L2) logistic regression with
#' cross-validation to identify the optimal regularization parameter
#' \eqn{\lambda} for a binary response.
#'
#' @param spe \linkS4class{SpatialExperiment}
#'   A SpatialExperiment containing at least:
#'   \itemize{
#'     \item assay(s) with nonzero \code{total} counts,
#'     \item \code{colData(spe)} columns including \code{log2CountArea},
#'     \code{Area_um}, \code{log2Ctrl_total_ratio}, etc.
#'   }
#'
#' @param train_df  \[data.frame\]
#'   A data frame for training that must include:
#'   \describe{
#'     \item{Predictor columns}{All columns referenced in the formula returned
#'     by \code{getModelFormula()}.}
#'     \item{\code{qcscore_train}}{A binary (0/1) response vector to be modeled.}
#'   }
#'
#' @param model_formula \[character\]
#'   A character string representing the model formula
#'   \describe{
#'    "\code{~ log2CountArea + ...}"), as returned by
#'   \code{getModelFormula()}.
#'   }
#'
#' @return
#' \[numeric\]
#'   The value of \eqn{\lambda} (i.e., \code{lambda.min}) from
#'   \code{\link[glmnet]{cv.glmnet}} that minimizes the cross-validation error.
#'
#' @details
#' Internally, the function:
#' \enumerate{
#'   \item Constructs the design matrix via \code{model.matrix()},
#'   \item Runs ridge logistic regression cross-validation using
#'         \code{\link[glmnet]{cv.glmnet}} with \code{alpha = 0},
#'   \item Extracts and returns \code{ridge_cv$lambda.min}.
#' }
#'
#' @examples
#' example(spatialPerCellQC)
#' withr::with_seed(1998, train_df <- computeTrainDF(spe))
#' best_lambda <- computeLambda(spe, train_df, model_formula)
#' print(best_lambda)
#'
#' @seealso
#' \code{\link[glmnet]{cv.glmnet}}
#'
#' @export

computeLambda <- function(spe, train_df, model_formula) {
    model_matrix <- model.matrix(as.formula(model_formula), data=train_df)
    ridge_cv <- cv.glmnet(model_matrix, train_df$qcscore_train,
                          family="binomial", alpha=0, lambda=NULL)
    best_lambda <- ridge_cv$lambda.min
    return(best_lambda)
}

#' computeQCScore
#' @name computeQCScore
#' @rdname computeQCScore
#' @description
#' Compute QC score and automatically define weights for QC score
#' through glm training. This function computes QC score with a formula
#' that is defined based on the metrics specified in metric_list and on the
#' number of available outliers for each metric.
#'
#' @details
#' For CosMx datasets, also CosMx Protein, the QC Score formula is
#' defined as follows:
#'
#' QC score ~ count density - aspect ratio - control-total ratio
#'
#' count density is total counts-to-area ratio, aspect ratio represents
#' FOV border effect typical of CosMx datasets and control-total ratio is
#' the aspecific signal. For each couple of variables interaction terms are
#' computed.
#'
#' For Xenium and Merscope datasets, QC score cannot depend on aspect ratio
#' as no FOV border effect was captured through this metric.
#'
#' Inclusion of metrics in the formula depends also on the number of available
#' outliers. If the number of outliers for each metric is < 0.1% out of the
#' entire dataset, the metric will be excluded from the QC score formula.
#'
#' To automatically define the formula coefficient weights, model training
#' is performed through ridge regression.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param metric_list A character vector containing the list of metrics to compute
#' QC score on. log2CountArea must be always included.
#' @param verbose logical for having a verbose output. Default is FALSE.
#' @param best_lambda the best lambda typically computed using `computeLambda`.
#'
#' @return The `SpatialExperiment` object with added QC score in `colData`.
#' @export
#' @importFrom dplyr case_when filter mutate distinct pull
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats as.formula model.matrix quantile predict
#' @examples
#' example(spatialPerCellQC)
#' set.seed(1998)
#' spe <- computeQCScore(spe)
#' summary(spe$training_status)
#' summary(spe$QC_score)

computeQCScore <- function(spe, metric_list = c("log2CountArea", "Area_um",
                                                "log2AspectRatio", "log2Ctrl_total_ratio"), best_lambda=NULL, verbose=TRUE) {
    stopifnot(is(spe, "SpatialExperiment"))
    if(dim(spe[,spe$total==0])[2]!=0){
        warning(paste0(dim(spe[,spe$total==0])[2],
        " cells with 0 counts were found. These cells will be removed."))
        spe <- spe[,spe$total>0]
    }
    spe <- computeOutliersQCScore(spe, metric_list = metric_list)
    spe <- checkOutliers(spe, verbose)
    train_df <- computeTrainDF(spe, verbose)
    model_formula <- getModelFormula(spe, verbose)
    model_matrix <- model.matrix(as.formula(model_formula), data=train_df)
    model <- trainModel(model_matrix, train_df)
    if(is.null(best_lambda)) {
        best_lambda <- computeLambda(spe, train_df, model_formula)
    }

    if (verbose){
        message("Model coefficients for every term used in the formula:")
        print(round(predict(model, s = best_lambda, type="coefficients"),2))
    }
    cd <- data.frame(colData(spe))
    full_matrix <- model.matrix(as.formula(model_formula), data = cd)
    cd$QC_score <- as.vector(predict(model, s=best_lambda,
                                     newx = full_matrix,
                                     type = "response"))
    spe$QC_score <- cd$QC_score
    train_identity <- rep("TEST", dim(spe)[2])
    train_bad <- train_df$cell_id[train_df$qcscore_train==0]
    train_good <- train_df$cell_id[train_df$qcscore_train==1]
    spe$training_status <- dplyr::case_when(
        spe$cell_id %in% train_bad ~ "BAD",
        spe$cell_id %in% train_good ~ "GOOD",
        TRUE ~ train_identity)
    return(spe)
}


#' computeQCScoreFlags
#' @name computeQCScoreFlags
#' @rdname computeQCScoreFlags
#' @description
#' Compute flagged cells based on a manually chosen threshold on QC score
#'
#' This function Compute flagged cells based on a manually chosen threshold on
#' QC score stored in `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param qs_threshold Numeric threshold or quantile for QC score. Default
#'   `0.5`.
#' @param use_qs_quantiles Logical; if `TRUE`, treat `qs_threshold` as a
#'   percentile.
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' example(computeQCScore)
#' spe <- computeQCScoreFlags(spe)
#' table(spe$is_qcscore_outlier)
#' # if fixed filters are defined we have an additional column
#' spe <- computeThresholdFlags(spe)
#' spe <- computeQCScoreFlags(spe)
#' table(spe$threshold_qcscore_flags)
computeQCScoreFlags <- function(spe, qs_threshold=0.5, use_qs_quantiles=FALSE) {
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("QC_score" %in% names(colData(spe)))

    if(use_qs_quantiles) {
        spe$is_qcscore_flags <- ifelse(
            spe$QC_score < quantile(spe$QC_score, probs=qs_threshold),
            TRUE, FALSE)
    } else {
        spe$is_qcscore_flags <- spe$QC_score < qs_threshold

    }

    if("threshold_flags" %in% names(colData(spe))) {
        spe$threshold_qcscore_flags <- (spe$is_qcscore_flags &
                                            spe$threshold_flags)
    }
    return(spe)
}
