#' spatialPerCellQC
#' @name spatialPerCellQC
#' @rdname spatialPerCellQC
#' @description
#' Computes quality‐control metrics for each cell and adds them to `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial data.
#' @param micronConvFact Numeric factor to convert pixels to microns. Default
#'    `0.12`.
#' @param rmZeros logical for removing zero counts cells (default is TRUE).
#' @param negProbList Character vector of patterns to identify negative probes.
#' Defaults include:
#'   Nanostring CosMx: `"NegPrb"`, `"Negative"`, `"SystemControl"`
#'   Xenium: `"NegControlProbe"`, `"NegControlCodeword"`, `"UnassignedCodeword"`
#'   MERFISH: `"Blank"`
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
    # using eps avoids NaNs and 0 in further computations
    eps <- 0.0001
    spe$log2Ctrl_total_ratio <- log2(spe$ctrl_total_ratio+eps)
    if(metadata(spe)$technology == "Nanostring_CosMx_Protein") {
        # Only for proteins will be included in QCScore
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
    cdf <- left_join(as.data.frame(cd), metadata(spe)$fov_positions, by="fov")
    spcn <- spatialCoordsNames(spe)
    fovpn <- colnames(metadata(spe)$fov_positions)[colnames(
        metadata(spe)$fov_positions) %in% c("x_global_px", "y_global_px")]
    cd$dist_border_x <- pmin(cdf[,spcn[1]] - cdf[,fovpn[1]],
                            (cdf[,fovpn[1]] + xwindim) - cdf[,spcn[1]])
    cd$dist_border_y <- pmin(cdf[,spcn[2]] - cdf[,fovpn[2]],
                            (cdf[,fovpn[2]] + ywindim) - cdf[,spcn[2]])
    cd$dist_border <- pmin(cd$dist_border_x, cd$dist_border_y)
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


#' computeLambda
#' @description
#' Compute Optimal Ridge Regularization Parameter \eqn{\lambda} via
#' Cross-Validation
#'
#' \code{computeLambda} performs ridge (L2) logistic regression with
#' cross-validation to identify the optimal regularization parameter
#' \eqn{\lambda} for a binary response.
#'
#' @param trainDF  `data.frame`
#'   A data frame for training that must include:
#'     Predictor columns: All columns referenced in the formula returned by `getModelFormula()`.
#'     `qscore_train` A binary (0/1) response vector to be modeled.
#' @param modelFormula `character`
#'   A character string representing the model formula
#'    `~ log2CountArea + ...`, as returned by `getModelFormula()`.
#'
#' @return
#' `numeric`
#'   The value of \eqn{\lambda} (i.e., `lambda.min`) from
#'   `cv.glmnet` that minimizes the cross-validation error.
#'
#' @details
#' Internally, the function:
#'   Constructs the design matrix via \code{model.matrix()},
#'   Runs ridge logistic regression cross-validation using `cv.glmnet` with `alpha = 0`,
#'   Extracts and returns `ridge_cv$lambda.min`.
#'
#' @examples
#' example(computeTrainDF)
#' modform <- getModelFormula(metadata(spe)$formula_variables)
#' best_lambda <- computeLambda(df_train, modform)
#' print(best_lambda)
#'
#'
#' @export
computeLambda <- function(trainDF, modelFormula) {
    # model_formula <- getModelFormula(technology)
    model_matrix <- model.matrix(as.formula(modelFormula), data=trainDF)
    ridge_cv <- cv.glmnet(model_matrix, trainDF$qcscore_train,
                        family="binomial", alpha=0, lambda=NULL)
    bestLambda <- ridge_cv$lambda.min
    return(bestLambda)
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
#' @param verbose logical for having a verbose output. Default is FALSE.
#' @param bestLambda the best lambda typically computed using `computeLambda`.
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
#' summary(spe$QC_score)
computeQCScore <- function(spe, bestLambda=NULL, verbose=FALSE) {
    stopifnot(is(spe, "SpatialExperiment"))
    if (dim(spe[,spe$total == 0])[2] != 0) {
        warning(paste0(dim(spe[,spe$total == 0])[2],
            " cells with 0 counts were found. These cells will be removed."))
        spe <- spe[,spe$total > 0]
    }
    metricList <- c("log2CountArea", "Area_um",
                    "log2AspectRatio", "log2Ctrl_total_ratio")

    stopifnot("Not all required metrics in the colData.\nPlease run spatialPerCellQC first." = all(metricList %in% names(colData(spe))))
    ctx <- .prepQCContext(spe, metricList, verbose)
    df <- ctx$df; out_var <- ctx$out_var; tech <- ctx$tech

    train_df <- computeTrainDF(df, out_var, tech, verbose)

    model_formula <- getModelFormula(out_var, verbose)

    model_matrix <- model.matrix(as.formula(model_formula), data=train_df)
    model <- trainModel(model_matrix, train_df)
    if(is.null(bestLambda)) {
        bestLambda <- computeLambda(train_df, model_formula)
    }

    if (verbose) {
        message("Model coefficients for every term used in the formula:")
        message(round(predict(model, s=bestLambda, type="coefficients"),2))
    }

    full_matrix <- model.matrix(as.formula(model_formula), data=df)
    cd <- colData(spe)
    cd$QC_score <- as.vector(predict(model, s=bestLambda,
                                    newx = full_matrix,
                                    type = "response"))
    spe$QC_score <- cd$QC_score
    # train_identity <- rep("TEST", dim(spe)[2])
    # train_bad <- train_df$cell_id[train_df$qcscore_train==0]
    # train_good <- train_df$cell_id[train_df$qcscore_train==1]
    # spe$training_status <- dplyr::case_when(
    #     spe$cell_id %in% train_bad ~ "BAD",
    #     spe$cell_id %in% train_good ~ "GOOD",
    #     TRUE ~ train_identity)
    return(spe)
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
#' @param trainDF `data.frame`
#'   A data frame containing at least the response column
#'   `qscore_train`, coded as 0/1.
#' @param modelMatrix a matrix describing the model variables, tipically created
#' with `getModelFormula` and `model.matrix` functions.
#'
#' @return
#' A \code{\link[glmnet]{glmnet}} model object fitted with
#' \code{family="binomial"}, \code{alpha=0} (ridge), and a sequence of
#' \eqn{\lambda} values.
#'
#' @export
#' @examples
#' example(computeTrainDF)
#' model_formula <- getModelFormula(metadata(spe)$formula_variables)
#' model_matrix <- model.matrix(as.formula(model_formula), data=df_train)
#' fit <- trainModel(model_matrix, df_train)
#' coef(fit, s = 0.01)
trainModel <- function(modelMatrix, trainDF)
{
    model <- glmnet(x=modelMatrix, y=trainDF$qcscore_train,
                    family="binomial", lambda=NULL, alpha=0)
    return(model)
}

#' computeTrainDF
#'
#' @rdname computeTrainDF
#'
#' @description
#' Build a Balanced Training Data Frame from a SpatialExperiment
#'
#' \code{computeTrainDF} takes a \code{SpatialExperiment} object
#' and assembles a balanced training set of “good” vs “bad” cells for
#' subsequent model fitting.
#' @param colData A per-cell metadata table. Typically
#'   `as.data.frame(colData(spe))`. Must include at least:
#'   `cell_id`, raw metric columns named in `formulaVars` (e.g.
#'   `log2CountArea`, `Area_um`, `log2Ctrl_total_ratio`,
#'   optionally `log2AspectRatio`), and the corresponding outlier-label
#'   columns referenced by `formulaVars`.
#' @param formulaVars A named character vector mapping variable name to
#'   its outlier label column name, e.g.
#'   `c(log2CountArea="log2CountArea_outlier_train", ...)`.
#' @param tech Character string with the acquisition technology. Used to
#'   enable CosMx-specific handling for `log2AspectRatio`. Expected
#'   values include `"Nanostring_CosMx"` or `"Nanostring_CosMx_Protein"`.
#' @param verbose Logical; print progress messages.
#'
#' @param verbose `logical(1)` (default \code{FALSE})
#'   If \code{TRUE}, prints the number of “bad” and “good” cells selected.
#'
#' @return
#' A \code{data.frame} with one row per cell, including:
#'   \code{qcscore_train} (0/1) indicating “bad” vs “good”,
#'   relevant \code{colData} columns used for modeling.
#'   Deduplicates and down-samples “good” cells to match the number of “bad” cells.
#'
#' @details The function builds a training set using the variables specified
#' in the `metadata` of the `SpatialExperiment` object.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter mutate distinct pull
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats as.formula model.matrix quantile predict
#'
#' @examples
#' example(spatialPerCellQC)
#' spe <- computeOutliersQCScore(spe)
#' spe <- checkOutliers(spe)
#' df_train <- computeTrainDF(colData(spe), metadata(spe)$formula_variables,
#'     metadata(spe)$technology)
#' table(df_train$qcscore_train)
#'
#' @export
computeTrainDF <- function(colData, formulaVars, tech, verbose=FALSE) {
    out_var <- formulaVars
    df <- as.data.frame(colData)
    train_bad_var <- character()
    train_good_var <- character()

    stopifnot("log2CountArea is not included in the QC score formula.\nQC score cannot be computed"="log2CountArea" %in% names(out_var))

    cfg <- list(
        log2CountArea=list(bad="LOW", good=c(0.90,0.99)),
        Area_um=list(bad="HIGH", good=c(0.25,0.75)),
        log2Ctrl_total_ratio=list(bad="HIGH", good=NULL)
    )

    vars <- intersect(names(cfg), names(out_var))

    for (v in vars) {
        r <- cfg[[v]]
        out_col <- out_var[names(out_var) == v]
        if (length(out_col) == 0L) next

        # BAD ids
        bad_ids <- df |>
            dplyr::filter(.data[[out_col]] == r$bad) |>
            dplyr::pull(cell_id)
        train_bad_var <- unique(c(train_bad_var, bad_ids))

        # GOOD ids (if a band is defined)
        if (!is.null(r$good)) {
            qs <- stats::quantile(df[[v]], probs=r$good)
            good_ids <- df |>
                dplyr::filter(.data[[v]] > qs[1] & .data[[v]] < qs[2]) |>
                dplyr::pull(cell_id)
            train_good_var <- unique(c(train_good_var, good_ids))
        }
    }

    # CosMx-specific handling for log2AspectRatio
    is_cosmx <- tech %in% c("Nanostring_CosMx", "Nanostring_CosMx_Protein")

    if (is_cosmx && "log2AspectRatio" %in% names(out_var)) {

        if (!("dist_border" %in% colnames(df))) {
            stop("dist_border column is required for CosMx handling")
        }

        out_col <- out_var[names(out_var) == "log2AspectRatio"]

        bad_ids <- df |>
            dplyr::filter(.data[[out_col]] %in% c("HIGH","LOW") &
                            dist_border < 50) |>
            dplyr::pull(cell_id)
        train_bad_var <- unique(c(train_bad_var, bad_ids))

        qs <- stats::quantile(df$log2AspectRatio, probs=c(0.25,0.75))
        good_ids <- df |>
            dplyr::filter(log2AspectRatio > qs[1] &
                            log2AspectRatio < qs[2] &
                            dist_border > 50) |>
            dplyr::pull(cell_id)
        train_good_var <- unique(c(train_good_var, good_ids))

        idx <- grep(pattern="log2AspectRatio_outlier", x=out_var)
        if (length(idx)) {
            names(out_var)[idx] <-
                "I(abs(log2AspectRatio) * as.numeric(dist_border < 50))"
        }
    }

    train_bad <- df |>
        dplyr::filter(cell_id %in% train_bad_var) |>
        dplyr::mutate(qcscore_train=0)

    train_good <- df |>
        dplyr::filter(cell_id %in% train_good_var) |>
        dplyr::mutate(
            qcscore_train=1,
            is_a_bad_boy=cell_id %in% train_bad$cell_id
        )

    train_bad <- train_bad |>
        dplyr::distinct(cell_id, .keep_all=TRUE)

    if (verbose) message(paste0("Chosen low qual examples: ", nrow(train_bad)))

    # good example duplicates removal without any warning to the user
    train_good <- train_good |>
        dplyr::distinct(cell_id, .keep_all=TRUE)

    train_good <- train_good[!train_good$is_a_bad_boy, ]

    n_bad <- nrow(train_bad)
    stopifnot("Not enough good examples to match bad examples"=
                nrow(train_good) > n_bad)
    idx <- sample(seq_len(nrow(train_good)), n_bad, replace=FALSE)
    train_good <- train_good[idx, ]

    if (verbose) {
        message(paste0(
            "Chosen good quality examples (should match bad): ",
            nrow(train_good)
        ))
    }

    train_good$is_a_bad_boy <- NULL
    train_df <- rbind(train_bad, train_good)
    train_df <- train_df |>
        dplyr::distinct(cell_id, .keep_all=TRUE)
    return(train_df)
}



#' getModelFormula
#' @name getModelFormula
#' @rdname getModelFormula
#' @description
#' Returns the right‐hand side of a model formula string based on formula
#' variables found in the `metadata` of a `SpatialExperiment` object.
#' @param formulaVars A named character vector mapping variable names
#'   (e.g. `"log2CountArea"`, `"Area_um"`, etc.) to their corresponding
#'   outlier label columns, typically from
#'   `metadata(spe)$formula_variables`.
#' @param verbose Logical. If `TRUE`, prints the final formula used for QC score
#' @return `character`
#'   A one‐sided formula as a string (e.g. "~ log2CountArea + ...").
#' @export
#' @examples
#' example(checkOutliers)
#' getModelFormula(metadata(spe)$formula_variables)
getModelFormula <- function(formulaVars, verbose=FALSE)
{
    out_var <- formulaVars
    if ("log2AspectRatio" %in% names(out_var)) {
        names(out_var)[grep("log2AspectRatio_outlier", out_var)] <-
            "I(abs(log2AspectRatio) * as.numeric(dist_border<50))"
    }
    model_formula <- paste0("~(", paste(names(out_var), collapse = " + "),
                        ")^2", sep = "")

    if (verbose) {
        message("Final formula used for QC score computation:")
        message(model_formula)
    }

    return(model_formula)
}

#' .computeXenMerTrainSet
#' @name dot-computeXenMerTrainSet
#' @rdname dot-computeXenMerTrainSet
#' @description
#' Internal: Build Training Set for Xenium & MERFISH
#' Splits a SpatialExperiment into “bad” vs “good” cells based on
#' pre-computed outlier labels on log2CountArea.
#' @param spe \code{SpatialExperiment}
#' @return
#' A list with elements \code{bad} and \code{good}, each a data.frame
#' with \code{qscore_train} and (for “good”) an \code{is_a_bad_boy} flag.
#' @keywords internal
.computeXenMerTrainSet <- function(spe)
{
    train_bad <- data.frame(colData(spe)) |>
        filter(log2CountArea_outlier_train=="LOW") |> mutate(qscore_train=0)
    train_good <- data.frame(colData(spe)) |>
        filter((log2CountArea > quantile(log2CountArea, probs = 0.90) &
                    log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
        mutate(qscore_train=1, is_a_bad_boy=cell_id %in% train_bad$cell_id)
    return(list(bad=train_bad, good=train_good))
}

#' .computeCosmxTrainSet
#' @name dot-computeCosmxTrainSet
#' @rdname dot-computeCosmxTrainSet
#' @description
#' Internal: Build Training Set for CosMx
#' Splits a SpatialExperiment into “bad” vs “good” cells based on
#' outliers in aspect ratio near tissue border or low count area.
#' @param spe \code{SpatialExperiment}
#' @return
#' A list with elements \code{bad} and \code{good}, each a data.frame
#' with \code{qscore_train} and (for “good”) an \code{is_a_bad_boy} flag.
#' @keywords internal
.computeCosmxTrainSet <- function(spe)
{
    spe <- computeSpatialOutlier(spe, "log2AspectRatio", "scuttle")
    train_bad <- data.frame(colData(spe)) |>
        filter((log2AspectRatio_outlier_sc == "HIGH" & dist_border < 50) |
                (log2AspectRatio_outlier_sc == "LOW" & dist_border < 50) |
                log2CountArea_outlier_train == "LOW") |>
        mutate(qscore_train = 0)

    train_good <- data.frame(colData(spe)) |>
        filter((log2AspectRatio > quantile(log2AspectRatio, probs = 0.25) &
                    log2AspectRatio < quantile(log2AspectRatio, probs = 0.75) &
                    dist_border > 50) |
                    (log2CountArea > quantile(log2CountArea, probs = 0.90) &
                    log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
        mutate(qscore_train=1, is_a_bad_boy=cell_id %in% train_bad$cell_id)
    return(list(bad=train_bad, good=train_good))
}

#' .computeCosmxProteinTrainSet
#' @name dot-computeCosmxProteinTrainSet
#' @rdname dot-computeCosmxProteinTrainSet
#' @description
#' Internal: Build Training Set for CosMx-Protein
#' Splits a SpatialExperiment into “bad” vs “good” cells based on
#' outliers in aspect ratio near tissue border or low count area.
#' @param spe \code{SpatialExperiment}
#' @return
#' A list with elements \code{bad} and \code{good}, each a data.frame
#' with \code{qscore_train} and (for “good”) an \code{is_a_bad_boy} flag.
#' @keywords internal
.computeCosmxProteinTrainSet <- function(spe)
{
    spe <- computeSpatialOutlier(spe, "log2AspectRatio", method="both")
    spe <- computeSpatialOutlier(spe, "log2Ctrl_total_ratio", method="both")
    train_bad <- data.frame(colData(spe)) |>
        filter((log2AspectRatio_outlier_mc == "HIGH" & dist_border < 50) |
                (log2AspectRatio_outlier_mc == "LOW" & dist_border < 50) |
                log2CountArea_outlier_train == "LOW" |
                log2Ctrl_total_ratio_outlier_sc == "HIGH") |>
        mutate(qscore_train = 0)

    train_good <- data.frame(colData(spe)) |>
        filter((log2AspectRatio > quantile(log2AspectRatio, probs = 0.25) &
                log2AspectRatio < quantile(log2AspectRatio, probs = 0.75) &
                dist_border > 50) |
                (log2CountArea > quantile(log2CountArea, probs = 0.90) &
                log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
        mutate(qscore_train=1, is_a_bad_boy=cell_id %in% train_bad$cell_id)
    return(list(bad=train_bad, good=train_good))
}


#' computeQCScoreFlags
#' @name computeQCScoreFlags
#' @rdname computeQCScoreFlags
#' @description
#' Compute flagged cells based on a manually chosen threshold on quality score
#'
#' This function Compute flagged cells based on a manually chosen threshold on
#' quality score stored in `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param qsThreshold Numeric threshold or quantile for quality score. Default
#'   `0.5`.
#' @param useQSQuantiles Logical; if `TRUE`, treat `qsThreshold` as a
#'   percentile.
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' example(computeQCScore)
#' spe <- computeQCScoreFlags(spe)
#' table(spe$low_qcscore)
#' # if fixed filters are defined we have an additional column
#' spe <- computeThresholdFlags(spe)
#' spe <- computeQCScoreFlags(spe)
#' table(spe$low_threshold_qcscore)
computeQCScoreFlags <- function(spe, qsThreshold=0.5, useQSQuantiles=FALSE) {
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("QC_score" %in% names(colData(spe)))

    if(useQSQuantiles) {
        spe$low_qcscore <- ifelse(
            spe$QC_score < quantile(spe$QC_score, probs=qsThreshold),
            TRUE, FALSE)
    } else {
        spe$low_qcscore <- spe$QC_score < qsThreshold
    }

    if("threshold_flags" %in% names(colData(spe))) {
        spe$low_threshold_qcscore <- (spe$low_qcscore &
                                        spe$threshold_flags)
    }
    return(spe)
}

#' .checkSkw
#' @name dot-checkSkw
#' @rdname dot-checkSkw
#' @description
#' Check skewness of metrics to choose outlier detection method.
#'
#' @param cd colData of `SpatialExperiment` object.
#' @param metricList A character vector specifying the metrics to include in
#' the QC score formula. Defaults are "log2CountArea", "Area_um",
#' "log2AspectRatio", "log2Ctrl_total_ratio".
#'
#' @return
#' A vector containing the list of chosen outlier detection method for each
#' metric.
#' @importFrom e1071 skewness
#' @keywords internal
.checkSkw <- function(cd, metricList=c("log2CountArea", "Area_um",
    "log2AspectRatio", "log2Ctrl_total_ratio")) {
    stopifnot(is(cd, "DataFrame"))
    stopifnot(all(metricList %in% names(cd)))
    method <- c()
    for (i in metricList) {
        skw <- e1071::skewness(cd[[i]], na.rm = TRUE)
        method[i] <- ifelse((skw>-1 & skw<1), "sc", "mc")
    }
    for (submeth in c("log2CountArea", "log2Ctrl_total_ratio"))
    {
        if (!submeth %in% names(method)) next
        idx <- switch(submeth,
                    "log2CountArea"        = cd[["total"]] > 0,
                    "log2Ctrl_total_ratio" = cd[["ctrl_total_ratio"]] != 0,
                    rep(TRUE, nrow(cd))
        )
        if (!any(idx, na.rm = TRUE)) next
        skw <- e1071::skewness(cd[[submeth]][idx], na.rm=TRUE)
        newm <- ifelse((skw>-1 & skw<1), "sc", "mc")
        if (!is.na(newm)) method[[submeth]] <- newm
    }
    return(method)
}

#' computeOutliersQCScore
#' @name computeOutliersQCScore
#' @rdname computeOutliersQCScore
#' @description
#' Compute outlier cells for each metric that can be used in QC score formula
#' for SpatialExperiment.
#'
#' This function calculates outlier cells for each variable specified
#' in `metricList` for a `SpatialExperiment`. Log2CountArea must be present in
#' the `colData` of the `SpatialExperiment` object as a minimum requirement.
#' The user can choose which metrics to include among the following: Area_um,
#' log2Ctrl_total_ratio, log2AspectRatio. For Xenium and Merfish datasets,
#' log2AspectRatio is automatically removed from the formula.
#'
#' @param spe A `SpatialExperiment` object with spatial omics data.
#' @param metricList A character vector specifying the metrics to include in
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
computeOutliersQCScore <- function(spe, metricList=c("log2CountArea","Area_um",
    "log2AspectRatio", "log2Ctrl_total_ratio")) {

    stopifnot(is(spe, "SpatialExperiment"))

    cd <- colData(spe)
    method <- .checkSkw(cd, metricList)

    spec <- list(
        log2CountArea = list(
            idx  = spe$total > 0,
            zero = spe$total == 0,
            tweak_lower = TRUE
        ),
        log2Ctrl_total_ratio = list(
            idx  = spe$ctrl_total_ratio != 0,
            zero = spe$ctrl_total_ratio == 0,
            tweak_lower = FALSE
        )
    )
    if(any(is.na(method)))
    {
        idx <- which(is.na(method))
        for (i in idx) {
            message("Metric ", names(method)[i], " produced NA skewness values.",
                "It will be removed from the QC score formula.")
            method <- method[-i]
        }
    }
    for (var in intersect(names(method), names(spec))) {
        s <- spec[[var]]
        spe_temp <- computeSpatialOutlier(spe[, s$idx], computeBy=var,
                                        method=method[[var]])
        out_var <- paste0(var, "_outlier_", method[[var]])
        low_thr  <- getFencesOutlier(spe_temp, out_var, "lower")
        high_thr <- getFencesOutlier(spe_temp, out_var, "higher")
        if (isTRUE(s$tweak_lower)) {
            if (low_thr < min(spe_temp[[var]], na.rm = TRUE)) {
                low_thr <- stats::quantile(spe[[var]], probs=0.01)
            }
        }
        tgt <- paste0(var, "_outlier_train")
        spe[[tgt]] <- dplyr::case_when(
            s$zero ~ "NO",
            spe[[var]] < low_thr  ~ "LOW",
            spe[[var]] > high_thr ~ "HIGH",
            TRUE ~ "NO"
        )
        spe[[tgt]] <- scuttle::outlier.filter(spe[[tgt]])#is this really needed?
        thr <- getFencesOutlier(spe_temp, out_var)
        if (isTRUE(s$tweak_lower)) thr[1] <- low_thr
        attr(spe[[tgt]], "thresholds") <- thr
    }

    submethod <- method[!names(method) %in%
        c("log2CountArea", "log2Ctrl_total_ratio")]

    for(j in names(submethod)){
        spe <- computeSpatialOutlier(spe, computeBy=j, method=submethod[j])
    }

    out_var <- paste0(names(method), "_outlier_", method)
    names(out_var) <- names(method)
    # gives warning if one of the variables is missing, but still works!
    # out_var[names(out_var) %in% c("log2CountArea", "log2Ctrl_total_ratio")] <-
    #     c("log2CountArea_outlier_train", "log2Ctrl_total_ratio_outlier_train")
    mapping <- c(
        log2CountArea = "log2CountArea_outlier_train",
        log2Ctrl_total_ratio = "log2Ctrl_total_ratio_outlier_train"
    )
    present <- intersect(names(out_var), names(mapping))
    if (length(present) > 0L) {
        for (nm in present) out_var[nm] <- mapping[nm]
    }
    metadata(spe)$formula_variables <- out_var
    return(spe)
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
#' spe <- checkOutliers(spe, verbose = TRUE)
#' metadata(spe)$formula_variables
checkOutliers <- function(spe, verbose=FALSE) {
    warnstopmsg <- function(var, warnstop=c("w","s")) {
        warnstop <- match.arg(warnstop)
        m1 <- paste0("Not enough outlier cells for ", var, ".\n")
        m2 <- switch(warnstop,
                s="QC score computation cannot be performed",
                w="This variable will not be used in the final formula")
        return(paste0(m1, m2))
    }
    out_var <- metadata(spe)$formula_variables
    cd <- colData(spe)
    if (verbose) {
        for (i in names(out_var)) {
            message("Outliers found for ", i, ":")
            message(table(cd[[out_var[i]]]))
        }
    }
    stopifnot(
        "log2CountArea is not included in the QC score formula.\n
        QC score cannot be computed"=
            "log2CountArea" %in% names(out_var)
    )
    cfg <- list(
        log2CountArea=list(pattern="log2CountArea_outlier",
            remove="log2CountArea_outlier_train", label="LOW", act=stop,
            code="s"),
        Area_um=list(pattern="Area_um_outlier", remove="Area_um_outlier",
            label="HIGH", act=warning, code="w"),
        log2Ctrl_total_ratio=list(pattern="log2Ctrl_total_ratio_outlier",
            remove="log2Ctrl_total_ratio_outlier_train", label="HIGH",
            act=warning, code="w")
    )
    for (v in intersect(names(cfg), names(out_var))) {
        r <- cfg[[v]]
        gi <- grep(r$pattern, out_var)
        if (length(gi)==0L) next
        col <- out_var[gi]
        labs <- factor(cd[[col]], levels=c("LOW","HIGH","NO"))
        tab <- table(labs)
        cnt <- tab[r$label]
        if (is.null(cnt) || is.na(cnt)) cnt <- 0L
        if (cnt < ncol(spe)*0.001) {
            r$act(warnstopmsg(v, r$code))
            out_var <- out_var[-grep(r$remove, out_var)]
        }
    }
    var <- "log2AspectRatio"
    pat <- paste0(var, "_outlier")
    is_cosmx <- metadata(spe)$technology %in%
                c("Nanostring_CosMx", "Nanostring_CosMx_Protein")
    if (is_cosmx && (var %in% names(out_var))) {
        idx <- grep(pat, out_var)
        if (length(idx)) {
            col <- out_var[idx]
            labs <- factor(cd[[col]], levels=c("LOW","HIGH","NO"))
            tab <- table(labs)
            nmin <- ncol(spe) * 0.001

            low  <- tab[["LOW"]];  if (is.na(low))  low  <- 0L
            high <- tab[["HIGH"]]; if (is.na(high)) high <- 0L

            if (low < nmin && high < nmin) {
                warning(warnstopmsg(var, "w"))
                out_var <- out_var[-grep(pat, out_var)]
            }
        }
    } else {
        out_var <- out_var[-grep(pat, out_var)]
    }
    metadata(spe)$formula_variables <- out_var
    return(spe)
}


.prepQCContext <- function(spe, metricList=c("log2CountArea", "Area_um",
    "log2AspectRatio", "log2Ctrl_total_ratio"), verbose=FALSE) {

    stopifnot(is(spe,"SpatialExperiment"))

    # remove zero-count cells once
    zerocells <- spe$total==0
    if (sum(zerocells) > 0) {
        warning(paste0(sum(zerocells),
            " cells with 0 counts were found. These cells will be removed."))
        spe <- spe[, sum(zerocells)]
    }

    spe1 <- computeOutliersQCScore(spe, metricList)
    spe1 <- checkOutliers(spe1, verbose)

    out_var <- metadata(spe1)$formula_variables
    df <- as.data.frame(colData(spe1))
    tech <- metadata(spe1)$technology

    return(list(df=df, out_var=out_var, tech=tech))
}

