library(SpaceTrooper)
library(SummarizedExperiment)
library(S4Vectors)

# Read a dataset with the appropriate SpaceTrooper reader according to the
# selected technology.
#
# Supported values are:
# - "cosmx"
# - "cosmx_protein"
# - "xenium"
read_dataset_for_transfer <- function(
    dir_name,
    technology=c("cosmx", "cosmx_protein", "xenium"),
    sample_name="sample01"
) {
    technology <- match.arg(technology)

    switch(
        technology,
        cosmx=readCosmxSPE(dir_name, sampleName=sample_name),
        cosmx_protein=readCosmxProteinSPE(dir_name, sampleName=sample_name),
        xenium=readXeniumSPE(
            dir_name,
            sampleName=sample_name,
            computeMissingMetrics=TRUE,
            keepPolygons=FALSE
        )
    )
}

# Compute the per-cell QC metrics required before training or applying QS.
#
# The same preprocessing is applied to both source and target datasets so that
# transferred coefficients operate on harmonised metric columns.
prepare_dataset_for_transfer <- function(
    spe,
    rm_zeros=TRUE,
    neg_prob_list=c("NegPrb", "Negative", "SystemControl")
) {
    spatialPerCellQC(
        spe,
        rmZeros=rm_zeros,
        negProbList=neg_prob_list
    )
}

# Identify the metric set that can be transferred from source to target.
#
# For heterogeneous transfers, the safest common metrics are:
# - log2SignalDensity
# - Area_um
# - log2Ctrl_total_ratio
#
# The border-dependent aspect-ratio term is included only if both datasets have
# `log2AspectRatio` and `dist_border`, which is typically the case for CosMx to
# CosMx transfers.
get_transferable_metrics <- function(source_spe, target_spe) {
    source_cols <- colnames(SummarizedExperiment::colData(source_spe))
    target_cols <- colnames(SummarizedExperiment::colData(target_spe))
    shared_cols <- intersect(source_cols, target_cols)

    metrics <- intersect(
        c("log2SignalDensity", "Area_um", "log2Ctrl_total_ratio"),
        shared_cols
    )

    if (all(c("log2AspectRatio", "dist_border") %in% source_cols) &&
            all(c("log2AspectRatio", "dist_border") %in% target_cols)) {
        metrics <- c(metrics, "log2AspectRatio")
    }

    unique(metrics)
}

# Extract a readable coefficient table from a fitted glmnet model.
#
# Coefficients are returned at the selected lambda, so they are exactly the ones
# transferred to the target dataset.
extract_qs_coefficients <- function(fit, lambda) {
    coef_mat <- as.matrix(stats::predict(fit, s=lambda, type="coefficients"))
    data.frame(
        term=rownames(coef_mat),
        coefficient=as.numeric(coef_mat[, 1]),
        row.names=NULL
    )
}

# Apply a fitted QS model to a new SpatialExperiment object.
#
# The design matrix is rebuilt from the target `colData` using the compatible
# transfer formula trained on the source dataset.
score_target_dataset <- function(spe, fit, lambda, model_formula) {
    df <- as.data.frame(SummarizedExperiment::colData(spe))
    x <- stats::model.matrix(stats::as.formula(model_formula), data=df)
    spe$QC_score <- as.vector(
        stats::predict(fit, s=lambda, newx=x, type="response")
    )
    spe
}

# Compute AUROC with a rank-based formula, avoiding extra package dependencies.
auc_rank <- function(score, label_positive) {
    label_positive <- as.integer(label_positive)
    n_pos <- sum(label_positive == 1L)
    n_neg <- sum(label_positive == 0L)

    if (n_pos == 0L || n_neg == 0L) {
        return(NA_real_)
    }

    ranks <- rank(score, ties.method="average")
    (sum(ranks[label_positive == 1L]) - n_pos * (n_pos + 1) / 2) /
        (n_pos * n_neg)
}

# Build pseudo-reference labels on a dataset using the same metric subset used
# for transfer.
#
# This is useful for evaluating how well transferred coefficients separate
# pseudo-bad from pseudo-good cells on the target dataset.
build_reference_labels_transfer <- function(spe, metric_list, verbose=FALSE) {
    spe_ref <- computeOutliersQCScore(spe, metricList=metric_list)
    spe_ref <- checkOutliers(spe_ref, verbose=verbose)

    df_ref <- computeTrainDF(
        colData=SummarizedExperiment::colData(spe_ref),
        formulaVars=S4Vectors::metadata(spe_ref)$formula_variables,
        tech=S4Vectors::metadata(spe_ref)$technology,
        verbose=verbose
    )

    df_ref[, c("cell_id", "qcscore_train")]
}

# Safe wrapper for pseudo-reference label generation.
safe_reference_labels_transfer <- function(spe, metric_list, verbose=FALSE) {
    tryCatch(
        list(
            data=build_reference_labels_transfer(
                spe,
                metric_list=metric_list,
                verbose=verbose
            ),
            error=NULL
        ),
        error=function(e) list(data=NULL, error=conditionMessage(e))
    )
}

# Summarise the separation between pseudo-bad and pseudo-good cells.
#
# The summary is computed on either source or target after scoring with the
# transferred model.
summarise_transfer_eval <- function(df_eval, threshold=0.5) {
    if (is.null(df_eval) || nrow(df_eval) == 0L) {
        return(data.frame())
    }

    is_bad <- df_eval$qcscore_train == 0
    is_good <- df_eval$qcscore_train == 1
    pred_bad <- df_eval$QC_score < threshold

    data.frame(
        n_labeled=nrow(df_eval),
        n_bad=sum(is_bad),
        n_good=sum(is_good),
        median_qs_bad=if (sum(is_bad) > 0) {
            stats::median(df_eval$QC_score[is_bad])
        } else {
            NA_real_
        },
        median_qs_good=if (sum(is_good) > 0) {
            stats::median(df_eval$QC_score[is_good])
        } else {
            NA_real_
        },
        auc_bad_vs_good=auc_rank(1 - df_eval$QC_score, as.integer(is_bad)),
        sensitivity_bad_qs_lt_threshold=if (sum(is_bad) > 0) {
            mean(pred_bad[is_bad])
        } else {
            NA_real_
        },
        specificity_good_qs_ge_threshold=if (sum(is_good) > 0) {
            mean(!pred_bad[is_good])
        } else {
            NA_real_
        }
    )
}

# Add good/bad/unlabelled labels to a SpatialExperiment object for plotting.
#
# Labels are derived from the evaluation table built for either the source or
# the target dataset after scoring with transferred coefficients.
append_transfer_labels <- function(
    spe,
    df_eval,
    label_col="transfer_label",
    colour_col="transfer_label_color"
) {
    labels <- rep("unlabelled", ncol(spe))

    if (!is.null(df_eval) && nrow(df_eval) > 0L) {
        idx <- match(spe$cell_id, df_eval$cell_id)
        keep <- !is.na(idx)
        labels[keep] <- ifelse(df_eval$qcscore_train[idx[keep]] == 0,
            "bad", "good")
    }

    spe[[label_col]] <- factor(labels, levels=c("bad", "good", "unlabelled"))
    palette <- c(bad="#D55E00", good="#009E73", unlabelled="#BDBDBD")
    spe[[colour_col]] <- unname(palette[as.character(spe[[label_col]])])
    spe
}

# Build a ROC data.frame from pseudo-reference labels and transferred QS.
build_transfer_roc_df <- function(df_eval) {
    if (is.null(df_eval) || nrow(df_eval) == 0L) {
        return(data.frame())
    }

    label_positive <- as.integer(df_eval$qcscore_train == 0)
    n_pos <- sum(label_positive == 1L)
    n_neg <- sum(label_positive == 0L)

    if (n_pos == 0L || n_neg == 0L) {
        return(data.frame())
    }

    score_bad <- 1 - df_eval$QC_score
    ord <- order(score_bad, decreasing=TRUE)
    score_bad <- score_bad[ord]
    label_positive <- label_positive[ord]

    tp <- cumsum(label_positive == 1L)
    fp <- cumsum(label_positive == 0L)
    change <- c(score_bad[-1] != score_bad[-length(score_bad)], TRUE)

    data.frame(
        fpr=c(0, fp[change] / n_neg),
        tpr=c(0, tp[change] / n_pos)
    )
}

# Plot source or target cells coloured as pseudo-good / pseudo-bad after
# scoring with transferred coefficients.
plot_transfer_labelled_cells <- function(
    result,
    dataset=c("source", "target"),
    size=0.05,
    alpha=0.8,
    scaleBar=TRUE
) {
    dataset <- match.arg(dataset)

    spe <- switch(dataset,
        source=result$source_spe,
        target=result$target_spe
    )
    df_eval <- switch(dataset,
        source=result$source_eval,
        target=result$target_eval
    )

    spe <- append_transfer_labels(spe, df_eval)

    plotCentroids(
        spe,
        colourBy="transfer_label",
        palette="transfer_label_color",
        size=size,
        alpha=alpha,
        scaleBar=scaleBar
    ) +
        ggplot2::labs(
            title=paste0("Transferred QS labels on ", dataset, " dataset"),
            colour="Reference label",
            fill="Reference label"
        )
}

# Plot the spatial distribution of transferred QS on source or target.
plot_transfer_qs <- function(
    result,
    dataset=c("source", "target"),
    size=0.05,
    alpha=0.8,
    scaleBar=TRUE
) {
    dataset <- match.arg(dataset)

    spe <- switch(dataset,
        source=result$source_spe,
        target=result$target_spe
    )

    plotCentroids(
        spe,
        colourBy="QC_score",
        size=size,
        alpha=alpha,
        scaleBar=scaleBar
    ) +
        ggplot2::labs(title=paste0("Transferred QC score on ", dataset, " dataset"))
}

# Plot ROC curves for the source or target transfer evaluation.
plot_transfer_roc <- function(result, dataset=c("source", "target")) {
    dataset <- match.arg(dataset)

    df_eval <- switch(dataset,
        source=result$source_eval,
        target=result$target_eval
    )
    roc_df <- build_transfer_roc_df(df_eval)

    if (nrow(roc_df) == 0L) {
        stop("No ROC curve available for the selected dataset")
    }

    auc_value <- switch(dataset,
        source=result$source_summary$auc_bad_vs_good,
        target=result$target_summary$auc_bad_vs_good
    )

    ggplot2::ggplot(roc_df, ggplot2::aes(x=fpr, y=tpr)) +
        ggplot2::geom_abline(intercept=0, slope=1, linetype=2,
            colour="grey60") +
        ggplot2::geom_path(linewidth=0.9, colour="#0072B2") +
        ggplot2::coord_fixed() +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title=paste0("ROC curve for transferred QS (", dataset, ")"),
            subtitle=paste0("AUC = ", round(auc_value, 4)),
            x="False positive rate",
            y="True positive rate"
        )
}

# Compare AUC between source and target after transfer.
plot_transfer_auc_summary <- function(result) {
    auc_df <- data.frame(
        dataset=c("source", "target"),
        auc_bad_vs_good=c(
            result$source_summary$auc_bad_vs_good,
            if (nrow(result$target_summary) == 0L) NA_real_ else result$target_summary$auc_bad_vs_good
        )
    )

    ggplot2::ggplot(auc_df, ggplot2::aes(x=dataset, y=auc_bad_vs_good, fill=dataset)) +
        ggplot2::geom_col(width=0.65, alpha=0.85, show.legend=FALSE) +
        ggplot2::geom_text(ggplot2::aes(label=round(auc_bad_vs_good, 4)),
            vjust=-0.4, na.rm=TRUE) +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title="AUC summary for transferred QS",
            x=NULL,
            y="AUC"
        )
}

# Train a native QS model directly on the target dataset, so it can be compared
# against the transferred QS.
#
# By default it uses the full SpaceTrooper native workflow on the target.
# If `use_transfer_metrics = TRUE`, it restricts the native training to the same
# metric subset used by the transferred model.
compute_native_target_qs <- function(
    result,
    use_transfer_metrics=FALSE,
    verbose=TRUE
) {
    target_native <- result$target_spe
    target_native$QC_score_transfer <- target_native$QC_score

    if (isTRUE(use_transfer_metrics)) {
        target_native <- computeOutliersQCScore(
            target_native,
            metricList=result$transfer_metrics
        )
        target_native <- checkOutliers(target_native, verbose=verbose)

        formula_vars <- S4Vectors::metadata(target_native)$formula_variables
        model_formula <- getModelFormula(formula_vars, verbose=verbose)
        train_df <- computeTrainDF(
            colData=SummarizedExperiment::colData(target_native),
            formulaVars=formula_vars,
            tech=S4Vectors::metadata(target_native)$technology,
            verbose=verbose
        )
        x_target <- stats::model.matrix(stats::as.formula(model_formula), data=train_df)
        best_lambda <- computeLambda(train_df, model_formula)
        fit_native <- trainModel(x_target, train_df)
        target_native$QC_score_native <- as.vector(
            stats::predict(
                fit_native,
                s=best_lambda,
                newx=stats::model.matrix(
                    stats::as.formula(model_formula),
                    data=as.data.frame(SummarizedExperiment::colData(target_native))
                ),
                type="response"
            )
        )
        return(target_native)
    }

    target_native <- computeQCScore(target_native, verbose=verbose)
    target_native$QC_score_native <- target_native$QC_score
    target_native$QC_score <- target_native$QC_score_transfer
    target_native
}

# Compare transferred QS with a native QS trained directly on the target.
#
# The function reports correlations and agreement on low-quality calls.
compare_transferred_vs_native <- function(
    result,
    threshold=0.5,
    low_fraction=0.1,
    use_transfer_metrics=FALSE,
    verbose=TRUE
) {
    stopifnot(length(low_fraction) == 1L, low_fraction > 0, low_fraction < 1)

    target_native <- compute_native_target_qs(
        result,
        use_transfer_metrics=use_transfer_metrics,
        verbose=verbose
    )

    df_compare <- data.frame(
        cell_id=target_native$cell_id,
        qc_transfer=target_native$QC_score_transfer,
        qc_native=target_native$QC_score_native
    )

    transfer_cut <- stats::quantile(df_compare$qc_transfer, probs=low_fraction, na.rm=TRUE)
    native_cut <- stats::quantile(df_compare$qc_native, probs=low_fraction, na.rm=TRUE)

    comparison_summary <- data.frame(
        pearson=stats::cor(df_compare$qc_transfer, df_compare$qc_native,
            method="pearson", use="complete.obs"),
        spearman=stats::cor(df_compare$qc_transfer, df_compare$qc_native,
            method="spearman", use="complete.obs"),
        low_qc_agreement_at_threshold=mean(
            (df_compare$qc_transfer < threshold) == (df_compare$qc_native < threshold),
            na.rm=TRUE
        ),
        overlap_lowest_fraction=mean(
            (df_compare$qc_transfer <= transfer_cut) &
                (df_compare$qc_native <= native_cut),
            na.rm=TRUE
        ) / low_fraction
    )

    list(
        summary=comparison_summary,
        df_compare=df_compare,
        target_native_spe=target_native,
        threshold=threshold,
        low_fraction=low_fraction,
        use_transfer_metrics=use_transfer_metrics
    )
}

# Scatter plot of transferred QS vs native QS on the target dataset.
plot_transfer_vs_native_scatter <- function(comparison_result) {
    ggplot2::ggplot(
        comparison_result$df_compare,
        ggplot2::aes(x=qc_native, y=qc_transfer)
    ) +
        ggplot2::geom_point(alpha=0.4, size=0.6, colour="#0072B2") +
        ggplot2::geom_abline(intercept=0, slope=1, linetype=2,
            colour="grey60") +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title="Transferred QS vs native QS on target",
            subtitle=paste0(
                "Spearman = ",
                round(comparison_result$summary$spearman, 4),
                ", agreement@",
                comparison_result$threshold,
                " = ",
                round(comparison_result$summary$low_qc_agreement_at_threshold, 4)
            ),
            x="Native target QS",
            y="Transferred QS"
        )
}

# Train a transferable QS model on one dataset and apply it to another.
#
# Workflow:
# 1. read source and target datasets
# 2. compute QC metrics on both
# 3. identify the metric subset that exists in both datasets
# 4. train the ridge-logistic QS model on the source dataset using only those metrics
# 5. transfer the fitted coefficients to the target dataset
# 6. optionally evaluate the transfer on pseudo-reference labels in the target
transfer_qs_coefficients <- function(
    source_dir,
    target_dir,
    source_technology=c("cosmx", "cosmx_protein", "xenium"),
    target_technology=c("cosmx", "cosmx_protein", "xenium"),
    source_sample_name="source_sample",
    target_sample_name="target_sample",
    threshold=0.5,
    evaluate_target=TRUE,
    verbose=TRUE
) {
    source_technology <- match.arg(source_technology)
    target_technology <- match.arg(target_technology)

    # Read the source dataset, from which coefficients will be learned.
    source_spe <- read_dataset_for_transfer(
        dir_name=source_dir,
        technology=source_technology,
        sample_name=source_sample_name
    )

    # Read the target dataset, on which coefficients will be transferred.
    target_spe <- read_dataset_for_transfer(
        dir_name=target_dir,
        technology=target_technology,
        sample_name=target_sample_name
    )

    # Harmonise per-cell QC metrics in both datasets.
    source_spe <- prepare_dataset_for_transfer(source_spe)
    target_spe <- prepare_dataset_for_transfer(target_spe)

    # Keep only the subset of QS metrics that is available in both datasets.
    transfer_metrics <- get_transferable_metrics(source_spe, target_spe)

    if (!"log2SignalDensity" %in% transfer_metrics) {
        stop("log2SignalDensity is required for transferable QS training")
    }
    if (length(transfer_metrics) < 2L) {
        stop("At least two shared metrics are recommended for coefficient transfer")
    }

    if (verbose) {
        message("Transferable metrics: ", paste(transfer_metrics, collapse=", "))
    }

    # Compute source outliers only on the transferable metric subset.
    source_spe <- computeOutliersQCScore(source_spe, metricList=transfer_metrics)
    source_spe <- checkOutliers(source_spe, verbose=verbose)

    formula_vars <- S4Vectors::metadata(source_spe)$formula_variables
    if (!"log2SignalDensity" %in% names(formula_vars)) {
        stop("The transferable source formula lost log2SignalDensity after outlier filtering")
    }

    # Build the compatible model formula and the balanced source training table.
    model_formula <- getModelFormula(formula_vars, verbose=verbose)
    source_train_df <- computeTrainDF(
        colData=SummarizedExperiment::colData(source_spe),
        formulaVars=formula_vars,
        tech=S4Vectors::metadata(source_spe)$technology,
        verbose=verbose
    )

    # Fit the ridge-logistic QS model on the source dataset.
    x_source <- stats::model.matrix(stats::as.formula(model_formula),
        data=source_train_df)
    best_lambda <- computeLambda(source_train_df, model_formula)
    fit <- trainModel(x_source, source_train_df)
    coefficient_table <- extract_qs_coefficients(fit, best_lambda)

    # Score the full source and target datasets with the same transferred model.
    source_spe <- score_target_dataset(source_spe, fit, best_lambda, model_formula)
    target_spe <- score_target_dataset(target_spe, fit, best_lambda, model_formula)

    # Build the labelled source evaluation set from the source training examples.
    source_eval <- merge(
        data.frame(cell_id=source_spe$cell_id, QC_score=source_spe$QC_score),
        source_train_df[, c("cell_id", "qcscore_train")],
        by="cell_id"
    )
    source_summary <- summarise_transfer_eval(source_eval, threshold=threshold)

    target_eval <- NULL
    target_summary <- data.frame()
    target_label_error <- NA_character_

    if (isTRUE(evaluate_target)) {
        # Build pseudo-reference labels independently on the target dataset,
        # using the same transferable metric subset.
        target_ref <- safe_reference_labels_transfer(
            target_spe,
            metric_list=transfer_metrics,
            verbose=verbose
        )
        target_label_error <- target_ref$error

        if (!is.null(target_ref$data)) {
            target_eval <- merge(
                data.frame(cell_id=target_spe$cell_id, QC_score=target_spe$QC_score),
                target_ref$data,
                by="cell_id"
            )
            target_summary <- summarise_transfer_eval(target_eval,
                threshold=threshold)
        }
    }

    list(
        source_technology=source_technology,
        target_technology=target_technology,
        transfer_metrics=transfer_metrics,
        model_formula=model_formula,
        lambda=best_lambda,
        coefficient_table=coefficient_table,
        source_summary=source_summary,
        target_summary=target_summary,
        target_label_error=target_label_error,
        source_train_df=source_train_df,
        source_eval=source_eval,
        target_eval=target_eval,
        source_spe=source_spe,
        target_spe=target_spe,
        fit=fit
    )
}

# Esempi di utilizzo:
# source("inst/scripts/transfer_qs_coefficients.R")
#
# Train su CosMx e applica a Xenium
# res_transfer_cx <- transfer_qs_coefficients(
#     source_dir="/path/to/cosmx_dataset",
#     source_technology="cosmx",
#     target_dir="/path/to/xenium_dataset",
#     target_technology="xenium",
#     source_sample_name="CosMx_source",
#     target_sample_name="Xenium_target"
# )
# res_transfer_cx$coefficient_table
# res_transfer_cx$target_summary
# plot_transfer_labelled_cells(res_transfer_cx, dataset="target", size=0.08)
# plot_transfer_qs(res_transfer_cx, dataset="target", size=0.08)
# plot_transfer_roc(res_transfer_cx, dataset="target")
# plot_transfer_auc_summary(res_transfer_cx)
#
# cmp_transfer_cx <- compare_transferred_vs_native(
#     res_transfer_cx,
#     threshold=0.5,
#     low_fraction=0.1
# )
# cmp_transfer_cx$summary
# plot_transfer_vs_native_scatter(cmp_transfer_cx)
#
# Train su Xenium e applica a CosMx
# res_transfer_xc <- transfer_qs_coefficients(
#     source_dir="/path/to/xenium_dataset",
#     source_technology="xenium",
#     target_dir="/path/to/cosmx_dataset",
#     target_technology="cosmx"
# )
# res_transfer_xc$coefficient_table
# res_transfer_xc$target_summary
# res_transfer_xc$target_summary
