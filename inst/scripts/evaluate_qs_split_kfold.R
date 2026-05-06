library(SpaceTrooper)
library(SummarizedExperiment)
library(S4Vectors)

# Apply a fitted ridge-logistic QS model to a new SpatialExperiment subset.
#
# The function rebuilds the design matrix from `colData`, using the same
# formula employed during training, and stores predicted probabilities in the
# `QC_score` column.
score_subset <- function(spe_subset, fit, lambda, model_formula) {
    # Convert per-cell metadata to a plain data.frame so model.matrix can use it.
    df_subset <- as.data.frame(SummarizedExperiment::colData(spe_subset))

    # Recreate the exact predictor matrix used by glmnet.
    x_subset <- stats::model.matrix(stats::as.formula(model_formula),
        data=df_subset)

    # Predict the probability of being a good-quality cell.
    spe_subset$QC_score <- as.vector(
        stats::predict(fit, s=lambda, newx=x_subset, type="response")
    )
    spe_subset
}

# Compute AUROC from ranks, without requiring extra packages.
#
# `label_positive` must be a 0/1 vector, where 1 indicates the positive class.
# In this script we use it both for bad-vs-good and for derived summaries.
auc_rank <- function(score, label_positive) {
    # Coerce labels to integer to avoid issues with logical/factor inputs.
    label_positive <- as.integer(label_positive)
    n_pos <- sum(label_positive == 1L)
    n_neg <- sum(label_positive == 0L)

    # AUROC is undefined if one of the two classes is missing.
    if (n_pos == 0L || n_neg == 0L) {
        return(NA_real_)
    }

    # Mann-Whitney / rank-based computation of AUROC.
    ranks <- rank(score, ties.method="average")
    (sum(ranks[label_positive == 1L]) - n_pos * (n_pos + 1) / 2) /
        (n_pos * n_neg)
}

# Build pseudo-reference labels on a subset, following the same SpaceTrooper
# logic used internally to define training examples.
#
# Returned labels are:
# - `qcscore_train = 0`: low-quality examples
# - `qcscore_train = 1`: good-quality examples
build_reference_labels <- function(spe_subset, verbose=FALSE) {
    # Detect outliers for the metrics used by QS.
    spe_ref <- computeOutliersQCScore(spe_subset)

    # Remove variables that do not have enough outliers to be informative.
    spe_ref <- checkOutliers(spe_ref, verbose=verbose)

    # Build a balanced table of bad and good examples.
    df_ref <- computeTrainDF(
        colData=SummarizedExperiment::colData(spe_ref),
        formulaVars=S4Vectors::metadata(spe_ref)$formula_variables,
        tech=S4Vectors::metadata(spe_ref)$technology,
        verbose=verbose
    )

    # Keep only the columns needed for downstream evaluation.
    df_ref[, c("cell_id", "qcscore_train")]
}

# Safe wrapper around `build_reference_labels()`.
#
# On small or unlucky splits/folds, SpaceTrooper may not find enough outliers
# to define a reference set. Instead of stopping the whole script, this returns
# the error message so the caller can inspect it.
safe_reference_labels <- function(spe_subset, verbose=FALSE) {
    tryCatch(
        list(data=build_reference_labels(spe_subset, verbose=verbose),
            error=NULL),
        error=function(e) list(data=NULL, error=conditionMessage(e))
    )
}

# Safely extract an optional scalar character field from a list.
#
# This is used when collecting fold-level error messages at the end of k-fold
# evaluation. Missing, NULL, or NA entries are converted to `NA_character_`.
extract_optional_chr1 <- function(x, name) {
    value <- x[[name]]

    if (is.null(value) || length(value) == 0L || all(is.na(value))) {
        return(NA_character_)
    }

    as.character(value[[1]])
}

# Add good/bad/unlabelled labels to a SpatialExperiment subset for plotting.
#
# Labels are derived from the evaluation table:
# - `qcscore_train = 0` -> bad
# - `qcscore_train = 1` -> good
# - missing label -> unlabelled
append_eval_labels <- function(
    spe_subset,
    df_eval,
    label_col="eval_label",
    colour_col="eval_label_color"
) {
    labels <- rep("unlabelled", ncol(spe_subset))

    if (!is.null(df_eval) && nrow(df_eval) > 0L) {
        idx <- match(spe_subset$cell_id, df_eval$cell_id)
        keep <- !is.na(idx)
        labels[keep] <- ifelse(df_eval$qcscore_train[idx[keep]] == 0,
            "bad", "good")
    }

    spe_subset[[label_col]] <- factor(labels,
        levels=c("bad", "good", "unlabelled"))

    palette <- c(
        bad="#D55E00",
        good="#009E73",
        unlabelled="#BDBDBD"
    )
    spe_subset[[colour_col]] <- unname(palette[as.character(spe_subset[[label_col]])])
    spe_subset
}

# Select either the whole result object (single split) or one fold result.
#
# This helper keeps the plotting functions compact and makes the API uniform
# across single-split and k-fold outputs.
get_result_level <- function(result, fold=NULL) {
    if (identical(result$split_type, "k_fold_cv")) {
        if (is.null(fold)) {
            stop("Please provide `fold` when plotting a k-fold result")
        }
        if (length(fold) != 1L || is.na(fold) || fold < 1 || fold > length(result$folds)) {
            stop("`fold` must be an integer between 1 and the number of folds")
        }
        return(result$folds[[fold]])
    }

    result
}

# Plot cells coloured as good, bad, or unlabelled.
#
# By default the function uses centroids because they are always available in
# the current workflow. For k-fold results, choose which fold to display.
plot_labelled_cells <- function(
    result,
    dataset=c("train", "test"),
    fold=NULL,
    size=0.05,
    alpha=0.8,
    scaleBar=TRUE
) {
    dataset <- match.arg(dataset)
    level_result <- get_result_level(result, fold=fold)

    spe_subset <- switch(dataset,
        train=level_result$spe_train,
        test=level_result$spe_test
    )
    df_eval <- switch(dataset,
        train=level_result$train_eval,
        test=level_result$test_eval
    )

    spe_subset <- append_eval_labels(spe_subset, df_eval)

    title_suffix <- if (identical(result$split_type, "k_fold_cv")) {
        paste0(" - fold ", fold)
    } else {
        ""
    }

    plotCentroids(
        spe_subset,
        colourBy="eval_label",
        palette="eval_label_color",
        size=size,
        alpha=alpha,
        scaleBar=scaleBar
    ) +
        ggplot2::labs(
            title=paste0("Labelled cells (", dataset, ")", title_suffix),
            colour="Reference label",
            fill="Reference label"
        )
}

# Compute the points of a ROC curve from evaluation labels and QS values.
#
# The positive class is the bad-quality group, so the curve is built using
# `1 - QC_score` as the decision score.
build_roc_df <- function(df_eval) {
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

# Plot ROC curves for a single split or across folds.
#
# In k-fold mode the function overlays one ROC curve per fold and reports the
# average AUC in the subtitle.
plot_roc_eval <- function(result, dataset=c("train", "test")) {
    dataset <- match.arg(dataset)

    if (identical(result$split_type, "k_fold_cv")) {
        roc_list <- lapply(seq_along(result$folds), function(i) {
            df_eval <- switch(dataset,
                train=result$folds[[i]]$train_eval,
                test=result$folds[[i]]$test_eval
            )
            roc_df <- build_roc_df(df_eval)
            if (nrow(roc_df) == 0L) {
                return(NULL)
            }
            roc_df$fold <- i
            roc_df
        })
        roc_list <- Filter(Negate(is.null), roc_list)

        if (length(roc_list) == 0L) {
            stop("No ROC curves available for the selected dataset")
        }

        roc_df <- do.call(rbind, roc_list)
        auc_mean <- switch(dataset,
            train=result$train_summary_average$auc_bad_vs_good,
            test=result$test_summary_average$auc_bad_vs_good
        )

        return(
            ggplot2::ggplot(roc_df,
                ggplot2::aes(x=fpr, y=tpr, colour=factor(fold), group=fold)) +
                ggplot2::geom_abline(intercept=0, slope=1, linetype=2,
                    colour="grey60") +
                ggplot2::geom_path(linewidth=0.7, alpha=0.8) +
                ggplot2::coord_fixed() +
                ggplot2::theme_bw() +
                ggplot2::labs(
                    title=paste0("ROC curves across folds (", dataset, ")"),
                    subtitle=paste0("Mean AUC = ", round(auc_mean, 4)),
                    x="False positive rate",
                    y="True positive rate",
                    colour="Fold"
                )
        )
    }

    df_eval <- switch(dataset,
        train=result$train_eval,
        test=result$test_eval
    )
    roc_df <- build_roc_df(df_eval)

    if (nrow(roc_df) == 0L) {
        stop("No ROC curve available for the selected dataset")
    }

    auc_value <- switch(dataset,
        train=result$train_summary$auc_bad_vs_good,
        test=result$test_summary$auc_bad_vs_good
    )

    ggplot2::ggplot(roc_df, ggplot2::aes(x=fpr, y=tpr)) +
        ggplot2::geom_abline(intercept=0, slope=1, linetype=2,
            colour="grey60") +
        ggplot2::geom_path(linewidth=0.9, colour="#0072B2") +
        ggplot2::coord_fixed() +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title=paste0("ROC curve (", dataset, ")"),
            subtitle=paste0("AUC = ", round(auc_value, 4)),
            x="False positive rate",
            y="True positive rate"
        )
}

# Plot AUC values for each fold.
#
# In single-split mode the function returns a one-point summary. In k-fold mode
# it shows one point per fold and a dashed line for the average AUC.
plot_auc_across_folds <- function(result, dataset=c("train", "test")) {
    dataset <- match.arg(dataset)

    if (identical(result$split_type, "k_fold_cv")) {
        auc_df <- switch(dataset,
            train=result$train_summary_per_fold,
            test=result$test_summary_per_fold
        )
        auc_mean <- switch(dataset,
            train=result$train_summary_average$auc_bad_vs_good,
            test=result$test_summary_average$auc_bad_vs_good
        )

        return(
            ggplot2::ggplot(auc_df,
                ggplot2::aes(x=fold, y=auc_bad_vs_good)) +
                ggplot2::geom_hline(yintercept=auc_mean, linetype=2,
                    colour="#D55E00") +
                ggplot2::geom_line(colour="#0072B2") +
                ggplot2::geom_point(size=2.2, colour="#0072B2") +
                ggplot2::theme_bw() +
                ggplot2::labs(
                    title=paste0("AUC across folds (", dataset, ")"),
                    subtitle=paste0("Dashed line = mean AUC = ", round(auc_mean, 4)),
                    x="Fold",
                    y="AUC"
                )
        )
    }

    auc_value <- switch(dataset,
        train=result$train_summary$auc_bad_vs_good,
        test=result$test_summary$auc_bad_vs_good
    )
    auc_df <- data.frame(split=1, auc_bad_vs_good=auc_value)

    ggplot2::ggplot(auc_df, ggplot2::aes(x=split, y=auc_bad_vs_good)) +
        ggplot2::geom_point(size=3, colour="#0072B2") +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(breaks=1, labels="split") +
        ggplot2::labs(
            title=paste0("AUC summary (", dataset, ")"),
            subtitle=paste0("AUC = ", round(auc_value, 4)),
            x=NULL,
            y="AUC"
        )
}

# Summarise QS performance against pseudo-reference labels.
#
# The summary reports:
# - number of labelled cells
# - number of bad/good examples
# - median QS in bad and good examples
# - AUROC for separating bad from good
# - sensitivity/specificity at the chosen threshold
summarise_eval <- function(df_eval, threshold=0.5) {
    # Return an empty data.frame if the evaluation table is missing.
    if (is.null(df_eval) || nrow(df_eval) == 0L) {
        return(data.frame())
    }

    # Derive boolean masks for the pseudo-reference classes.
    is_bad <- df_eval$qcscore_train == 0
    is_good <- df_eval$qcscore_train == 1

    # A cell is predicted as bad if its QS is below the chosen threshold.
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
        # `1 - QC_score` is used so that larger values correspond to worse cells.
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

# Create train/test indices for a single random split.
#
# `train_fraction` is only used when `k_folds = 1`. It controls the proportion
# of cells assigned to the training set, allowing e.g. 50/50, 70/30, 80/20.
make_random_split <- function(n_cells, train_fraction, seed) {
    stopifnot(length(n_cells) == 1L, n_cells > 1L)
    stopifnot(length(train_fraction) == 1L, !is.na(train_fraction))

    if (train_fraction <= 0 || train_fraction >= 1) {
        stop("train_fraction must be strictly between 0 and 1")
    }

    set.seed(seed)

    # Sample the training cells at random.
    n_train <- floor(n_cells * train_fraction)
    if (n_train < 1L || n_train >= n_cells) {
        stop("train_fraction produces an empty training or test set")
    }

    train_idx <- sample(seq_len(n_cells), size=n_train, replace=FALSE)
    test_idx <- setdiff(seq_len(n_cells), train_idx)

    list(train_idx=train_idx, test_idx=test_idx)
}

# Assign each cell to one of the k folds used in cross-validation.
#
# Each fold acts as test set once, while the remaining folds are used for
# training. This function is only used when `k_folds > 1`.
make_fold_ids <- function(n_cells, k_folds, seed) {
    stopifnot(length(n_cells) == 1L, n_cells > 1L)
    stopifnot(length(k_folds) == 1L, !is.na(k_folds), k_folds >= 1)

    k_folds <- as.integer(k_folds)
    if (k_folds < 1L) {
        stop("k_folds must be >= 1")
    }

    set.seed(seed)

    if (k_folds > n_cells) {
        stop("k_folds cannot be greater than the number of cells")
    }

    # Shuffle cells and then distribute them approximately evenly across folds.
    shuffled <- sample(seq_len(n_cells))
    fold_id <- integer(n_cells)
    fold_id[shuffled] <- rep(seq_len(k_folds), length.out=n_cells)
    fold_id
}

# Train the QS model on one training split/fold and evaluate it on both train
# and test cells.
#
# This is the core routine used by both single-split validation and k-fold CV.
run_one_fold <- function(spe, test_idx, threshold=0.5, verbose=TRUE) {
    # Partition the input SpatialExperiment into training and test subsets.
    spe_train <- spe[, -test_idx]
    spe_test <- spe[, test_idx]

    # Compute outliers only on training cells, so the model is fitted without
    # information leakage from the test set.
    spe_train <- computeOutliersQCScore(spe_train)
    spe_train <- checkOutliers(spe_train, verbose=verbose)

    # Extract the training formula selected by SpaceTrooper.
    formula_vars <- S4Vectors::metadata(spe_train)$formula_variables
    model_formula <- getModelFormula(formula_vars, verbose=verbose)

    # Create the balanced training table of pseudo-bad and pseudo-good cells.
    df_train <- computeTrainDF(
        colData=SummarizedExperiment::colData(spe_train),
        formulaVars=formula_vars,
        tech=S4Vectors::metadata(spe_train)$technology,
        verbose=verbose
    )

    # Build the design matrix and estimate the ridge penalty.
    x_train <- stats::model.matrix(stats::as.formula(model_formula), data=df_train)
    best_lambda <- computeLambda(df_train, model_formula)

    # Fit the ridge logistic model.
    fit <- trainModel(x_train, df_train)

    # Score both the training and test subsets using the trained model.
    spe_train <- score_subset(spe_train, fit, best_lambda, model_formula)
    spe_test <- score_subset(spe_test, fit, best_lambda, model_formula)

    # Training evaluation uses the very same labelled examples used for fitting.
    train_eval <- merge(
        data.frame(cell_id=spe_train$cell_id, QC_score=spe_train$QC_score),
        df_train[, c("cell_id", "qcscore_train")],
        by="cell_id"
    )

    # On the test set, build pseudo-reference labels independently.
    test_ref <- safe_reference_labels(spe_test, verbose=verbose)
    test_eval <- if (is.null(test_ref$data)) {
        NULL
    } else {
        merge(
            data.frame(cell_id=spe_test$cell_id, QC_score=spe_test$QC_score),
            test_ref$data,
            by="cell_id"
        )
    }

    # Return both summaries and raw objects/tables for further inspection.
    list(
        formula_variables=formula_vars,
        model_formula=model_formula,
        lambda=best_lambda,
        train_class_balance=table(df_train$qcscore_train),
        train_label_error=NA_character_,
        train_summary=summarise_eval(train_eval, threshold=threshold),
        test_summary=summarise_eval(test_eval, threshold=threshold),
        test_label_error=test_ref$error,
        spe_train=spe_train,
        spe_test=spe_test,
        train_eval=train_eval,
        test_eval=test_eval
    )
}

# Combine per-fold summaries into a single table and compute their mean.
#
# This is only used in k-fold mode. The output contains:
# - `per_fold`: one row per fold
# - `average`: average across folds for each metric
combine_fold_summaries <- function(fold_results, summary_name) {
    summaries <- lapply(seq_along(fold_results), function(i) {
        x <- fold_results[[i]][[summary_name]]
        if (is.null(x) || nrow(x) == 0L) {
            return(NULL)
        }
        x$fold <- i
        x
    })
    summaries <- Filter(Negate(is.null), summaries)

    if (length(summaries) == 0L) {
        return(list(per_fold=data.frame(), average=data.frame()))
    }

    per_fold <- do.call(rbind, summaries)
    metric_cols <- setdiff(colnames(per_fold), "fold")

    # Average each numeric summary metric across folds.
    average <- as.data.frame(
        lapply(per_fold[, metric_cols, drop=FALSE], function(x) mean(x, na.rm=TRUE))
    )

    list(per_fold=per_fold, average=average)
}

# Main entry point for evaluating SpaceTrooper QS on a CosMx dataset.
#
# Modes:
# - `k_folds = 1`: single random split using `train_fraction`
# - `k_folds > 1`: k-fold cross-validation
#
# Important notes:
# - `train_fraction` is ignored when `k_folds > 1`
# - `threshold` is used only for classification summaries, not for model fitting
evaluate_qs_split <- function(
    data_dir="/Users/inzirio/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2",
    sample_name="CosMx_Case2",
    seed=1713,
    threshold=0.5,
    k_folds=1,
    train_fraction=0.5,
    verbose=TRUE
) {
    # Validate the high-level input parameters.
    stopifnot(length(k_folds) == 1L, !is.na(k_folds), k_folds >= 1)
    k_folds <- as.integer(k_folds)
    stopifnot(length(train_fraction) == 1L, !is.na(train_fraction))

    if (train_fraction <= 0 || train_fraction >= 1) {
        stop("train_fraction must be strictly between 0 and 1")
    }

    # Read the CosMx dataset from disk.
    spe <- readCosmxSPE(data_dir, sampleName=sample_name)

    # Compute the per-cell QC metrics required by the QS model.
    spe <- spatialPerCellQC(
        spe,
        rmZeros=TRUE,
        negProbList=c("NegPrb", "Negative", "SystemControl")
    )

    # Count the number of cells remaining after QC preprocessing.
    n_cells <- ncol(spe)

    if (k_folds == 1L) {
        # Single random split: build train/test indices according to the
        # requested training fraction.
        split_idx <- make_random_split(
            n_cells=n_cells,
            train_fraction=train_fraction,
            seed=seed
        )

        res <- run_one_fold(
            spe=spe,
            test_idx=split_idx$test_idx,
            threshold=threshold,
            verbose=verbose
        )

        return(c(
            list(
                split_type="random_split",
                k_folds=k_folds,
                seed=seed,
                threshold=threshold,
                train_fraction=train_fraction,
                test_fraction=1 - train_fraction,
                n_cells=n_cells,
                n_train=length(split_idx$train_idx),
                n_test=length(split_idx$test_idx)
            ),
            res
        ))
    }

    # K-fold mode: assign each cell to a fold once.
    fold_id <- make_fold_ids(n_cells=n_cells, k_folds=k_folds, seed=seed)

    fold_results <- lapply(seq_len(k_folds), function(i) {
        if (verbose) {
            message("Running fold ", i, " of ", k_folds)
        }

        # Fold i is the test set; all other folds are the training set.
        run_one_fold(
            spe=spe,
            test_idx=which(fold_id == i),
            threshold=threshold,
            verbose=verbose
        )
    })

    train_combined <- combine_fold_summaries(fold_results, "train_summary")
    test_combined <- combine_fold_summaries(fold_results, "test_summary")
    test_label_errors <- vapply(
        fold_results,
        extract_optional_chr1,
        character(1),
        name = "test_label_error"
    )

    train_label_errors <- vapply(
        fold_results,
        extract_optional_chr1,
        character(1),
        name = "train_label_error"
    )

    # if you only want real errors later, drop missing values explicitly
    test_label_errors <- stats::na.omit(test_label_errors)
    train_label_errors <- stats::na.omit(train_label_errors)

    list(
        split_type="k_fold_cv",
        k_folds=k_folds,
        seed=seed,
        threshold=threshold,
        n_cells=n_cells,
        fold_assignment=fold_id,
        train_summary_per_fold=train_combined$per_fold,
        train_summary_average=train_combined$average,
        test_summary_per_fold=test_combined$per_fold,
        test_summary_average=test_combined$average,
        test_label_errors=test_label_errors,
        train_label_errors=train_label_errors,
        folds=fold_results
    )
}

# Esempi di utilizzo:
# source("inst/scripts/evaluate_qs_split_kfold.R")
#
# Split singolo 50/50
# res_split_50_50 <- evaluate_qs_split(k_folds=1, train_fraction=0.50)
# res_split_50_50$test_summary
#
# Split singolo 70/30
# res_split_70_30 <- evaluate_qs_split(k_folds=1, train_fraction=0.70)
# res_split_70_30$test_summary
res_split <- evaluate_qs_split(k_folds=1, train_fraction=0.70)

plot_labelled_cells(res_split, dataset="train", size=0.08)
plot_labelled_cells(res_split, dataset="test", size=0.08)

plot_roc_eval(res_split, dataset="test")
plot_auc_across_folds(res_split, dataset="test")
#
# Split singolo 80/20
# res_split_80_20 <- evaluate_qs_split(k_folds=1, train_fraction=0.80)
# res_split_80_20$test_summary
#
# 5-fold cross-validation
res_kfold_5 <- evaluate_qs_split(k_folds=5)
res_kfold_5$test_summary_per_fold
res_kfold_5$test_summary_average
res_kfold_5$test_label_errors
#
# Celle colorate good/bad sul test set del fold 1
plot_labelled_cells(res_kfold_5, dataset="test", fold=1, size=0.1)
#
# ROC curve across folds
plot_roc_eval(res_kfold_5, dataset="test")
#
# AUC across folds
plot_auc_across_folds(res_kfold_5, dataset="test")
