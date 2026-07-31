library(testthat)
library(SpaceTrooper)

qscore_formula_spe <- function() {
    path <- system.file(
        "extdata",
        "CosMx_DBKero_Tiny",
        package="SpaceTrooper"
    )
    spatialPerCellQC(
        readCosmxSPE(path, sampleName="DBKero_Tiny")
    )
}

model_terms <- function(spe) {
    attr(
        stats::terms(stats::as.formula(
            metadata(spe)$QScore_model$model_formula
        )),
        "term.labels"
    )
}

captured_training_metrics <- function(spe, model_formula, seed=206) {
    seen <- new.env(parent=emptyenv())
    original_prep <- SpaceTrooper:::.prepQCContext
    testthat::local_mocked_bindings(
        .prepQCContext=function(spe, metricList, verbose=FALSE) {
            seen$metricList <- metricList
            original_prep(spe, metricList, verbose)
        },
        .package="SpaceTrooper"
    )
    set.seed(seed)
    scored <- computeQScore(
        spe,
        bestLambda=0.01,
        modelFormula=model_formula
    )
    list(metrics=seen$metricList, scored=scored)
}

test_that("supported Quality Score predictors are explicit and stable", {
    expect_identical(
        SpaceTrooper:::.qscoreSupportedPredictors(),
        c(
            "log2SignalDensity",
            "Area_um",
            "log2AspectRatio",
            "log2Ctrl_total_ratio"
        )
    )
})

test_that("default formula keeps the established pairwise interactions", {
    spe <- qscore_formula_spe()
    set.seed(201)
    scored <- computeQScore(spe)
    terms <- model_terms(scored)

    expect_true("log2SignalDensity" %in% terms)
    expect_true("Area_um" %in% terms)
    expect_true(any(grepl(":", terms, fixed=TRUE)))
})

test_that("additive custom formula remains additive", {
    spe <- qscore_formula_spe()
    supplied <- "~ log2SignalDensity + Area_um"
    set.seed(202)
    scored <- computeQScore(spe, modelFormula=supplied)

    expect_identical(
        metadata(scored)$QScore_model$model_formula,
        supplied
    )
    expect_identical(
        model_terms(scored),
        c("log2SignalDensity", "Area_um")
    )
    expect_identical(
        metadata(scored)$QScore_model$model_matrix_colnames,
        c("(Intercept)", "log2SignalDensity", "Area_um")
    )
})

test_that("custom formula preserves selected interactions and column order", {
    spe <- qscore_formula_spe()
    supplied <- paste(
        "~ log2SignalDensity + Area_um +",
        "log2SignalDensity:log2Ctrl_total_ratio"
    )
    set.seed(203)
    scored <- computeQScore(spe, modelFormula=supplied)

    expect_identical(
        model_terms(scored),
        c(
            "log2SignalDensity",
            "Area_um",
            "log2SignalDensity:log2Ctrl_total_ratio"
        )
    )
    expect_identical(
        metadata(scored)$QScore_model$model_matrix_colnames,
        c(
            "(Intercept)",
            "log2SignalDensity",
            "Area_um",
            "log2SignalDensity:log2Ctrl_total_ratio"
        )
    )
})

test_that("custom formula is not rebuilt by getModelFormula", {
    spe <- qscore_formula_spe()
    testthat::local_mocked_bindings(
        getModelFormula=function(...) stop("getModelFormula was called"),
        .package="SpaceTrooper"
    )

    set.seed(204)
    expect_no_error(computeQScore(
        spe,
        modelFormula=~ log2SignalDensity + Area_um
    ))
})

test_that("custom formulas define metrics passed to QScore preparation", {
    spe <- qscore_formula_spe()

    single <- captured_training_metrics(
        spe,
        ~ log2SignalDensity,
        seed=206
    )
    expect_identical(single$metrics, "log2SignalDensity")

    additive <- captured_training_metrics(
        spe,
        ~ log2SignalDensity + Area_um,
        seed=207
    )
    expect_identical(
        additive$metrics,
        c("log2SignalDensity", "Area_um")
    )

    interaction <- captured_training_metrics(
        spe,
        ~ log2SignalDensity * Area_um,
        seed=208
    )
    expect_identical(
        interaction$metrics,
        c("log2SignalDensity", "Area_um")
    )
    expect_identical(
        model_terms(interaction$scored),
        c(
            "log2SignalDensity",
            "Area_um",
            "log2SignalDensity:Area_um"
        )
    )

    selected <- captured_training_metrics(
        spe,
        ~ log2SignalDensity +
            log2SignalDensity:log2Ctrl_total_ratio,
        seed=209
    )
    expect_identical(
        selected$metrics,
        c("log2SignalDensity", "log2Ctrl_total_ratio")
    )
})

test_that("default formula keeps the established training metrics", {
    spe <- qscore_formula_spe()
    default <- captured_training_metrics(spe, NULL, seed=210)

    expect_identical(
        default$metrics,
        SpaceTrooper:::.qscoreSupportedPredictors()
    )
})

test_that("supported CosMx border expression is accepted", {
    spe <- qscore_formula_spe()
    supplied <- paste0(
        "~ log2SignalDensity + ",
        "I(abs(log2AspectRatio) * as.numeric(dist_border < 50))"
    )
    set.seed(205)
    scored <- computeQScore(spe, modelFormula=supplied)

    expect_identical(
        metadata(scored)$QScore_model$model_formula,
        supplied
    )
    formula_info <- SpaceTrooper:::.validateQScoreFormula(
        modelFormula=supplied,
        dataNames=names(colData(spe)),
        technology=metadata(spe)$technology
    )
    expect_identical(
        formula_info$training_metrics,
        c("log2SignalDensity", "log2AspectRatio")
    )
    expect_false("dist_border" %in% formula_info$training_metrics)
})

test_that("custom metric subsets change deterministic training labels", {
    n_cells <- 100L
    training_data <- data.frame(
        cell_id=paste0("cell", seq_len(n_cells)),
        log2SignalDensity=seq_len(n_cells),
        Area_um=seq_len(n_cells),
        log2SignalDensity_outlier_train=c(
            "LOW",
            rep("NO", n_cells - 1L)
        ),
        Area_um_outlier_sc=c(
            "NO",
            "HIGH",
            rep("NO", n_cells - 2L)
        )
    )
    outlier_columns <- c(
        log2SignalDensity="log2SignalDensity_outlier_train",
        Area_um="Area_um_outlier_sc"
    )
    formula_metrics <- function(model_formula) {
        SpaceTrooper:::.validateQScoreFormula(
            modelFormula=model_formula,
            dataNames=names(training_data),
            technology="10X_Xenium"
        )$training_metrics
    }

    signal_metrics <- formula_metrics(~ log2SignalDensity)
    subset_metrics <- formula_metrics(
        ~ log2SignalDensity + Area_um
    )

    set.seed(211)
    signal_train <- computeTrainDF(
        training_data,
        outlier_columns[signal_metrics],
        tech="10X_Xenium"
    )
    set.seed(211)
    subset_train <- computeTrainDF(
        training_data,
        outlier_columns[subset_metrics],
        tech="10X_Xenium"
    )

    expect_identical(
        sort(signal_train$cell_id[signal_train$QScore_train == 0]),
        "cell1"
    )
    expect_identical(
        sort(subset_train$cell_id[subset_train$QScore_train == 0]),
        c("cell1", "cell2")
    )
})

test_that("unsupported predictors and transformations are rejected", {
    spe <- qscore_formula_spe()

    expect_error(
        computeQScore(
            spe,
            modelFormula=~ log2SignalDensity + customMetric
        ),
        "Unsupported Quality Score predictor.*customMetric"
    )
    expect_error(
        computeQScore(
            spe,
            modelFormula=~ log2SignalDensity + log1p(Area_um)
        ),
        "Unsupported Quality Score formula term.*log1p"
    )

    spe$customMetric <- seq_len(ncol(spe))
    expect_error(
        computeOutliersQScore(
            spe,
            metricList=c("log2SignalDensity", "customMetric")
        ),
        "Unsupported Quality Score metric.*customMetric"
    )
})

test_that("custom formulas retain the signal-density requirement", {
    spe <- qscore_formula_spe()

    expect_error(
        computeQScore(spe, modelFormula=~ Area_um),
        "'log2SignalDensity' is required.*training labels"
    )
})

test_that("missing and technology-incompatible variables are reported", {
    spe <- qscore_formula_spe()
    spe$Area_um <- NULL
    expect_error(
        computeQScore(
            spe,
            modelFormula=~ log2SignalDensity + Area_um
        ),
        "Missing variables required.*Area_um"
    )

    xenium_like <- qscore_formula_spe()
    metadata(xenium_like)$technology <- "10X_Xenium"
    expect_error(
        computeQScore(
            xenium_like,
            modelFormula=~ log2SignalDensity + log2AspectRatio
        ),
        "supported only for Nanostring CosMx"
    )
})
