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

test_that("custom fit formula does not redefine training-label metrics", {
    spe <- qscore_formula_spe()
    seen <- new.env(parent=emptyenv())
    original_prep <- SpaceTrooper:::.prepQCContext
    testthat::local_mocked_bindings(
        .prepQCContext=function(spe, metricList, verbose=FALSE) {
            seen$metricList <- metricList
            original_prep(spe, metricList, verbose)
        },
        .package="SpaceTrooper"
    )

    set.seed(206)
    computeQScore(
        spe,
        modelFormula=~ log2SignalDensity
    )

    expect_identical(
        seen$metricList,
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
