library(testthat)
library(SpaceTrooper)

# Load the example SpatialExperiment directly so tests do not depend on an
# installed help database.
cosmx_path <- system.file(
    "extdata",
    "CosMx_DBKero_Tiny",
    package="SpaceTrooper"
)
spe0 <- readCosmxSPE(cosmx_path, sampleName="DBKero_Tiny")

test_that("QC functions are exported", {
    expect_true(exists("spatialPerCellQC",   mode = "function"))
    expect_true(exists("computeQScore",       mode = "function"))
    expect_true(exists("computeSpatialOutlier", mode = "function"))
    expect_true(exists("computeQScoreFlags",  mode = "function"))
    expect_true(exists("computeThresholdFlags",  mode = "function"))
    ## deprecated wrappers should still be exported
    expect_true(exists("computeQCScore",      mode = "function"))
    expect_true(exists("computeQCScoreFlags", mode = "function"))
    expect_true(exists("computeOutliersQCScore", mode = "function"))
    expect_true(exists("applyQScoreModel", mode = "function"))
    expect_true(exists("applyQCScoreModel", mode = "function"))
})


test_that("spatialPerCellQC adds per-cell metrics to colData", {
    spe <- spatialPerCellQC(spe0, micronConvFact = 0.15)
    expect_s4_class(spe, "SpatialExperiment")
    cd <- colData(spe)
    # at least these should now be in colData()
    required <- c("sum", "detected", "total", "control_sum", "target_sum",
                  "Area_um", "log2AspectRatio", "log2SignalDensity")
    expect_true(all(required %in% colnames(cd)))
})

test_that("computeQScore adds a QScore between 0 and 1", {
    spe <- spatialPerCellQC(spe0)
    set.seed(42)
    spe2 <- computeQScore(spe)
    cd2 <- colData(spe2)
    expect_true("QScore" %in% colnames(cd2))
    fs <- cd2$QScore
    expect_true(is.numeric(fs))
    expect_true(all(fs[!is.na(fs)] >= 0 & fs[!is.na(fs)] <= 1))
})

test_that("computeQCScore preserves legacy score and model names", {
    spe <- spatialPerCellQC(spe0)
    set.seed(42)
    expect_warning(
        spe2 <- computeQCScore(spe),
        "deprecated"
    )
    expect_true("QC_score" %in% colnames(colData(spe2)))
    expect_false("QScore" %in% colnames(colData(spe2)))
    expect_true("QCScore_model" %in% names(metadata(spe2)))
    expect_false("QScore_model" %in% names(metadata(spe2)))
})

test_that("computeQCScoreFlags preserves legacy flag names", {
    spe <- spe0
    spe$QC_score <- seq(0, 1, length.out=ncol(spe))
    spe$threshold_flags <- rep(c(TRUE, FALSE), length.out=ncol(spe))

    expect_warning(
        flagged <- computeQCScoreFlags(spe, qsThreshold=0.5),
        "deprecated"
    )

    expect_true(all(c(
        "low_qcscore",
        "low_threshold_qcscore"
    ) %in% colnames(colData(flagged))))
    expect_false("low_QScore" %in% colnames(colData(flagged)))
})

test_that("computeSpatialOutlier flags outliers for a chosen metric", {
    spe <- spatialPerCellQC(spe0)
    # use 'Area_um' and both methods
    out <- computeSpatialOutlier(spe, computeBy="Area_um", method="scuttle")
    cd_out <- colData(out)
    expect_true(is.logical(cd_out$Area_um_outlier_sc) ||
                    is(cd_out$Area_um_outlier_sc, "outlier.filter"))
})


test_that("computeQScoreFlags combines filters and returns filter_out", {
    spe <- spatialPerCellQC(spe0)
    set.seed(42)
    spe <- computeQScore(spe)
    ff <- computeThresholdFlags(spe,
                            totalThreshold = 10,
                            ctrlTotRatioThreshold = 0.2)
    cd_ff <- colData(ff)
    # expect logical columns and combined filter_out
    expect_true(all(c("is_zero_counts", "is_ctrl_tot_outlier",
                      "threshold_flags") %in% colnames(cd_ff)))
    expect_true(is.logical(cd_ff$is_zero_counts))
    expect_true(is.logical(cd_ff$is_ctrl_tot_outlier))
    expect_true(is.logical(cd_ff$threshold_flags))
    # filter_out should only be TRUE where all three flags are TRUE
    combined <- (ff$is_zero_counts & ff$is_ctrl_tot_outlier)
    expect_identical(combined, cd_ff$threshold_flags)
})
