library(testthat)
library(SpaceTrooper)

qscore_transfer_data <- function() {
    path <- system.file(
        "extdata",
        "CosMx_DBKero_Tiny",
        package="SpaceTrooper"
    )
    spe <- spatialPerCellQC(
        readCosmxSPE(path, sampleName="DBKero_Tiny")
    )
    split <- floor(ncol(spe) * 0.7)
    list(
        train=spe[, seq_len(split)],
        query=spe[, seq.int(split + 1L, ncol(spe))]
    )
}

test_that("applyQScoreModel transfers a supported custom model", {
    data <- qscore_transfer_data()
    supplied <- paste(
        "~ log2SignalDensity + Area_um +",
        "log2SignalDensity:log2Ctrl_total_ratio"
    )
    set.seed(301)
    trained <- computeQScore(
        data$train,
        modelFormula=supplied
    )
    qscore_model <- metadata(trained)$QScore_model

    applied <- applyQScoreModel(
        data$query,
        qsModel=qscore_model
    )

    expect_true("QScore" %in% colnames(colData(applied)))
    expect_true(all(
        applied$QScore[!is.na(applied$QScore)] >= 0 &
            applied$QScore[!is.na(applied$QScore)] <= 1
    ))
    expect_identical(
        metadata(applied)$QScore_model_applied$model_matrix_colnames,
        qscore_model$model_matrix_colnames
    )
    expect_identical(
        qscore_model$model_matrix_colnames,
        c(
            "(Intercept)",
            "log2SignalDensity",
            "Area_um",
            "log2SignalDensity:log2Ctrl_total_ratio"
        )
    )
})

test_that("deprecated applyQCScoreModel preserves legacy output names", {
    data <- qscore_transfer_data()
    set.seed(302)
    trained <- computeQScore(
        data$train,
        modelFormula=~ log2SignalDensity + Area_um
    )
    qscore_model <- metadata(trained)$QScore_model

    expect_warning(
        applied <- applyQCScoreModel(
            data$query,
            qcModel=qscore_model
        ),
        "deprecated"
    )

    expect_true("QC_score" %in% colnames(colData(applied)))
    expect_false("QScore" %in% colnames(colData(applied)))
    expect_true("QCScore_model_applied" %in% names(metadata(applied)))
    expect_false("QScore_model_applied" %in% names(metadata(applied)))
})

test_that("model transfer rejects technology-incompatible border terms", {
    data <- qscore_transfer_data()
    border_formula <- paste0(
        "~ log2SignalDensity + ",
        "I(abs(log2AspectRatio) * as.numeric(dist_border < 50))"
    )
    set.seed(303)
    trained <- computeQScore(
        data$train,
        modelFormula=border_formula
    )

    metadata(data$query)$technology <- "10X_Xenium"
    expect_error(
        applyQScoreModel(
            data$query,
            qsModel=metadata(trained)$QScore_model
        ),
        "supported only for Nanostring CosMx"
    )
})

test_that("computeTrainDF and trainModel accept legacy training response", {
    data <- qscore_transfer_data()
    outliers <- computeOutliersQScore(data$train)
    outliers <- checkOutliers(outliers)
    set.seed(304)
    train_df <- computeTrainDF(
        colData(outliers),
        metadata(outliers)$formula_variables,
        metadata(outliers)$technology
    )

    expect_identical(
        train_df$QScore_train,
        train_df$qcscore_train
    )

    matrix <- stats::model.matrix(
        ~ log2SignalDensity + Area_um,
        data=train_df
    )
    legacy_only <- train_df
    legacy_only$QScore_train <- NULL
    expect_s3_class(
        trainModel(matrix, legacy_only),
        "glmnet"
    )

    conflicting <- train_df
    conflicting$qcscore_train[1] <-
        1L - conflicting$QScore_train[1]
    expect_error(
        trainModel(matrix, conflicting),
        "different values"
    )
})

test_that("stored coefficient metadata contains the complete vector", {
    data <- qscore_transfer_data()
    set.seed(305)
    trained <- computeQScore(
        data$train,
        modelFormula=~ log2SignalDensity + Area_um
    )
    model <- metadata(trained)$QScore_model

    expect_equal(
        nrow(model$coefficients),
        length(model$model_matrix_colnames) + 1L
    )
    expect_equal(
        nrow(model$coefficients_table),
        nrow(model$coefficients)
    )
})
