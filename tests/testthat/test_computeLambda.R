library(testthat)
library(SpaceTrooper)

lambda_train_data <- function(response_name="QScore_train") {
    set.seed(101)
    train_df <- data.frame(
        x1=rnorm(80),
        x2=rnorm(80),
        x3=rnorm(80)
    )
    response <- as.integer(
        train_df$x1 - 0.5 * train_df$x2 + rnorm(80) > 0
    )
    train_df[[response_name]] <- response
    train_df
}

test_that("computeLambda preserves historical positional and named calls", {
    train_df <- lambda_train_data("qcscore_train")
    model_formula <- "~ x1 + x2"

    set.seed(11)
    positional <- computeLambda(train_df, model_formula)
    set.seed(11)
    named <- computeLambda(
        trainDF=train_df,
        modelFormula=model_formula
    )

    expect_type(positional, "double")
    expect_length(positional, 1L)
    expect_true(is.finite(positional))
    expect_identical(positional, named)
})

test_that("computeLambda accepts canonical and legacy responses", {
    canonical <- lambda_train_data("QScore_train")
    legacy <- canonical
    names(legacy)[names(legacy) == "QScore_train"] <- "qcscore_train"

    set.seed(22)
    canonical_lambda <- computeLambda(canonical, ~ x1 + x2)
    set.seed(22)
    legacy_lambda <- computeLambda(legacy, ~ x1 + x2)
    expect_identical(canonical_lambda, legacy_lambda)

    both <- canonical
    both$qcscore_train <- both$QScore_train
    set.seed(22)
    expect_identical(
        computeLambda(both, ~ x1 + x2),
        canonical_lambda
    )

    both$qcscore_train[1] <- 1L - both$QScore_train[1]
    expect_error(
        computeLambda(both, ~ x1 + x2),
        "different values"
    )
})

test_that("computeLambda filters complete cases", {
    train_df <- lambda_train_data()
    train_df$x2[c(2, 7)] <- NA_real_
    complete_df <- train_df[complete.cases(
        train_df[, c("x1", "x2", "QScore_train")]
    ), ]

    set.seed(33)
    expect_warning(
        with_missing <- computeLambda(train_df, ~ x1 + x2),
        "2 training cells"
    )
    set.seed(33)
    without_missing <- computeLambda(complete_df, ~ x1 + x2)

    expect_identical(with_missing, without_missing)
})

test_that("computeLambda preserves additive and selected interaction formulas", {
    train_df <- lambda_train_data()

    set.seed(44)
    additive <- computeLambda(train_df, ~ x1 + x2 + x3)
    set.seed(44)
    selected_interaction <- computeLambda(
        train_df,
        ~ x1 + x2 + x1:x3
    )

    expect_true(is.finite(additive))
    expect_true(is.finite(selected_interaction))
})

test_that("public computeLambda agrees with matrix implementation and glmnet", {
    train_df <- lambda_train_data()
    model_formula <- ~ x1 + x2 + x1:x2
    model_matrix <- stats::model.matrix(model_formula, data=train_df)
    response <- train_df$QScore_train

    set.seed(55)
    public <- computeLambda(train_df, model_formula)
    set.seed(55)
    internal <- SpaceTrooper:::.computeLambda(
        modelMatrix=model_matrix,
        response=response
    )
    set.seed(55)
    reference <- glmnet::cv.glmnet(
        model_matrix,
        response,
        family="binomial",
        alpha=0,
        lambda=NULL
    )$lambda.min

    expect_identical(public, internal)
    expect_identical(public, reference)
})
