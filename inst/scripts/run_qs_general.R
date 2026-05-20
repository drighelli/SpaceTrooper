
package_path <- "/Users/inzirio/My Drive/works/coding/SpaceTrooper"
input_rds <- "~/SpaceTrooper_qs_check/cosmx_case2_spe_raw.rds"
output_rds <- "~/SpaceTrooper_qs_check/cosmx_case2_qs_result_interaction.rds"


message("Package path: ", package_path)
message("Input RDS: ", input_rds)
message("Output RDS: ", output_rds)

library(devtools)
library(SpatialExperiment)

devtools::load_all(package_path, quiet=FALSE)


obj <- readRDS(input_rds)
metadata(obj)$formula_variables <- c(log2CountArea="log2CountArea",
                                        log2AspectRatio="log2AspectRatio")
obj <- spatialPerCellQC(obj)

metadata(obj)
obj <- computeQCScore(obj)

qc_nointeract <- obj$QC_score
qc_interact <- obj$QC_score
saveRDS(qc_nointeract, "~/SpaceTrooper_qs_check/qc_nointeract.rds")
saveRDS(qc_interact, "~/SpaceTrooper_qs_check/qc_interact.rds")

library(ggplot2)

# assume qc_nointeract and qc_interact exist in the environment
df <- data.frame(qc_nointeract = qc_nointeract, qc_interact = qc_interact)

# correlation
cor_val <- cor(df$qc_nointeract, df$qc_interact, use = "complete.obs", method = "pearson")

# plot: use 2D binning for large datasets, points otherwise; add y=x line and annotation
plot_qc <- ggplot(df, aes(x = qc_nointeract, y = qc_interact)) +
    geom_point(alpha = 0.35, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "QC score (no interaction)",
    y = "QC score (with interaction)",
    title = "QC: no-interaction vs interaction"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0.75) +
  geom_vline(xintercept = 0.75) +
  annotate(
    "text",
    x = quantile(df$qc_nointeract, 0.99, na.rm = TRUE),
    y = quantile(df$qc_interact, 0.01, na.rm = TRUE),
    label = paste0("r = ", formatC(cor_val, digits = 3, format = "f")),
    hjust = 1, vjust = 0
  )

# return the plot object
plot_qc


library(ggplot2)

df <- data.frame(
    qc_nointeract=qc_nointeract,
    qc_interact=qc_interact
)

threshold <- 0.75

df$quadrant <- with(
    df,
    ifelse(
        qc_nointeract < threshold & qc_interact < threshold,
        "low / low",
        ifelse(
            qc_nointeract >= threshold & qc_interact < threshold,
            "high no-interaction / low interaction",
            ifelse(
                qc_nointeract < threshold & qc_interact >= threshold,
                "low no-interaction / high interaction",
                "high / high"
            )
        )
    )
)

table(df$quadrant)

cor_val <- cor(
    df$qc_nointeract,
    df$qc_interact,
    use="complete.obs",
    method="pearson"
)

plot_qc <- ggplot(df, aes(x=qc_nointeract, y=qc_interact, color=quadrant)) +
    geom_point(alpha=0.35, size=0.6) +
    geom_abline(
        slope=1,
        intercept=0,
        linetype="dashed",
        color="red"
    ) +
    geom_hline(yintercept=threshold) +
    geom_vline(xintercept=threshold) +
    labs(
        x="QC score (no interaction)",
        y="QC score (with interaction)",
        title="QC: no-interaction vs interaction",
        color="Quadrant"
    ) +
    theme_minimal() +
    annotate(
        "text",
        x=quantile(df$qc_nointeract, 0.99, na.rm=TRUE),
        y=quantile(df$qc_interact, 0.01, na.rm=TRUE),
        label=paste0("r = ", formatC(cor_val, digits=3, format="f")),
        hjust=1,
        vjust=0,
        color="black"
    )

plot_qc

plot_qc <- plot_qc +
    scale_color_manual(
        values=c(
            "low / low"="grey60",
            "high no-interaction / low interaction"="orange",
            "low no-interaction / high interaction"="dodgerblue",
            "high / high"="black"
        )
    )

ggsave(
    filename="~/SpaceTrooper_qs_check/qc_nointeraction_vs_interaction_A4_landscape.pdf",
    plot=plot_qc,
    device="pdf",
    width=11.69,
    height=8.27,
    units="in"
)
# get_coldata_df <- function(x) {
#     cd <- SummarizedExperiment::colData(x)
#     as.data.frame(cd)
# }

# get_qscore_from_coldata <- function(x) {
#     cd <- get_coldata_df(x)

#     candidate_names <- c(
#         "quality_score",
#         "QualityScore",
#         "qualityScore",
#         "qscore",
#         "QScore",
#         "Q_score",
#         "q_score",
#         "score",
#         "flag_score",
#         "qc_score",
#         "QC_score",
#         "QCScore",
#         "cell_quality_score",
#         "spatial_quality_score"
#     )

#     exact_hits <- intersect(candidate_names, colnames(cd))

#     if (length(exact_hits) > 0) {
#         selected <- exact_hits[1]
#         return(list(
#             values=cd[[selected]],
#             column=selected,
#             method="exact_name_match",
#             all_candidates=exact_hits,
#             coldata=cd
#         ))
#     }

#     pattern_hits <- grep(
#         "quality|qscore|q_score|qc_score|flag_score|score",
#         colnames(cd),
#         ignore.case=TRUE,
#         value=TRUE
#     )

#     numeric_pattern_hits <- pattern_hits[
#         vapply(cd[pattern_hits], is.numeric, logical(1))
#     ]

#     if (length(numeric_pattern_hits) > 0) {
#         selected <- numeric_pattern_hits[1]
#         return(list(
#             values=cd[[selected]],
#             column=selected,
#             method="pattern_numeric_match",
#             all_candidates=numeric_pattern_hits,
#             coldata=cd
#         ))
#     }

#     stop(
#         "Could not automatically identify Quality Score column.\n",
#         "Available colData columns are:\n",
#         paste(colnames(cd), collapse=", ")
#     )
# }

# call_function_if_exists <- function(fun_name, x) {
#     if (!exists(fun_name, mode="function")) {
#         stop("Function not found in loaded SpaceTrooper branch: ", fun_name)
#     }

#     fun <- get(fun_name, mode="function")
#     fun(x)
# }

# obj_before <- obj

# if (run_spatial_qc) {
#     message("Running spatial QC function...")
#     obj <- call_function_if_exists(spatial_qc_fun, obj)
# }

# message("Running Quality Score function...")
# obj_after <- call_function_if_exists(qscore_fun, obj)

# message("Extracting Quality Score...")
# qs_info <- get_qscore_from_coldata(obj_after)

# qs <- qs_info$values

# if (!is.numeric(qs)) {
#     warning("Detected Quality Score column is not numeric: ", qs_info$column)
# }

# out <- list(
#     git_branch=git_branch,
#     git_hash=git_hash,
#     package_path=normalizePath(package_path),
#     input_rds=normalizePath(input_rds),
#     spatial_qc_fun=spatial_qc_fun,
#     qscore_fun=qscore_fun,
#     run_spatial_qc=run_spatial_qc,
#     qscore_column=qs_info$column,
#     qscore_detection_method=qs_info$method,
#     qscore_candidates=qs_info$all_candidates,
#     n=length(qs),
#     object_class=class(obj_after),
#     object_dim=tryCatch(dim(obj_after), error=function(e) NULL),
#     object_colnames=tryCatch(colnames(obj_after), error=function(e) NULL),
#     quality_score=qs,
#     quality_score_summary=summary(qs),
#     colData=qs_info$coldata,
#     session_info=sessionInfo()
# )

# saveRDS(out, output_rds)

# message("Saved result to: ", output_rds)
# message("Detected Quality Score column: ", qs_info$column)
# message("Detection method: ", qs_info$method)
# message("Number of values: ", length(qs))
# print(summary(qs))


# computeQCScore(obj)

# # model_formula <- paste0("~ log2CountArea + I(abs(log2AspectRatio) ",
# #                                 "* as.numeric(dist_border<50)) + ",
# #                                 " log2CountArea:I(abs(log2AspectRatio)",
# #                                 "* as.numeric(dist_border<50))") #for cosmx
