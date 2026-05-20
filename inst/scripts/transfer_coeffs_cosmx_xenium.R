library(SpaceTrooper)
# library(SummarizedExperiment)
# library(S4Vectors)

set.seed(1998)

common_formula <- "~(log2SignalDensity + Area_um + log2AspectRatio)^2"

## -----------------------------
## Read Xenium data
## -----------------------------

dbkx_path <- "~/Downloads/Xenium_data/db_kero_xen"

spexen <- readXeniumSPE(dbkx_path)

spexen <- spatialPerCellQC(spexen)

spexen <- computeQCScore(
    spexen,
    modelFormula=common_formula,
    verbose=TRUE
)

sum(is.na(spexen$QC_score))
colSums(is.na(as.data.frame(colData(spexen))[
    , c("log2SignalDensity", "Area_um", "log2AspectRatio")
]))

length(spexen$QC_score) == ncol(spexen)

sum(is.na(spexen$QC_score))
summary(spexen$QC_score)

metadata(spexen)$QCScore_model$model_formula
metadata(spexen)$QCScore_model$model_matrix_colnames
metadata(spexen)$QCScore_model$coefficients_table


## -----------------------------
## Read CosMx data
## -----------------------------

cosm_path <- "~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/"

specosm <- readCosmxSPE(cosm_path)

specosm <- spatialPerCellQC(specosm)

specosm <- computeQCScore(
    specosm,
    modelFormula=common_formula,
    verbose=TRUE
)

length(specosm$QC_score) == ncol(specosm)

sum(is.na(specosm$QC_score))
summary(specosm$QC_score)

metadata(specosm)$QCScore_model$model_formula
metadata(specosm)$QCScore_model$model_matrix_colnames
metadata(specosm)$QCScore_model$coefficients_table

## -----------------------------
## Save native scores before transfer
## -----------------------------

spexen$QC_score_xenium_native <- spexen$QC_score
specosm$QC_score_cosmx_native <- specosm$QC_score

xen_model <- metadata(spexen)$QCScore_model
cosmx_model <- metadata(specosm)$QCScore_model

## -----------------------------
## Apply transferred models
## -----------------------------

specosm <- applyQCScoreModel(
    spe=specosm,
    qcModel=xen_model,
    scoreName="QC_score_xenium_model"
)

spexen <- applyQCScoreModel(
    spe=spexen,
    qcModel=cosmx_model,
    scoreName="QC_score_cosmx_model"
)


## -----------------------------
## Summaries
## -----------------------------

message("Xenium native score:")
print(summary(spexen$QC_score_xenium_native))

message("Xenium scored with CosMx model:")
print(summary(spexen$QC_score_cosmx_model))

message("CosMx native score:")
print(summary(specosm$QC_score_cosmx_native))

message("CosMx scored with Xenium model:")
print(summary(specosm$QC_score_xenium_model))


## -----------------------------
## Correlations
## -----------------------------

message("Correlation in Xenium: native vs CosMx model")
print(
    cor(
        spexen$QC_score_xenium_native,
        spexen$QC_score_cosmx_model,
        use="complete.obs"
    )
)

message("Correlation in CosMx: native vs Xenium model")
print(
    cor(
        specosm$QC_score_cosmx_native,
        specosm$QC_score_xenium_model,
        use="complete.obs"
    )
)


## -----------------------------
## Model details
## -----------------------------

message("Xenium model formula:")
print(metadata(spexen)$QCScore_model$model_formula)

message("CosMx model formula:")
print(metadata(specosm)$QCScore_model$model_formula)

message("Xenium model coefficients:")
print(metadata(spexen)$QCScore_model$coefficients_table)

message("CosMx model coefficients:")
print(metadata(specosm)$QCScore_model$coefficients_table)


## -----------------------------
## Correlations of scores with model variables
## -----------------------------
cor(
    spexen$QC_score_xenium_native,
    spexen$QC_score_cosmx_model,
    use="complete.obs",
    method="spearman"
)

cor(
    specosm$QC_score_cosmx_native,
    specosm$QC_score_xenium_model,
    use="complete.obs",
    method="spearman"
)

range(spexen$QC_score_cosmx_model, na.rm=TRUE)
sum(is.na(spexen$QC_score_cosmx_model))

range(specosm$QC_score_xenium_model, na.rm=TRUE)
sum(is.na(specosm$QC_score_xenium_model))

res_xen <- plot_qc_score_comparison(
    x=spexen$QC_score_xenium_native,
    y=spexen$QC_score_cosmx_model,
    x_label="Native Xenium QC score",
    y_label="CosMx-trained QC score",
    title="Xenium: native vs CosMx-trained QC score",
    cor_method = "pearson",
    threshold=0.25,
    output_file="~/SpaceTrooper_qs_check/xenium_native_vs_cosmx_model.pdf"
)

res_xen$plot
res_xen$correlation
res_xen$quadrant_table

res_cosmx <- plot_qc_score_comparison(
    x=specosm$QC_score_cosmx_native,
    y=specosm$QC_score_xenium_model,
    x_label="Native CosMx QC score",
    y_label="Xenium-trained QC score",
    title="CosMx: native vs Xenium-trained QC score",
    cor_method = "pearson",
    threshold=0.75,
    output_file="~/SpaceTrooper_qs_check/cosmx_native_vs_xenium_model.pdf"
)

res_cosmx$plot
res_cosmx$correlation
res_cosmx$quadrant_table
