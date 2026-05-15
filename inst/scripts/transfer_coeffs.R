library(SpaceTrooper)
library(SummarizedExperiment)
library(S4Vectors)

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

xen_model <- metadata(spexen)$QCScore_model


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

cosmx_model <- metadata(specosm)$QCScore_model


## -----------------------------
## Save native scores before transfer
## -----------------------------

spexen$QC_score_xenium_native <- spexen$QC_score
specosm$QC_score_cosmx_native <- specosm$QC_score


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
