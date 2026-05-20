spexen_bc1 <- readXeniumSPE("/Users/inzirio/Downloads/Xenium_data/Xenium_HumanBreast1_Janesick/")
spexen_bc2 <- readXeniumSPE("/Users/inzirio/Downloads/Xenium_data/Xenium_HumanBreast2_Janesick")


spexen_bc1 <- spatialPerCellQC(spexen_bc1)
spexen_bc2 <- spatialPerCellQC(spexen_bc2)

xenium_formula <- paste0(
    "~(",
    "log2SignalDensity + ",
    "Area_um + ",
    "log2AspectRatio + ",
    "log2Ctrl_total_ratio",
    ")^2"
)

spexen_bc1 <- computeQCScore(
    spe=spexen_bc1,
    modelFormula=xenium_formula,
    verbose=TRUE
)

spexen_bc2 <- computeQCScore(
    spe=spexen_bc2,
    modelFormula=xenium_formula,
    verbose=TRUE
)

spexen_bc1$QC_score_bc1_native <- spexen_bc1$QC_score
spexen_bc2$QC_score_bc2_native <- spexen_bc2$QC_score

bc1_model <- metadata(spexen_bc1)$QCScore_model
bc2_model <- metadata(spexen_bc2)$QCScore_model

spexen_bc1 <- applyQCScoreModel(
    spe=spexen_bc1,
    qcModel=bc2_model,
    scoreName="QC_score_bc2_model"
)

spexen_bc2 <- applyQCScoreModel(
    spe=spexen_bc2,
    qcModel=bc1_model,
    scoreName="QC_score_bc1_model"
)

cor_bc1_spearman <- cor(
    spexen_bc1$QC_score_bc1_native,
    spexen_bc1$QC_score_bc2_model,
    use="complete.obs",
    method="spearman"
)

cor_bc1_pearson <- cor(
    spexen_bc1$QC_score_bc1_native,
    spexen_bc1$QC_score_bc2_model,
    use="complete.obs",
    method="pearson"
)

cor_bc2_spearman <- cor(
    spexen_bc2$QC_score_bc2_native,
    spexen_bc2$QC_score_bc1_model,
    use="complete.obs",
    method="spearman"
)

cor_bc2_pearson <- cor(
    spexen_bc2$QC_score_bc2_native,
    spexen_bc2$QC_score_bc1_model,
    use="complete.obs",
    method="pearson"
)

res_bc1 <- plot_qc_score_comparison(
    x=spexen_bc1$QC_score_bc1_native,
    y=spexen_bc1$QC_score_bc2_model,
    x_label="Breast Cancer 1 native QS",
    y_label="Breast Cancer 2-trained QS",
    title="Xenium Breast Cancer 1: native QS vs Breast Cancer 2-trained QS",
    threshold=0.75,
    cor_method="spearman",
    output_file=file.path(
        out_dir,
        "xenium_breast1_native_vs_breast2_model_A4_landscape.pdf"
    )
)

res_bc2 <- plot_qc_score_comparison(
    x=spexen_bc2$QC_score_bc2_native,
    y=spexen_bc2$QC_score_bc1_model,
    x_label="Breast Cancer 2 native QS",
    y_label="Breast Cancer 1-trained QS",
    title="Xenium Breast Cancer 2: native QS vs Breast Cancer 1-trained QS",
    threshold=0.75,
    cor_method="spearman",
    output_file=file.path(
        out_dir,
        "xenium_breast2_native_vs_breast1_model_A4_landscape.pdf"
    )
)

res_bc1$plot
res_bc2$plot

color_vars <- c(
    "log2SignalDensity",
    "Area_um",
    "log2AspectRatio",
    "log2Ctrl_total_ratio"
)

bc1_colored_plots <- lapply(color_vars, function(v) {
    plot_qc_score_coldata(
        x=spexen_bc1$QC_score_bc1_native,
        y=spexen_bc1$QC_score_bc2_model,
        spe=spexen_bc1,
        color_col=v,
        x_label="Breast Cancer 1 native QS",
        y_label="Breast Cancer 2-trained QS",
        title=paste0("Breast Cancer 1 native vs Breast Cancer 2-trained QS: ", v),
        threshold=0.75,
        cor_method="spearman",
        output_file=file.path(
            out_dir,
            paste0("xenium_breast1_native_vs_breast2_model_", v, ".pdf")
        )
    )
})
names(bc1_colored_plots) <- color_vars

bc2_colored_plots <- lapply(color_vars, function(v) {
    plot_qc_score_coldata(
        x=spexen_bc2$QC_score_bc2_native,
        y=spexen_bc2$QC_score_bc1_model,
        spe=spexen_bc2,
        color_col=v,
        x_label="Breast Cancer 2 native QS",
        y_label="Breast Cancer 1-trained QS",
        title=paste0("Breast Cancer 2 native vs Breast Cancer 1-trained QS: ", v),
        threshold=0.75,
        cor_method="spearman",
        output_file=file.path(
            out_dir,
            paste0("xenium_breast2_native_vs_breast1_model_", v, ".pdf")
        )
    )
})
names(bc2_colored_plots) <- color_vars

xenium_breast_transfer_summary <- data.frame(
    comparison=c(
        "Breast Cancer 1 native vs Breast Cancer 2-trained",
        "Breast Cancer 2 native vs Breast Cancer 1-trained"
    ),
    spearman=c(cor_bc1_spearman, cor_bc2_spearman),
    pearson=c(cor_bc1_pearson, cor_bc2_pearson),
    stringsAsFactors=FALSE
)

xenium_breast_transfer_summary

write.csv(
    xenium_breast_transfer_summary,
    file=file.path(out_dir, "xenium_breast1_breast2_transfer_summary.csv"),
    row.names=FALSE
)

metadata(spexen_bc1)$QCScore_model$model_matrix_colnames
metadata(spexen_bc2)$QCScore_model$model_matrix_colnames

########################

dbkx_path <- "~/Downloads/Xenium_data/db_kero_xen"

spexen <- readXeniumSPE(dbkx_path)
spexen_bc1 <- readXeniumSPE("/Users/inzirio/Downloads/Xenium_data/Xenium_HumanBreast1_Janesick/")



spexen <- spatialPerCellQC(spexen)
spexen_bc1 <- spatialPerCellQC(spexen_bc1)

xenium_formula <- paste0(
    "~(",
    "log2SignalDensity + ",
    "Area_um + ",
    "log2AspectRatio + ",
    "log2Ctrl_total_ratio",
    ")^2"
)

spexen <- computeQCScore(
    spe=spexen,
    modelFormula=xenium_formula,
    verbose=TRUE
)

spexen_bc1 <- computeQCScore(
    spe=spexen_bc1,
    modelFormula=xenium_formula,
    verbose=TRUE
)

spexen$QC_score_dbk_native <- spexen$QC_score
spexen_bc1$QC_score_bc1_native <- spexen_bc1$QC_score

dbk_model <- metadata(spexen)$QCScore_model
bc1_model <- metadata(spexen_bc1)$QCScore_model

spexen <- applyQCScoreModel(
    spe=spexen,
    qcModel=bc1_model,
    scoreName="QC_score_bc1_model"
)

spexen_bc1 <- applyQCScoreModel(
    spe=spexen_bc1,
    qcModel=dbk_model,
    scoreName="QC_score_dbk_model"
)

cor_dbk_spearman <- cor(
    spexen$QC_score_dbk_native,
    spexen$QC_score_bc1_model,
    use="complete.obs",
    method="spearman"
)

cor_dbk_pearson <- cor(
    spexen$QC_score_dbk_native,
    spexen$QC_score_bc1_model,
    use="complete.obs",
    method="pearson"
)

cor_bc1_spearman <- cor(
    spexen_bc1$QC_score_bc1_native,
    spexen_bc1$QC_score_dbk_model,
    use="complete.obs",
    method="spearman"
)

cor_bc1_pearson <- cor(
    spexen_bc1$QC_score_bc1_native,
    spexen_bc1$QC_score_dbk_model,
    use="complete.obs",
    method="pearson"
)

res_dbk <- plot_qc_score_comparison(
    x=spexen$QC_score_dbk_native,
    y=spexen$QC_score_bc1_model,
    x_label="DBK Xenium native QS",
    y_label="Breast 1-trained QS",
    title="Xenium DBK: native QS vs Breast 1-trained QS",
    threshold=0.75,
    cor_method="spearman",
    output_file=file.path(
        out_dir,
        "xenium_dbk_native_vs_breast1_model_A4_landscape.pdf"
    )
)

res_bc1 <- plot_qc_score_comparison(
    x=spexen_bc1$QC_score_bc1_native,
    y=spexen_bc1$QC_score_dbk_model,
    x_label="Breast 1 native QS",
    y_label="DBK-trained QS",
    title="Xenium Breast 1: native QS vs DBK-trained QS",
    threshold=0.75,
    cor_method="spearman",
    output_file=file.path(
        out_dir,
        "xenium_breast1_native_vs_dbk_model_A4_landscape.pdf"
    )
)

res_dbk$plot
res_bc1$plot

color_vars <- c(
    "log2SignalDensity",
    "Area_um",
    "log2AspectRatio",
    "log2Ctrl_total_ratio"
)

dbk_colored_plots <- lapply(color_vars, function(v) {
    plot_qc_score_coldata(
        x=spexen$QC_score_dbk_native,
        y=spexen$QC_score_bc1_model,
        spe=spexen,
        color_col=v,
        x_label="DBK Xenium native QS",
        y_label="Breast 1-trained QS",
        title=paste0("DBK Xenium native vs Breast 1-trained QS: ", v),
        threshold=0.75,
        cor_method="spearman",
        output_file=file.path(
            out_dir,
            paste0("xenium_dbk_native_vs_breast1_model_", v, ".pdf")
        )
    )
})
names(dbk_colored_plots) <- color_vars

bc1_colored_plots <- lapply(color_vars, function(v) {
    plot_qc_score_coldata(
        x=spexen_bc1$QC_score_bc1_native,
        y=spexen_bc1$QC_score_dbk_model,
        spe=spexen_bc1,
        color_col=v,
        x_label="Breast 1 native QS",
        y_label="DBK-trained QS",
        title=paste0("Breast 1 native vs DBK-trained QS: ", v),
        threshold=0.75,
        cor_method="spearman",
        output_file=file.path(
            out_dir,
            paste0("xenium_breast1_native_vs_dbk_model_", v, ".pdf")
        )
    )
})
names(bc1_colored_plots) <- color_vars

xenium_dbk_breast1_transfer_summary <- data.frame(
    comparison=c(
        "DBK Xenium native vs Breast 1-trained",
        "Breast 1 native vs DBK-trained"
    ),
    spearman=c(cor_dbk_spearman, cor_bc1_spearman),
    pearson=c(cor_dbk_pearson, cor_bc1_pearson),
    stringsAsFactors=FALSE
)

xenium_dbk_breast1_transfer_summary

write.csv(
    xenium_dbk_breast1_transfer_summary,
    file=file.path(out_dir, "xenium_dbk_breast1_transfer_summary.csv"),
    row.names=FALSE
)

metadata(spexen)$QCScore_model$model_matrix_colnames
metadata(spexen_bc1)$QCScore_model$model_matrix_colnames

#########################
