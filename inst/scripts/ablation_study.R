## ============================================================
## 0. Libraries and output directory
## ============================================================

library(SpaceTrooper)
library(ggplot2)

out_dir <- "~/SpaceTrooper_qs_check"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive=TRUE)
}


## ============================================================
## 1. Define ablation formulas
## ============================================================
## The full model contains all Quality Score (QS) components.
## Each ablated model removes one component at a time.

formulas <- list(
    full="~(log2SignalDensity + Area_um + log2AspectRatio + log2Ctrl_total_ratio)^2",

    no_log2SignalDensity=
        "~(Area_um + log2AspectRatio + log2Ctrl_total_ratio)^2",

    no_Area_um=
        "~(log2SignalDensity + log2AspectRatio + log2Ctrl_total_ratio)^2",

    no_log2AspectRatio=
        "~(log2SignalDensity + Area_um + log2Ctrl_total_ratio)^2",

    no_log2Ctrl_total_ratio=
        "~(log2SignalDensity + Area_um + log2AspectRatio)^2"
)


## ============================================================
## 2. Function to run QS ablation models
## ============================================================
## For each formula, computeQCScore() is run with a forced model formula.
## The output stores:
## - the computed QS vector
## - the trained QS model stored in metadata(spe)$QCScore_model

run_qs_ablation <- function(spe, formulas, verbose=FALSE) {
    out <- list()

    for (nm in names(formulas)) {
        message("Running model: ", nm)

        spe_i <- computeQCScore(
            spe=spe,
            modelFormula=formulas[[nm]],
            verbose=verbose
        )

        out[[nm]] <- list(
            spe=spe_i,
            score=spe_i$QC_score,
            model=metadata(spe_i)$QCScore_model
        )
    }

    return(out)
}


## ============================================================
## 3. Read and preprocess CosMx data
## ============================================================
## The object is read, then per-cell QC metrics are computed.
## These metrics are required by computeQCScore().

cosm_path <- "~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/"

specosm <- readCosmxSPE(cosm_path)

specosm <- spatialPerCellQC(specosm)


## ============================================================
## 4. Run ablation study on CosMx
## ============================================================
## This computes the full QS and one QS for each ablated formula.

set.seed(1998)

abl_cosmx <- run_qs_ablation(
    spe=specosm,
    formulas=formulas,
    verbose=TRUE
)


## ============================================================
## 5. Functions to summarize ablation results
## ============================================================
## bottom_overlap() compares the low-QS tails of two score vectors.
## By default q=0.10 means the bottom 10% of cells.
##
## summarize_qs_ablation() compares every ablated QS against the full QS
## using:
## - Pearson correlation
## - Spearman correlation
## - absolute score differences
## - overlap of bottom-tail low-quality cells

bottom_overlap <- function(x, y, q=0.10) {
    x_bad <- x <= quantile(x, probs=q, na.rm=TRUE)
    y_bad <- y <= quantile(y, probs=q, na.rm=TRUE)

    overlap <- sum(x_bad & y_bad, na.rm=TRUE)
    union <- sum(x_bad | y_bad, na.rm=TRUE)

    data.frame(
        quantile=q,
        overlap=overlap,
        union=union,
        jaccard=overlap / union,
        prop_ref_recovered=overlap / sum(x_bad, na.rm=TRUE)
    )
}


summarize_qs_ablation <- function(ablation_results, reference="full", q=0.10) {
    ref_score <- ablation_results[[reference]]$score

    res <- lapply(names(ablation_results), function(nm) {
        score <- ablation_results[[nm]]$score

        ov <- bottom_overlap(ref_score, score, q=q)

        data.frame(
            model=nm,
            pearson=cor(
                ref_score,
                score,
                use="complete.obs",
                method="pearson"
            ),
            spearman=cor(
                ref_score,
                score,
                use="complete.obs",
                method="spearman"
            ),
            mean_abs_diff=mean(abs(ref_score - score), na.rm=TRUE),
            median_abs_diff=median(abs(ref_score - score), na.rm=TRUE),
            max_abs_diff=max(abs(ref_score - score), na.rm=TRUE),
            bottom_jaccard=ov$jaccard,
            bottom_ref_recovered=ov$prop_ref_recovered,
            n_na=sum(is.na(score)),
            stringsAsFactors=FALSE
        )
    })

    do.call(rbind, res)
}


## ============================================================
## 6. Build ablation summary table
## ============================================================

tab_cosmx <- summarize_qs_ablation(
    ablation_results=abl_cosmx,
    reference="full",
    q=0.10
)

tab_cosmx


## ============================================================
## 7. Pairwise plots: full QS vs each ablated QS
## ============================================================
## This assumes plot_qc_score_comparison() is already available.
## Each plot is also saved as an A4 landscape PDF.

ablation_names <- setdiff(names(abl_cosmx), "full")

plots_cosmx <- lapply(ablation_names, function(nm) {
    plot_qc_score_comparison(
        x=abl_cosmx$full$score,
        y=abl_cosmx[[nm]]$score,
        x_label="Full QS",
        y_label=paste0("QS: ", nm),
        title=paste0("CosMx: full QS vs ", nm),
        threshold=0.75,
        cor_method="spearman",
        output_file=file.path(
            out_dir,
            paste0("cosmx_full_vs_", nm, "_A4_landscape.pdf")
        )
    )
})

names(plots_cosmx) <- ablation_names

## Print plots
lapply(plots_cosmx, function(x) print(x$plot))


## ============================================================
## 8. Heatmap: complete ablation summary
## ============================================================
## This heatmap includes correlations, score differences, and bottom-tail
## overlap metrics.

heat_df <- tab_cosmx[tab_cosmx$model != "full", ]

heat_df <- heat_df[, c(
    "model",
    "pearson",
    "spearman",
    "mean_abs_diff",
    "median_abs_diff",
    "bottom_jaccard",
    "bottom_ref_recovered"
)]

heat_long <- reshape(
    heat_df,
    varying=names(heat_df)[names(heat_df) != "model"],
    v.names="value",
    timevar="metric",
    times=names(heat_df)[names(heat_df) != "model"],
    direction="long"
)

rownames(heat_long) <- NULL

p_heat <- ggplot(
    heat_long,
    aes(x=metric, y=model, fill=value)
) +
    geom_tile() +
    geom_text(
        aes(label=formatC(value, digits=3, format="f")),
        color="white",
        size=3
    ) +
    labs(
        title="CosMx QS ablation summary",
        x="Metric",
        y="Ablation model",
        fill="Value"
    ) +
    theme_minimal() +
    theme(
        axis.text.x=element_text(angle=45, hjust=1)
    )

p_heat

ggsave(
    filename=file.path(
        out_dir,
        "cosmx_ablation_summary_heatmap_A4_landscape.pdf"
    ),
    plot=p_heat,
    device="pdf",
    width=11.69,
    height=8.27,
    units="in"
)


## ============================================================
## 9. Heatmap: correlation-only summary
## ============================================================
## Cleaner figure showing only Pearson and Spearman correlations
## between the full QS and each ablated QS.

cor_heat_df <- tab_cosmx[
    tab_cosmx$model != "full",
    c("model", "pearson", "spearman")
]

cor_heat_long <- reshape(
    cor_heat_df,
    varying=c("pearson", "spearman"),
    v.names="correlation",
    timevar="metric",
    times=c("pearson", "spearman"),
    direction="long"
)

rownames(cor_heat_long) <- NULL

p_cor_heat <- ggplot(
    cor_heat_long,
    aes(x=metric, y=model, fill=correlation)
) +
    geom_tile() +
    geom_text(
        aes(label=formatC(correlation, digits=3, format="f")),
        color="white",
        size=4
    ) +
    labs(
        title="CosMx QS ablation correlation with full model",
        x="Correlation metric",
        y="Ablation model",
        fill="Correlation"
    ) +
    theme_minimal()

p_cor_heat

ggsave(
    filename=file.path(
        out_dir,
        "cosmx_ablation_correlation_heatmap.pdf"
    ),
    plot=p_cor_heat,
    device="pdf",
    width=8,
    height=5,
    units="in"
)


## ============================================================
## 10. Optional: save tabular results
## ============================================================

write.csv(
    tab_cosmx,
    file=file.path(out_dir, "cosmx_ablation_summary_table.csv"),
    row.names=FALSE
)


plot_ablation_colored <- function(ablation_results, ablation_name,
    reference="full", color_col=NULL,
    x_label="Full QS", y_label=NULL, title=NULL,
    cor_method="spearman", output_file=NULL) {

    stopifnot(reference %in% names(ablation_results))
    stopifnot(ablation_name %in% names(ablation_results))

    ref_score <- ablation_results[[reference]]$score
    abl_score <- ablation_results[[ablation_name]]$score

    stopifnot(length(ref_score) == length(abl_score))

    spe_ref <- ablation_results[[reference]]$spe
    cd <- as.data.frame(colData(spe_ref))

    df <- data.frame(
        full_qs=ref_score,
        ablated_qs=abl_score
    )

    if (!is.null(color_col)) {
        stopifnot(color_col %in% colnames(cd))
        df$color_var <- cd[[color_col]]
    }

    cor_val <- cor(
        df$full_qs,
        df$ablated_qs,
        use="complete.obs",
        method=cor_method
    )

    if (is.null(y_label)) {
        y_label <- paste0("QS: ", ablation_name)
    }

    if (is.null(title)) {
        title <- paste0("Full QS vs ", ablation_name)
    }

    if (is.null(color_col)) {
        p <- ggplot(df, aes(x=full_qs, y=ablated_qs)) +
            geom_point(alpha=0.35, size=0.6)
    } else {
        p <- ggplot(df, aes(x=full_qs, y=ablated_qs, color=color_var)) +
            geom_point(alpha=0.35, size=0.6)
    }

    p <- p +
        geom_abline(
            slope=1,
            intercept=0,
            linetype="dashed",
            color="red"
        ) +
        labs(
            x=x_label,
            y=y_label,
            title=title,
            color=color_col
        ) +
        theme_minimal() +
        annotate(
            "text",
            x=quantile(df$full_qs, 0.99, na.rm=TRUE),
            y=quantile(df$ablated_qs, 0.01, na.rm=TRUE),
            label=paste0(
                cor_method,
                " r = ",
                formatC(cor_val, digits=3, format="f")
            ),
            hjust=1,
            vjust=0,
            color="black"
        )

    if (!is.null(output_file)) {
        ggsave(
            filename=output_file,
            plot=p,
            device="pdf",
            width=11.69,
            height=8.27,
            units="in"
        )
    }

    return(list(
        plot=p,
        correlation=cor_val,
        data=df
    ))
}

p_no_signal_col <- plot_ablation_colored(
    ablation_results=abl_cosmx,
    ablation_name="no_log2SignalDensity",
    color_col="log2Ctrl_total_ratio",
    title="CosMx: full QS vs QS without log2SignalDensity",
    y_label="QS without log2SignalDensity",
    output_file="~/SpaceTrooper_qs_check/cosmx_full_vs_no_log2SignalDensity_log2Ctrl_total_ratiocolored.pdf"
)
p_no_signal_col_area <- plot_ablation_colored(
    ablation_results=abl_cosmx,
    ablation_name="no_log2SignalDensity",
    color_col="Area_um",
    title="CosMx: full QS vs QS without log2SignalDensity",
    y_label="QS without log2SignalDensity",
    output_file="~/SpaceTrooper_qs_check/cosmx_full_vs_no_log2SignalDensity_Area_um_colored.pdf"
)

p_no_signal_col_aspect <- plot_ablation_colored(
    ablation_results=abl_cosmx,
    ablation_name="no_log2SignalDensity",
    color_col="log2AspectRatio",
    title="CosMx: full QS vs QS without log2SignalDensity",
    y_label="QS without log2SignalDensity",
    output_file="~/SpaceTrooper_qs_check/cosmx_full_vs_no_log2SignalDensity_log2AspectRatio_colored.pdf"
)

p_no_signal_col_aspect$plot
p_no_signal_col$plot

plot_ablation_colored(
    ablation_results=abl_cosmx,
    ablation_name="no_log2Ctrl_total_ratio",
    color_col="log2SignalDensity",
    title="CosMx: full QS vs QS without log2Ctrl_total_ratio",
    y_label="QS without log2Ctrl_total_ratio",
    output_file="~/SpaceTrooper_qs_check/cosmx_full_vs_no_log2ctrl_total_ratio_log2SignalDensitycolored.pdf"
)$plot

plot_ablation_colored(
    ablation_results=abl_cosmx,
    ablation_name="no_log2Ctrl_total_ratio",
    color_col="Area_um",
    title="CosMx: full QS vs QS without log2Ctrl_total_ratio",
    y_label="QS without log2Ctrl_total_ratio",
    output_file="~/SpaceTrooper_qs_check/cosmx_full_vs_no_log2Ctrl_total_ratio_Area_um_colored.pdf"
)$plot
plot_ablation_colored(
    ablation_results=abl_cosmx,
    ablation_name="no_log2Ctrl_total_ratio",
    color_col="log2AspectRatio",
    title="CosMx: full QS vs QS without log2Ctrl_total_ratio",
    y_label="QS without log2Ctrl_total_ratio",
    output_file="~/SpaceTrooper_qs_check/cosmx_full_vs_no_log2Ctrl_total_ratio_log2AspectRatio_colored.pdf"
)$plot
