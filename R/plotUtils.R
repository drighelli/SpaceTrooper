#' .fov_image_theme
#' @name .fov_image_theme
#' @rdname dot-fov_image_theme
#' @description
#' internal function to setup the theme for the fov background on the whole
#' image
#'
#' @param backColor not used
#' @param backBorder color for the borders of the background (default=NA)
#' @param titleCol character indicating the color of the title
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect element_text
#' @keywords internal
.fov_image_theme <- function(backColor="black", backBorder=NA,
                            titleCol="white")
{
    theme(panel.border=element_blank(),
        legend.key=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.title=element_text(color=titleCol, hjust=0.5, face="bold"),
        plot.background=element_rect(fill="transparent", colour=backBorder))
}

#' .centroid_image_theme
#' @name .centroid_image_theme
#' @rdname dot-centroid_image_theme
#' @description
#' internal function to setup the theme for the centroid plot background
#'
#' @param backBorder color for the borders of the background (default=NA)
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect
#' @keywords internal
.centroid_image_theme <- function(backBorder=NA) {
    theme(panel.border=element_rect(color="black"),
        legend.key=element_blank(),
        axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.title=element_blank(),
        plot.background=element_rect(fill="transparent", colour=backBorder))
}


#' .negative_image_theme
#' @name .negative_image_theme
#' @rdname dot-negative_image_theme
#' @description
#' internal function to setup the theme for the negative background for
#' negative plots
#'
#' @param fillColor color to fill the element_rect (default is "black")
#' @param foreColor color for all the other elements (default is "white")
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_line element_rect element_text
#' element_blank
#' @keywords internal
.negative_image_theme <- function(fillColor="black", foreColor="white") {
    theme(panel.border = element_rect(color = foreColor),
        panel.background=element_rect(fill=fillColor, color=NA),
        plot.background=element_rect(fill=fillColor, color=NA),
        title=element_blank(),
        axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        legend.background=element_rect(fill=fillColor, color=NA),
        legend.text=element_text(color=foreColor),
        legend.title=element_text(color=foreColor),
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
}

#' .light_theme
#' @name .light_theme
#' @rdname dot-light_theme
#' @description
#' internal function to setup the white background theme for the First Filter
#' plot
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect
#' @keywords internal
.light_theme <- function(fillColor="white", foreColor="black") {
    theme(panel.border=element_rect(color = foreColor, fill = NA,
        linewidth = 0.1),
        panel.background=element_rect(fill=fillColor, color=NA),
        plot.background=element_rect(fill=fillColor, color=NA),
        axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
}

#' .dark_theme
#' @name .dark_theme
#' @rdname dot-dark_theme
#' @description
#' internal function to setup the black background theme for the First Filter
#' plot
#' @param fillColor color to fill the element_rect (default is "black")
#' @param foreColor color for all the other elements (default is "white")
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect element_text
#' @keywords internal
.dark_theme <- function(fillColor="black", foreColor="white") {
    theme(panel.background=element_rect(fill=fillColor, color=NA),
        plot.background=element_rect(fill=fillColor, color=NA),
        panel.border=element_rect(color=foreColor, fill=NA, linewidth=0.1),
        title=element_text(color=foreColor),
        axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        legend.background=element_rect(fill=fillColor, color=NA),
        legend.text=element_text(color=foreColor),
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
}

#' createPaletteFromColData
#' @name createPaletteFromColData
#' @rdname createPaletteFromColData
#' @description
#' Create a Palette from colData in a SpatialExperiment Object
#'
#' This function generates a palette mapping based on specified columns in
#' the `colData` of a `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param paletteNames A character string specifying the column in
#' `colData(spe)` to be used for the names in the palette.
#' @param paletteColors A character string specifying the column in
#' `colData(spe)` to be used for the colors in the palette.
#'
#' @return A character vector representing the palette mapping, where each
#' element is a string in the format `"name=color"`.
#'
#' @details The function creates a new palette based on the unique combinations
#' of values in the specified `paletteNames` and `paletteColors` columns in
#' `colData(spe)`.
#'
#' @importFrom SummarizedExperiment colData
#' @keywords internal
createPaletteFromColData <- function(spe, paletteNames, paletteColors)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(all(c(paletteNames, paletteColors) %in% names(colData(spe))))

    tb <- table(spe[[paletteNames]], spe[[paletteColors]])
    newpal <- colnames(tb)[which(tb!=0, arr.ind=TRUE)[,2]]
    names(newpal) <- rownames(tb)[which(tb!=0, arr.ind=TRUE)[,1]]
    return(newpal)
}

#' firstFlagPalette
#' @name firstFlagPalette
#' @rdname firstFlagPalette
#' @description
#' neon color palette for firstFlagPlot
#'
#' @return a palette for firstFlagPlot
#' @keywords internal
firstFlagPalette <- c(
    "unflagged" = "#7f7f7f",
    "ctrl/total ratio > 0.1"     = "magenta",
    "< area um lower thr."       = "darkturquoise",
    "> area um higher thr."      = "red",
    "< logged aspect ratio lower thr."  = "purple",
    "> logged aspect ratio higher thr." = "greenyellow"
)

# To be decided if this should be moved to a separate file or included in the main package
# plot_qc_score_comparison <- function(x, y, x_label="x", y_label="y",
#     title="QC score comparison", threshold=0.75,
#     cor_method="spearman", output_file=NULL) {

#     stopifnot(length(x) == length(y))

#     df <- data.frame(x=x, y=y)

#     df$quadrant <- with(df, ifelse(
#         x < threshold & y < threshold, "low / low",
#         ifelse(x >= threshold & y < threshold, "high x / low y",
#         ifelse(x < threshold & y >= threshold, "low x / high y",
#             "high / high"))
#     ))

#     cor_val <- cor(df$x, df$y, use="complete.obs", method=cor_method)

#     p <- ggplot(df, aes(x=x, y=y, color=quadrant)) +
#         geom_point(alpha=0.35, size=0.6) +
#         geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
#         geom_hline(yintercept=threshold) +
#         geom_vline(xintercept=threshold) +
#         scale_color_manual(values=c(
#             "low / low"="grey60",
#             "high x / low y"="orange",
#             "low x / high y"="dodgerblue",
#             "high / high"="black"
#         )) +
#         labs(x=x_label, y=y_label, title=title, color="Quadrant") +
#         theme_minimal() +
#         annotate(
#             "text",
#             x=quantile(df$x, 0.99, na.rm=TRUE),
#             y=quantile(df$y, 0.01, na.rm=TRUE),
#             label=paste0(cor_method, " r = ",
#                 formatC(cor_val, digits=3, format="f")),
#             hjust=1,
#             vjust=0,
#             color="black"
#         )

#     if (!is.null(output_file)) {
#         ggsave(output_file, p, device="pdf",
#             width=11.69, height=8.27, units="in")
#     }

#     return(list(
#         plot=p,
#         correlation=cor_val,
#         quadrant_table=table(df$quadrant)
#     ))
# }

