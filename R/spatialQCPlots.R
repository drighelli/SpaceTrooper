#' plotCellsFovs
#' @name plotCellsFovs
#' @rdname plotCellsFovs
#'
#' @description
#' Plot cell centroids in FoVs, creating a map of the whole experiment, where
#' cells are plotted as points and FoV boundaries and numbers are overlaid.
#'
#' @param spe A `SpatialExperiment` object with `fov` in `colData`.
#' @param sampleId Character string identifying which sample to plot.
#'   Default: `unique(spe$sample_id)`.
#' @param pointCol Color for the cell centroids. Default: `"firebrick"`.
#' @param pointSize Numeric point size for the cell centroids. Default: `0.05`.
#' @param pointAlpha Numeric transparency for the cell centroids. Default: `0.8`.
#' @param numbersCol Color for the FoV labels. Default: `"black"`.
#' @param numberSize Numeric size for the FoV labels. Default: `1`.
#' @param numbersAlpha Numeric transparency for FoV labels. Default: `0.8`.
#' @param fovDim numeric with two named dimensions xdim, ydim. (Default is
#' metadata(spe)$fov_dim)
#' @param scaleBar A logical value indicating whether to add a scale bar to the
#' plot. (Default is `TRUE`)
#' @param micronConvFact Numeric conversion factor from pixels to microns.
#' Default is `0.12`.
#'
#' @return A `ggplot` object showing cell centroids and FoV boundaries.
#'
#' @details
#' The function expects spe (a SpatialExperiment) to include field-of-view
#' metadata in metadata(spe):
#' - `metadata(spe)$fov_positions`: a matrix or data.frame
#' (or list with named elements) containing at minimum `x_global_px`, `y_global_px`,
#' and `fov`. Values `x_global_px`/`y_global_px` are in pixels and represent
#' the origin (bottom-left) of each FoV.
#' - `metadata(spe)$fov_dim` (or the `fovDim` argument): a named numeric with
#' `xdim` and `ydim` giving FoV width/height in pixels.
#'
#' For each FoV the rectangle is drawn as:
#' - `xmin` = `x_global_px`,
#' - `xmax` = `x_global_px` + `xdim`,
#' - `ymin` = `y_global_px`,
#' - `ymax` = `y_global_px` + `ydim`.
#'
#' `sampleId`: if `NULL` no plot title is added; if a vector of length > 1 is
#' provided the first element is used. Prefer passing a single sample identifier
#' string.
#' `micronConvFact` is the conversion factor from pixels to micrometers (µm/pixel).
#' It is forwarded to `.plotScaleBar()` and controls how the scale bar is labeled
#' (useful for technologies where coordinates are in pixels but the scale should
#' display µm).
#' The function uses the coordinates returned by `spatialCoords(spe)` and the
#' names from `spatialCoordsNames(spe)`. It performs basic input checks but
#' assumes the described metadata structures are present; if they are missing
#' or malformed the function will raise an informative error.
#'
#' @importFrom ggplot2 ggplot aes annotate geom_point geom_text ggtitle
#' coord_fixed
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' example(readCosmxSPE)
#' g <- plotCellsFovs(spe)
#' print(g)
plotCellsFovs <- function(spe, sampleId=unique(spe$sample_id),
                        pointCol="firebrick", pointSize=0.05,
                        pointAlpha=0.8, numbersCol="black",
                        numberSize= 1, numbersAlpha=0.8,
                        fovDim=metadata(spe)$fov_dim,
                        scaleBar=TRUE, micronConvFact = 0.12)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("fov" %in% names(colData(spe)))
    stopifnot( all(names(fovDim) %in% c("xdim","ydim")) )

    spd <- as.data.frame(spatialCoords(spe))
    x_coord <- spatialCoordsNames(spe)[1]
    y_coord <- spatialCoordsNames(spe)[2]
    ggp <- ggplot() +
        geom_point(data=spd, mapping=aes(x=.data[[x_coord]],
                                        y=.data[[y_coord]]),
                    colour=pointCol,
                    fill=pointCol,
                    size=pointSize, alpha=pointAlpha) +
        annotate("rect",
            xmin=metadata(spe)$fov_positions["x_global_px"][ , , drop=TRUE],
            xmax=metadata(spe)$fov_positions["x_global_px"][ , , drop=TRUE] +
                fovDim[["xdim"]],
            ymin=metadata(spe)$fov_positions["y_global_px"][ , , drop=TRUE],
            ymax=metadata(spe)$fov_positions["y_global_px"][ , , drop=TRUE] +
                fovDim[["ydim"]],
            alpha=.2, color="black", linewidth=0.2) +
        geom_text(aes(x=metadata(spe)$fov_positions["x_global_px"][,,drop=TRUE]+
                        fovDim[["xdim"]]/2,
                    y=metadata(spe)$fov_positions["y_global_px"][,,drop=TRUE]+
                        fovDim[["ydim"]]/2,
                    label=metadata(spe)$fov_positions["fov"][,,drop=TRUE]),
                    color=numbersCol, size = numberSize, alpha=numbersAlpha,
                    fontface="bold") +
        ggtitle(sampleId) +
        .fov_image_theme(backColor="white", backBorder="white",
                        titleCol="black") + ggplot2::coord_fixed()
    if(scaleBar) {
        ggp <- .plotScaleBar(ggp, spe, micronConvFact)
    }
    return(ggp)
}

#' plotCentroids
#' @name plotCentroids
#' @rdname plotCentroids
#' @description
#' Plot Spatial Coordinates for a SpatialExperiment Object
#' This function generates a ggplot of spatial coordinates from a
#' `SpatialExperiment` object, optionally coloring the points by a specified
#' column in `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial
#' transcriptomics data.
#' @param colourBy An optional character string specifying the column in
#' `colData(spe)` to use for coloring the points. If `NULL`, all points will be
#' colored the same.
#' @param colourLog Logical to log-transform the data to enhance visualization
#' (Default is FALSE).
#' @param sampleId A character string specifying the sample identifier to be
#' used as the plot title. (Default is the unique sample ID from `spe`)
#' @param isNegativeProbe A logical value indicating whether to apply a custom
#' color gradient for negative probe data. (Default is `FALSE`)
#' @param palette A vector of colors to be used as a custom palette. For
#' categorical data, this should be a vector of colors with the same length as
#' the number of levels in `colourBy`. For continuous data, this should be a
#' vector of colors used to create a gradient.
#' @param pointCol A character string specifying the color of the points when
#' `colourBy` is `NULL`. (Default is `"darkmagenta"`)
#' @param size A numeric value specifying the size of the points. (Default is
#' `0.05`)
#' @param alpha A numeric value specifying the transparency level of the points.
#' (Default is `0.8`)
#' @param aspectRatio A numeric value specifying the aspect ratio of the plot.
#' (Default is `1`)
#' @param scaleBar A logical value indicating whether to add a scale bar to the
#' plot. (Default is `TRUE`)
#' @param micronConvFact Numeric conversion factor from pixels to microns.
#' DEFAULT is `0.12`.
#'
#' @return A `ggplot` object representing the spatial coordinates plot of
#' polygon centroids.
#'
#' @details
#' This function plots cell centroids from a `SpatialExperiment` object using
#' coordinates returned by `spatialCoords(spe)` and the coordinate names from
#' `spatialCoordsNames(spe)` (the function expects at least two spatial
#' coordinate dimensions). If `colourBy` is `NULL` all points are drawn using
#' `pointCol`; otherwise the function colors points using the specified
#' `colData` column or a provided palette.
#'
#' Requirements and inputs:
#' - `spe` must be a `SpatialExperiment` with valid spatial coordinates.
#' - If `colourBy` is provided and is the name of a `colData` column, that
#'   column must exist and be suitable for plotting (factor/logical for
#'   discrete palettes, numeric for continuous palettes).
#' - `palette` may be either (a) a vector of color values (used directly) or
#'   (b) the name of a column in `colData(spe)` from which a palette will be
#'   built via `createPaletteFromColData()`. The function detects which mode
#'   to use at runtime.
#'
#' @import SpatialExperiment
#' @importFrom ggplot2 geom_point aes_string theme_bw theme ggtitle guides
#' guide_legend scale_color_gradient scale_color_manual scale_color_gradientn
#' @importFrom scater plotColData
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot geom_point aes ggtitle theme_bw coord_fixed
#' @importFrom ggplot2 scale_color_manual scale_color_gradientn
#' @export
#'
#' @examples
#' example(readCosmxSPE)
#' g <- plotCentroids(spe, colourBy="Mean.DAPI")
#' print(g)
plotCentroids <- function(spe, colourBy=NULL, colourLog=FALSE,
                        sampleId=unique(spe$sample_id),
                        isNegativeProbe=FALSE, palette=NULL,
                        pointCol="darkmagenta", size=0.05, alpha=0.8,
                        aspectRatio=1, scaleBar=TRUE, micronConvFact=0.12) {
    stopifnot(is(spe, "SpatialExperiment"))
    if(is.null(colourBy)) {
        ggp <- ggplot(data.frame(spatialCoords(spe)),
                    aes(x=.data[[spatialCoordsNames(spe)[1]]],
                        y=.data[[spatialCoordsNames(spe)[2]]])) +
            geom_point(colour=pointCol, fill=pointCol,
                        size=size, alpha=alpha) +
            ggplot2::theme_bw() + ggplot2::coord_fixed()
    } else {
        if(colourLog) {
            stopifnot(colourBy %in% names(colData(spe)))
            colourByo <- colourBy
            colourBy <- paste0("log(", colourByo, ")")
            colData(spe)[[colourBy]] <- log1p(colData(spe)[[colourByo]])
        }
        ggp <- ggplot(data.frame(colData(spe), spatialCoords(spe)),
                    aes(x=.data[[spatialCoordsNames(spe)[1]]],
                        y=.data[[spatialCoordsNames(spe)[2]]],
                        colour=.data[[colourBy]],
                        fill=.data[[colourBy]])) +
            geom_point(size=size, alpha=alpha) + coord_fixed()
        if( all(!is(spe[[colourBy]], "factor"),
                !is(spe[[colourBy]], "logical"))) {
            ggp <- ggp + ggplot2::scale_colour_viridis_c() +
                ggplot2::scale_fill_viridis_c() + coord_fixed()
        }
        if(isNegativeProbe) {
            ggp <- ggp + scale_color_gradient(low="white", high="red",
                name=colourBy) + .negative_image_theme()
        } else if(all(!is.null(palette), (palette %in% names(colData(spe))))) {
            palette <- createPaletteFromColData(spe, paletteNames=colourBy,
                                                    paletteColors=palette)
            if(is.factor(colData(spe)[[colourBy]])) {
                ggp <- ggp + scale_color_manual(values=palette)
            } else if(is.numeric(colData(spe)[[colourBy]])) {
                ggp <- ggp + scale_color_gradientn(colors=palette)
            }
        }
    }
    ggp <- ggp + ggtitle(sampleId) +
        theme(aspect.ratio=aspectRatio, plot.title=element_text(hjust=0.5))
    if(!isNegativeProbe) ggp <- ggp + theme_bw()
    if(scaleBar) {
        ggp <- .plotScaleBar(ggp, spe, micronConvFact)
    }
    return(ggp)
}


#' plotMetricHist
#' @name plotMetricHist
#' @rdname plotMetricHist
#' @description Plot a Histogram for a Given Metric in a SpatialExperiment
#' Object
#'
#' This function generates a histogram for a specified metric in a
#' `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param metric A character string specifying the name of the metric
#' (column in `colData(spe)`) to plot.
#' @param fillColor A character string specifying the fill color of the
#' histogram bars. (Default is `"#69b3a2"`)
#' @param useFences A character string specifying the name of the column in
#' `colData(spe)` that contains the fence thresholds (typically from an outlier
#' filter). If `NULL`, no fences will be plotted. (Default is `NULL`)
#' @param fencesColors A named character vector specifying the colors to use
#' for the lower and higher fences. The names should be `"lower"` and `"higher"`
#'. (Default is `c("lower"="purple4", "higher"="tomato")`)
#' @param bins An integer specifying the number of bins to use in the histogram.
#' (Default is `30`)
#' @param binWidth A numeric value specifying the width of the bins. If `NULL`,
#' the bin width will be automatically determined based on the `bins` parameter.
#' (Default is `NULL`)
#'
#' @return A `ggplot` object representing the histogram of the specified metric.
#'
#' @importFrom ggplot2 ggplot geom_histogram aes ggtitle theme_bw geom_vline
#' labs scale_colour_manual
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#' @export
#' @examples
#' example(readCosmxSPE)
#' g <- plotMetricHist(spe, metric="Mean.DAPI")
#' print(g)
plotMetricHist <- function(spe, metric, fillColor="#c0c8cf",
        useFences=NULL, fencesColors=c("lower"="purple4", "higher"="tomato"),
        bins=30, binWidth=NULL)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(metric %in% names(colData(spe)))

    ggp <- ggplot(data=as.data.frame(colData(spe))) +
            geom_histogram(aes(x=.data[[metric]]), fill=fillColor,
                    bins=bins, binwidth=binWidth)
    if (!is.null(useFences))
    {
        stopifnot(useFences %in% names(colData(spe)))
        fences <- getFencesOutlier(spe, useFences, "both", 2)
        fences_labs <- paste0(names(fences), ": ", fences)
        names(fencesColors) <- fences_labs
        ggp <- ggp +
            geom_vline(aes(xintercept=fences[1], color=fences_labs[1])) +
            geom_vline(aes(xintercept=fences[2], color=fences_labs[2])) +
            labs(color=paste0("Fences ", useFences)) +
            scale_colour_manual(values=fencesColors)
    }
    ggp <- ggp + ggtitle(metric) + theme_bw()

    return(ggp)
}


#' plotPolygons
#' @name plotPolygons
#' @rdname plotPolygons
#'
#' @description Plot polygons from a `SpatialExperiment` object using ggplot2.
#'
#' @param spe A `SpatialExperiment` object with polygon data as an `sf` object.
#' @param colourBy A column in `colData(spe)` for coloring the polygons or a
#' string color in colors(). (Default is "darkgrey")
#' @param colourLog Logical to log-transform the data to enhance visualization
#' (Default is FALSE).
#' @param sampleId Sample ID for plot title. Default is the unique sample ID.
#' @param fillAlpha Transparency level for polygon fill. Default is `1`.
#' @param palette Colors to use if `colourBy` is a factor. Default is `NULL`.
#' @param borderCol Color of polygon borders. Default is `"black"`.
#' @param borderAlpha Transparency level for borders. Default is `1`.
#' @param borderLineWidth Width of polygon borders. Default is `0.1`.
#' @param drawBorders Logical; whether to draw borders. Default is `TRUE`.
#' @param polyColumn character for the name of the column where the polygons sf
#' are stored (default is "polygons.global")
#' @param bgColor character indicating color for the background
#' (default is "white")
#' @param scaleBar A logical value indicating whether to add a scale bar to the
#' plot. (Default is `TRUE`)
#' @param micronConvFact Numeric conversion factor from pixels to microns.
#' DEFAULT is `0.12`.
#'
#' @return A `ggplot` object representing the polygon plot of the spatial data.
#'
#' @details
#' Renders polygon geometries stored in `colData` of a `SpatialExperiment`
#' object using `ggplot2::geom_sf()` and provides options for fill, borders
#' and an optional scale bar.
#'
#' Input expectations:
#' - `spe` must be a `SpatialExperiment` and `polyColumn` must name a column
#'   in `colData(spe)` that contains polygon geometries (`sfc`/`sf` or a list
#'   of geometry objects) suitable for `geom_sf()`.
#' - `colourBy` supports two modes: (a) the name of a `colData` column to
#'   map fill to data values, or (b) a literal colour name (e.g. `"darkgrey"`)
#'   used as a constant fill.
#'
#' Colour / palette behaviour:
#' - If `colourBy` refers to a `colData` column and `colourLog = TRUE`, the
#'   function computes `log1p()` of that column in the local plotting data
#'   frame (it does not mutate `spe`).
#' - `palette` may be a vector of colour values (used for discrete scales).
#'   Numeric `colourBy` values use a continuous Viridis-like scale by default.
#'
#' Scale bar and units:
#' - If `scaleBar = TRUE` the function calls `.plotScaleBar()` to add a length
#'   scale. `micronConvFact` is the conversion factor (µm per pixel) used by
#'   `.plotScaleBar()` to label the bar in micrometres when appropriate.
#'
#' Robustness notes:
#' - If `colourBy` is neither a `colData` column nor a valid colour name the
#'   function issues a warning and falls back to the default fill.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_sf aes scale_fill_manual scale_fill_viridis_c
#' scale_fill_identity theme_minimal theme element_text margin labs
#' @importFrom sf st_as_sf
#' @importFrom SummarizedExperiment colData
#' @importFrom grDevices colors
#'
#' @examples
#' example(readAndAddPolygonsToSPE)
#' plotPolygons(spe, colourBy="Mean.DAPI")
plotPolygons <- function(spe, colourBy="darkgrey", colourLog=FALSE,
    polyColumn="polygons.global", sampleId=unique(spe$sample_id),
    bgColor="white", fillAlpha=1, palette=NULL, borderCol=NA,
    borderAlpha=1, borderLineWidth=0.1, drawBorders=TRUE,
    scaleBar=TRUE, micronConvFact=0.12) {

    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("polygons" %in% names(colData(spe)))
    df <- data.frame(colData(spe))
    polflag <- FALSE
    if(!is.null(colourBy)) {
        if(colourBy %in% names(colData(spe))) {
            if(colourLog) {
                colourByo <- colourBy
                colourBy <- paste0("log(", colourByo, ")")
                df[[colourBy]] <- log1p(df[[colourByo]])
            }
            polflag <- TRUE
        } else {
            if(!(colourBy %in% grDevices::colors())) {
                warning(colourBy, " not in known colors nor in colData")
                colourBy <- "darkgrey"
            }}}
    border_params <- if(drawBorders) {
        list(color=borderCol, size=borderLineWidth)
    } else { list(color=NA, size=0) }
    if(polflag) {
        p <- ggplot(df, aes(geometry=.data[[polyColumn]],
            fill=.data[[colourBy]])) + geom_sf(alpha=fillAlpha,
            color=border_params$color, size=border_params$size)
    } else {
        p <- ggplot(df, aes(geometry = .data[[polyColumn]])) +
            geom_sf(fill=colourBy, alpha=fillAlpha,
                color=border_params$color, size=border_params$size)
    }
    if(!is.null(colourBy) && (is.factor(df[[colourBy]]) ||
                                is.logical(df[[colourBy]]))) {
        if(!is.null(palette)) { p <- p + scale_fill_manual(values=palette) }
    } else if(!is.null(colourBy) && is.numeric(df[[colourBy]])) {
        p <- p + scale_fill_viridis_c(option="D")
    } else { p <- p + scale_fill_identity() }
    p <- p + theme_minimal() + theme(legend.position="right",
        plot.title.position="plot",plot.title=element_text(face="bold",size=14),
            plot.margin=margin(0, 0, 0, 0),
            panel.background=element_rect(fill=bgColor, color=bgColor, linewidth=1),
            panel.border = element_rect(color="black", fill=NA, linewidth=0.1),
            panel.grid.minor=element_blank()
        ) + labs(title=sampleId, fill=colourBy)
    if(scaleBar) {
        p <- .plotScaleBar(p, spe, micronConvFact)
    }
    return(p)
}


#' plotZoomFovsMap
#' @name plotZoomFovsMap
#' @rdname plotZoomFovsMap
#' @description
#'
#' Plot Zoomed-in FOVs with Map and Polygons
#'
#' This function generates a plot that shows a map of all fields of view (FOVs)
#' within a `SpatialExperiment` object, alongside a zoomed-in view of the
#' specified FOVs with an overlay of polygons and optional coloring.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param fovs A character vector specifying the FOVs to be zoomed in and
#' plotted. Must match values in the `fov` column of `colData(spe)`.
#' @param title An optional character string specifying the title of the final
#' plot. If `NULL`, no title is added. Default is `NULL`.
#' @param mapPointCol A character string specifying the color of the points
#' in the map. Default is `"darkmagenta"`.
#' @param mapPointSize Numeric size for points in the map. Default: `0.5`.
#' @param mapPointAlpha Numeric transparency for points in the map. Default: `0.8`.
#' @param fovNumbersCol A character string specifying the color of the
#' numbers on the FoV zoom-in. Default is `"black"`.
#' @param fovNumberSize Numeric size for the FoV labels. Default: `1`.
#' @param fovNumbersAlpha Numeric transparency for FoV labels. Default: `0.8`.
#' @param scaleBars Logical or NULL. Default is `NULL`.
#' Master switch controlling the presence of scale bars in both panels.
#' If \code{TRUE}, scale bars are shown in both the map and polygon panels.
#' If \code{FALSE}, scale bars are hidden in both panels.
#' If \code{NULL}, individual settings defined by \code{scaleBarMap} and
#' \code{scaleBarPol} are used.
#' @param scaleBarMap,scaleBarPol Logical. Default is \code{TRUE}.
#' Control the presence of the scale bar in the map (cells/FOV overview)
#' and polygon (segmentation) panels, respectively.
#' These parameters are only used when \code{scaleBars} is \code{NULL};
#' otherwise they are overridden by \code{scaleBars}.
#' @param ... Additional arguments passed to `plotPolygons`.
#'
#' @return A combined plot showing a map of all FOVs with zoomed-in views of
#' the specified FOVs and their associated polygons.
#'
#' @details The function first filters the `SpatialExperiment` object to the
#' specified FOVs, generates a plot of the cells for the entire map, then
#' creates a detailed polygon plot of the selected FOVs, and finally combines
#' these into a single side-by-side visualization. If `title` is not `NULL`, it
#' adds a title to the combined plot.
#'
#' @importFrom ggpubr ggarrange annotate_figure
#' @export
#'
#' @examples
#' example(readAndAddPolygonsToSPE)
#' plotZoomFovsMap(spe, fovs=16, title="FOV 16")
plotZoomFovsMap <- function(spe, fovs=NULL, title=NULL,
                            mapPointCol="darkmagenta",
                            mapPointSize=0.5, mapPointAlpha=0.8,
                            fovNumbersCol="black", fovNumberSize=1,
                            fovNumbersAlpha=0.8,
                            scaleBars=NULL,
                            scaleBarMap=TRUE,
                            scaleBarPol=TRUE,
                            ...) {
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("fov" %in% names(colData(spe)))
    stopifnot(all(fovs %in% spe$fov))
    spefovs <- spe[, spe$fov %in% fovs]
    if (!is.null(scaleBars)) {
        stopifnot(is.logical(scaleBars), length(scaleBars) == 1L)
        scaleBarMap <- scaleBars
        scaleBarPol <- scaleBars
    }
    map <- plotCellsFovs(spefovs, pointCol=mapPointCol,
        pointSize=mapPointSize, pointAlpha=mapPointAlpha,
        numbersCol=fovNumbersCol, numberSize=fovNumberSize,
        numbersAlpha=fovNumbersAlpha, sampleId=NULL, scaleBar=scaleBarMap)
    g2 <- plotPolygons(spefovs, sampleId=NULL, scaleBar=scaleBarPol, ...)
    final_plot <- ggpubr::ggarrange(map, g2, ncol=2)
    if (!is.null(title)) {
        final_plot <- ggpubr::annotate_figure(final_plot,
            top=ggpubr::text_grob(title, face="bold", size=14))
    }
    return(final_plot)
}

#' plotQScoreTerms
#' @name plotQScoreTerms
#' @rdname plotQScoreTerms
#'
#' @description
#' Plots the individual terms that combine into the quality score formula,
#' allowing assessment of each term's impact on the final score.
#'
#' @param spe A `SpatialExperiment` object with `quality_score` and term
#'   columns in `colData`.
#' @param sampleId Character string for plot title. Must match values in the
#'   `fov` column of `colData(spe)`. Default: `unique(spe$sample_id)`.
#' @param size Numeric point size for the scatter plots. Default: `0.05`.
#' @param alpha Numeric transparency for the scatter plots. Default: `0.2`.
#' @param aspectRatio Numeric aspect ratio of the plots. Default: `1`.
#' @param custom Logical; if `TRUE`, use custom polygon-derived metrics.
#'
#' @return A combined plot (via `cowplot::plot_grid`) showing spatial maps
#'   of each QS term.
#'
#' @importFrom scater plotColData
#' @importFrom ggplot2 ggtitle coord_fixed
#' @importFrom cowplot plot_grid
#' @export
#'
#' @examples
#' example(readAndAddPolygonsToSPE)
#' example(spatialPerCellQC)
#' p <- plotQScoreTerms(spe)
#' print(p)
plotQScoreTerms <- function(spe, sampleId=unique(spe$sample_id), size=0.05,
                            alpha=0.8, aspectRatio=1, custom = FALSE) {
    if(metadata(spe)$technology=="Nanostring_CosMx") {
        if(custom==TRUE) {
            ggp <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                        y=spatialCoordsNames(spe)[2],
                                        colour_by="cust_log2SignalDensity",
                                        point_size=size, point_alpha=alpha)+
                ggplot2::ggtitle(sampleId)+ .centroid_image_theme() +
                ggplot2::coord_fixed()
            ggp2 <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                        y=spatialCoordsNames(spe)[2],
                                        colour_by="cust_log2AspectRatio",
                                        point_size=size, point_alpha=alpha)+
                ggplot2::ggtitle(sampleId) + .centroid_image_theme() +
                ggplot2::coord_fixed()
        } else {
            ggp <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                        y=spatialCoordsNames(spe)[2],
                                        colour_by="log2SignalDensity",
                                        point_size=size, point_alpha=alpha)+
                ggtitle(sampleId)+ .centroid_image_theme() + coord_fixed()
            ggp2 <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                        y=spatialCoordsNames(spe)[2],
                                        colour_by="log2AspectRatio",
                                        point_size=size, point_alpha=alpha)+
                ggplot2::ggtitle(sampleId) + .centroid_image_theme() +
                ggplot2::coord_fixed()
        }
        ggp3 <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                    y=spatialCoordsNames(spe)[2],
                                    colour_by="dist_border",
                                    point_size=size, point_alpha=alpha)+
            ggplot2::ggtitle(sampleId) + .centroid_image_theme() +
            ggplot2::coord_fixed()
        ggp <- cowplot::plot_grid(ggp, ggp2, ggp3, ncol = 2)
    } else {
        ## check if column variable is logical to impose our colors
        ggp <- ggp1 <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                            y=spatialCoordsNames(spe)[2],
                                            colour_by="log2SignalDensity",
                                            point_size=size, point_alpha=alpha)+
            ggplot2::ggtitle(sampleId)+ .centroid_image_theme() +
            ggplot2::coord_fixed()
    }
    return(ggp)
}



.assign_outlier_color <- function(vals, outlier_mc, high_label, low_label) {
    thr <- round(attr(outlier_mc, "thresholds"), 2)
    dplyr::case_when(
        vals > thr[2] ~ high_label,
        vals < thr[1] ~ low_label,
        TRUE ~ "unflagged"
    )
}

.assign_collapsed_color <- function(is_ctrl, area_col, dapi_col) {
    color <- rep("unflagged", length(is_ctrl))
    color[is_ctrl] <- "ctrl/total ratio > 0.1"
    color[area_col == "> area um higher thr."] <- "> area um higher thr."
    color[area_col == "< area um lower thr."] <- "< area um lower thr."
    color[dapi_col == "> DAPI higher thr."] <- "> DAPI higher thr."
    color[dapi_col == "< DAPI lower thr."] <- "< DAPI lower thr."
    color
}

.make_outlier_plot <- function(polygons, fov, fillvar, pal, title = NULL,
                                leg = FALSE) {
    filtered_polygons <- polygons[polygons$fov %in% fov, ]

    if (nrow(filtered_polygons) == 0 || !fillvar %in% names(filtered_polygons)){
        # Create an empty plot with appropriate theme
        p <- ggplot2::ggplot() +
            ggplot2::theme_void()
        if (!is.null(title)) {
            p <- p + ggplot2::ggtitle(title)
        }
        return(p)
    }
    newpal <- pal[which(names(pal)%in%unique(filtered_polygons[[fillvar]]))]
    p <- ggplot2::ggplot() +
        ggplot2::geom_sf(
            data = filtered_polygons,
            mapping = ggplot2::aes(
                fill  = .data[[fillvar]],
                color = .data[[fillvar]]
            ),
            lwd = 0,
            show.legend = TRUE
        ) +
        ggplot2::scale_fill_manual(values = newpal) +
        ggplot2::scale_color_manual(values = newpal)

    if (!leg) {
        p <- p + ggplot2::theme(legend.position = "none")
    }
    if (!is.null(title)) {
        p <- p + ggplot2::ggtitle(title)
    }
    return(p)
}



#' qcFlagPlots
#' @name qcFlagPlots
#' @rdname qcFlagPlots
#' @description
#' Plots the flagged cells identified with first filter, based on control count
#' on total count ratio, area in um and DAPI signal.
#'
#' This function generates a plot that shows selected (FOVs)
#' within a `SpatialExperiment` object, with cells flagged in different colors
#' over a light or dark layout chosen by the user.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param fov An integer or numeric vector specifying the FOVs to be plotted
#' Must match values in the `fov` column of `colData(spe)`.
#' @param theme A character string among "light" or "dark".
#' @param custom A boolean value. If TRUE, custom polygons derived metrics will
#' be used.
#'
#' @return A panel with multiple plots showing flagged cells for different
#' variables.
#'
#' @importFrom dplyr case_when
#' @importFrom ggplot2 ggplot geom_sf aes scale_fill_manual scale_color_manual
#' @importFrom ggplot2 ggtitle theme theme_void
#' @importFrom cowplot get_legend plot_grid
#' @export
#'
#' @examples
#'
#' example(readAndAddPolygonsToSPE)
#' spe <- spatialPerCellQC(spe)
#' spe <- computeThresholdFlags(spe)
#' p <- qcFlagPlots(spe, fov=16, theme="dark")
#' print(p)
qcFlagPlots <- function(spe, fov=unique(spe$fov),
                            theme=c("light","dark"), custom=FALSE) {
    if (!all(c("is_zero_counts", "is_ctrl_tot_outlier") %in%
            names(colData(spe))))
        stop("Fixed thresholds flag cells not found. Run computeFixedFlags().")
    spe <- computeSpatialOutlier(spe, computeBy="Area_um", method="both")
    spe <- computeSpatialOutlier(spe, computeBy="Mean.DAPI", method="both")
    spe$polygons$fixed_flags_color <- dplyr::case_when(
        spe$is_zero_counts ~ "0 counts",
        spe$is_ctrl_tot_outlier ~ "ctrl/total ratio > 0.1",
        TRUE ~ "unflagged"
    )
    if(all(spe$polygons$fixed_flags_color=="unflagged"))
        warning("No 0 counts or control/total ratio > 0.1 found")
    dapi_col <- .assign_outlier_color(spe$Mean.DAPI, spe$Mean.DAPI_outlier_mc,
                                    "> DAPI higher thr.", "< DAPI lower thr.")
    spe$polygons$dapi_outlier_color <- dapi_col
    area_vals <- if (custom) spe$cust_Area_um else spe$Area_um
    area_mc   <- if (custom) spe$cust_Area_um_outlier_mc else
        spe$Area_um_outlier_mc
    area_col  <- .assign_outlier_color(area_vals, area_mc,
        "> area um higher thr.", "< area um lower thr.")
    spe$polygons$area_outlier_color <- area_col
    spe$polygons$collapsed_color <-
        .assign_collapsed_color(spe$is_ctrl_tot_outlier, area_col, dapi_col)
    outlier_palette <- c(
        "unflagged"="#c0c8cf", "ctrl/total ratio > 0.1"="magenta",
        "< area um lower thr."="darkturquoise", "> area um higher thr."="red",
        "< DAPI lower thr."="purple", "> DAPI higher thr."="greenyellow"
    )
    plot_func <- if(theme[1]=="light") .light_theme else .dark_theme
    ggp1 <- .make_outlier_plot(spe$polygons, fov, "fixed_flags_color",
                            outlier_palette, "Control counts ratio") +
        plot_func()
    ggp2 <- .make_outlier_plot(spe$polygons, fov, "area_outlier_color",
                            outlier_palette, "Area in um") + plot_func()
    ggp3 <- .make_outlier_plot(spe$polygons, fov, "dapi_outlier_color",
                            outlier_palette, "Mean DAPI") + plot_func()
    # legp <- .make_outlier_plot(spe$polygons, fov, "collapsed_color",
    #                         outlier_palette, leg=TRUE) + plot_func() +
    #     ggplot2::theme(legend.title=element_blank())
    # ggp4 <- cowplot::get_legend(legp)
    # final <- cowplot::plot_grid(ggp1, ggp2, ggp3, ggp4, ncol=2)
    # if(theme[1]=="dark")
    # final <- final +
    #     ggplot2::theme(panel.background=element_rect(fill="black"))
    # final
    legp <- .make_outlier_plot(spe$polygons, fov, "collapsed_color",
                                outlier_palette, leg = TRUE) +
        plot_func() +
        ggplot2::theme(legend.title = ggplot2::element_blank()) +
        # scale "robuste" per evitare il drop delle classi assenti
        ggplot2::scale_fill_manual(
            values = outlier_palette,
            limits = names(outlier_palette),
            drop   = FALSE,
            na.value = "grey80"
        ) +
        ggplot2::scale_color_manual(
            values = outlier_palette,
            limits = names(outlier_palette),
            drop   = FALSE,
            na.value = "grey80"
        )

    # Estrazione sicura della legenda
    ggp4 <- tryCatch(
        cowplot::get_legend(legp),
        error = function(e) NULL
    )

    # Composizione finale
    if (!is.null(ggp4)) {
        final <- cowplot::plot_grid(ggp1, ggp2, ggp3, ggp4, ncol = 2)
    } else {
        final <- cowplot::plot_grid(ggp1, ggp2, ggp3, ncol = 2)
    }

    if (theme[1] == "dark") {
        final <- final +
            ggplot2::theme(panel.background = ggplot2::element_rect(fill="black"))
    }
    final


}


#' @importFrom ggplot2 geom_rect aes geom_text
#' @author Benedetta Banzi
.plotScaleBar <- function(p, spe, micronConvFact=0.12) {
    stopifnot(is(spe, "SpatialExperiment"))
    if (!inherits(p, "ggplot")) {
        stop("`p` must be a ggplot object (e.g. created by ggplot()).")
    }
    if (!is.numeric(micronConvFact) || length(micronConvFact) != 1 || micronConvFact <= 0) {
        stop("`micronConvFact` must be a single positive numeric value.")
    }
    x_vals <- spatialCoords(spe)[,1]
    y_vals <- spatialCoords(spe)[,2]
    x_range <- range(x_vals, na.rm = TRUE)
    x_length <- diff(x_range)
    if(metadata(spe)$technology %in% c("Nanostring_CosMx",
                                        "Nanostring_CosMx_Protein")){
        target <- x_length / 4 * micronConvFact
    } else {
        target <- x_length / 4
    }
    order <- 10^floor(log10(target))
    # guard rails: if target is NA or non-positive, fallback to 1*10^0
    if (!is.finite(order) || order <= 0) order <- 1
    mantissa <- target / order
    nice_mantissa <- if (mantissa < 1.5) {
        1
    } else if (mantissa < 3.5) {
        2
    } else if (mantissa < 7.5) {
        5
    } else {
        10
    }
    scale_length <- nice_mantissa * order
    scale_label <- if (scale_length >= 1000) {
        c(as.character(scale_length/2000), paste0(scale_length / 2000, " mm"))
    } else {
        micro <- "\u00B5"
        c(as.character(scale_length/2), paste0(scale_length/2, " ", micro, "m"))
    }
    if(metadata(spe)$technology %in% c("Nanostring_CosMx",
                                       "Nanostring_CosMx_Protein")){
        scale_length <- scale_length / micronConvFact
    }
    # Set coordinates for the scale bar
    x_start <- x_range[2] - scale_length * 1.1
    x_end <- x_range[2] - scale_length * 0.1
    y_pos <- min(y_vals, na.rm = TRUE) + scale_length * 0.1
    box_data <- data.frame(
        xmin = x_start - scale_length * 0.05,
        xmax = x_end + scale_length * 0.05,
        ymin = y_pos - scale_length * 0.05,
        ymax = y_pos + scale_length * 0.3
    )

    scale_bar_data <- data.frame(
        xmin = c(x_start, x_start + (scale_length/2)),
        xmax = c(x_start + (scale_length/2), x_start + (scale_length)),
        ymin = c(y_pos, y_pos),
        ymax = c(y_pos + (scale_length * 0.05), y_pos + (scale_length * 0.05))
    )

    label_data <- data.frame(
        x = c(x_start + (scale_length/4), x_start + 0.7*scale_length),
        y = c(y_pos + (scale_length * 0.2), y_pos + (scale_length * 0.2)),
        label = scale_label)

    p <- p+geom_rect(data = box_data, aes(xmin = xmin, xmax = xmax,
                                            ymin = ymin, ymax = ymax),
                       fill = "white", color = "grey", alpha = 0.8,
                    inherit.aes = FALSE) +
        geom_rect(data = scale_bar_data, aes(xmin = xmin, xmax = xmax,
                                             ymin = ymin, ymax = ymax),
                  fill = c("black", "white"), color = "grey",
                inherit.aes = FALSE) +
        geom_text(data = label_data,
                aes(x = x, y = y, label = label), color = "grey40",
                inherit.aes = FALSE, size = 3, hjust = 0.2)
    return(p)
}
