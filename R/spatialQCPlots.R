
#' plotCellsFovs
#'
#' @description
#' Creates a scatter plot with all the centroids of the cells distributed inside
#' their own FoVs as an overlapping grid.
#'
#' @param spe A SpatialExperiment object
#' @param point_col the color to use for the dots in the plot
#' (default is "darkmagenta")
#' @param sample_id string indicating the sample to plot (default is NULL,
#' implying a single sample in the spe)
#'
#' @return a ggplot2 object
#' @importFrom ggplot2 ggplot aes aes_string annotate geom_point geom_text
#' ggtitle
#' @import SpatialExperiment
#' @export
#'
#' @examples
#' # TBD
plotCellsFovs <- function(spe, sample_id=unique(spe$sample_id),
                point_col="darkmagenta", numbers_col="black", alpha_numbers=0.8)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("fov" %in% names(colData(spe)))
    # fov_positions <- data.table::fread(fovpos_file, header = T)
    # fov_positions <- fov_positions[fov_positions$fov%in%unique(metadata$fov),]
    # if( !is.null(sample_id) ) spe <- spe[,spe$sample_id]
    spd <- as.data.frame(spatialCoords(spe))
    x_coord <- spatialCoordsNames(spe)[1]
    y_coord <- spatialCoordsNames(spe)[2]
    ggp <- ggplot() +
        geom_point(data=spd,
                mapping=aes(x=.data[[x_coord]],
                            y=.data[[y_coord]]),
                            colour=point_col,
                            fill=point_col,
                            size=0.05, alpha=0.2) +
        annotate("rect",
                 xmin=metadata(spe)$fov_positions[2][ , , drop=TRUE],
                 xmax=metadata(spe)$fov_positions[2][ , , drop=TRUE] +
                     metadata(spe)$fov_dim[["xdim"]],
                 ymin=metadata(spe)$fov_positions[3][ , , drop=TRUE],
                 ymax=metadata(spe)$fov_positions[3][ , , drop=TRUE] +
                     metadata(spe)$fov_dim[["ydim"]],
                 alpha=.2, color="black", linewidth=0.2) +
        geom_text(aes(x=metadata(spe)$fov_positions[2][ , , drop=TRUE] +
                          metadata(spe)$fov_dim[["xdim"]]/2,
                      y=metadata(spe)$fov_positions[3][ , , drop=TRUE] +
                          metadata(spe)$fov_dim[["ydim"]]/2,
                      label=metadata(spe)$fov_positions[1][ , , drop=TRUE]),
                            color=numbers_col, fontface="bold",
                            alpha=alpha_numbers) +
        ggtitle(sample_id) +
        .fov_image_theme(back.color="white", back.border="white",
                         title.col="black")
    return(ggp)
}

#' plotCentroidsSPE
#'
#' @description
#' Plot Spatial Coordinates for a SpatialExperiment Object
#' This function generates a ggplot of spatial coordinates from a
#' `SpatialExperiment` object, optionally coloring the points by a specified
#' column in `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial
#' transcriptomics data.
#' @param colour_by An optional character string specifying the column in
#' `colData(spe)` to use for coloring the points. If `NULL`, all points will be
#' colored the same.
#' @param order_by An optional character string specifying the column in
#' `colData(spe)` to use for ordering the points.
#' @param sample_id A character string specifying the sample identifier to be
#' used as the plot title. (Default is the unique sample ID from `spe`)
#' @param isNegativeProbe A logical value indicating whether to apply a custom
#' color gradient for negative probe data. (Default is `FALSE`)
#' @param palette A vector of colors to be used as a custom palette. For
#' categorical data, this should be a vector of colors with the same length as
#' the number of levels in `colour_by`. For continuous data, this should be a
#' vector of colors used to create a gradient.
#' @param point_col A character string specifying the color of the points when
#' `colour_by` is `NULL`. (Default is `"darkmagenta"`)
#' @param size A numeric value specifying the size of the points. (Default is
#' `0.05`)
#' @param alpha A numeric value specifying the transparency level of the points.
#' (Default is `0.2`)
#' @param aspect_ratio A numeric value specifying the aspect ratio of the plot.
#' (Default is `1`)
#' @param legend_point_size A numeric value specifying the size of the points
#' in the legend. (Default is `2`)
#' @param legend_point_alpha A numeric value specifying the transparency level
#' of the points in the legend. (Default is `0.8`)
#'
#' @return A `ggplot` object representing the spatial coordinates plot of
#' polygon centroids.
#'
#' @importFrom ggplot2 geom_point aes_string theme_bw theme ggtitle guides
#' guide_legend scale_color_gradient scale_color_manual scale_color_gradientn
#' @importFrom scater plotColData
#' @importFrom SummarizedExperiment colData
#' @export
#'
#' @examples
#' #TBD
plotCentroidsSPE <- function(spe, colour_by=NULL, order_by=NULL,
                        sample_id=unique(spe$sample_id),
                        isNegativeProbe=FALSE, palette=NULL,
                        point_col="darkmagenta", size=0.05, alpha=0.2,
                        aspect_ratio=1,
                        legend_point_size=2, legend_point_alpha=0.8)
{

    stopifnot( all( is(spe, "SpatialExperiment"),
                    (colour_by %in% names(colData(spe)))))
    if(is.null(colour_by))
    {
        ggp <- ggplot() +
            geom_point(data=as.data.frame(spatialCoords(spe)),
                       mapping=aes_string(x=spatialCoordsNames(spe)[1],
                                          y=spatialCoordsNames(spe)[2]),
                       colour=point_col,
                       fill=point_col,
                       size=size, alpha=alpha)
    } else {
        ## check if column variable is logical to impose our colors
        ggp <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                    y=spatialCoordsNames(spe)[2],
                    colour_by=colour_by, order_by=order_by,
                    point_size=size, point_alpha=alpha)
        if(isNegativeProbe)
        {
            ggp <- ggp + scale_color_gradient(low="white", high="red",
                                              name=colour_by) +
                .negative_image_theme()
        } else if(all(!is.null(palette), (palette %in% names(colData(spe))))) {
            palette <- createPaletteFromColData(spe, palette_names=colour_by,
                                                    palette_colors=palette)
            if(is.factor(colData(spe)[[colour_by]]))
            {
                ggp <- ggp + scale_color_manual(values=palette)
            } else if(is.numeric(colData(spe)[[colour_by]])) {
                ggp <- ggp + scale_color_gradientn(colors=palette)
            }
        }
    }


    ggp <- ggp + ggtitle(sample_id) +
        theme(aspect.ratio=aspect_ratio, plot.title=element_text(hjust=0.5))+
        guides(colour=guide_legend(override.aes=list(size=legend_point_size,
                                                    alpha=legend_point_alpha)))

    if(!isNegativeProbe) ggp <- ggp + theme_bw()

    return(ggp)
}


#' plotMetricHists
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
#' @param fill_color A character string specifying the fill color of the
#' histogram bars. (Default is `"#69b3a2"`)
#' @param use_fences A character string specifying the name of the column in
#' `colData(spe)` that contains the fence thresholds (typically from an outlier
#' filter). If `NULL`, no fences will be plotted. (Default is `NULL`)
#' @param fences_colors A named character vector specifying the colors to use
#' for the lower and higher fences. The names should be `"lower"` and `"higher"`
#'. (Default is `c("lower"="purple4", "higher"="tomato")`)
#' @param bins An integer specifying the number of bins to use in the histogram.
#' (Default is `30`)
#' @param bin_width A numeric value specifying the width of the bins. If `NULL`,
#' the bin width will be automatically determined based on the `bins` parameter.
#' (Default is `NULL`)
#'
#' @return A `ggplot` object representing the histogram of the specified metric.
#'
#' @importFrom ggplot2 ggplot geom_histogram aes ggtitle theme_bw geom_vline
#' labs scale_colour_manual
#' @importFrom SummarizedExperiment colData
#' @export
#' @example
#' #TBD
plotMetricHist <- function(spe, metric, fill_color="#69b3a2",
        use_fences=NULL, fences_colors=c("lower"="purple4", "higher"="tomato"),
        bins=30, bin_width=NULL)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(metric %in% names(colData(spe)))

    ggp <- ggplot(data=as.data.frame(colData(spe))) +
            geom_histogram(aes(x=.data[[metric]]), fill=fill_color,
                    bins=bins, binwidth=bin_width)
    if (!is.null(use_fences))
    {
        stopifnot(use_fences %in% names(colData(spe)))
        fences <- getFencesOutlier(spe, use_fences, "both", 2)
        fences_labs <- paste0(names(fences), ": ", fences)
        names(fences_colors) <- fences_labs
        ggp <- ggp +
            geom_vline(aes(xintercept=fences[1], color=fences_labs[1])) +
            geom_vline(aes(xintercept=fences[2], color=fences_labs[2])) +
            labs(color=paste0("Fences ", use_fences)) +
            scale_colour_manual(values=fences_colors)
    }
    ggp <- ggp + ggtitle(metric) + theme_bw()

    return(ggp)
}


#' plotPolygonsSPE
#'
#' @description
#' Plot Polygons from a SpatialExperiment Object.
#' This function generates a plot of polygons stored in a `SpatialExperiment`
#' object.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data, including polygon data as an `sf` object.
#' @param colour_by An optional character string specifying the column in
#' `colData(spe)` to use for coloring the polygons. If `NULL`, all polygons
#' will be colored the same. Default is `NULL`.
#' @param sample_id A character string specifying the title of the plot.
#' Default is the unique sample ID from `spe`.
#' @param fill_alpha A numeric value specifying the transparency level of the
#' polygon fill. Default is `NA`.
#' @param palette A character vector specifying the colors to use for filling
#' the polygons when `colour_by` is a factor. Default is `NULL`.
#' @param border_col A character string specifying the color of the polygon
#' borders. Default is `NA`.
#' @param border_alpha A numeric value specifying the transparency level of the
#' polygon borders. Default is `NA`.
#' @param border_line_width A numeric value specifying the width of the polygon
#' borders. Default is `0.1`.
#'
#' @return A `tmap` plot object displaying the polygons.
#'
#' @importFrom tmap tm_shape tm_borders tm_layout tm_fill tm_polygons
#' @importFrom SummarizedExperiment colData
#' @export
#'
#' @examples
#' # Assuming `spe` is a SpatialExperiment object with polygon data:
#' # plotPolygonsSPE(spe, colour_by="gene_expression")
plotPolygonsSPE <- function(spe, colour_by=NULL,sample_id=unique(spe$sample_id),
                    fill_alpha=NA, palette=NULL, border_col=NA, border_alpha=NA,
                    border_line_width=0.1, bg_color="black")
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("polygons" %in% names(colData(spe)))

    pols <- spe$polygons
    if (!is.null(colour_by))
    {
        stopifnot(colour_by %in% names(colData(spe)))
        if (is(spe[[colour_by]], "logical"))
        {
            sums <- sum(colData(spe)[[colour_by]])
            pols[[colour_by]] <- ifelse(colData(spe)[[colour_by]]==TRUE,
                                paste0("TRUE (", sums,")"),
                                paste0("FALSE (", dim(pols)[1]-sums,")"))
        } else {
            pols[[colour_by]] <- colData(spe)[[colour_by]]
        }
    }
    tmm <- tm_shape(pols)

    if (is.null(colour_by))
    {
        colour_by="grey50"
        border_line_width=0.1
    }
    if(!is.null(palette)) if (palette %in% names(colData(spe)))
    {
        palette <- createPaletteFromColData(spe, palette_names=colour_by,
                                            palette_colors=palette)
        palette <- setNames(palette, NULL)
    }
    tmm <- tmm + tm_polygons(col=colour_by, alpha=fill_alpha, palette=palette,
                             border.col=border_col, border.alpha=border_alpha,
                             lwd=border_line_width)

    tmm <- tmm + tm_layout(legend.outside=TRUE,
                            main.title.position=c("left", "top"),
                            main.title = sample_id,
                            main.title.fontface = 2,
                            main.title.size = 1,
                            inner.margins = c(0, 0, 0, 0),
                            outer.margins = c(0, 0, 0, 0),
                            bg.color=bg_color)
    return(tmm)
}


#' plotPolygonsSPE_ggplot
#'
#' @description Plot polygons from a `SpatialExperiment` object using ggplot2.
#'
#' @param spe A `SpatialExperiment` object with polygon data as an `sf` object.
#' @param colour_by A column in `colData(spe)` for coloring the polygons.
#' @param sample_id Sample ID for plot title. Default is the unique sample ID.
#' @param fill_alpha Transparency level for polygon fill. Default is `1`.
#' @param palette Colors to use if `colour_by` is a factor. Default is `NULL`.
#' @param border_col Color of polygon borders. Default is `"black"`.
#' @param border_alpha Transparency level for borders. Default is `1`.
#' @param border_line_width Width of polygon borders. Default is `0.1`.
#' @param draw_borders Logical; whether to draw borders. Default is `TRUE`.
#'
#' @return A `ggplot` object representing the polygon plot of the spatial data.
#' @export
#'
#' @importFrom ggplot2 ggplot geom_sf aes scale_fill_manual scale_fill_viridis_c
#' scale_fill_identity theme_minimal theme element_text margin labs
#' @importFrom sf st_as_sf
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' # Assuming `spe` is a SpatialExperiment object with polygon data:
#' # plotPolygonsSPE_ggplot(spe, colour_by="gene_expression")
plotPolygonsSPE_ggplot <- function(spe, colour_by=NULL,
                                   sample_id=unique(spe$sample_id),
                                   fill_alpha=1, palette=NULL,
                                   border_col=NA,
                                   border_alpha=1,
                                   border_line_width=0.1,
                                   draw_borders=TRUE) {
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("polygons" %in% names(colData(spe)))
    stopifnot(!is.null(colour_by))
    pols <- spe$polygons

    if(!is.null(colour_by)) {
        stopifnot(colour_by %in% names(colData(spe)))
        pols[[colour_by]] <- colData(spe)[[colour_by]]
    }

    border_params <- if(draw_borders) {
        list(color=border_col, size=border_line_width)
    } else {
        list(color=NA, size=0)
    }

    p <- ggplot(pols) + geom_sf(aes(fill=.data[[colour_by]]),
                     alpha=fill_alpha,
                     color=border_params$color,
                     size=border_params$size)

    if(!is.null(colour_by) && is.factor(pols[[colour_by]])) {
        if(!is.null(palette)) {
            p <- p + scale_fill_manual(values=palette)
        }
    } else if(!is.null(colour_by) && is.numeric(pols[[colour_by]])) {
        p <- p + scale_fill_viridis_c(option="D")
    } else {
        p <- p + scale_fill_identity()
    }

    p <- p + theme_minimal() +
        theme(
            legend.position="right",
            plot.title.position="plot",
            plot.title=element_text(face="bold", size=14),
            plot.margin=margin(0, 0, 0, 0)
        ) +
        labs(title=sample_id, fill=colour_by)

    return(p)
}


#' plotZoomFovsMap
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
#' @param colour_by An optional character string specifying the column in
#' `colData(spe)` to use for coloring the polygons. Default is `NULL`.
#' @param map_point_col A character string specifying the color of the points
#' in the map. Default is `"darkmagenta"`.
#' @param map_numbers_col A character string specifying the color of the
#' numbers on the map. Default is `"black"`.
#' @param map_alpha_numbers A numeric value specifying the transparency of the
#' numbers on the map. Default is `0.8`.
#' @param title An optional character string specifying the title of the final
#' plot. If `NULL`, no title is added. Default is `NULL`.
#' @param ... Additional arguments passed to `plotPolygonsSPE`.
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
#' @importFrom tmap tmap_grob
#' @export
#'
#' @examples
#' # Assuming 'spe' is a SpatialExperiment object with FOVs and polygon data:
#' # plotZoomFovsMap(spe, fovs = c("FOV1", "FOV2"), colour_by = "cell_type",
#' #                title = "Zoomed FOVs with Polygons")
plotZoomFovsMap <- function(spe, fovs = NULL, colour_by = NULL,
                            map_point_col = "darkmagenta",
                            map_numbers_col = "black",
                            map_alpha_numbers = 0.8,
                            title = NULL, ..., useggplot=TRUE)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("fov" %in% names(colData(spe)))
    stopifnot(all(fovs %in% spe$fov))
    stopifnot(!is.null(colour_by))

    spefovs <- spe[, spe$fov %in% fovs]

    map <- plotCellsFovs(spefovs, point_col=map_point_col,
                         numbers_col=map_numbers_col,
                         alpha_numbers=map_alpha_numbers,
                         sample_id=NULL)

    if(useggplot)
    {
        g2 <- plotPolygonsSPE_ggplot(spefovs, colour_by=colour_by,
                                     sample_id=NULL, ...)
    } else {
        g2 <- plotPolygonsSPE(spefovs, colour_by=colour_by,
                                     sample_id=NULL, ...)
        g2 <- tmap_grob(g2)
    }

    final_plot <- ggpubr::ggarrange(map, g2, ncol=2)

    if (!is.null(title))
    {
        final_plot <- ggpubr::annotate_figure(final_plot,
                            top=ggpubr::text_grob(title, face="bold", size=14))
    }

    return(final_plot)
}



