
#' plotCellsFovs
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
#' TBD
plotCellsFovs <- function(spe, point_col="darkmagenta", sample_id=NULL)
{
    stopifnot(is(spe, "SpatialExperiment"))
    # fov_positions <- data.table::fread(fovpos_file, header = T)
    # fov_positions <- fov_positions[fov_positions$fov%in%unique(metadata$fov),]
    if( !is.null(sample_id) ) spe <- spe[,spe$sample_id]
    ggp <- ggplot() +
        geom_point(data=as.data.frame(spatialCoords(spe)),
                mapping=aes_string(x=spatialCoordsNames(spe)[1],
                            y=spatialCoordsNames(spe)[2]),
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
                            color="black", fontface="bold", alpha=.7) +
        ggtitle(unique(spe$sample_id)) +
        .fov_image_theme(back.color="white", back.border="white",
                         title.col="black")
    return(ggp)
}

#' plotCentroidsSpe
#'
#' @param spe
#' @param colour_by
#' @param sample_id
#' @param point_col
#' @param size
#' @param alpha
#'
#' @return
#' @export
#' @importFrom ggplot2 geom_point aes_string theme_bw theme ggtitle
#' @importFrom scater plotColData
#'
#' @examples
plotCentroidsSpe <- function(spe, colour_by=NULL, order_by=NULL,
                        sample_id=unique(spe$sample_id),
                        isNegativeProbe=FALSE,
                        point_col="darkmagenta", size=0.05, alpha=0.2,
                        aspect_ratio=1)
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
            ggp <- ggp + scale_color_gradient(low="white", high="red",
                                              name=colour_by) +
                .negative_image_theme()
    }
    ggp <- ggp + ggtitle(sample_id) + theme(aspect.ratio=aspect_ratio,
                                            plot.title=element_text(hjust=0.5))

    if(!isNegativeProbe) ggp <- ggp + theme_bw()

    return(ggp)
}


plotMetricHist <- function(spe, metric, fill_color="#69b3a2",
                            bins=30, bin_width=NULL)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(metric %in% names(colData(spe)))

    ggp <- ggplot(data=as.data.frame(colData(spe))) +
            geom_histogram(aes(x=.data[[metric]]), fill=fill_color,
                    bins=bins, binwidth=bin_width) +
            ggtitle(metric) + theme_bw()

    return(ggp)
}
#' plotPolygonsSPE
#'
#' @param spe
#' @param title
#'
#' @return
#' @export
#' @importFrom tmap tm_shape tm_borders tm_layout
#'
#' @examples
plotPolygonsSPE <- function(spe, title=unique(spe$sample_id))
{
    stopifnot(all(is(spe, "SpatialExperiment"),
                  ("polygons" %in% colnames(colData(spe)))))

    tm_shape(spe$polygons) +
        tm_borders(lwd = 0.1, col = "grey50") +
        tm_layout(legend.outside = TRUE,
                  main.title.position = c("center", "top"),
                  main.title = title,
                  main.title.size = 0.5,
                  inner.margins = c(0, 0, 0, 0),
                  outer.margins = c(0, 0, 0, 0))
}
