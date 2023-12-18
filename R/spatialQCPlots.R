#' .fov_image_theme
#' @description
#' internal function to setup the theme for the fov background on the whole
#' image
#'
#' @param back.color not used
#' @param back.border color for the borders of the background (default=NA)
#' @param title.col character indicating the color of the title
#'
#' @return a ggplot2 theme object
#' @importFrom ggplot2 theme element_blank element_rect element_text
#' @keywords internal
.fov_image_theme <- function(back.color="black", back.border=NA,
                            title.col="white")
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
          panel.background=element_rect(fill = "transparent", colour = NA),
          plot.title=element_text(color=title.col, hjust=0.5, face="bold"),
          plot.background=element_rect(fill="transparent", colour=back.border))
}


#' plotCellsFovs
#' @description
#' Creates a scatter plot with all the centroids of the cells distributed inside
#' their own FoVs as an overlapping grid.
#'
#' @param spe A SpatialExperiment object
#' @param point_col the color to use for the dots in the plot (default is "darkmagenta")
#' @param sample_id string indicating the sample to plot (default is NULL,
#' implying a single sample in the spe)
#'
#' @return a ggplot2 object
#' @importFrom ggplot2 ggplot aes aes_string annotate geom_point geom_text ggtitle
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
                      label=metadata(spe)$fov_positions[1][ , , drop=TRUE],
                            color="black", fontface="bold")) +
        ggtitle(unique(spe$sample_id)) +
        .fov_image_theme(back.color="white", back.border="white", title.col="black")
    return(ggp)
}
