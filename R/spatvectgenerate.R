####################################################################
### spatVectGenerate reads ST images and generate a terra        ###
### SpatVector                                                   ###
### Wrote with terra  version 1.7-29                             ###
### Written by Silvio Bicciato, adapted by Dario Righelli        ###
####################################################################
#' spatVectGenerate
#' @description
#' Generates/loads the shape file for each FoV image present in the
#' `main.path\image.type`.
#' See details for additional information.
#'
#' @details
#' Generates/loads the shape file for each FoV image present in the
#' `main.path\image.type`.
#' In case the shape file already exists, the function loads the already
#' existing shape file, returning the loaded object.
#' If the `force` param is set to `TRUE` it computes the polygons even when the
#' file already exists and also overwrites it.
#' @param main.path character for the main path of the CosMx sample
#' @param image.type character for the images folder name (default is
#' `CellLabels`)
#' @param fov.number numeric for the FoV id to process
#' @param flip.vertical logical indicating if the image has to be vertically
#' flipped
#' @param force logical, when TRUE it computes the fov polygons even if the file
#' already exists, it also overwrites it.
#'
#' @return a terra polygon object or an sf object in case of reading the object
#' @export
#' @importFrom tools file_ext
#' @importFrom terra rast fillHoles is.valid flip shift subset write.vector
#' @importFrom sf read_sf
#' @examples
spatVectGenerate <- function(main.path, image.type="CellLabels", fov.number,
                             flip.vertical=FALSE, force=FALSE)
{
    spatVect.dir <- file.path(main.path, "spatVect")
    if (!file.exists(spatVect.dir)) { dir.create(spatVect.dir) }
    # set the image directory and the image file name
    data.dir <- file.path(main.path, image.type)
    prefix <- ifelse(fov.number < 10, "F00", "F0")
    final.filename <- file.path(spatVect.dir, paste0(prefix, fov.number))
    if ( !file.exists(final.filename) )
    {
        terra.polygon <- .spatVectGenerate(data.dir=data.dir, prefix=prefix,
            fov.number=fov.number, flip.vertical=flip.vertical)
        # write the SpatVector file
        terra::writeVector(terra.polygon, final.filename)
    } else {
        if (force)
        {
            terra.polygon <- .spatVectGenerate(data.dir=data.dir, prefix=prefix,
                fov.number=fov.number, flip.vertical=flip.vertical)
            # write the SpatVector file
            terra::writeVector(terra.polygon, final.filename)
        } else {
            ## read the terra.polygon from the existing file
            terra.polygon <- sf::read_sf(final.filename)
            names(terra.polygon)[1] <- "cell.id"
        }
    }
    return(terra.polygon)
}

.spatVectGenerate <- function(data.dir, prefix, fov.number, flip.vertical)
{
    # Read the image file and generate the SpatRaster object and polygons
    fov.tif <- list.files(data.dir)
    file.ext <- tools::file_ext(fov.tif[grep(paste0(prefix, fov.number),
                                             fov.tif, ignore.case=FALSE)])
    # create raster for the image
    terra.rast <- terra::rast(file.path(data.dir,
            paste0(image.type, "_", prefix, fov.number, ".", file.ext)))
    rast.dimensions <- dim(terra.rast)
    terra.polygon <- terra::as.polygons(terra.rast, value = TRUE)
    names(terra.polygon) <- "cellID"
    # Remove the holes in SpatVector polygons
    terra.polygon <- terra::fillHoles(terra.polygon)
    # Check the validity of polygons and remove invalid polygons
    valid.index <- terra::is.valid(terra.polygon)
    terra.polygon <- terra.polygon[valid.index]
    # flip and shift vertical
    if (flip.vertical == TRUE)
    {
        terra.polygon <- terra::flip(terra.polygon, direction = "vertical")
        shift_horizontal_step <- 0
        shift_vertical_step <- dim(terra.rast)[1]
        terra.polygon = terra::shift(terra.polygon,
                                     dx = shift_horizontal_step,
                                     dy = shift_vertical_step)
    }
    # check and eventually remove an the background polygon
    if (terra::values(terra.polygon[1,])==0)
    {
        terra.polygon <- terra::subset(terra.polygon, terra.polygon$cellID!=0)
    }
    return(terra.polygon)
}
