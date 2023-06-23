####################################################################
### Function that read ST images and generate a terra SpatVector ###
### Depends on terra  version 1.6-47                             ###
####################################################################

SpatVect_generate <- function(main.path,
                             image.type = image.type,
                             fov.number = fov.number,
                             flip.vertical = FALSE) {
  
  spatVect.dir <- file.path(main.path,"spatVect")
  if (!file.exists(spatVect.dir)){dir.create(spatVect.dir)}
# set the image directory and the image file name
  data.dir <- file.path(main.path,image.type)
  prefix <- "F0"
  if (fov.number < 10) {
    prefix <- "F00"
  }
# Read the image file file and generate the SpatRaster object and polygons
  file.ext <- tools::file_ext(dir(data.dir)[grep(paste0(prefix,fov.number), dir(data.dir), ignore.case=FALSE)])
  terra.rast <- terra::rast(file.path(data.dir,paste0(image.type,"_",prefix,fov.number,".",file.ext)))
  rast.dimensions <- dim(terra.rast)
  terra.polygon <- terra::as.polygons(terra.rast, value = TRUE)
  names(terra.polygon)<-"cellID"
# Remove the holes in SpatVector polygons
  terra.polygon <- terra::fillHoles(terra.polygon)
# Check the validity of polygons and remove invalid polygons
  valid.index <- terra::is.valid(terra.polygon)
  terra.polygon <- terra.polygon[valid.index]
# flip and shift vertical
  if (flip.vertical == TRUE) {
    terra.polygon <- terra::flip(terra.polygon, direction = "vertical")
    shift_horizontal_step <- 0
    shift_vertical_step <- dim(terra.rast)[1]
    terra.polygon = terra::shift(terra.polygon,
                                 dx = shift_horizontal_step, 
                                 dy = shift_vertical_step)
  }
# check and eventually remove an the background polygon 
  if (terra::values(terra.polygon[1,])==0){
    terra.polygon <- terra::subset(terra.polygon,terra.polygon$cellID!=0)
  }
# write the SpatVector file
  if (!file.exists(file.path(spatVect.dir,paste0(prefix,fov.number)))){
    terra::writeVector(terra.polygon,file.path(spatVect.dir,paste0(prefix,fov.number)))
  }
  return(terra.polygon)
}


