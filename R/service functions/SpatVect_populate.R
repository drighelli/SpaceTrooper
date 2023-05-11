#######################################################################
### Function that add expression levels and meta-data to a SpatVect ###
### Depends on terra  version 1.6-47 ##################################
#######################################################################

SpatVect_populate <- function(main.path,
                              sample.name = sample.name,
                              fov.number = fov.number,
                              meta.to.add = NULL,
                              feat.to.add = NULL,
                              assay = "SCT") {
  
  data.dir <- file.path(main.path,"data",sample.name)
  spatVect.dir <- file.path(data.dir,"spatVect")
  
############################################
### Load the SpatVector for the given FOV  
############################################
  prefix <- "F0"
  if (fov.number < 10) {
    prefix <- "F00"
  }
  Spat.obj <- sf::read_sf(file.path(spatVect.dir,paste0(prefix,fov.number)))
  names(Spat.obj)[1] <- "cell.id"
  
####################################################
### Get the un-shifted centroids for the given FOV  
###################################################
  centr.geom <- terra::centroids(terra::vect(Spat.obj))
  # Convert the SpatVector of centroids to a data table
  CT.geom <- data.table::as.data.table(terra::geom(centr.geom))
  CT.geom <- CT.geom[,c(1,3,4)]
  names(CT.geom) <- c("cell.id","sdimx", "sdimy")
  
####################################################
### Add metadata to the SpatVect  
####################################################
  if (!is.null(meta.to.add)){
    meta.data <- readRDS(file = file.path(data.dir,paste0(sample.name,"_metadata.rds")))
    meta.data.fov <- meta.data[meta.data$fov==fov.number,c("cell.id",meta.to.add)]
    MD.to.add <- merge(CT.geom, meta.data.fov, by="cell.id")
    Spat.obj <- merge(Spat.obj, MD.to.add, by = "cell.id")
  }
  
####################################################
### Add features to the SpatVect  
####################################################
  if (!is.null(feat.to.add)){
    CosMx.obj <- readRDS(file = file.path(data.dir,paste0(sample.name,"_data.rds")))
    feature.data <- Seurat::GetAssayData(object = CosMx.obj, assay = assay)
    if (!is.null(meta.to.add)){
      FT.to.add <- feature.data[feat.to.add, rownames(meta.data.fov)]
      if (length(feat.to.add) == 1){
        FT.to.add <- data.frame(cell.id = meta.data.fov$cell.id, FT.to.add)
        colnames(FT.to.add)[2] <- feat.to.add
      } else {
        FT.to.add <- data.frame(cell.id = meta.data.fov$cell.id, t(FT.to.add))
      }
      Spat.obj <- merge(Spat.obj, FT.to.add, by = "cell.id")
    } else {
      meta.data <- CosMx.obj[[]]
      meta.data.fov <- meta.data[meta.data$fov==fov.number,]
      FT.to.add <- feature.data[feat.to.add, rownames(meta.data.fov)]
      if (length(feat.to.add) == 1){
        FT.to.add <- data.frame(cell.id = meta.data.fov$cell.id, FT.to.add)
        colnames(FT.to.add)[2] <- feat.to.add
      } else {
        FT.to.add <- data.frame(cell.id = meta.data.fov$cell.id, t(FT.to.add))
      }
      FT.to.add <- merge(CT.geom, FT.to.add, by="cell.id")
      Spat.obj <- merge(Spat.obj, FT.to.add, by = "cell.id")
    }
  }
  
################################################################
### Construct and return the sf object  
################################################################
  Spat.obj$cell.id<-paste(Spat.obj$cell.id,fov.number,sep="_")
  poly.geom <- sf::st_as_sf(Spat.obj)
  return(poly.geom)
}

   