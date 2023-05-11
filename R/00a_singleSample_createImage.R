####################################################################################
### Script to generate the polygons and centroids file for Seurat LoadNanostring ###
####################################################################################

library(sf)
library(tmap)

project <- "CosMx_Lung"
sample.name <- "Lung5_Rep1"
type.of.image <- "CellLabels"

### directory settings ###
  main.path <- file.path("/Users","silviobicciato","Documents","GeneChip_Data",project)
  data.dir <- file.path(main.path,"data",sample.name)
  meta.data.all <- read.csv(file.path(data.dir,paste0(sample.name,"_metadata_file.csv")),header = T)
  fov.positions <- read.csv(file.path(data.dir,paste0(sample.name,"_fov_positions_file.csv")),header = T)
  
### merge and plot all analyzed FOVs ###
  fov.positions <- fov.positions[fov.positions$fov%in%unique(meta.data.all$fov),]
  Spat.obj <- NULL
  for (k in 1: dim(fov.positions)[1]){
    if (fov.positions[k,"fov"] < 10) {prefix <- "F00"} else {prefix <- "F0"}
    print(paste0("Adding ",prefix,fov.positions[k,"fov"]))
    temp.obj <- read_sf(file.path(data.dir,"spatVect",paste0(prefix,fov.positions[k,"fov"])))
    # shift local FOV coordinates to global
    st_geometry(temp.obj) <- temp.obj$geometry + c(fov.positions$x_global_px[k],fov.positions$y_global_px[k]) 
    Spat.obj <- rbind(Spat.obj, temp.obj)
  }

### plot global image ###  
  tm_shape(Spat.obj) + 
    tm_fill(col="black") + 
    tm_borders(lwd = 0.3, alpha=0.5, col = "grey")

### save global spatVect ###
  global.image.dir <- file.path(data.dir,"spatVect","global")
  if (!file.exists(global.image.dir)){dir.create(global.image.dir)}
  st_write(Spat.obj, file.path(global.image.dir,paste0(sample.name,".shp")))
  
  # temp.obj <- read_sf(global.image.dir)
  