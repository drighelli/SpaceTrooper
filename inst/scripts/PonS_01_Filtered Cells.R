###############################################################################
# Plot on the slides the cells with low number of features and counts
###############################################################################

library(Seurat)
library(sf)
library(tmap)
library(raster)

project <- "S2-S3"
sample.name <- "R5363_S3_188-B"
image.type <- "CellOverlay"
image.extn <- ".jpg"
q.thr <- 0.05

### settings ###
  main.path <- file.path("/Users","silviobicciato","Documents",
                         "GeneChip_Data","AIRC_Piccolo","Nanostring_CosMx",project)
  main.data.dir <- file.path(main.path,"data",sample.name)
  spatVect.dir <- file.path(main.data.dir,"spatVect")
  image.dir <- file.path(main.data.dir,image.type)
  result.dir <- file.path(main.path,"results","QC and Cell phenotyping",sample.name)
  image.res.dir <- file.path(result.dir,"Filter_Overlays")
  if (!file.exists(image.res.dir)){dir.create(image.res.dir)}
  fov.positions <- read.csv(file.path(main.data.dir,paste0(sample.name,"_fov_positions_file.csv")),header = T)

############################################
### Creation of the unfiltered Seurat object
############################################
  CosMx.obj <- LoadNanostring(data.dir = main.data.dir, fov = "fov")
  meta.data <- CosMx.obj[[]]
  fov.subnames <- strsplit(rownames(meta.data), "_")
  fov.ids <- matrix(unlist(fov.subnames), ncol = 2, byrow = TRUE)
  meta.data$fov <- as.numeric(fov.ids[,2])
  meta.data$cell.id <- as.numeric(fov.ids[,1])

############################################
### Filter on nCount and nFeatures
############################################
  thr.count <- as.numeric(quantile(meta.data$nCount_Nanostring,q.thr))
  thr.feat <- as.numeric(quantile(meta.data$nFeature_Nanostring,q.thr))
  meta.filt <- meta.data[meta.data$nCount_Nanostring<thr.count&meta.data$nFeature_Nanostring<thr.feat,]

############################################
### Plotting and overlay with slides
############################################
  tmap_options(max.raster=c(plot = 3e6, view = 3e5))
  for (k in 1:length(unique(meta.filt$fov))){
    cells.to.plot <- data.frame(cell.id=meta.filt[meta.filt$fov==unique(meta.filt$fov)[k],"cell.id"])
    if (unique(meta.filt$fov)[k] < 10) {prefix <- "F00"} else {prefix <- "F0"}
    Spat.obj <- read_sf(file.path(spatVect.dir,paste0(prefix,unique(meta.filt$fov)[k])))
    names(Spat.obj)[1] <- "cell.id"
    Spat.obj <- merge(Spat.obj,cells.to.plot,by="cell.id")
    bg.image.file <- paste0(image.type,"_",prefix,unique(meta.data$fov)[k],".jpg")
    bg.image <- brick(file.path(image.dir,bg.image.file))
    tm.plot <- tm_shape(bg.image, raster.downsample = T) +
      tm_rgb(r=1, g=2, b=3) +
      tm_shape(Spat.obj) + 
      tm_fill(col="orange", alpha=1) +
      tm_borders(alpha=0.5, col = "white") +
      tm_layout(frame = FALSE, bg.color="black")
    tmap_save(tm.plot, file.path(image.res.dir,paste0("FilterOverlay_F0",unique(meta.filt$fov)[k],".jpg")),
              width = as.numeric(st_bbox(Spat.obj)[3]), 
              height = as.numeric(st_bbox(Spat.obj)[4]),
              units = "px", 
              outer.margins = F)
  }
  
###############################################################
### create the TileConfiguration file for Fiji Grip/stitching
### need to resize 50% first and then move the TileConfiguration
### file into the output directory and run Fiji Grip/stitching
###############################################################
  fov.positions$y_global_px <- -1*fov.positions$y_global_px
  filename <- file.path(image.res.dir,"TileConfiguration.txt")
  cat("# Define the number of dimensions we are working on","\n",file = filename)
  cat("dim = 2","\n",file=filename,append=T)
  cat("","\n",file=filename,append=T)
  cat("# Define the image coordinates","\n",file=filename,append=T)
  for (i in 1: dim(fov.positions)[1]){
    cat(paste0("FilterOverlay_F0",fov.positions$fov[i],".jpg","; ; (",
               fov.positions$x_global_px[i],", ",fov.positions$y_global_px[i],")"),"\n",file=filename,append=T)
  }
  
  
  