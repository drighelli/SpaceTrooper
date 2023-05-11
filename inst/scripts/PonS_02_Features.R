###############################################################################
### Plot gene expression levels on images
###############################################################################

library(Seurat)
library(sf)
library(tmap)
library(raster)
library(RColorBrewer)

project <- "S2-S3"
sample.name <- "R5363_S2_187-B"
image.type <- "CellComposite"
image.extn <- ".jpg"

feat.to.plot <- "KRT19"
overlay <- TRUE

### settings ###
  main.path <- file.path("/Users","silviobicciato","Documents",
                         "GeneChip_Data","AIRC_Piccolo","Nanostring_CosMx",project)
  main.data.dir <- file.path(main.path,"data",sample.name)
  spatVect.dir <- file.path(main.data.dir,"spatVect")
  image.dir <- file.path(main.data.dir,image.type)
  result.dir <- file.path(main.path,"results","QC and Cell phenotyping",sample.name)
  image.res.dir <- file.path(result.dir,paste0(feat.to.plot,"_Overlays"))
  if (!file.exists(image.res.dir)){dir.create(image.res.dir)}
  fov.positions <- read.csv(file.path(main.data.dir,paste0(sample.name,"_fov_positions_file.csv")),header = T)

####################################
### Data loading and subsetting
####################################
# load Seurat CosMx object #
  CosMx.obj <- readRDS(file = file.path(main.data.dir,paste0(sample.name,"_data.rds")))
  feature.data <- GetAssayData(object = CosMx.obj, assay = "SCT")
  meta.data <- CosMx.obj[[]]
  fov.subnames <- strsplit(rownames(meta.data), "_")
  fov.ids <- matrix(unlist(fov.subnames), ncol = 2, byrow = TRUE)
  meta.data$fov <- as.numeric(fov.ids[,2])
  meta.data$cell.id <- as.numeric(fov.ids[,1])
  
############################################
### Plotting and overlay with slides
############################################
  tmap_options(max.raster=c(plot = 3e6, view = 3e5))
  cols <- brewer.pal(7, "Reds")
  for (k in 1:length(unique(meta.data$fov))){
    feat.fov <- data.frame(cell.id=meta.data[meta.data$fov==unique(meta.data$fov)[k],"cell.id"],
                           feature.data[feat.to.plot,rownames(meta.data[meta.data$fov==unique(meta.data$fov)[k],])])
    colnames(feat.fov)[2] <- feat.to.plot
    if (unique(meta.data$fov)[k] < 10) {prefix <- "F00"} else {prefix <- "F0"}
    Spat.obj <- read_sf(file.path(spatVect.dir,paste0(prefix,unique(meta.data$fov)[k])))
    names(Spat.obj)[1] <- "cell.id"
    Spat.obj <- merge(Spat.obj,feat.fov,by="cell.id")
    # select only cells with an expression level greater than 0 for the given feature
    Spat.obj.filt <- Spat.obj %>%
      dplyr::filter(Spat.obj[[feat.to.plot]]!=0)
    if (overlay == FALSE){
       tm.plot <- tm_shape(Spat.obj) + 
         tm_fill(col="black") + 
         tm_borders(lwd = 0.3, alpha=0.5, col = "grey")
       file.name <- file.path(image.res.dir,paste0(feat.to.plot,"_F0",unique(meta.data$fov)[k],".jpg"))
    }
    # load the background image and overlay the expression level 
    if (overlay == TRUE){
      bg.image.file <- paste0(image.type,"_",prefix,unique(meta.data$fov)[k],".jpg")
      bg.image <- brick(file.path(image.dir,bg.image.file))
      tm.plot <- tm_shape(bg.image, raster.downsample = T) +
        tm_rgb(r=1, g=2, b=3, alpha = 0.5)
      file.name <- file.path(image.res.dir,paste0(feat.to.plot,"_Overlay_F0",unique(meta.data$fov)[k],".jpg"))
    }
    tm.plot <- tm.plot +
      tm_shape(Spat.obj.filt) + 
      tm_fill(col = feat.to.plot,
              style = "pretty",
              n = 7,
              palette = cols,
              legend.hist = FALSE,
              title = feat.to.plot) +
      tm_borders(alpha=0.5, col = "black") +
      tm_layout(frame = FALSE,
                bg.color="black",
                legend.text.color= "white", 
                outer.margins = c(0,0,0,0),
                inner.margins = c(0,0,0,0),
                legend.outside =F)
    tmap_save(tm.plot, file.name,
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
  if (overlay == TRUE){
    middle.text <- "_Overlay_F0"
    filename <- file.path(image.res.dir,"Overlay_TileConfiguration.txt")
  } else {
    middle.text <- "_F0"
    filename <- file.path(image.res.dir,"TileConfiguration.txt")
  }
  cat("# Define the number of dimensions we are working on","\n",file = filename)
  cat("dim = 2","\n",file=filename,append=T)
  cat("","\n",file=filename,append=T)
  cat("# Define the image coordinates","\n",file=filename,append=T)
  for (i in 1: dim(fov.positions)[1]){
    cat(paste0(feat.to.plot,middle.text,fov.positions[i,1],".jpg","; ; (",
               fov.positions$x_global_px[i],", ",fov.positions$y_global_px[i],")"),"\n",file=filename,append=T)
  }
  
  
  