##################################################################
### Eventually add single images for each FOV in the Seuat object
##################################################################

library(Seurat)
library(data.table)
  
project <- "S2-S3"

### sample settings ###
  sample.name <- "R5363_S2_187-B"
  field.of.view <- "S2"

### directories ###
  main.path <- file.path("/Users","silviobicciato","Documents","GeneChip_Data",
                         "AIRC_Piccolo","Nanostring_CosMx",project)
  data.dir <- file.path(main.path,"data")
  
########################################
### Load Seurat object and polygon file 
########################################
  CosMx.obj <- readRDS(file = file.path(data.dir,sample.name,paste0(sample.name,"_data.rds")))
  coords <- fread(file = file.path(data.dir,sample.name,paste0(sample.name,"-polygons.csv")))

########################################
### Construct an image for each FOV 
########################################
  for (i in 1: length(unique(coords$fov))){
    print(paste0("adding FOV ",unique(coords$fov)[i]))
    # extract centroid coordinates
    centr.coords<-coords[coords$fov==unique(coords$fov)[i],c("cellID","sdimx","sdimy")]
    colnames(centr.coords) <- c("cell","x", "y")
    centr.coords$cell <- paste0(centr.coords$cell,"_",unique(coords$fov)[i])
    centr.coords <- data.frame(cell= unique(centr.coords$cell),
                               x = unique(centr.coords$x),
                               y = unique(centr.coords$y))
    # extract segment coordinates
    poly.coords<-coords[coords$fov==unique(coords$fov)[i],c("cellID","x_global_px","y_global_px")]
    colnames(poly.coords) <- c("cell","x", "y")
    poly.coords$cell <- paste0(poly.coords$cell,"_",unique(coords$fov)[i])
    # get segment and centroid coordinates for cells in the Seurat object
    centr.coords <- centr.coords[centr.coords$cell%in%intersect(centr.coords$cell,colnames(CosMx.obj)),]
    poly.coords <- poly.coords[poly.coords$cell%in%intersect(poly.coords$cell,colnames(CosMx.obj)),]
    # create and add the FOV
    cents <- CreateCentroids(centr.coords)
    segs <- CreateSegmentation(poly.coords)
    segmentations.data <- list(centroids = cents, segmentation = segs)
    fov.coords <- CreateFOV(coords = segmentations.data,
                            type = c("centroids","segmentation"), 
                            assay = "SCT",
                            key = paste0("fov",unique(coords$fov)[i],"_"))
    CosMx.obj[[paste0("fov",unique(coords$fov)[i])]] <- fov.coords
  }
  
#########################################################
### Remove the complete image and save the Seurat object 
#########################################################
  CosMx.obj@images[field.of.view] <- NULL
  saveRDS(CosMx.obj, file = file.path(data.dir,sample.name,paste0(sample.name,"_data.rds")))

  