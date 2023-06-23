####################################################################################
### Script to generate the polygons and centroids file for Seurat LoadNanostring ###
####################################################################################

library(data.table)

# project <- "CosMx_Lung"
sample.name <- "Lung5_Rep1"
type.of.image <- "CellLabels"
# force=FALSE
### directory settings ###
  # main.path <- file.path("/Users","silviobicciato","Documents","GeneChip_Data",project)
  main.path <- file.path("~/Downloads","CosMx_data")
  data.dir <- file.path(main.path,sample.name, "Lung5_Rep1-Flat_files_and_images")
  meta.data.all <- read.csv(file.path(data.dir,paste0(sample.name,"_metadata_file.csv")),header = T)
  fov.positions <- read.csv(file.path(data.dir,paste0(sample.name,"_fov_positions_file.csv")),header = T)
  source(file.path("inst/scripts/SpatVect_generate.R"))

###############################################################
### Generate the polygon and centroids .csv file for Seurat ###
###############################################################
  fov.positions <- fov.positions[fov.positions$fov%in%unique(meta.data.all$fov),]
  poly.mat <- NULL
  for (k in 1:dim(fov.positions)[1])
  {
    print(paste0("Processing FOV ", fov.positions[k,"fov"]))
    Spat.obj <- spatVectGenerate(main.path = data.dir,
                                 image.type = type.of.image,
                                 fov.number = fov.positions[k,"fov"],
                                 flip.vertical = TRUE)
        # SpatVect_generate(main.path = data.dir,
        #                           image.type = type.of.image,
        #                           fov.number = fov.positions[k,"fov"],
        #                           flip.vertical = TRUE) #force=force
    # Convert the SpatVector to a data table
    DT.geom <- as.data.table(terra::geom(Spat.obj))
    DT.values <- as.data.table(terra::values(Spat.obj))
    DT.values[, `:=`(geom, 1:nrow(DT.values))]
    spatVecDT <- merge.data.table(DT.geom, DT.values,by = "geom")
    spatVecDT <- spatVecDT[,c("cellID","x","y")]
    names(spatVecDT) <- c("cellID","x_global_px", "y_global_px")
    # calculate centroids
    centr.geom <- terra::centroids(Spat.obj)
    # Convert the SpatVector of centroids to a data table
    CT.geom <- as.data.table(terra::geom(centr.geom))
    CT.values <- as.data.table(terra::values(centr.geom))
    CT.values[, `:=`(geom, 1:nrow(CT.values))]
    spatVecCT <- merge.data.table(CT.geom, CT.values,by = "geom")
    spatVecCT <- spatVecCT[,c("cellID","x","y")]
    names(spatVecCT) <- c("cellID","sdimx", "sdimy")
    # merge polygons and centroids
    sample.poly.mat <- merge(spatVecDT,spatVecCT, by = "cellID")
    # shift coordinates base on the underlying tissue geometry
    sample.poly.mat[,c("x_global_px","sdimx")] <- sample.poly.mat[,c("x_global_px","sdimx")] + fov.positions[k,"x_global_px"]
    sample.poly.mat[,c("y_global_px","sdimy")] <- sample.poly.mat[,c("y_global_px","sdimy")] + fov.positions[k,"y_global_px"]
    # add a column with the fov and transform into a data.table
    sample.poly.mat <- data.frame(fov = fov.positions[k,"fov"], as.data.frame(sample.poly.mat))
    poly.mat <- rbind(poly.mat, sample.poly.mat)
  }
  fwrite(poly.mat, file.path(data.dir,paste0(sample.name,"-polygons.csv")))


