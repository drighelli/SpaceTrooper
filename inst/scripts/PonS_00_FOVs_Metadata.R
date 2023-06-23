##################################################
### Script to plot metadata into a whole image ###
##################################################

library(sf)
library(tmap)
library(RColorBrewer)

project <- "CosMx_Lung"
sample.name <- "Lung5_Rep1"
type.of.image <- "CellLabels"

### directory settings ###
  main.path <- file.path("/Users","silviobicciato","Documents","GeneChip_Data",project)
  data.dir <- file.path(main.path,"data",sample.name)
  global.image.dir <- file.path(data.dir,"spatVect","global")
  result.dir <- file.path(main.path,"results")
  if (!file.exists(result.dir)){dir.create(result.dir)}
  result.dir <- file.path(result.dir,"QC and Cell phenotyping")
  if (!file.exists(result.dir)){dir.create(result.dir)}
  result.dir <- file.path(result.dir,sample.name)
  if (!file.exists(result.dir)){dir.create(result.dir)}
  meta.data.all <- read.csv(file.path(data.dir,paste0(sample.name,"_metadata_file.csv")),header = T)
  fov.positions <- read.csv(file.path(data.dir,paste0(sample.name,"_fov_positions_file.csv")),header = T)
  fov.positions <- fov.positions[fov.positions$fov%in%unique(meta.data.all$fov),]
  
####################################################
### Load entire image and plot FOVs and metadata  
###################################################
  # load entire image
  Spat.obj <- read_sf(global.image.dir)
# create the FOVs borders
  fov.xdim <- 5472
  fov.ydim <- 3648
  sfc = st_sfc(lapply(data.frame(t(fov.positions[,2:3])),st_point))
  # create the bounding boxes for each FOV
  rect_around_point <- function(x,xsize,ysize){
    bbox <- st_bbox(x)
    bbox <- bbox + c(xsize,ysize,-0,0)
    y = st_as_sfc(bbox)
    return(y)
  }
  FOV.bbox <- lapply(sfc, rect_around_point, xsize = fov.xdim, ysize = fov.ydim)
  FOV.list = do.call(c, FOV.bbox)
  FOV.geolist <- st_as_sf(FOV.list)
  st_crs(FOV.geolist)<-st_crs(Spat.obj)
  FOV.geolist <- cbind(FOV.geolist,1:dim(FOV.geolist)[1])
  names(FOV.geolist)[1] <- "FOVs"
  # plot FOVs borders and numbers on entire image
  tm.plot <- tm_shape(Spat.obj) + 
    tm_fill(col="grey90") +
    tm_borders(lwd = 0.1, alpha=0.5, col = "grey50") +
    tm_shape(FOV.geolist) + 
    tm_borders(col="black", lwd = 2) +
    tm_text("FOVs", fontface="bold", size =1.3)
  tmap_save(tm.plot, file.path(result.dir,paste0("01_",sample.name,"_FOVs.pdf")),
            outer.margins = F)
# add and plot fluorescence markers
  cell.markers <- c("PanCK","CD45","CD3","DAPI")
  cols <- c("Greens","Oranges","Reds", "Blues")
  meta.data.all$cellID <- paste0(meta.data.all$cell_ID,"_",meta.data.all$fov)
  colnames(meta.data.all)[2] <- "fov_cellID"
  meta.to.add <- meta.data.all[,c("cellID",paste0("Mean.",cell.markers))]
  Spat.obj <- merge(Spat.obj,meta.to.add,by="cellID")
  for (k in 1: length(cell.markers)){
    Spat.obj.filt <- Spat.obj %>%
      dplyr::filter(Spat.obj[[paste0("Mean.",cell.markers[k])]]>=mean(Spat.obj[[paste0("Mean.",cell.markers[k])]]))
    tm.plot <- tm_shape(Spat.obj.filt) + 
      tm_fill(col = paste0("Mean.",cell.markers[k]),
              style = "pretty",
              n = 7,
              palette = brewer.pal(7, cols[k]),
              legend.hist = FALSE) +
      tm_borders(lwd = 0.1, alpha=0.5, col = "grey90") +
      tm_layout(frame = FALSE,
                bg.color="black",
                legend.text.color= "white", 
                legend.outside = F)
    tmap_save(tm.plot,
              file.path(result.dir,paste0("02_",sample.name,"_",cell.markers[k],".pdf")),
              outer.margins = F)
  }
  
  
