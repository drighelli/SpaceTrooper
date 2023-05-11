# library(remotes)
# remotes::install_github("satijalab/seurat", "feat/imaging", quiet = TRUE)
# library(future)
# plan("multisession", workers = 10)
# trace("SingleImagePlot", edit=TRUE)

library(Seurat)
library(ggplot2)
library(cowplot)
library(openxlsx)

project <- "S2-S3"
sample.name <- "R5363_S3_188-B"
field.of.view <-"S3"

### settings ###
  main.path <- file.path("/Users","silviobicciato","Documents","GeneChip_Data",
                         "AIRC_Piccolo","Nanostring_CosMx",project)
  main.data.dir <- file.path(main.path,"data")
  data.dir <- file.path(main.data.dir,sample.name)
  main.result.dir<-file.path(main.path,"results")
  if (!file.exists(main.result.dir)){dir.create(main.result.dir)}
  main.result.dir<-file.path(main.path,"results","QC and Cell phenotyping")
  if (!file.exists(main.result.dir)){dir.create(main.result.dir)}
  result.dir <- file.path(main.result.dir,sample.name)
  if (!file.exists(result.dir)){dir.create(result.dir)}
  meta.data.all <- read.csv(file.path(data.dir,paste0(sample.name,"_metadata_file.csv")),header = T)
  fov.positions <- read.csv(file.path(data.dir,paste0(sample.name,"_fov_positions_file.csv")),header = T)
  ### define my image theme ###
  my_image_theme <- function(back.color="black",back.border=NA,title.col="white") {
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.title = element_text(color = title.col,hjust = 0.5,face = "bold"),
          plot.background = element_rect(fill = "transparent",colour = back.border))
  }

####################################
### Data loading and feature plots 
####################################
  CosMx.obj <- LoadNanostring(data.dir = data.dir, fov = field.of.view)
  Idents(CosMx.obj) <- sample.name
  drop.out <- dim(meta.data.all)[1]-dim(CosMx.obj)[2]
  print(paste0(drop.out, " cells have been removed because have no expressed genes"))
  rownames(meta.data.all) <- paste(meta.data.all$cell_ID,meta.data.all$fov,sep="_")
  # filter the original metadata table
  meta.data.all <- meta.data.all[colnames(CosMx.obj),]
  # populate the meta.data slot of the Seurat object
  CosMx.obj$orig.ident <- sample.name
  CosMx.obj$image.id <- rep(Images(CosMx.obj),dim(CosMx.obj)[2])
  CosMx.obj$fov <- meta.data.all$fov
  CosMx.obj$cell.id <- meta.data.all$cell_ID
  CosMx.obj$x <- GetTissueCoordinates(CosMx.obj)$x
  CosMx.obj$y <- GetTissueCoordinates(CosMx.obj)$y
  mean.antibody <- meta.data.all[,grep("Mean.",colnames(meta.data.all))]
  colnames(mean.antibody) <- c("Membrane", "PanCK", "CD45", "CD3", "DAPI")
  CosMx.obj@meta.data <- data.frame(CosMx.obj@meta.data,mean.antibody)
  # change the HLA- gene symbols into HLA.
  rownames(CosMx.obj@assays$Nanostring@counts) <-gsub("HLA-", "HLA.", rownames(CosMx.obj))
  rownames(CosMx.obj@assays$Nanostring@data) <-gsub("HLA-", "HLA.", rownames(CosMx.obj))
  # plot tissue slide with FOVs #
  fov.positions <- fov.positions[fov.positions$fov%in%unique(meta.data.all$fov),]
  fov.xdim <- 5472
  fov.ydim <- 3648
  pdf(file.path(result.dir,paste0("01_",sample.name,"_FOVs.pdf")), useDingbats=FALSE)
  print(ggplot()+geom_point(data=meta.data.all,mapping = aes(x=CenterX_global_px,y=CenterY_global_px),
                            colour="darkorange2",
                            fill="darkorange2",
                            size = 0.05, alpha = 0.2) +
          annotate("rect", xmin = fov.positions$x_global_px,
                   xmax = fov.positions$x_global_px+fov.xdim, 
                   ymin = fov.positions$y_global_px, 
                   ymax = fov.positions$y_global_px+fov.ydim,
                   alpha = .2, color = "black",linewidth = 0.2) +
          geom_text(aes(x=fov.positions$x_global_px+fov.xdim/2,
                        y=fov.positions$y_global_px+fov.ydim/2,
                        label = fov.positions$fov), color="black", fontface = "bold") +
          ggtitle(sample.name) +
          my_image_theme(back.color="white",back.border="white",title.col="black"))
  dev.off()
  # plot markers for cell detection #
  cell.markers <- c("PanCK","CD45","CD3")
  upper.lim <- max(matrixStats::colQuantiles(as.matrix(meta.data.all[,paste0("Mean.",cell.markers)]),probs=0.75))
  pdf(file.path(result.dir, paste0("02_", sample.name, "_CellMarkersIF.pdf")), 
      width=8, height=8, useDingbats=FALSE)
  for (i in 1:length(cell.markers)) {
    print(ggplot(data=meta.data.all, aes(CenterX_global_px,y=CenterY_global_px))+
      geom_point(aes(color=meta.data.all[,paste0("Mean.",cell.markers[i])]),
                size = 0.1) +
      scale_color_gradient(low = "white", high="orange",na.value = "darkorange2",
                           limits=c(0,upper.lim))+
      annotate("rect", xmin = fov.positions$x_global_px,
               xmax = fov.positions$x_global_px+fov.xdim, 
               ymin = fov.positions$y_global_px, 
               ymax = fov.positions$y_global_px+fov.ydim,
               alpha = 0, color = "black",linewidth = 0.1) +
      ggtitle(cell.markers[i])+
      NoLegend() +
      my_image_theme(back.color="white",back.border="black",title.col="black"))
  }
 
  #print(plot_grid(plotlist=plot.list, ncol=2, nrow=2))
  invisible(dev.off())

################################################################################
### QC metrics
################################################################################
  pdf(file.path(result.dir, paste0("03_", sample.name, "_QCmetrics.pdf")), width=14, height=14, useDingbats=FALSE)
  # VlnPlot of nCount and nFeature in raw and normalized
  plot.list <- NULL
  plot.list[[1]] <- VlnPlot(CosMx.obj, features = "nCount_Nanostring",cols = 4, pt.size=0) +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 0,hjust = 0.5)) +ggplot2::xlab("") +
    geom_boxplot(width = 0.1, outlier.shape = NA,fill = "yellow",alpha = 0.3) +
    ylab("Counts per cells") +
    ggtitle("nCount (raw)")
  plot.list[[2]] <- VlnPlot(CosMx.obj, features = "nFeature_Nanostring",cols = 2, pt.size=0) +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 0,hjust = 0.5)) +ggplot2::xlab("") +
    geom_boxplot(width = 0.1, outlier.shape = NA,fill = "yellow",alpha = 0.3)+
    ylab("Features per cells") +
    ggtitle("nFeature (raw)")
  print(plot_grid(plotlist=plot.list, ncol=2, nrow=1))
  # Distributions of nCount and nFeature in raw
  plot.list <- NULL
  plot.list[[1]] <- ggplot(CosMx.obj[[]], aes(x = nCount_Nanostring)) + 
    geom_histogram(binwidth=10, color = 4,fill = 4, alpha = 0.25) +
    xlab("Counts per cell") +
    ylab("Number of cells") +
    ggtitle("nCounts (raw)")
  plot.list[[2]] <- ggplot(CosMx.obj[[]], aes(x = nFeature_Nanostring)) + 
    geom_histogram(binwidth=5, color = 2,fill = 2, alpha = 0.25) +
    xlab("Features per cell") +
    ylab("Number of cells") +
    ggtitle("nFeatures (raw)")
  print(plot_grid(plotlist=plot.list, ncol=2, nrow=1))
  print(FeatureScatter(CosMx.obj, feature1 = "nFeature_Nanostring", 
                                   feature2 = "nCount_Nanostring", pt.size = 0.5, cols = 4, raster=FALSE) +
    ggtitle("Features vs Counts (raw)") + NoLegend())
  # Position of cells with low nCount and nFeatures ≤ of a given quantile
  temp<-CosMx.obj[[]]
  thresholds <- c(0.01,0.05,0.1,0.25)
  plot.list <- NULL
  for (i in 1:length(thresholds)){
    thr.count <- as.numeric(quantile(temp$nCount_Nanostring,thresholds[i]))
    thr.feat <- as.numeric(quantile(temp$nFeature_Nanostring,thresholds[i]))
    high.label <- paste0("nCounts > ",thr.count,
                         "; nFeat > ",thr.feat)
    low.label <- paste0("nCounts < ",thr.count,
                        "; nFeat < ",thr.feat)
    temp$filter <- rep(low.label,dim(temp)[1])
    temp[temp$nCount_Nanostring>thr.count&temp$nFeature_Nanostring>thr.feat,"filter"] <- high.label
    plot.list[[i]] <- ggplot(data=temp,aes(x=x,y=y,
                                           colour = factor(filter),
                                           size= factor(filter)))+
      geom_point()+
      scale_size_manual(values=c(0.3,0.05), guide = "none")+
      scale_colour_manual(values=c("black","grey90"))+
      xlab("")+
      ylab("")+
      theme(legend.position="top")+
      theme(legend.title = element_blank())+
      ggtitle(paste0(dim(temp[temp$filter==low.label,])[1], 
                     " cells with nCount and nFeature lower than the ",
                     thresholds[i]*100,"th quantile"))+
      my_image_theme(back.color="white",back.border="white",title.col="black")
  }
  print(plot_grid(plotlist=plot.list, ncol=2, nrow=2))
  # Position of cells with small area
  # temp<-CosMx.obj[[]]
  # thresholds <- c(0.01,0.05,0.1,0.25)
  # plot.list <- NULL
  # for (i in 1:length(thresholds)){
  #   area.thr <- as.numeric(quantile(meta.data.all$Area,thresholds[i]))
  #   high.label <- paste0("cell Area > ",area.thr)
  #   low.label <- paste0("cell Area < ",area.thr)
  #   temp$filter <- rep(low.label,dim(temp)[1])
  #   temp[meta.data.all$Area>area.thr,"filter"] <- high.label
  #   plot.list[[i]] <- ggplot(data=temp,aes(x=x,y=y,
  #                                          colour = factor(filter),
  #                                          size= factor(filter)))+
  #     geom_point()+
  #     scale_size_manual(values=c(0.3,0.05), guide = "none")+
  #     scale_colour_manual(values=c("black","grey90"))+
  #     xlab("")+
  #     ylab("")+
  #     theme(legend.position="top")+
  #     theme(legend.title = element_blank())+
  #     ggtitle(paste0(dim(temp[temp$filter==low.label,])[1], 
  #                    " cells with Area lower than the ",
  #                    thresholds[i]*100,"th quantile"))+
  #     my_image_theme(back.color="white",back.border="white",title.col="black")
  # }
  # print(plot_grid(plotlist=plot.list, ncol=2, nrow=2))
  dev.off()
  
################################################################################
### Filtering cell with nCounts and nFeatures lower than a given quantile
################################################################################
  qtl.thr <- 0.05
  nCountLow <- as.numeric(quantile(CosMx.obj$nCount_Nanostring,qtl.thr))
  nFeatLow <- as.numeric(quantile(CosMx.obj$nFeature_Nanostring,qtl.thr))
  CosMx.obj.filt <- subset(CosMx.obj, 
                      subset = nCount_Nanostring > nCountLow & nFeature_Nanostring > nFeatLow)
  
######################################################
### Normalization and scaling
######################################################
  # Normalization and scaling using SCTransform #
  # The default assay is now "SCT" with counts being (corrected) counts,
  # data being log1p(counts), scale.data being Pearson residuals #
  CosMx.obj.norm <- SCTransform(CosMx.obj.filt, assay = "Nanostring",
                                clip.range = c(-10, 10), verbose = T)
  # identification and removal of negative probes
  all.probes <- rownames(CosMx.obj.norm)
  neg.probes <- all.probes[grep("Neg",all.probes)]
  not.neg.probes <- all.probes[-grep("Neg",all.probes)]
  CosMx.obj.norm.neg <- subset(CosMx.obj.norm, features = neg.probes)
  CosMx.obj <- subset(CosMx.obj, features = not.neg.probes)
  CosMx.obj.norm <- subset(CosMx.obj.norm, features = not.neg.probes)
  # plot of nCounts and nFeatures after normalization #
  pdf(file.path(result.dir,paste0("04_",sample.name,"_Normalization.pdf")), 
      width = 14, height = 12, useDingbats=FALSE)
  p1 <- ggplot(CosMx.obj[[]], aes(x = nCount_Nanostring)) + 
    geom_histogram(binwidth=10, color = 4,fill = 4, alpha = 0.25) +
    xlab("Total counts per cell before normalization") +
    ylab("Number of cells")
  p2 <- ggplot(CosMx.obj[[]], aes(x = nFeature_Nanostring)) + 
    geom_histogram(binwidth=5, color = 2,fill = 2, alpha = 0.25) +
    xlab("Total feature per cell before normalization") +
    ylab("Number of cells")
  p3 <- ggplot(CosMx.obj.norm[[]], aes(x = nCount_SCT)) + 
    geom_histogram(binwidth=1, color = 4,fill = 4, alpha = 0.25) +
    xlab("Total counts per cell after normalization") +
    ylab("Number of cells")
  p4 <- ggplot(CosMx.obj.norm[[]], aes(x = nFeature_SCT)) + 
    geom_histogram(binwidth=2, color = 2,fill = 2, alpha = 0.25) +
    xlab("Total features per cell after normalization") +
    ylab("Number of cells")
  plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
  dev.off()

######################################################################
### Write expression levels, meta.data and filtering information
######################################################################
  # log1p data file with coordinates 
  wb<-createWorkbook(title="results")
  bold.style <- createStyle(textDecoration = "Bold")
  out.data <- data.frame(cell_id = colnames(CosMx.obj.norm),
                         fov = CosMx.obj.norm$fov,
                         GetTissueCoordinates(CosMx.obj.norm)[,1:2],
                         round(t(as.matrix(GetAssayData(CosMx.obj.norm))),2))
  addWorksheet(wb,"Expression matrix")
  modifyBaseFont(wb, fontSize = 14, fontName = "Arial")
  writeData(wb, sheet = "Expression matrix", out.data ,
            rowNames = F, headerStyle = bold.style)
  saveWorkbook(wb, file.path(result.dir,paste0(sample.name,"_expression_data.xlsx")), overwrite = TRUE)    
  # original metadata 
  wb<-createWorkbook(title="results")
  bold.style <- createStyle(textDecoration = "Bold")
  cellIDs.all <- paste(meta.data.all$cell_ID,meta.data.all$fov,sep="_")
  rownames(meta.data.all) <- cellIDs.all
  meta.data.out <- data.frame(cell_id = rownames(meta.data.all),
                              meta.data.all)
  colnames(meta.data.out)[3] <- "cell_id.original"
  addWorksheet(wb,"meta data")
  modifyBaseFont(wb, fontSize = 14, fontName = "Arial")
  writeData(wb, sheet = "meta data", meta.data.out,
            rowNames = F, headerStyle = bold.style)
  # filtering information
  filter.column <- rep("NO",dim(meta.data.all)[1])
  names(filter.column)<-rownames(meta.data.all)
  filtered.cells <- setdiff(rownames(meta.data.all),colnames(CosMx.obj.norm))
  filter.column[filtered.cells] <- "YES"
  filter.data <- data.frame(cell_id = colnames(CosMx.obj),
                            fov = CosMx.obj$fov,
                            GetTissueCoordinates(CosMx.obj)[,1:2],
                            CosMx.obj@meta.data[,grep("Nanostring", colnames(CosMx.obj[[]]))],
                            filtered = filter.column)
  addWorksheet(wb,"filter infos")
  modifyBaseFont(wb, fontSize = 14, fontName = "Arial")
  writeData(wb, sheet = "filter infos", filter.data,
            rowNames = F, headerStyle = bold.style)
  saveWorkbook(wb, file.path(data.dir,paste0(sample.name,"_meta_data.xlsx")), overwrite = TRUE)    
  
################################################################################
### Save Seurat object
################################################################################
  saveRDS(CosMx.obj.norm, file = file.path(data.dir,paste0(sample.name,"_data.rds")))
  saveRDS(CosMx.obj.norm[[]], file = file.path(data.dir,paste0(sample.name,"_metadata.rds")))
  
