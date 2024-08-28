cosmxpolfile <- "~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/Run5810_Case2-polygons.csv"
xenium_polfile <- "~/Downloads/Xenium_data/pancreas/outs/cell_boundaries.parquet"
merscopepolfile <- "~/Downloads/Merfish_data/human_uterine_cancer_patient2/"

(polcosm <- readPolygonsCosmx(cosmxpolfile))
sum(polcosm$is_multi)
sum(polcosm$multi_n) == dim(polcosm)[1]
(polxen <- readPolygonsXenium(xenium_polfile))
sum(polxen$is_multi)
sum(polxen$multi_n) == dim(polxen)[1]
polxen[polxen$is_multi,]
(polmersc <- readPolygonsMerfish(merscopepolfile, type="parquet"))
sum(polxen$is_multi)
sum(polxen$multi_n) == dim(polxen)[1]
polxen[polxen$is_multi,]


xeniumFolder <- "~/Downloads/Xenium_data/pancreas/"
spe <- readXeniumSPE(xeniumFolder)

merscopeFolder <- "~/Downloads/Merfish_data/human_uterine_cancer_patient2/"
spe <- readMerfishSPE(merscopeFolder, boundaries_type="parquet")

cosmxFolder <- "~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/"
spe <- readCosmxSPE(cosmxFolder, loadPolygons=TRUE)


#############################
library(devtools)
load_all()
cosmxFolder <- "~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/"
spe <- readCosmxSPE(cosmxFolder, sample_name="DBKero_CosMx", loadPolygons=FALSE)
spe <- spatialPerCellQC(spe)
colData(spe)
spe <- computeQCScore(spe)
colData(spe)
spe <- computeSpatialOutlier(spe, compute_by="Area_um", method="both")
colData(spe)
spe <- computeSpatialOutlier(spe, compute_by="Mean.DAPI", method="both")
colData(spe)
spe <- computeFilterFlags(spe)
colData(spe)
plotCentroidsSpe(spe)
plotCentroidsSpe(spe, colour_by="total")
plotCentroidsSpe(spe, colour_by="is_fscore_outlier")

plotMetricHist(spe, metric="Area_um")
plotMetricHist(spe, metric="Area_um", use_fences="Area_um_outlier_sc")

pols <- readPolygonsCosmx(metadata(spe)$polygons)
spe <- addPolygonsToSPE(spe, pols)
plotCellsFovs(spe)
spe1 <- spe[,spe$fov %in% c(16:19)]
colData(spe1)
plotPolygonsSPE(spe1)
plotPolygonsSPEold(spe1)
plotPolygonsSPE(spe1, colour_by="Area_um")
plotPolygonsSPEold(spe1, color_by="Area_um")
plotPolygonsSPE(spe1, colour_by="Area_um_outlier_mc")
plotPolygonsSPEold(spe1, color_by="Area_um_outlier_mc")
plotPolygonsSPE(spe1, colour_by="is_fscore_outlier")
plotPolygonsSPEold(spe1, color_by="is_fscore_outlier")
labels <- data.table::fread("/Users/inzirio/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/cosmx_dbkero_IST_labels_complete_simple.txt")
id <- strsplit(labels$orig.ident, "_")
labels$fov <- lapply(id, function(i) return(i[1]))
labels$cellID <- lapply(id, function(i) return(i[2]))
labels$cell_id <- paste0("f",labels$fov, "_c", labels$cellID)
spe$labels <- NA
spe$labels[match(labels$cell_id, spe$cell_id)] <- labels$InSituType_simple

colData(spe)

spe$polygons$outlier_fence_color <- dplyr::case_when(spe$is_zero_counts == TRUE ~ "darkturquoise",
                                              spe$is_ctrl_tot_outlier == TRUE ~ "magenta",
                                              spe$Mean.DAPI > round(metadata(spe)$dapi_fence[2,1], 2) ~ "greenyellow",
                                              spe$Mean.DAPI < round(metadata(spe)$dapi_fence[1,1], 2) ~ "purple",
                                              spe$Area_um > round(metadata(spe)$area_fence[2,1], 2) ~ "red",
                                              spe$Area_um < round(metadata(spe)$area_fence[1,1], 2) ~ "white",
                                              TRUE ~ "black")




