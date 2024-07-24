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
spe <- readCosmxSPE(cosmxFolder)


#############################

spe <- readCosmxSPE(cosmxFolder)
spe <- spatialPerCellQC(spe)
spe <- computeQCScore(spe)#, threshold=0.01)
# spe[, spe$flag_score>quantile(spe$flag_score, probs=0.15)]
# spe <- computeSpatialOutlier(spe)
colData(spe)
spe <- computeFilterFlags(spe)
spe[,!spe$is_fscore_outlier]
# spe[, spe$to_keep]
# (spefiltered <- spe[,!spe$fscore_outlier])

