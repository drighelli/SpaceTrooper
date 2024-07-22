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
