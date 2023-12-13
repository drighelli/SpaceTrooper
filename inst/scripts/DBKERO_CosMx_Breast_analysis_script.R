#### DBKERO COSMX BREAST DATASET PREPROCESSING ####
### SPATIAL TRANSCRIPTOMICS ###

library(SpatialExperiment)
library(data.table)
library(scater) # it imports scuttle
library(cowplot)
library(ggplot2)
theme_set(theme_bw())
library(matrixStats)
library(dplyr)
library(tidyr)
library(tibble)
library(arrow)
library(scales)
library(sf)
library(tmap)

# read-in function

readCosmxSPE <- function(dirname = dirname, 
                         countmatfpattern = "exprMat_file.csv", 
                         metadatafpattern = "metadata_file.csv", 
                         coord_names = c("CenterX_global_px",
                                         "CenterY_global_px")){
  
  countmat_file <- file.path(dirname, list.files(dirname, countmatfpattern))
  metadata_file <- file.path(dirname, list.files(dirname, metadatafpattern))
  
  # Read in 
  countmat <- data.table::fread(countmat_file)
  metadata <- data.table::fread(metadata_file)
  
  # Count matrix   
  counts <- merge(countmat, metadata[, c("fov", "cell_ID")])
  counts <- subset(counts, select = -c(fov, cell_ID))
  cell_codes <- rownames(counts)
  features <- colnames(counts)
  counts <- t(counts)
  
  rownames(counts) <- features
  colnames(counts) <- cell_codes
  
  # rowData (does not exist)
  
  # colData
  colData <- merge(metadata, countmat[, c("fov", "cell_ID")])
  
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    # rowData = rowData,
    colData = colData,
    spatialCoordsNames = coord_names
  )
  
  return(spe)
}

dirname <- "/NAS06/work/Spatial_omics/DBKero/CosMx_Breast/CosMX_data_Case2"
countmatfpattern <- "Run5810_Case2_exprMat_file.csv"
metadatafpattern <- "Run5810_Case2_metadata_file.csv"
cell_dir <- getwd()

spe <- readCosmxSPE(dirname = dirname, countmatfpattern = countmatfpattern, metadatafpattern = metadatafpattern)

# Preprocessing: QC analysis

spe <- scuttle::addPerCellQC(spe, subsets=list(NegProb = grep("^NegPrb", rownames(spe))))
names(colData(spe)@listData)

#[1] "fov"                      "cell_ID"                 
#[3] "Area"                     "AspectRatio"             
#[5] "CenterX_local_px"         "CenterY_local_px"        
#[7] "Width"                    "Height"                  
#[9] "Mean.PanCK"               "Max.PanCK"               
#[11] "Mean.CD68"                "Max.CD68"                
#[13] "Mean.MembraneStain_B2M"   "Max.MembraneStain_B2M"   
#[15] "Mean.CD45"                "Max.CD45"                
#[17] "Mean.DAPI"                "Max.DAPI"                
#[19] "sample_id"                "sum"                     
#[21] "detected"                 "subsets_NegProb_sum"     
#[23] "subsets_NegProb_detected" "subsets_NegProb_percent" 
#[25] "total"  

dim(spe)

# selecting only numeric and integer variables from metadata
temp <- unlist(lapply(spe@colData@listData, function(x) is.numeric(x) | is.integer(x)))
num_variables <- c()
for (i in 1:length(temp)){
  num_variables[i] <- ifelse(temp[i] == TRUE, names(temp[i]), NA)
  num_variables <- num_variables[!is.na(num_variables)]
}

num_variables

# selecting only some variables of interest
variables <- c("total", "detected", "Mean.PanCK", "Mean.CD68", "Mean.MembraneStain_B2M", "Mean.CD45", "Mean.DAPI", "subsets_NegProb_percent")

# Large-scale to small-scale
# numeric/integer dependence on global position

# importing global coordinates since they are stored only in spatialCoords slot

colData(spe)$CenterX_global_px <- spatialCoords(spe)[,1]
colData(spe)$CenterY_global_px <- spatialCoords(spe)[,2]

# variables ~ global coordinates scatterplot, could highlight spatial patterns along continuous global coordinates
# total = sum = nCounts_RNA in Seurat
# detected = nFeatures_RNA in Seurat

cell_dir <- file.path(getwd(), "cell_plots")
if(!dir.exists(cell_dir)) dir.create(cell_dir)

# transcriptomics-specific variables ~ global coordinates

png(file.path(cell_dir,"total-detected~global_coord_scatterplots.png"), units = "cm", height = 20, width = 30, res = 180)
p1 <- plotColData(spe, "total", "CenterX_global_px")
p2 <- plotColData(spe, "total", "CenterY_global_px")
p3 <- plotColData(spe, "detected", "CenterX_global_px")
p4 <- plotColData(spe, "detected", "CenterX_global_px")
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

# marker signal detection variables ~ global coordinates

png(file.path(cell_dir,"Mean.CD68-subsets_NegProb_percent~global_coord_scatterplots.png"), units = "cm", height = 20, width = 30, res = 180)
p1 <- plotColData(spe, "Mean.CD68", "CenterX_global_px")
p2 <- plotColData(spe, "Mean.CD68", "CenterY_global_px")
p3 <- plotColData(spe, "subsets_NegProb_percent", "CenterX_global_px")
p4 <- plotColData(spe, "subsets_NegProb_percent", "CenterX_global_px")
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

png(file.path(cell_dir,"PanCK-memb~global_coord_scatterplots.png"), units = "cm", height = 20, width = 30, res = 180)
p1 <- plotColData(spe, "Mean.PanCK", "CenterX_global_px")
p2 <- plotColData(spe, "Mean.PanCK", "CenterY_global_px")
p3 <- plotColData(spe, "Mean.MembraneStain_B2M", "CenterX_global_px")
p4 <- plotColData(spe, "Mean.MembraneStain_B2M", "CenterX_global_px")
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

png(file.path(cell_dir,"CD45-DAPI~global_coord_scatterplots.png"), units = "cm", height = 20, width = 30, res = 180)
p1 <- plotColData(spe, "Mean.CD45", "CenterX_global_px")
p2 <- plotColData(spe, "Mean.CD45", "CenterY_global_px")
p3 <- plotColData(spe, "Mean.DAPI", "CenterX_global_px")
p4 <- plotColData(spe, "Mean.DAPI", "CenterX_global_px")
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

# cell geometry variables ~ global coordinates

png(file.path(cell_dir,"Area-AspectRatio~global_coord_scatterplots.png"), units = "cm", height = 20, width = 30, res = 180)
p1 <- plotColData(spe, "Area", "CenterX_global_px")
p2 <- plotColData(spe, "Area", "CenterY_global_px")
p3 <- plotColData(spe, "AspectRatio", "CenterX_global_px")
p4 <- plotColData(spe, "AspectRatio", "CenterX_global_px")
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

png(file.path(cell_dir,"Height-width~global_coord_scatterplots.png"), units = "cm", height = 20, width = 30, res = 180)
p1 <- plotColData(spe, "Height", "CenterX_global_px")
p2 <- plotColData(spe, "Height", "CenterY_global_px")
p3 <- plotColData(spe, "Width", "CenterX_global_px")
p4 <- plotColData(spe, "Width", "CenterX_global_px")
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()


# binning global coordinates to perform boxplots of numeric/integer variables along global coordinate bins

spe$X_glo_bin <- cut(spe$CenterX_global_px, breaks=quantile(spe$CenterX_global_px, seq(0, 1, by=0.05)))
spe$Y_glo_bin <- cut(spe$CenterY_global_px, breaks=quantile(spe$CenterY_global_px, seq(0, 1, by=0.05)))

# variable distribution by global coordinates bins, spatial patterns are better visible on scatterplot

pdf(file.path(cell_dir,"num_variable~global_bins_boxplots.pdf"),height = 8, width = 6, useDingbats = T)
for (var in 1:length(variables)){
  n <- grep(paste0("^", variables[var], "$"), colnames(spe@colData))
  par(mfrow = c(2,1))
  boxplot(spe@colData@listData[[n]] ~ spe$X_glo_bin, ylab = colnames(spe@colData)[n])
  boxplot(spe@colData@listData[[n]] ~ spe$Y_glo_bin, ylab = colnames(spe@colData)[n])
}
dev.off()

# global coordinate scatterplot colored by variables showing spatial pattern in boxplots

spe$logtot <- log10(spe$total)
spe$logdet <- log10(spe$detected)

png(file.path(cell_dir,"total-detected_tissue_slide_scatterplots.png"), units = "cm", height = 20, width = 30, res = 180)
p1 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "total", point_size = 0.1)
p2 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "logtot", point_size = 0.1)
p3 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "detected", point_size = 0.1)
p4 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "logdet", point_size = 0.1)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

png(file.path(cell_dir,"marker_signals_tissue_slide_scatterplots.png"), units = "cm", height = 30, width = 30, res = 180)
p1 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "Mean.PanCK", point_size = 0.1)
p2 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "Mean.MembraneStain_B2M", point_size = 0.1)
p3 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "Mean.CD68", point_size = 0.1)
p4 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "Mean.CD45", point_size = 0.1)
p5 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "Mean.DAPI", point_size = 0.1)
plot_grid(p1, p2, p3, p4, p5, ncol = 2)
dev.off()

png(file.path(cell_dir,"area-aspect_ratio_tissue_slide_scatterplots.png"), units = "cm", height = 10, width = 30, res = 180)
p1 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "Area", point_size = 0.1)
p2 <- plotColData(spe, "CenterX_global_px", "CenterY_global_px", colour_by = "AspectRatio", point_size = 0.1)
plot_grid(p1, p2, ncol = 2)
dev.off()

# check dependence on local coordinates

# variables ~ local coordinates scatterplot

spe$fov <- as.factor(spe$fov)

png(file.path(cell_dir,"total-detected~local_coord_scatterplots.png"), units = "cm", height = 20, width = 15, res = 180) 
p1 <- plotColData(spe, "total", "CenterX_local_px")
p2 <- plotColData(spe, "total", "CenterY_local_px")
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

# binning local coordinates to obtain a discrete variable as classes to a boxplot

spe$X_loc_bin <- cut(spe$CenterX_local_px, breaks=quantile(spe$CenterX_local_px, seq(0, 1, by=0.05)))
spe$Y_loc_bin <- cut(spe$CenterY_local_px, breaks=quantile(spe$CenterY_local_px, seq(0, 1, by=0.05)))

# boxplot of each of the numeric/integer variable by local coordinate bin

pdf(file.path(cell_dir,"num_variable~local_bins_exploratory_plots.pdf"), height = 8, width = 6, useDingbats = T)
for (var in 1:length(variables)){
  n <- grep(paste0("^", variables[var], "$"), colnames(spe@colData))
  par(mfrow = c(2,1))
  boxplot(spe@colData@listData[[n]] ~ spe$X_loc_bin, ylab = colnames(spe@colData)[n])
  boxplot(spe@colData@listData[[n]] ~ spe$Y_loc_bin, ylab = colnames(spe@colData)[n])
}
dev.off()

# check dependence on fov

fov_dir <- file.path(getwd(), "fov_plots")
if (!dir.exists(fov_dir)) dir.create(fov_dir)

# boxplot of each of the numeric/integer variable distribution by fov, useful to see spatial patterns in data

pdf(file.path(fov_dir, "num_variables~FOV_boxplots.pdf"), height = 12, width = 20, useDingbats = F)
for (var in 1:length(variables)){
  n <- grep(paste0("^", variables[var], "$"), colnames(spe@colData))
  boxplot(spe@colData@listData[[n]] ~ spe$fov, ylab = colnames(spe@colData)[n], cex = 0.2)
  title(paste0(colnames(spe@colData)[n], "~fov"))
}
dev.off()

# Check relation between total/detected and area (stratified by space or FOV)

# Area~total/detected scatterplot faced by fov

cell_df <- as.data.frame(spe@colData@listData)

png(file.path(fov_dir,"Area~total_by_fov_exploratory_plots.png"), units = "cm", height = 40, width = 30, res = 180)
ggplot(data = cell_df, aes(x = Area, y=total)) +
  geom_point() +
  facet_wrap(~fov)
dev.off()

png(file.path(fov_dir,"Area~detected_by_fov_exploratory_plots.png"), units = "cm", height = 40, width = 30, res = 180)
ggplot(data=as.data.frame(colData(spe)), aes(x=Area, y=detected)) +
  geom_point() +
  facet_wrap(~fov)
dev.off()

# Area~total colored by local X-Y coordinates faceted by fov

png(file.path(fov_dir,"Area~total_localX_by_fov_exploratory_plots.png"), units = "cm", height = 40, width = 30, res = 180)
ggplot(data=as.data.frame(colData(spe)), aes(x=Area, y=total, color=CenterX_local_px)) +
  geom_point() +
  facet_wrap(~fov)
dev.off()

png(file.path(fov_dir,"Area~total_localY_by_fov_exploratory_plots.png"), units = "cm", height = 40, width = 30, res = 180)
ggplot(data=as.data.frame(colData(spe)), aes(x=Area, y=total, color=CenterY_local_px)) +
  geom_point() +
  facet_wrap(~fov)
dev.off()

## Goodness-of-fit

# Let's see what could be a good distribution for the data.

# linear model fitting on mean~variance scatterplot

means <- rowMeans(counts(spe))
vars <- rowVars(counts(spe))

png(file.path(cell_dir,"linear_model_fit.png"), units = "cm", height = 10, width = 15, res = 180)
ggplot(data.frame(mean=means, var=vars), aes(mean, var)) +
  geom_point() + 
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10") +
  geom_abline(slope=1, intercept=0, colour="red") +
  geom_smooth()
dev.off()

# linear model fitting on mean~variance scatterplot faced by fov

means_fov <- t(apply(counts(spe), 1, tapply, spe$fov, mean))
means_fov <- rownames_to_column(as.data.frame(means_fov), var="gene")
vars_fov <- t(apply(counts(spe), 1, tapply, spe$fov, var))
vars_fov <- rownames_to_column(as.data.frame(vars_fov), var="gene")

meanvar_fov <- pivot_longer(means_fov, cols=2:31, cols_vary = "slowest", names_to="fov", values_to="mean")
tmp_var <- pivot_longer(vars_fov, cols=2:31, cols_vary = "slowest", names_to="fov", values_to="var")

meanvar_fov$var <- tmp_var$var
meanvar_fov$fov <- as.factor(as.numeric(meanvar_fov$fov))


png(file.path(fov_dir,"by_fov_linear_fit.png"), units = "cm", height = 40, width = 30, res = 180)
ggplot(meanvar_fov, aes(mean, var)) +
  geom_point() + 
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10") +
  geom_abline(slope=1, intercept=0, colour="red") +
  geom_smooth() +
  facet_wrap(~fov)
dev.off()

# linear model fitting only of negative probe signal on mean~variance scatterplot

png(file.path(cell_dir,"subset_neg_probes_linear_model_fit.png"), units = "cm", height = 10, width = 15, res = 180)
ggplot(data.frame(mean=means[grep("^NegPrb", rownames(spe))],
                  var=vars[grep("^NegPrb", rownames(spe))]),
       aes(mean, var)) +
  geom_point() + 
  #   scale_x_continuous(trans="log10") +
  #   scale_y_continuous(trans="log10") +
  geom_abline(slope=1, intercept=0, colour="red") +
  geom_smooth()
dev.off()

# linear model fitting only of negative probe signal on mean~variance scatterplot faced by fov

png(file.path(fov_dir,"by_fov_subset_neg_probes_linear_fit.png"), units = "cm", height = 40, width = 30, res = 180)
ggplot(meanvar_fov[grep("^NegPrb", meanvar_fov$gene),], aes(mean, var)) +
  geom_point() + 
  #    scale_x_continuous(trans="log10") +
  #    scale_y_continuous(trans="log10") +
  geom_abline(slope=1, intercept=0, colour="red") +
  geom_smooth() +
  facet_wrap(~fov)
dev.off()

# transcript file

tx <- fread("/NAS06/work/Spatial_omics/DBKero/CosMx_Breast/CosMX_data_Case2/Run5810_Case2_tx_file.csv")
head(tx)
# fov cell_ID x_global_px y_global_px x_local_px y_local_px z target CellComp
#1:   1       0    135100.8    40123.79    1581.82    3838.00 0   PSAP     None
#2:   1       0    137053.4    37709.06    3534.38    1423.27 0  HDAC3     None
#3:   1       0    133866.5    37268.25     347.50     982.46 0   BMP2     None
#4:   1       0    135309.6    36878.97    1790.57     593.18 0 LGALS1     None
#5:   1       0    135297.0    36838.74    1778.05     552.95 0   MSR1     None
#6:   1       0    135370.0    36732.39    1851.00     446.60 0   CD80     None

# Converting to parquet format as Davide did

write_parquet(tx, file.path(getwd(), "/Run5810_Case2_tx_file.parquet"))

tx <- read_parquet(file.path(getwd(), "/Run5810_Case2_tx_file.parquet"))
head(tx)
# str(tx)
#Classes ‘data.table’ and 'data.frame':  35343787 obs. of  9 variables:
#$ fov        : int  1 1 1 1 1 1 1 1 1 1 ...
#$ cell_ID    : int  0 0 0 0 0 0 0 0 0 0 ...
#$ x_global_px: num  135101 137053 133866 135310 135297 ...
#$ y_global_px: num  40124 37709 37268 36879 36839 ...
#$ x_local_px : num  1582 3534 348 1791 1778 ...
#$ y_local_px : num  3838 1423 982 593 553 ...
#$ z          : int  0 0 0 0 0 0 0 0 0 0 ...
#$ target     : chr  "PSAP" "HDAC3" "BMP2" "LGALS1" ...
#$ CellComp   : chr  "None" "None" "None" "None" ...
#- attr(*, ".internal.selfref")=<externalptr> 

#Transcript distribution by compartment.

table(tx$CellComp)/nrow(tx)

tx_dir <- file.path(getwd(), "tx_plots")
if (!dir.exists(tx_dir)) dir.create(tx_dir)

png(file.path(tx_dir,"cell_comp_barplot.png"), units = "cm", height = 10, width = 15, res = 180)
ggplot(tx, aes(x = CellComp, fill = CellComp)) +
  geom_bar()
dev.off()

png(file.path(tx_dir,"stats_cell_comp_barplot.png"), units = "cm", height = 10, width = 15, res = 180)
ggplot(tx, aes(x = CellComp, fill = CellComp)) +
  geom_bar(aes(y = (after_stat(count))/sum(after_stat(count)))) +
  scale_y_continuous(labels = percent)
dev.off()

png(file.path(tx_dir,"cell_comp_barplot_by_fov.png"), units = "cm", height = 30, width = 40, res = 180)
ggplot(tx, aes(x = CellComp, fill = CellComp)) +
  geom_bar() + 
  facet_wrap(~fov)
dev.off()

#Repeat with negative controls.

negcon <- filter(tx, grepl("^Neg", tx$target))
table(negcon$CellComp)/nrow(negcon)

png(file.path(tx_dir,"neg_ctrl_cell_comp_barplot.png"), units = "cm", height = 10, width = 15, res = 180)
ggplot(negcon, aes(x = CellComp, fill = CellComp)) +
  geom_bar()
dev.off()

png(file.path(tx_dir,"stats_neg_ctrl_cell_comp_barplot.png"), units = "cm", height = 10, width = 15, res = 180)
ggplot(negcon, aes(x = CellComp, fill = CellComp)) +
  geom_bar(aes(y = (after_stat(count))/sum(after_stat(count)))) +
  scale_y_continuous(labels = percent)

png(file.path(tx_dir,"neg_ctrl_cell_comp_barplot_by_fov.png"), units = "cm", height = 30, width = 40, res = 180)
ggplot(negcon, aes(x = CellComp, fill = CellComp)) +
  geom_bar() + 
  facet_wrap(~fov)
dev.off()

# Let's look at the position of transcripts by cell compartment within a FOV. Warning, the next plot is very heavy.
# observing the position of all fov1 transcripts split by subcellular compartment. In red there are the fov1 cell coordinates 

cell_coords <- as.data.frame(colData(spe)) |>
  filter(fov==1)

png(file.path(tx_dir,"transcript_position_in_fov1.png"), units = "cm", height = 30, width = 40, res = 180)
filter(tx, fov==1) |>
  ggplot(aes(x=x_local_px, y=y_local_px)) +
  geom_point() +
  geom_point(data=cell_coords, aes(x=CenterX_local_px, y=CenterY_local_px, color="red")) +
  facet_wrap(~CellComp)
dev.off()

# Let's look at the density of transcripts.

png(file.path(tx_dir,"transcript_global_density_by_cellcomp.png"), units = "cm", height = 30, width = 40, res = 180)
filter(tx) |>
  ggplot(aes(x=x_local_px, y=y_local_px)) +
  geom_hex(bins=100) +
  facet_wrap(~CellComp)
dev.off()

# Can we use the transcripts outside of cells to check the spatial patterns in the data?
# first with local coordinates: not a very good idea because the data tend to be overlapped by fov

tx$X_loc_bin <- cut(tx$x_local_px, breaks=quantile(tx$x_local_px, seq(0, 1, by=0.01)))
tx$Y_loc_bin <- cut(tx$y_local_px, breaks=quantile(tx$y_local_px, seq(0, 1, by=0.01)))

png(file.path(tx_dir,"per_compartment_local_tiles.png"), units = "cm", height = 30, width = 40, res = 180)
p1 <- filter(tx, CellComp=="None") |>
  count(X_loc_bin, Y_loc_bin) |>
  ggplot(aes(X_loc_bin, Y_loc_bin, fill=log10(n))) +
  geom_tile()

p2 <- filter(tx, CellComp=="Nuclear") |>
  count(X_loc_bin, Y_loc_bin) |>
  ggplot(aes(X_loc_bin, Y_loc_bin, fill=log10(n))) +
  geom_tile()

p3 <- filter(tx, CellComp=="Cytoplasm") |>
  count(X_loc_bin, Y_loc_bin) |>
  ggplot(aes(X_loc_bin, Y_loc_bin, fill=log10(n))) +
  geom_tile()

p4 <- filter(tx, CellComp=="Membrane") |>
  count(X_loc_bin, Y_loc_bin) |>
  ggplot(aes(X_loc_bin, Y_loc_bin, fill=log10(n))) +
  geom_tile()
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

# now with global coordinates, there are clear differential spatial patterns per each cell compartment

tx$X_glob_bin <- cut(tx$x_global_px, breaks=quantile(tx$x_global_px, seq(0, 1, by=0.01)))
tx$Y_glob_bin <- cut(tx$y_global_px, breaks=quantile(tx$y_global_px, seq(0, 1, by=0.01)))

png(file.path(tx_dir,"per_compartment_global_tiles.png"), units = "cm", height = 30, width = 40, res = 180)
p1 <- filter(tx, CellComp=="None") |>
  count(X_glob_bin, Y_glob_bin) |>
  ggplot(aes(X_glob_bin, Y_glob_bin, fill=log10(n))) +
  geom_tile() +
  ggtitle("None")

p2 <- filter(tx, CellComp=="Nuclear") |>
  count(X_glob_bin, Y_glob_bin) |>
  ggplot(aes(X_glob_bin, Y_glob_bin, fill=log10(n))) +
  geom_tile() +
  ggtitle("Nuclear")

p3 <- filter(tx, CellComp=="Cytoplasm") |>
  count(X_glob_bin, Y_glob_bin) |>
  ggplot(aes(X_glob_bin, Y_glob_bin, fill=log10(n))) +
  geom_tile() +
  ggtitle("Cytoplasm")

p4 <- filter(tx, CellComp=="Membrane") |>
  count(X_glob_bin, Y_glob_bin) |>
  ggplot(aes(X_glob_bin, Y_glob_bin, fill=log10(n))) +
  geom_tile() +
  ggtitle("Membrane")

plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

# cell boundary analysis
# loading the polygons. I hope the file I found contains the polygons, starting with fov 1

polygon_dir <- file.path(getwd(), "polygon_plots")
if (!dir.exists(polygon_dir)) dir.create(polygon_dir)

spat_obj <- fread("/NAS06/work/Spatial_omics/DBKero/CosMx_Breast/CosMX_data_Case2/Run5810_Case2-polygons.csv")
str(spat_obj)
sf_data <- st_as_sf(spat_obj, coords = c("x_local_px", "y_local_px"))

polygons <- sf_data %>%
  group_by(fov, cellID) %>%
  summarize(geometry = st_combine(geometry))#st_combine = combine several feature geometries into one, without unioning or resolving internal boundaries

polygons <- st_cast(polygons, "POLYGON")

# plot it for fov 1
# ?st_as_sf: convert foreign object into sf object
polygons1 <- filter(polygons, fov==1) 
tx1 <- filter(tx, fov==1) |> st_as_sf(coords=c("x_local_px", "y_local_px"))

png(file.path(polygon_dir,"fov1_polygon.png"), units = "cm", height = 30, width = 40, res = 180)
tm_shape(tx1) +
  tm_dots(col = "CellComp") +
  tm_shape(polygons1) + 
  tm_fill(col="grey90", alpha = 0.5) +
  tm_borders(lwd = 0.1, col = "grey50") +
  tm_layout(legend.outside = TRUE)
dev.off()

# add image overlay

f001_composite <- file.path("/NAS06/work/Spatial_omics/DBKero/CosMx_Breast/CosMX_data_Case2/CellComposite","CellComposite_F001.jpg")
comp_rast <- terra::rast(f001_composite)

f001_overlay <- file.path("/NAS06/work/Spatial_omics/DBKero/CosMx_Breast/CosMX_data_Case2/CellOverlay","CellOverlay_F001.jpg")
over_rast <- terra::rast(f001_overlay)

png(file.path(polygon_dir,"fov1_polygon_with_overlay.png"), units = "cm", height = 30, width = 40, res = 180)
tm_shape(comp_rast) +
  tm_rgb() +
  tm_shape(polygons1) +
  tm_borders(lwd = 1, col = "grey50") # +
#    tm_shape(tx1) +
#    tm_dots(col = "CellComp", alpha=0.1) +
#    tm_layout(legend.outside = TRUE)
dev.off()

# ## Cell compartments

# It looks like there are a bunch of transcripts where there are no cells. Is this enriched by specific genes? Again, let's focus on one FOV.

tot_counts <- filter(tx, fov==1) |>
  count(target, .drop=FALSE)
ec_counts <- filter(tx, fov==1 & CellComp=="None") |>
  count(target, .drop=FALSE)
ec_prop <- ec_counts$n/tot_counts$n
names(ec_prop) <- ec_counts$target

png(file.path(tx_dir,"ec_tot_counts_prop.png"), units = "cm", height = 30, width = 40, res = 180)
hist(ec_prop)
dev.off()

head(sort(ec_prop))
tail(sort(ec_prop))

# Let's repeat for all FOVs and look at the distribution across FOVs.
tot_counts <- tx |>
  count(target, fov, .drop = FALSE)
ec_counts <- filter(tx, CellComp=="None") |>
  count(target, fov, .drop = FALSE)

colnames(tot_counts)[3] <- "count"
colnames(ec_counts)[3] <- "ec_count"

tot_counts <- inner_join(tot_counts, ec_counts)

tot_counts$ec_prop <- tot_counts$ec_count/tot_counts$count

png(file.path(tx_dir,"ec_prop_hist_boxplot_by_fov.png"), units = "cm", height = 50, width = 50, res = 180)
p1 <- ggplot(tot_counts, aes(x=ec_prop)) +
  geom_histogram()
p2 <- ggplot(tot_counts, aes(y=ec_prop, group=fov)) +
  geom_boxplot()
plot_grid(p1, p2, ncol = 1)
dev.off()

arrange(tot_counts, ec_prop)
arrange(tot_counts, desc(ec_prop))

# Focus on negative controls.
tot_counts$control <- grepl("^Neg", tot_counts$target)

png(file.path(tx_dir,"neg_ctrl_ec_prop_hist_boxplot_by_fov.png"), units = "cm", height = 50, width = 50, res = 180)
p1 <- ggplot(tot_counts, aes(x=ec_prop)) +
  geom_histogram((aes(y=after_stat(density)))) +
  facet_grid(rows = vars(control))

p2 <- ggplot(tot_counts, aes(y=ec_prop, x=factor(fov), color=control)) +
  geom_boxplot()
plot_grid(p1, p2, ncol = 1)
dev.off()

# same analysis but for nuclear vs cytoplasm vs membrane

nuc_counts <- filter(tx, CellComp=="Nuclear") |>
  count(target, fov, .drop = FALSE)
cyto_counts <- filter(tx, CellComp=="Cytoplasm") |>
  count(target, fov, .drop = FALSE)
mem_counts <- filter(tx, CellComp=="Membrane") |>
  count(target, fov, .drop = FALSE)

colnames(nuc_counts)[3] <- "nuc_count"
colnames(cyto_counts)[3] <- "cyto_count"
colnames(mem_counts)[3] <- "mem_count"

tot_counts <- inner_join(tot_counts, nuc_counts)
tot_counts <- inner_join(tot_counts, cyto_counts)
tot_counts <- inner_join(tot_counts, mem_counts)

tot_counts$nuc_prop <- tot_counts$nuc_count/tot_counts$count
tot_counts$cyto_prop <- tot_counts$cyto_count/tot_counts$count
tot_counts$mem_prop <- tot_counts$mem_count/tot_counts$count

png(file.path(tx_dir,"cell_comp_ec_prop_hist.png"), units = "cm", height = 20, width = 50, res = 180)
p1 <- ggplot(tot_counts, aes(x=nuc_prop)) +
  geom_histogram()
p2 <- ggplot(tot_counts, aes(x=cyto_prop)) +
  geom_histogram()
p3 <- ggplot(tot_counts, aes(x=mem_prop)) +
  geom_histogram()
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

png(file.path(tx_dir,"cell_comp_ec_prop_boxplot_by_fov.png"), units = "cm", height = 60, width = 60, res = 180)
p1 <- ggplot(tot_counts, aes(y=nuc_prop, group=fov)) +
  geom_boxplot()
p2 <- ggplot(tot_counts, aes(y=cyto_prop, group=fov)) +
  geom_boxplot()
p3 <- ggplot(tot_counts, aes(y=mem_prop, group=fov)) +
  geom_boxplot()
plot_grid(p1, p2, p3, ncol = 1)
dev.off()

arrange(tot_counts, desc(nuc_prop))
arrange(tot_counts, desc(cyto_prop))
arrange(tot_counts, desc(mem_prop))

# Compare nuclear vs cell summaries

cell_counts <- tx |>
  filter(cell_ID != 0) |>
  count(cell_ID, fov, .drop = FALSE)
cell_counts

nucleus_counts <- tx |>
  filter(cell_ID != 0 & CellComp == "Nuclear") |>
  count(cell_ID, fov, .drop = FALSE)
nucleus_counts

colnames(cell_counts)[3] <- "cell_count"
colnames(nucleus_counts)[3] <- "nucleus_count"

cell_counts <- inner_join(cell_counts, nucleus_counts)

png(file.path(tx_dir,"cell_count~nucleus_count.png"), units = "cm", height = 20, width = 20, res = 180)
cell_counts |>
  ggplot(aes(cell_count, nucleus_count)) +
  geom_point() +
  geom_abline(color="red") +
  geom_smooth()
dev.off()

png(file.path(tx_dir,"cell_count~nucleus_count_by_fov.png"), units = "cm", height = 30, width = 40, res = 180)
cell_counts |>
  ggplot(aes(cell_count, nucleus_count)) +
  geom_point() +
  geom_abline(color="red") +
  geom_smooth() +
  facet_wrap(~fov)
dev.off()












