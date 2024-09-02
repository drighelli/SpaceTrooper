addLabelsDBKero <- function(spe, filename)
{
    labels <- read.table(file="~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/cosmx_dbkero_IST_labels_complete_simple_colors.tsv", sep="\t", header=TRUE)
    spe$labels <- NA
    spe$labels_colors <- "black"
    spe$labels[match(labels$cell_id, spe$cell_id)] <- labels$InSituType_simple
    spe$labels_colors[match(labels$cell_id, spe$cell_id)] <- labels$IST_simple_colors
    spe
}
