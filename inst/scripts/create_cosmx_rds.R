suppressPackageStartupMessages({
    library(devtools)
    library(SpatialExperiment)
    library(SummarizedExperiment)
})


input_dir <- "/Users/inzirio/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2"

output_rds <- "~/SpaceTrooper_qs_check/cosmx_case2_spe_raw.rds"
output_rds <- path.expand(output_rds)


library(SpaceTrooper)

if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
}

message("Reading CosMx dataset from:")
message(input_dir)

spe <- readCosmxSPE(input_dir)

message("Object class:")
print(class(spe))

message("Object dimensions:")
print(dim(spe))

message("colData columns:")
print(colnames(SummarizedExperiment::colData(spe)))

message("Saving RDS to:")
message(output_rds)

saveRDS(spe, output_rds)

message("Done.")
