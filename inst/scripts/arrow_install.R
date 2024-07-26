Sys.setenv(LIBARROW_MINIMAL = "false")
Sys.setenv(ARROW_WITH_ZSTD = "ON")
# Sys.setenv(ARROW_S3="ON")
Sys.setenv(NOT_CRAN="true")
# install.packages("arrow", repos = "https://arrow-r-nightly.s3.amazonaws.com")
system("echo $ARROW_WITH_ZSTD")
system("echo $LIBARROW_MINIMAL")
install.packages("arrow")
library(arrow)
arrow::install_arrow(minimal=FALSE,nightly=TRUE,verbose=TRUE)
arrow::arrow_info()
pol_file <- "~/Downloads/Xenium_data/pancreas/outs/cell_boundaries.parquet"
devtools::load_all()
polygons <- readPolygonsXenium(pol_file, keepMultiPol=TRUE)

## to read xenium parquet it's needed the support for ZSTD in arrow package
## just install arrow from R not with external libraries
