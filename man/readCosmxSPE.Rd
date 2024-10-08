% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readCosMx.R
\name{readCosmxSPE}
\alias{readCosmxSPE}
\title{readCosmxSPE}
\usage{
readCosmxSPE(
  dirname,
  sample_name = "sample01",
  coord_names = c("CenterX_global_px", "CenterY_global_px"),
  countmatfpattern = "exprMat_file.csv",
  metadatafpattern = "metadata_file.csv",
  polygonsfpattern = "polygons.csv",
  keep_polygons = FALSE,
  fovposfpattern = "fov_positions_file.csv",
  fov_dims = c(xdim = 4256, ydim = 4256)
)
}
\arguments{
\item{dirname}{A character string specifying the directory containing the
CosMx data files.}

\item{sample_name}{A character string specifying the sample name. Default is
`"sample01"`.}

\item{coord_names}{A character vector specifying the names of the spatial
coordinate columns in the data. Default is `c("CenterX_global_px",
"CenterY_global_px")`.}

\item{countmatfpattern}{A character string specifying the pattern to match
the count matrix file. Default is `"exprMat_file.csv"`.}

\item{metadatafpattern}{A character string specifying the pattern to match
the metadata file. Default is `"metadata_file.csv"`.}

\item{polygonsfpattern}{A character string specifying the pattern to match
the polygons file. Default is `"polygons.csv"`.}

\item{keep_polygons}{A logical value indicating whether to include polygon
data in the resulting `SpatialExperiment` object. Default is `FALSE`.}

\item{fovposfpattern}{A character string specifying the pattern to match the
FOV positions file. Default is `"fov_positions_file.csv"`.}

\item{fov_dims}{A named numeric vector specifying the dimensions of the FOV
in pixels. Default is `c(xdim=4256, ydim=4256)`.}
}
\value{
A `SpatialExperiment` object containing the read CosMx data,
including count matrices, metadata, and optionally polygons.
}
\description{
Read and Construct a SpatialExperiment Object from CosMx Data

This function reads in data from Nanostring CosMx files and constructs a
`SpatialExperiment` object, optionally including polygon data.
}
\details{
The function reads in the specified files for count matrices,
metadata, and FOV positions, and constructs a `SpatialExperiment` object.
Optionally, polygon data can be read and added to the object.
}
\examples{
# Assuming the data files are located in "path/to/cosmx_data":
#spe <- readCosmxSPE(dirname = "path/to/cosmx_data")
}
