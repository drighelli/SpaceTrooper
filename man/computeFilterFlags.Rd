% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.R
\name{computeFilterFlags}
\alias{computeFilterFlags}
\title{computeFilterFlags}
\usage{
computeFilterFlags(
  spe,
  fs_threshold = 0.5,
  use_fs_quantiles = FALSE,
  total_threshold = 0,
  ctrl_tot_ratio_threshold = 0.1
)
}
\arguments{
\item{spe}{A `SpatialExperiment` object with spatial transcriptomics data.}

\item{fs_threshold}{A numeric value for the threshold of `flag_score` to
identify outliers. Default is `0.5`.}

\item{use_fs_quantiles}{A logical value indicating whether to use quantiles
for the `flag_score` threshold. Default is `FALSE`.}

\item{total_threshold}{A numeric value for the threshold of total counts to
identify cells with zero counts. Default is `0`.}

\item{ctrl_tot_ratio_threshold}{A numeric value for the threshold of
control-to-total ratio to flag outliers. Default is `0.1`.}
}
\value{
The `SpatialExperiment` object with added filter flags in `colData`.
}
\description{
Compute Filter Flags for SpatialExperiment

This function calculates various flags to identify outliers in a
`SpatialExperiment` object based on quality control metrics.
}
\details{
The function flags cells based on zero counts, control-to-total
ratio, and `flag_score` to identify potential outliers. It also combines
these flags into a single filter flag.
}
\examples{
# TBD
}
