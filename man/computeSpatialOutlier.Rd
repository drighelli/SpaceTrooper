% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.R
\name{computeSpatialOutlier}
\alias{computeSpatialOutlier}
\title{computeSpatialOutlier}
\usage{
computeSpatialOutlier(
  spe,
  compute_by = NULL,
  method = c("mc", "scuttle", "both"),
  mcDoScale = FALSE,
  scuttleType = c("both", "lower", "higher")
)
}
\arguments{
\item{spe}{a SpatialExperiment object with target_counts, area in micron
and log2 of the aspect ratio in the `colData`.}

\item{compute_by}{character indicating a `colData` column name on which compute
the outlier.}

\item{method}{one of `mc`, `scuttle`, `both`.
Use `mc` for medcouple, `scuttle` for median absolute deviations as computed
in `scuttle`, `both` for computing both of them.}

\item{mcDoScale}{logical indicating if the values to compute the medcouple
for the outlier detection should be scaled (default is FALSE, as suggested
by the original Medcouple authors.). See \link[robustbase]{mc} for further
readings.}
}
\value{
a SpatialExperiment object with additional column(s) (named as
the column name indicated in `column_by` followed by the outlier_sc/mc
nomenclature) with the outlier detection as `outlier.filter` logical class
object. This allows to store the thresholds as attributes of the column.
use attr(,"thresholds") to retrieve them.
}
\description{
Computes outliers based on the Area (in micron) of the experiment.
It gives the possibility to choose between the medcouple (mc method argument)
and the MADs (scuttle method argument).
}
\details{
The medcouple method is a measure for the skeweness of univariate distribution
as described in Hubert M. et al. (2008).
In particular, the computed medcouple value must be in a range between -0.6
and 0.6 to computed adjusted boxplots and perform the outlier detection.
For median absolute deviations (MADs) method we just wrap the isOutlier
function in the scuttle package. Please see McCarthy DJ et al (2017)
for further details.
}
\examples{
# TBD
}
