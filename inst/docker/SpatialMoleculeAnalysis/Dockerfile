# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Update apt-get
RUN apt-get update \
	## Install the python package tensorflow
	## && pip install tensorflow		\
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install required Bioconductor package
RUN R -e 'BiocManager::install(c("SpatialExperiment", "data.table", "scater", ' \
    "cowplot", "ggplot2", "matrixStats", "dplyr", "tidyr", "tibble", "arrow", ' \
    "scales", "sf", "tmap")'
