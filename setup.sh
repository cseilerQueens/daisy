module load  StdEnv/2023  gcc/12.3  udunits/2.2.28  hdf/4.2.16  gdal/3.7.2  r/4.3.1 meta-farm/1.0.2

# Install renv
R -e "install.packages("renv")"

# Initiate renv
R -e "renv::init()"

# Check local paths
.libPaths()

# If .libPaths() does not point to your local renv folder, then add is mannualy to your .Rprofile:
# .libPaths("/home/cseiler/daisy/renv")

# Install dependencies
export R_LIBS="/home/cseiler/daisy/renv"
R -e "install.packages(c('sf', 'terra'), repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('ncdf4', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('ncdf4.helpers', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('foreach', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('raster', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('devtools', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('roxygen2', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('GA', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('latex2exp', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages(c('shiny', 'foreach', 'slickR'), repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
R -e "install.packages('readr', repos='https://mirror.csclub.uwaterloo.ca/CRAN/', dep=TRUE)"
