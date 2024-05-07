setwd("/home/cseiler/daisy")

.libPaths("/home/cseiler/daisy/renv")

myDirectory <- "/home/cseiler/daisy/daisy"
myPackage <- "daisy"
library(devtools)
library(roxygen2)
setwd(myDirectory)

try(remove.packages("daisy"))

devtools::document(file.path(myDirectory))
devtools::install(file.path(myDirectory), dep=FALSE)

library(daisy)

# renv::snapshot()

# Check package
# setwd(myDirectory)
# devtools::check(args = c('--as-cran','--run-donttest'))
