################################################################################
#' Data frame of model with model output and reference data
#' @description This function extracts the data from (i) a model run conducted for a single grid cell and (ii) multiple globally gridded reference data sets.
#' The data frame only includes the time steps that all data sets have in common.
#' 
#' @param parameterFile A string that gives the path and name of the model parameter file
#' @param parameterName A string that gives the name of the parameter of interest
#' 
#' @return A matrix with the deired parameter values
#'
#' @examples
#' parameterFile <- '/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/simulations/transient_CRUJRAv2/run_parameters.txt'
#' parameterName <- 'vmax'
#' getParameterValues(parameterFile = parameterFile, parameterName = parameterName)
#'
#' @export

getParameterValues <- function(parameterFile, parameterName) {

    # Read in the run_parameters.txt file
    data <- readLines(parameterFile)

    # Get the section relevant to your parameter of choice
    data <- data[grep(parameterName, data)]

    # Clean data from unwanted information
    toBeExcluded <- grep(pattern = "!", x = data)
    if (length(toBeExcluded) > 0) {
        data <- data[-toBeExcluded]
    }

    if (length(data) > 4) {
        data <- data[1:4]
    }

    # Convert to matrix
    data <- strsplit(data, ",|=")
    data <- do.call(rbind, data)
    rowNames <- data[, 1]
    data <- data[, -1]

    numberOfColumns <- ncol(data)
    data <- as.numeric(data)

    if (is.integer(numberOfColumns)) {
        data <- matrix(data, ncol = numberOfColumns)
    } else {

        data <- matrix(data, ncol = 1)
    }
    rownames(data) <- rowNames
    return(data)

}
