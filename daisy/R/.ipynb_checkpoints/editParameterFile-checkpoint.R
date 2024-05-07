################################################################################
#' Edited text of parameter file
#' @description This function reads the content of the parameter file (run_parameters.txt) and replaces the parameter values of a selected parameter with new values.
#' 
#' @param parameterFile A string that gives the path and name of the model parameter file
#' @param parameterName A string that gives the parameter name
#' @param parameterValues An R object that contains the new parameter values
#' 
#' @return An R object that contains the text of the updated parameter file
#'
#' @examples
#' # Get the values of a particular parameter, e.g. vmax and change it
#' parameterFile <- '/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/simulations/transient_CRUJRAv2/run_parameters.txt'
#' parameterName <- 'vmax'
#' parameterValues <- getParameterValues(parameterFile = parameterFile, parameterName = parameterName)
#' parameterValues <- parameterValues * 2 # new parameter values
#' # Insert the new parameter values into the parameter file 
#' lines <- editParameterFile(parameterFile, parameterName, parameterValues)
#' # writeLines(lines, 'new_run_parameters.txt')
#'
#' @export

editParameterFile <- function(parameterFile, parameterName, parameterValues) {

    # Read in the parameter file run_parameters.txt
    lines <- readLines(parameterFile)

    # Find the lines containing the parameter that will be updated
    linesToReplace <- grep(pattern = parameterName, x = lines)

    # Replace the lines with the specified replacement string

    counter <- 0
    for (line in linesToReplace) {
        counter <- counter + 1
        rowName <- row.names(parameterValues)[counter]
        vector <- as.character(parameterValues[counter, ])
        csv <- paste(vector, collapse = ", ")
        updatedLine <- paste(rowName, csv, sep = "= ")
        lines[line] <- updatedLine
    }

    return(lines)
}
