#' Function to adjust code instance-specific runs.
#' @param Identifier The identifier for which code instance we are executing.
#' @param oldOutDir The old output directory
#' @param newOutDir The new output directory
#' @param argNamesFile File listing all argument names in use.
#' @param argValuesFile File listing the values of all arguments in use.
#' @param argsToChange A list of the code instance-specific arguments.
#' @param farmName The name of the farm.
setArgs <- function(identifier, oldOutDir, newOutDir, argNamesFile, argValuesFile, 
                    argsToChange, farmName) {
  
  # Names of all arguments in use.
  names <- readLines(argNamesFile)
  
  # Split up arguments into a list.
  argsToChange <- unlist(strsplit(argsToChange, ","))
  # Trim white space.
  argsToChange <- lapply(argsToChange, trimws)
  
  for (argument in argsToChange) {
    # .libPaths(c(.libPaths(), "/home/cseiler/projects/def-cseiler-ab/cseiler/AMBER/renv/"))
    library(readr)
    
    # The line where the argument appears
    line <- match(argument, names)
    
    # Adjust the argument if it is in use.
    if (!is.na(line)) {
      
      # Argument is parameterFile.
      if (argument == "parameterFile") {
        originalValue <- read_lines(argValuesFile, skip = line-1, n_max = 1)
        # Copy file.
        cmd <- paste("cp", originalValue, newOutDir)
        system(cmd)
        Sys.sleep(1)
        # Adjust argument value.
        cmd <- paste0("sed -i '", line, " s|", originalValue, "|", newOutDir, 
                      "/run_parameters.txt|' ", argValuesFile)
        system(cmd)
        Sys.sleep(1)
        
        # Argument is run_classic_file.  
      } else if (argument == "run_classic_file") {
        originalValue <- read_lines(argValuesFile, skip = line-1, n_max = 1)
        # Copy file.
        cmd <- paste("cp", originalValue, newOutDir)
        system(cmd)
        Sys.sleep(1)
        # Adjust parameter value.
        cmd <- paste0("sed -i '", line, "s|", originalValue, "|", newOutDir, 
                      "/run_classic.sh|' ", argValuesFile)
        system(cmd)
        Sys.sleep(1)
        
        # Adjust file path in run_classic_file.
        cmd <- paste0("sed -i 's|", oldOutDir, "/$simulationID|", newOutDir, 
                      "/$simulationID|g' ", newOutDir, "/run_classic.sh")
        system(cmd)
        Sys.sleep(1)
      }
      
      # Both mod.list and modelOutputFolder have the data-assimilation/simulations/daisyRun
      # file path in their values. Since each case executed by meta-farm has its own
      # RUN file, the file paths need to be changed from /simulations/ to /RUNX/.
      # As such, each case will get its own copy of argValues.txt with the correct file paths.
      cmd <- paste0("sed -i 's|", oldOutDir, "|", newOutDir, "|g' ", argValuesFile)
      system(cmd)
      
      Sys.sleep(1)
    }
  }
}

