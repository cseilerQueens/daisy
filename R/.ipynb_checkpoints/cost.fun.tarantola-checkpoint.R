################################################################################
#' Data frame of model with model output and reference data
#' @description This function is a cost function defined by Tarantola (2005) (p. 64).
#' @param parameterValues Parameter values
#' @param parameterNames Parameter names
#' @param parameterFile Parameter file (path and name)
#' @param run_classic_file Shell script that submits a CLASSIC run on an HPC
#' @param classicFinished_file path and name of the classic.out file
#' 
#' @return A textfile with the optimized parameter values
#'
#' @examples
#' \dontrun{
#' library(daisy)
#' library(foreach)
#' 
#' # Provide path to currently existing parameter file
#' # Provide cost function with new parameter values
#' # The cost function will first read the old values from the parameter file
#' # The cost function will then overwrite the parameter file with the new values and run the model with the new values
#' 
#' setwd('/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation')
#' # Copy default parameter run_parameters.txt file to current work directory
#' system('cp /home/cseiler/classic/classic_code/CLASSIC/configurationFiles/template_run_parameters.txt run_parameters.txt')
#' 
#' # Select parameters and obtain their parameter values. This is currently implemented for one PFT-specific parameter only 
#' parameterFile <- '/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/run_parameters.txt'
#' parameterName <- 'vmax'
#' vmax <- getParameterValues(parameterFile = parameterFile, parameterName = parameterName)
#' 
#' # For the purpose of this example, create some new parameter values that differ from the default values.
#' # The new values will be used as the posterior values, while the default values
#' # are the prior values.  This is not necessary when using the cost function
#' # inside the optim() function
#' vmax <- vmax * 1.02
#'
#' parameterNames <- c('vmax')
#' parameterValues <- c(vmax)

#' run_classic_file <- '/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/run_classic.sh'
#' classicFinished_file <- '/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/simulations/classic.out'
#' result <- cost.fun.tarantola(parameterValues, parameterNames, parameterFile, run_classic_file, classicFinished_file)
#' }
#' 
#' @export

cost.fun.tarantola <- function(normParameterValues, range, minValue, parameterNames,
    parameterFile, run_classic_file, classicFinished_file) {

    # un-normalize parameter values
    parameterValues <- intFun.unnormalize(normParameterValues, range, minValue)

    # Get the parameter names
    parameterName <- parameterNames[1]

    # Get the prior parameter values that are provided by the current parameter
    # file
    priorParameterValues <- getParameterValues(parameterFile, parameterName)

    # Get the posterior parameter values that are provided as an input to this
    # function.

    # CLASSIC has 12 PFTs, however, the default configuration only uses 9 PFTs.
    # The PFTs that are not being used have parameter values that are equal to
    # zero.  Those values need to be excluded, otherwise the matrix inversion
    # does not work as you can't divide by zero. Find the location where PFT
    # values are zero:

    non_zero_indices <- which(priorParameterValues != 0, arr.ind = TRUE)
    newValues <- priorParameterValues

    # Replace prior with new parameter values, skipping all instances where the
    # parameter values equal zero in the parameter file

    newValues[non_zero_indices] <- priorParameterValues[non_zero_indices] - priorParameterValues[non_zero_indices] +
        parameterValues

    # Overwrite the parameter file with the posterior parameter values
    lines <- editParameterFile(parameterFile, parameterName, parameterValues = newValues)
    writeLines(lines, parameterFile)

    # Run CLASSIC with new parameter file
    removeFile <- paste("rm", classicFinished_file, sep = " ")
    system(removeFile)
    system(run_classic_file)

    # Wait until CLASSIC finishes
    while (!file.exists(classicFinished_file)) {
        Sys.sleep(10)
    }

    print("CLASSIC run completed")

    # Obtain model and reference data. This is currently hard-coded and will be
    # replaced later
    nc.mod <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/simulations/transient_CRUJRAv2/outputFiles/gpp_monthly.nc"
    nc.ref01 <- "/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63/gpp_GOSIF_128x64.nc"
    nc.ref02 <- "/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63/GPP.ensembleMedian.CRUNCEPv6.monthly_128x64.nc"
    nc.ref.list <- list(nc.ref01, nc.ref02)
    ref.id.list <- list("GOSIF", "FLUXCOM")

    unit.conv.mod <- 86400 * 1000  # optional unit conversion for model data
    unit.conv.ref <- c(1, 1)
    variable.unit <- "gC m$^{-2}$ day$^{-1}$"  # unit after conversion (LaTeX notation)
    data <- getDataModRef(nc.mod, nc.ref.list, ref.id.list, unit.conv.mod, unit.conv.ref,
        variable.unit, timeInt = "month")

    # move outputs to a new folder so that they are not overwritten
    timeStamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    modelOutputFolder <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/simulations/transient_CRUJRAv2"
    NewModelOutputFolder <- paste(modelOutputFolder, timeStamp, sep = "_")
    moveModelOutputFolder <- paste("mv", modelOutputFolder, NewModelOutputFolder,
        sep = " ")
    system(moveModelOutputFolder)

    # get model data
    mod <- data$mod

    # get reference data
    ref01 <- subset(data, select = ref.id.list[[1]])
    ref02 <- subset(data, select = ref.id.list[[2]])

    # compute mean of two reference data for each time step
    ref <- (ref01 + ref02)/2
    colnames(ref) <- "ref"
    ref <- ref$ref

    ### TEST different cost function start

    N <- length(mod)
    E <- (1/N * sum((mod - ref)^2))^(1/2)

    # When using GA, the lagorithm maximizes the fir function. So, let's use
    # 1/E as a fit, assuming E will never equal zero

    E <- 1/E

    ### TEST different cost function end

    # compute the root-mean-square-error (RMSE) from both reference data sets
    n <- length(mod)
    rmse <- (sum((ref01 - ref02)^2)/n)^0.5  # this is a single number

    # Get the squared RMSE (single value)
    rmse2 <- rmse^2

    # Put the model output data into a single-column matrix
    mod <- matrix(mod, ncol = 1)

    # Put the reference data into a single-column matrix
    ref <- matrix(ref, ncol = 1)

    # Create the covariance operator that gives the Gaussian uncertainty of
    # your observed data (diagonal matrix)
    C_D <- diag(n) * rmse2

    # print(mod) print(ref) print(C_D)

    # Rename to be consistent with naming convention used by Tarantola, 2005
    # (p. 64)
    g <- mod
    d_obs <- ref

    # new parameter values
    m <- matrix(parameterValues, ncol = 1)

    # prior parameter values
    m_prior <- matrix(priorParameterValues[non_zero_indices], ncol = 1)

    # Next, we need to create the covariance operator that gives the Gaussian
    # uncertainty of the model parameter (diagonal matrix) Similarly to the
    # estimation of CD, we have set an 'automatic' approach to estimate CM. The
    # *standard deviation* (not variance) is set to X% of the range of
    # variation of a given parameter. In our earlier work, X was indeed 40, but
    # recent studies (Kuppel et al. and Bacour et al., 2023 - although the
    # final value is not mentionned) has reduced it to 15%.

    C_M_fraction <- 0.15

    # You need to read in the prescribed uncertainty ranges here Quick solution
    # for now, replace this for later:
    range <- 1.1 * m - 0.9 * m  # REPLACE THIS LINE SUCH THAT THE RANGES ARE NOT HARD-CODED!
    sd <- range * C_M_fraction

    n <- length(m)
    C_M_diag <- diag(n)
    m.var <- sd^2
    m.var <- rep(m.var, n)
    m.var <- matrix(m.var, ncol = n)
    C_M <- C_M_diag * m.var

    # Cost function by Tarantola (2005), equation 3.32, p. 64
    termI <- t(g - d_obs) %*% solve(C_D) %*% (g - d_obs)
    termII <- t(m - m_prior) %*% solve(C_M) %*% (m - m_prior)
    S <- 0.5 * (termI + termII)

    # Write current parameter values to text file for analysis purpose
    # write.table(x = matrix(c(S, parameterValues), nrow = 1), file =
    # 'allParameterValues', append = TRUE, row.names = FALSE, col.names =
    # FALSE)
    write.table(x = matrix(c(E, parameterValues), nrow = 1), file = "allParameterValues",
        append = TRUE, row.names = FALSE, col.names = FALSE)

    # return(S)
    return(E)
}
