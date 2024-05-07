################################################################################
#' Data frame of model with model output and reference data
#' @description This function is a cost function
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
#' 

cost.fun <- function(normParameterValues, range, minValue, parameterNames, parameterFile,
    run_classic_file, classicFinished_file) {

    # un-normalize parameter values
    parameterValues <- intFun.unnormalize(normParameterValues, range, minValue)

    # Get the parameter names
    parameterName <- parameterNames[1]

    # Get the default parameter values 
    
    defaultParameterValues <- getParameterValues(parameterFile, parameterName)

    # CLASSIC has 12 PFTs, however, the default configuration only uses 9 PFTs.
    # The PFTs that are not being used have parameter values that are equal to
    # zero.  Those values need to be excluded, otherwise the matrix inversion
    # does not work as you can't divide by zero. Find the location where PFT
    # values are zero:

    non_zero_indices <- which(defaultParameterValues != 0, arr.ind = TRUE)
    newValues <- defaultParameterValues

    # Replace prior with new parameter values, skipping all instances where the
    # parameter values equal zero in the parameter file

    newValues[non_zero_indices] <- defaultParameterValues[non_zero_indices] - defaultParameterValues[non_zero_indices] +
        parameterValues

    # Overwrite the parameter file with the new parameter values
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
    Sys.sleep(2)
    print("CLASSIC run completed")

    # Obtain model and reference data. This is currently hard-coded and will be
    # replaced later
    nc.mod <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/simulations/daisyRun/netcdf_files/gpp_monthly.nc"
    nc.ref01 <- "/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63/gpp_GOSIF_128x64.nc"
    nc.ref02 <- "/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63/GPP.ensembleMedian.CRUNCEPv6.monthly_128x64.nc"
    nc.ref.list <- list(nc.ref01, nc.ref02)
    ref.id.list <- list("GOSIF", "FLUXCOM")

    unit.conv.mod <- 86400 * 1000  # optional unit conversion for model data
    unit.conv.ref <- c(1, 1)
    variable.unit <- "gC m$^{-2}$ day$^{-1}$"  # unit after conversion (LaTeX notation)

    myLevel <- 1

    mod <- raster::brick(nc.mod, level = myLevel)
    suppressWarnings(mod <- raster::rotate(mod))


    ref <- raster::brick(nc.ref01)
    suppressWarnings(ref <- raster::rotate(ref))

    # model data dates
    dates.mod <- intFun.getZ(nc.mod)
    mod <- raster::setZ(mod, dates.mod)
    names(mod) <- dates.mod
    dates.mod <- format(as.Date(dates.mod), "%Y-%m")  # only year and month
    start.date.mod <- min(dates.mod)
    end.date.mod <- max(dates.mod)

    # reference data dates
    dates.ref <- intFun.getZ(nc.ref01)
    ref <- raster::setZ(ref, dates.ref)
    names(ref) <- dates.ref
    dates.ref <- format(as.Date(dates.ref), "%Y-%m")  # only year and month
    start.date.ref <- min(dates.ref)
    end.date.ref <- max(dates.ref)

    #---------------------------------------------------------------------------

    # 1.3 Remaining part applies to both regular and irregular gridded data

    #---------------------------------------------------------------------------

    # find common time period
    start.date <- max(start.date.mod, start.date.ref)
    end.date <- min(end.date.mod, end.date.ref)

    # subset common time period
    mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date &
        format(as.Date(raster::getZ(mod)), "%Y-%m") <= end.date)]]
    ref <- ref[[which(format(as.Date(raster::getZ(ref)), "%Y-%m") >= start.date &
        format(as.Date(raster::getZ(ref)), "%Y-%m") <= end.date)]]

    # get layer names
    mod.names <- base::names(mod)
    ref.names <- base::names(ref)

    # unit conversion if appropriate
    mod <- mod * unit.conv.mod
    ref <- ref * unit.conv.ref

    #---------------------------------------------------------------------------

    # II Statistical analysis

    #---------------------------------------------------------------------------

    # (1) Bias

    #---------------------------------------------------------------------------
    # create a mask to excludes all grid cells that the model and reference
    # data do not have in common.  This mask varies in time.
    mask <- (mod * ref)
    mask <- mask - mask + 1
    mod <- mod * mask
    names(mod) <- mod.names  # this adds the corresponding dates
    ref <- ref * mask
    names(ref) <- ref.names  # this adds the corresponding dates
    # now mod and ref are based on the same grid cells

    mod.mean <- raster::mean(mod, na.rm = TRUE)  # time mean
    ref.mean <- raster::mean(ref, na.rm = TRUE)  # time mean
    bias <- mod.mean - ref.mean  # time mean

    mod.sd <- raster::calc(mod, fun = sd, na.rm = TRUE)  # standard deviation of model data
    ref.sd <- raster::calc(ref, fun = sd, na.rm = TRUE)  # standard deviation of reference data
    epsilon_bias <- abs(bias)/ref.sd
    epsilon_bias[epsilon_bias == Inf] <- NA  # relative error
    bias.score <- exp(-1 * epsilon_bias)  # bias score as a function of space
    S_bias <- mean(raster::getValues(bias.score), na.rm = TRUE)  # scalar score

    #---------------------------------------------------------------------------

    # (2) root mean square error (rmse)

    #---------------------------------------------------------------------------
    ESM.mode <- FALSE
    # Set ESM.mode to FALSE when forcing model with quasi-observed data (e.g.
    # reanalysis)
    if (ESM.mode == FALSE) {
        rmse <- intFun.rmse(mod, ref)  # rmse
        mod.anom <- mod - mod.mean  # anomaly
        ref.anom <- ref - ref.mean  # anomaly
    }

    # Set ESM.mode to TRUE for online experiments or when forcing model with
    # earth system model data
    if (ESM.mode == TRUE) {
        index <- format(as.Date(names(ref), format = "X%Y.%m.%d"), format = "%m")
        index <- as.numeric(index)
        mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
        ref.clim.mly <- raster::stackApply(ref, index, fun = raster::mean)

        mod.cycle <- mod.clim.mly  # does not necessarily start in Jan or end in Dec
        ref.cycle <- ref.clim.mly  # does not necessarily start in Jan or end in Dec

        # Ensure that the order is correct
        JanToDec <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6",
            "index_7", "index_8", "index_9", "index_10", "index_11", "index_12")
        mod.clim.mly <- mod.clim.mly[[JanToDec]]
        ref.clim.mly <- ref.clim.mly[[JanToDec]]

        rmse <- intFun.rmse(mod.clim.mly, ref.clim.mly)  # rmse

        mod.anom <- mod.clim.mly - mod.mean  # anomaly
        ref.anom <- ref.clim.mly - ref.mean  # anomaly

    }

    crmse <- intFun.crmse(mod.anom, ref.anom)  # centralized rmse

    #-------------------------------------------------------------
    epsilon_rmse <- crmse/ref.sd
    epsilon_rmse[epsilon_rmse == Inf] <- NA  # relative error
    rmse.score <- exp(-1 * epsilon_rmse)  # rmse score as a function of space
    S_rmse <- mean(raster::getValues(rmse.score), na.rm = TRUE)

    #---------------------------------------------------------------------------

    # (3) phase shift

    #---------------------------------------------------------------------------

    index <- format(as.Date(names(ref), format = "X%Y.%m.%d"), format = "%m")
    index <- as.numeric(index)
    mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
    ref.clim.mly <- raster::stackApply(ref, index, fun = raster::mean)

    mod.cycle <- mod.clim.mly  # does not necessarily start in Jan or end in Dec
    ref.cycle <- ref.clim.mly  # does not necessarily start in Jan or end in Dec

    # Ensure that the order is correct
    JanToDec <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6",
        "index_7", "index_8", "index_9", "index_10", "index_11", "index_12")
    mod.clim.mly <- mod.clim.mly[[JanToDec]]
    ref.clim.mly <- ref.clim.mly[[JanToDec]]

    bias.clim.mly <- mod.clim.mly - ref.clim.mly

    # find month of seasonal peak

    # In most cases, we are interested in the timing of the seasonal maximum
    # value

    # In some cases, however, the seasonal peak is a minimum, e.g. NEE = RECO -
    # GPP

    phaseMinMax <- "phaseMax"
    #
    if (phaseMinMax == "phaseMax") {
        mod.max.month <- raster::which.max(mod.clim.mly)
        ref.max.month <- raster::which.max(ref.clim.mly)
    }

    if (phaseMinMax == "phaseMin") {
        mod.max.month <- raster::which.min(mod.clim.mly)
        ref.max.month <- raster::which.min(ref.clim.mly)
    }

    # get shortest time distance between these months
    abs.diff <- abs(mod.max.month - ref.max.month)  # absolute difference from 0 to 12 months
    phase <- raster::calc(abs.diff, intFun.theta)  # shortest distance from 0 to 6 months (theta)
    phase.score <- 0.5 * (1 + cos(2 * pi * phase/12))  # score from 0 (6 months) to 1 (0 months)
    S_phase <- mean(raster::getValues(phase.score), na.rm = TRUE)  # scalar score

    #---------------------------------------------------------------------------

    # (4) interannual variability

    #---------------------------------------------------------------------------

    years <- floor(raster::nlayers(mod)/12)  # total number of years
    months <- years * 12  # number of months considering complete years only
    mod.fullyear <- raster::subset(mod, 1:months)
    ref.fullyear <- raster::subset(ref, 1:months)
    c.mod <- raster::calc(mod.cycle, fun = function(x) {
        rep(x, years)
    })  # climatological cycle for all months (mod)
    c.ref <- raster::calc(ref.cycle, fun = function(x) {
        rep(x, years)
    })  # climatological cycle for all months (ref)
    mod.iav <- sqrt(raster::mean((mod.fullyear - c.mod)^2, na.rm = TRUE))  # interannual variability  (mod)
    ref.iav <- sqrt(raster::mean((ref.fullyear - c.ref)^2, na.rm = TRUE))  # interannual variability  (ref)
    # set values close to zero to NA
    ref.iav.na <- ref.iav
    ref.iav.na[ref.iav.na < 10^(-5)] <- NA

    epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav.na
    epsilon_iav[epsilon_iav == Inf] <- NA  # I changed Eq. 26 so that epsilon_iav >= 0
    iav.score <- exp(-1 * epsilon_iav)  # iav score as a function of space
    S_iav <- mean(raster::getValues(iav.score), na.rm = TRUE)  # scalar score (not weighted)

    #---------------------------------------------------------------------------

    # (5) dist

    #---------------------------------------------------------------------------

    mod.sigma.scalar <- stats::sd(raster::getValues(mod.mean), na.rm = TRUE)  # standard deviation of period mean data
    ref.sigma.scalar <- stats::sd(raster::getValues(ref.mean), na.rm = TRUE)  # standard deviation of period mean data
    sigma <- mod.sigma.scalar/ref.sigma.scalar
    y <- raster::getValues(mod.mean)
    x <- raster::getValues(ref.mean)
    reg <- stats::lm(y ~ x)
    R <- sqrt(summary(reg)$r.squared)
    S_dist <- 2 * (1 + R)/(sigma + 1/sigma)^2

    #---------------------------------------------------------------------------

    # overall scores

    #---------------------------------------------------------------------------

    S <- mean(c(S_bias, S_rmse, S_phase, S_iav, S_dist))

    timeStamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    # Write current scores and parameter values to text file for analysis purpose
    write.table(x = matrix(c(S, S_bias, S_rmse, S_phase, S_iav, S_dist, parameterValues, timeStamp), nrow = 1), file = "allParameterValues",
        append = TRUE, row.names = FALSE, col.names = FALSE)

    # Move outputs to a new folder so that they are not overwritten hard-coded
    # path will be replaced
    modelOutputFolder <- "/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/simulations/daisyRun"
    NewModelOutputFolder <- paste(modelOutputFolder, timeStamp, sep = "_")
    moveModelOutputFolder <- paste("mv", modelOutputFolder, NewModelOutputFolder,
        sep = " ")
    system(moveModelOutputFolder)

    return(S)
}
