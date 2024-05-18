################################################################################
#' Normalize parameter values
#' @description This function noamrilszes parameter values so that they range between zero and one.
#' 
#' @param x parameter Values
#' @return A list with the mornalized parameter values, the range, and the minimum value
#'
#' @examples
#' intFun.normalize(x=5,u=6,l=2)
#' @export

intFun.normalize <- function(x, u, l) {
  # x = value, u = upper bound, l = lower bound
  # normalize x, such that x ranges between 0 and 1:
  xn <- (x - l) / (u - l)
  return(xn)
}

#-----------------------------------------------------------
#' Convert normalized to original parameter values
#' @description This function reverses the normalization
#' 
#' @param normValues normalized values
#' @param range Range of original parameterValues
#' @param minValue Minimum value of original parameter values
#' @return original parameter values
#'
#' @examples
#' x <- -100:100
#' normalization <- int.fun.normalize(x)
#' intFun.unnormalize(normValues = normalization[[1]], range = normalization[[2]], minValue = normalization[[3]])
#' 
#' @export
intFun.unnormalize <- function(normValues, upperValues, lowerValues) {
    x <- normValues * (upperValues - lowerValues) + lowerValues
    return(x)
}

#-----------------------------------------------------------

#-------------------------------------------------------------------------------
# intFun.getZ
#' Get the time dimension of a NetCDF file
#' @description This function tries to obtain the dates of a NetCDf file with
#' monthly data. The third attempt used the nc.get.time.series function from the
#' ncdf4.helpers package. NetCDf files differ among modeling groups. This
#' function accounts for some of those differences.
#' @param nc.mod A string that gives the path and name of the netcdf file that contains the model output, e.g. '/home/model_gpp.nc'
#' @param min.date Earliest possible date, e.g. '1970-01-01'.
#' @param max.date Latest possible date, e.g. '1970-01-01'. Dates outside min.date and max.date are assumed to be wrong.
#' @return Dates, e.g. '1980-01-15'
#' @examples
#'
#' library(raster)
#' library(ncdf4)
#' library(ncdf4.helpers)
#'
#' nc.mod <- system.file('extdata/modelRegular', 'gpp_monthly.nc', package = 'amber')
#' intFun.getZ(nc.mod)
#'
#' @keywords internal
#' @export
intFun.getZ <- function(nc.mod, min.date = "1000-01-01", max.date = "9999-01-01") {

    # print('vobjtovarid4: error #F message may be safely ignored.')

    # function that verifies date format:
    is.Date <- function(x) {
        inherits(x, c("Date", "POSIXt"))
    }

    # function that evaluates date range:
    is.Date.reasonable <- function(x) {
        if (x < as.Date(min.date) | x > as.Date(max.date)) {
            return("BadDates")
        } else {
            "GoodDates"
        }
    }

    # Initial value:
    dates.reasonable <- "BadDates"

    # First try:
    f1 <- ncdf4::nc_open(nc.mod)
    dates.mod <- try(ncdf4.helpers::nc.get.time.series(f1), silent = TRUE)
    dates.mod <- try(format(dates.mod, "%Y-%m-15"), silent = TRUE)
    dates.mod <- try(as.Date(dates.mod), silent = TRUE)
    ncdf4::nc_close(f1)

    # check whether dates look good
    dates.TrueFalse <- is.Date(try(as.Date(dates.mod), silent = TRUE))

    # In some cases dates appear reasonable, but some dates occur twice.

    # Those cases can be detected by looking at the length of the dates vector,
    # as done below.

    if (dates.TrueFalse == TRUE & length(dates.mod) > 1) {
        n <- length(dates.mod)
        deltaTime <- dates.mod[2] - dates.mod[1]
        if (deltaTime > 300) {
            timeInterval <- "year"
        } else {
            timeInterval <- "month"
        }
        seqDate <- seq.Date(from = dates.mod[1], to = dates.mod[n], by = timeInterval)
        lengthOne <- length(dates.mod)
        lengthTwo <- length(seqDate)

        if (lengthOne == lengthTwo) {
            dates.TrueFalse = TRUE
        } else {
            dates.TrueFalse = FALSE
        }
    }


    if (dates.TrueFalse == TRUE) {
        dates.reasonable <- is.Date.reasonable(try(as.Date(dates.mod[1]), silent = TRUE))
    }

    # Second Try:
    if (dates.TrueFalse == FALSE | dates.reasonable == "BadDates") {
        f1 <- ncdf4::nc_open(nc.mod)
        try(time <- ncdf4::ncvar_get(f1, "time"), silent = TRUE)
        try(time <- ncdf4::ncvar_get(f1, "time_counter"), silent = TRUE)
        nTime <- base::length(time)
        try(tunits <- ncdf4::ncatt_get(f1, "time", attname = "units"), silent = TRUE)
        try(tunits <- ncdf4::ncatt_get(f1, "time_counter", attname = "units"), silent = TRUE)
        ncdf4::nc_close(f1)

        tunits <- strsplit(tunits$value, " ")
        origin <- unlist(tunits)[3]
        # Not all models provide YYYY-MM-DD. Some provide only YYYY or YYYY-MM.
        # To add the days or months and days:
        if (nchar(origin) == 7) {
            origin <- paste(origin, "15", sep = "-")
        }
        if (nchar(origin) == 4) {
            origin <- paste(origin, "01-15", sep = "-")
        }

        timeInterval <- unlist(tunits)[1]

        if (grepl("day", timeInterval) == TRUE) {
            timeConversion <- 1
        }
        if (grepl("month", timeInterval) == TRUE) {
            timeConversion <- 365.25/12
        }
        if (grepl("year", timeInterval) == TRUE) {
            timeConversion <- 365.25
        }

        time <- time * timeConversion

        start.date <- as.Date(origin) + time[1]
        start.date <- format(as.Date(start.date), "%Y-%m-15")  # center of month

        dates.mod <- base::seq(as.Date(start.date), by = timeInterval, length = nTime)

        # check whether that was succesful:
        dates.TrueFalse <- is.Date(try(as.Date(dates.mod), silent = TRUE))
        dates.reasonable <- is.Date.reasonable(try(as.Date(dates.mod[1]), silent = TRUE))
    }

    if (dates.TrueFalse == FALSE | dates.reasonable == "BadDates") {
        print(paste("ERROR in reading time axis in ", nc.mod, sep = ""))
    }

    if (dates.TrueFalse == TRUE | dates.reasonable == "GoodDates") {
        print(paste("Time axis in ", nc.mod, " OK", sep = ""))
        return(dates.mod)
    }
}

#-------------------------------------------------------------------------------
# intFun.isRaster
#' Reproject coastline
#' @description This function assesses whether an R object is a raster. The original
#' code was copied from intFun.isRaster (spatial.tools_1.6.0)
#' @param x An R object such as a raster or a number
#' @examples
#'
#' x <- 1
#' intFun.isRaster(x)
#' y <- raster::raster(matrix(seq(1,10), ncol = 2))
#' intFun.isRaster(y)
#'
#' @keywords internal
#' @export
intFun.isRaster <- function(x) {
    return((class(x) == "RasterLayer" || class(x) == "RasterBrick" || class(x) ==
        "RasterStack"))
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# intFun.rmse
#' Root mean square error (RMSE)
#' @description This function computes the root mean square error (RMSE), which is defined as:
#'
#' \eqn{$rmse(\lambda, \phi)=\sqrt{\frac{1}{t_{f}-t_{0}}\int_{t_{0}}^{t_{f}}(v_{mod}(t,\lambda, \phi)-v_{ref}(t,\lambda, \phi))^{2}dt}$}
#'
#' where \eqn{\lambda} is the longitude, \eqn{\phi} is the latitude, \eqn{t}
#' is the time, \eqn{t_0} is the initial time step, \eqn{t_f} is the final time
#' time step, \eqn{v_{mod}} is a modelled variable and \eqn{v_{ref}} is the
#' corresponding reference variable.
#'
#' @param mod An R object (model output data)
#' @param ref An R object (reference data)
#' @return An R object that gives the root mean square error when comparing
#' \code{mod} against \code{ref}.
#' @examples
#'
#' library(raster)
#' # create two raster stacks
#' for(i in 1:100)
#' {
#'  mod <- raster::raster(matrix(runif(100,-1,1), ncol=10))
#'  ref <- raster::raster(matrix(runif(100,-2,2), ncol=10))
#'  assign(paste('mod', i , sep='_'), mod)
#'  assign(paste('ref', i , sep='_'), ref)
#' }
#' my.list.mod <- lapply(ls(pattern='mod_'), get)
#' my.list.ref <- lapply(ls(pattern='ref_'), get)
#' mod <- do.call(stack, my.list.mod)
#' ref <- do.call(stack, my.list.ref)
#' # compute RMSE
#' rmse <- intFun.rmse(mod,ref)
#' plot(rmse); text(rmse, digits=2)
#'
#' @keywords internal
#' @export
intFun.rmse <- function(mod, ref) {
    if (intFun.isRaster(ref) == TRUE)
        sqrt(raster::mean((mod - ref)^2, na.rm = TRUE)) else sqrt(mean((mod - ref)^2, na.rm = TRUE))
}

#-------------------------------------------------------------------------------

# intFun.anom
#' Anomalies
#' @description This function computes anomalies.
#' @param x An R object
#' @return An R object that gives the anomalies
#' @examples
#'
#' # make some data
#' month <- seq(1,12,1)
#' month <- rep(month,10)
#' mod <- runif(length(month), 0,100)
#' mod <- data.frame(month, mod)
#' # compute anomalies
#' mod.anom <- data.frame(apply(mod[2], 2, intFun.anom))
#'
#' @keywords internal
#' @export
intFun.anom <- function(x) {
    if (intFun.isRaster(x) == TRUE)
        x - raster::mean(x, na.rm = TRUE) else x - mean(x, na.rm = TRUE)
}

#-------------------------------------------------------------------------------

# intFun.crmse
#' Centralized root mean square error (CRMSE)
#' @description This function computes the centralized root mean square error,
#' which is defined as:
#'
#' \eqn{$crmse(\lambda, \phi) = \sqrt{\frac{1}{t_{f}-t_{0}}\int_{t_{0}}^{t_{f}}[(v_{mod}(t,\lambda, \phi)-\overline{v_{mod}}(\lambda, \phi))-(v_{ref}(t,\lambda, \phi)-\overline{v_{ref}}(\lambda, \phi))]^{2}dt}$}
#'
#' where \eqn{\lambda} is the longitude, \eqn{\phi} is the latitude, \eqn{t}
#' is the time, \eqn{t_0} is the initial time step, \eqn{t_f} is the final time
#' time step, \eqn{v_{mod}} is a modelled variable, \eqn{v_{ref}} is the
#' corresponding reference variable, \eqn{\overline{v_{mod}}} is the time-mean
#' modelled variable, and \eqn{\overline{v_{ref}}} is the time-mean reference
#' variable.
#' @param mod.anom An R object (e.g. monthly anomalies from model output)
#' @param ref.anom An R object (e.g. monthly anomalies from reference data)
#' @return An R object that shows the centralized root mean square error
#' @examples
#'
#' library(raster)
#' # create two raster stacks
#' for(i in 1:100)
#' {
#'  mod <- raster::raster(matrix(runif(100,0,10), ncol=10))
#'  ref <- raster::raster(matrix(runif(100,0,10), ncol=10))
#'  assign(paste('mod', i , sep='_'), mod)
#'  assign(paste('ref', i , sep='_'), ref)
#' }
#' my.list.mod <- lapply(ls(pattern='mod_'), get)
#' my.list.ref <- lapply(ls(pattern='ref_'), get)
#' mod <- do.call(stack, my.list.mod)
#' ref <- do.call(stack, my.list.ref)
#' # compute anomalies
#' mod.anom <- intFun.anom(mod)
#' ref.anom <- intFun.anom(ref)
#' # compute CRMSE
#' crmse <- intFun.crmse(mod.anom, ref.anom)
#' plot(crmse); text(crmse, digits=2)
#'
#' @keywords internal
#' @export
intFun.crmse <- function(mod.anom, ref.anom) {
    if (intFun.isRaster(ref.anom) == TRUE)
        sqrt(raster::mean((mod.anom - ref.anom)^2, na.rm = TRUE)) else sqrt(mean((mod.anom - ref.anom)^2, na.rm = TRUE))
}

#-------------------------------------------------------------------------------

# intFun.theta
#' Time distance in months
#' @description Consider two rasters that show the month of the maximum value of
#' a variable during the climatological mean annual cycle. Calculating the
#' absolute difference between both rasters may yield values that range from
#' 0 to 11 (type a). A more meaningful range that considers the circular nature
#' of the annual cycle would range from 0 to 6 (type b). This function converts
#' a raster of type (a) to type (b).
#' @param x A raster object of type (a)
#' @return A raster object of type (b)
#' @examples
#'
#' library(raster)
#' # create a raster object with NA
#' data <- runif(100,0,11)
#' data <- matrix(data, ncol=10)
#' data <- raster::raster(data)
#' # replace NA with zero in raster object
#' data <- calc(data,intFun.theta)
#' plot(data); text(data)
#'
#' @keywords internal
#' @export
intFun.theta <- function(x) {
    ifelse(x <= 6, x, (12 - x))
}

#-------------------------------------------------------------------------------

# intFun.anom.mly
#' Monthly anomalies
#' @description This function computes monthly anomalies.
#' @param mly An R object with a monthly time series
#' @param clim An R object with climatological monthly mean data
#' @return A raster object with monthly anomalies
#' @examples
#'
#' # make some data
#' month <- seq(1,12,1)
#' month <- rep(month,10)
#' mod <- runif(length(month), 0,10)
#' mod <- data.frame(month, mod)
#' # make an index
#' index <- list(mod$month)
#' clim.mly <- apply(mod, 2, function(x) {tapply(x, index, mean, na.rm=TRUE)})
#' mod.anom <- intFun.anom.mly(mod, clim.mly)
#'
#' @keywords internal
#' @export
intFun.anom.mly <- function(mly, clim) {
    data <- merge(mly, clim, by = "month")
    n <- (ncol(data) - 1)/2 + 1
    mod <- data[2:n]
    clim <- data[(n + 1):ncol(data)]
    anom <- mod - clim
    return(anom)
}

#-------------------------------------------------------------------------------

# intFun.iav
#' Inter-annual variability
#' @description This function computes the inter-annual variability. All data must start in January. All months after the last Dec will be dropped if data does not end in December.
#' @param anom R object with monthly anomalies
#' @return R object with inter-annual variability
#' @examples
#'
#' # make some data
#' month <- seq(1,12,1)
#' month <- rep(month,10)
#' mod <- runif(length(month), 0,100)
#' mod <- data.frame(month, mod)
#' # compute climatological mean monthly values
#' index <- list(mod$month)
#' mod.clim.mly <- apply(mod, 2, function(x) {tapply(x, index, mean, na.rm=TRUE)})
# compute anomalies
#' mod.anom <- intFun.anom.mly(mod, mod.clim.mly)
# compute monthly inter-annual variability
#' mod.iav <- intFun.iav(mod.anom)
#'
#' @keywords internal
#' @export
intFun.iav <- function(anom) {
    sqrt(mean((anom)^2, na.rm = TRUE))
}

