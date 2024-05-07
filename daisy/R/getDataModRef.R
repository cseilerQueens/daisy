################################################################################
#' Data frame of model with model output and reference data
#' @description This function extracts the data from (i) a model run conducted for a single grid cell and (ii) multiple globally gridded reference data sets.
#' The data frame only includes the time steps that all data sets have in common.
#' 
#' @param nc.mod A string that gives the path and name of the netcdf file that contains the model output, e.g. '/home/model_gpp.nc'
#' @param nc.ref.list A list of paths and names that point to the reference data netcdf files e.g. list('/home/ref01_gpp.nc', '/home/ref02_gpp.nc')
#' @param ref.id.list A list of reference data set names, e.g. list('MODIS', 'FLUXCOM')
#' @param unit.conv.mod A number that is used as a factor to convert the unit of the model data, e.g. 86400
#' @param unit.conv.ref An R object that contains the conversion factors of the reference data e.g. c(1, 86400)
#' @param variable.unit A string that gives the final units using LaTeX notation, e.g. 'gC m$^{-2}$ day$^{-1}$'
#' @param timeInt A string that gives the time interval of the model data, e.g. 'month' or 'year'
#'
#' @return A data frame that shows the model and reference values, as well as the dates, variable name, and units. The data frame serves as an input for running the cost function.
#'
#' @examples

#' library(foreach)
#' 
#' nc.mod <- '/home/cseiler/projects/def-cseiler-ab/cseiler/data-assimilation/simulations/transient_CRUJRAv2/outputFiles/gpp_monthly.nc'
#' nc.ref01 <- '/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63/gpp_GOSIF_128x64.nc'
#' nc.ref02 <- '/home/cseiler/projects/def-cseiler-ab/cseiler/reference_T63/GPP.ensembleMedian.CRUNCEPv6.monthly_128x64.nc'
#' nc.ref.list <- list(nc.ref01, nc.ref02)
#' ref.id.list <- list('GOSIF', 'FLUXCOM')

#' unit.conv.mod <- 86400 * 1000  # optional unit conversion for model data
#' unit.conv.ref <- c(1, 1)
#' variable.unit <- 'gC m$^{-2}$ day$^{-1}$'  # unit after conversion (LaTeX notation)
#' getDataModRef(nc.mod, nc.ref.list, ref.id.list, unit.conv.mod, unit.conv.ref, variable.unit, timeInt = 'month')
#'
#' @export

getDataModRef <- function(nc.mod, nc.ref.list, ref.id.list, unit.conv.mod, unit.conv.ref,
    variable.unit, timeInt = "month") {

    # Model output from single grid cell

    nc <- nc.mod
    nc <- ncdf4::nc_open(nc)
    variable.name <- names(nc[["var"]])
    variable.name <- variable.name[length(variable.name)]  # take the last variable
    longName <- ncdf4::ncatt_get(nc, variable.name)$long_name
    data <- ncdf4::ncvar_get(nc)  # get values
    lon <- ncdf4::ncvar_get(nc, "longitude")
    lat <- ncdf4::ncvar_get(nc, "latitude")
    dates <- ncdf4.helpers::nc.get.time.series(nc)
    dates <- format(dates, "%Y-%m-15")
    start.date <- min(dates)
    end.date <- max(dates)
    start.date.mod <- start.date
    end.date.mod <- end.date

    mod <- data
    mod <- mod * unit.conv.mod
    names(mod) <- dates

    # Get the dates from all reference data sets

    dates.ref <- foreach::foreach(i = 1:length(nc.ref.list)) %do% {
        nc.ref <- nc.ref.list[[i]]
        nc <- ncdf4::nc_open(nc.ref)
        dates <- ncdf4.helpers::nc.get.time.series(nc)
        dates <- format(dates, "%Y-%m-15")
        start.date <- min(dates)
        end.date <- max(dates)
        start.end.dates <- data.frame(start.date, end.date)
        return(start.end.dates)
    }
    start.end.dates <- do.call(rbind, dates.ref)
    start.date.ref <- max(start.end.dates$start.date)
    end.date.ref <- min(start.end.dates$end.date)

    # find common time period for model output and all reference data sets
    start.date <- max(start.date.mod, start.date.ref)
    end.date <- min(end.date.mod, end.date.ref)

    # Subset model data for common time period
    mod <- mod[which(names(mod) >= start.date & names(mod) <= end.date)]

    # Convert longitude from 0/360 to -180/180:
    lon <- ifelse(lon > 180, lon - 360, lon)

    # Create spatial point to extract values from gridded reference data
    regular <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
    points <- data.frame(lon, lat, 1)
    sp::coordinates(points) <- ~lon + lat
    raster::projection(points) <- regular

    # Subset reference data for common time period and extract values

    ref.data <- foreach::foreach(i = 1:length(nc.ref.list)) %do% {
        nc.ref <- nc.ref.list[[i]]
        nc <- ncdf4::nc_open(nc.ref)
        dates <- ncdf4.helpers::nc.get.time.series(nc)
        dates <- format(dates, "%Y-%m-15")
        ref <- raster::brick(nc.ref)
        ref <- ref * unit.conv.ref[i]
        ref <- raster::setZ(ref, dates)
        names(ref) <- dates
        ref <- ref[[which(format(as.Date(raster::getZ(ref)), "%Y-%m-15") >= start.date &
            format(as.Date(raster::getZ(ref)), "%Y-%m-15") <= end.date)]]
        suppressWarnings(ref <- raster::rotate(ref))
        ref <- raster::extract(x = ref, y = points)
        ref <- t(ref)
        colnames(ref) <- ref.id.list[[i]]
        return(ref)
    }
    ref <- do.call(cbind, ref.data)

    dates <- seq.Date(from = as.Date(start.date), to = as.Date(end.date), by = timeInt)

    # mod and ref differ by one month
    df <- data.frame(dates, mod, ref, variable.name, variable.unit, longName)

    return(df)

}

utils::globalVariables(c("%do%", "i"))
