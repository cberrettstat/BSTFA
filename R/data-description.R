#' Utah Minimum Temperatures
#'
#' Zero-centered daily minimum temperatures averaged across 30-day windows from 1919 to 2014 across the US state of Utah; also includes dates, coordinates, and station names.
#'
#' @format A data set with 4 variables:
#' \describe{
#'   \item{TemperatureVals}{A 1251 by 146 matrix of zero-centered 30-day average daily minimum temperatures from 1912 through 2014. Missing observations are denoted using \code{NA}.}
#'   \item{Dates}{A vector of length 1251 and class \code{Date} providing the day of each observation in the same order as the rows of \code{TemperatureVals}.  Note that this package expects observation times to be regularly spaced.}
#'   \item{Coords}{A 142 by 2 data frame containing the longitude (first variable) and latitude (second variable) of measured locations.}
#'   \item{Locations}{A character vector of length 146 containing the station names for measured locations in the same order as the columns of \code{TemperatureVals} and the rows of \code{Coords}.}
#' }
#' @source \url{https://climate.usu.edu}
"utahDataList"


#' Output of BSTFA evaluated on a subset of utahDataList
#'
#' List object named \code{out.sm} containing the output from running the BSTFA function provided in the example code below using a subset of the \code{utahDataList}.
#'
#' @format See \code{help(BSTFA)} for details of what is included as output from the BSTFA function.
#' @examples
#' data(BSTFAoutput)
#' \dontrun{
#' #Code used to obtain this output
#' data("utahDataList")
#' attach(utahDataList)
#' set.seed(466)
#' dates.ind <- 1191:1251
#' locs.use <- c(4, 21, 22, 23, 32, 33, 36,
#'               40, 42, 48, 66, 78, 85, 89,
#'               94, 96, 114, 118, 124, 144)
#' temps.sm <- TemperatureVals[dates.ind, locs.use]
#' coords.sm <- Coords[locs.use,]
#' dates.sm <- Dates[dates.ind]
#' locsm.names <- Locations[locs.use]
#' out.sm <- BSTFA(ymat=temps.sm,
#'      dates=dates.sm,
#'      coords=coords.sm,
#'      iters=200000,
#'      burn=10001,
#'      thin=1826,
#'      save.missing=F,
#'      save.output=T,
#'      factors.fixed=c(20, 14, 11, 12))
#' }
"out.sm"
