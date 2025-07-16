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
#' data(out.sm)
#' \dontrun{
#' #Code used to obtain this output
#' data("utahDataList")
#' attach(utahDataList)
#' dates.ind <- 1151:1251
#' locs.use <- c(3, 8, 11, 16, 17,
#'               20, 23, 29, 30, 46,
#'               47, 49, 60, 62, 66, 73,
#'               75, 76, 77, 78, 85, 89, 94,
#'               96, 98, 100, 109, 112,
#'               115, 121, 124, 128, 133, 144)
#' temps.sm <- TemperatureVals[dates.ind, locs.use]
#' coords.sm <- Coords[locs.use,]
#' dates.sm <- Dates[dates.ind]
#' locsm.names <- Locations[locs.use]
#' set.seed(466)
#' out.sm <- BSTFA(ymat=temps.sm, 
#'                 dates=dates.sm, 
#'                 coords=coords.sm, 
#'                 iters=5000, 
#'                 burn=1000,
#'                 thin=40, 
#'                 factors.fixed=c(14, 22, 15, 20), 
#'                 n.temp.bases=45, 
#'                 save.missing=FALSE)
#' }
"out.sm"
