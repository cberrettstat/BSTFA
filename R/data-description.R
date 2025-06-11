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


#' Output of BSTFA evaluated on utahDataList
#'
#' List object named \code{out} containing the output from running the BSTFA function provided in the example code below.
#'
#' @format See \code{help(BSTFA)} for details of what is included as output from the BSTFA function.
#' @examples
#' data(BSTFAoutput)
#' \dontrun{
#' #Code used to obtain this output
#' data("utahDataList")
#' attach(utahDataList)
#' out <- BSTFA(ymat=TemperatureVals,
#'   dates=Dates,
#'   coords=Coords,
#'   iters=100000,
#'   burn=10001,
#'   thin=300,
#'   save.missing=F,
#'   save.output=TRUE,
#'   factors.fixed=c(144,89,129,78))
#' }
"out"
