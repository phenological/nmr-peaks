#' Interpolates a \code{\linkS4class{NMRSignal1D}} on the given set of points
#' 
#' @param signal, \code{\linkS4class{NMRSignal1D}}, signal to interpolate
#' @param ppm, chemical shift values to interpolate
#' @returns numeric vector of interpolated intensities
#' @export
signalToY <- function(signal,ppm){
  if("NMRSignal1D" %in% is(signal)){
    y <- 0
    fwhm <- signal@shape$params$fwhm
    mu <- signal@shape$params$mu
    for (p in signal@peaks)
      y <- y + fpseudoVoigt(ppm, fwhm=fwhm, mu=mu, mean=p@x, max=p@y)
    return(y)  
  }
  cat(crayon::red("nmrSpectraProcessing::signalToY >>",
                  "Argument is not a NMRSignal1D object\n"))
}

#' Calculate the domain of a \code{\linkS4class{NMRSignal1D}}
#' 
#' Estimates the domain of the signal as the coordinates of the outer peaks
#'  offset by \emph{n} x linewidth (fwhm)
#' @param signal, \code{\linkS4class{NMRSignal1D}}
#' @param n numeric, linewidth multiplier used to offset the outer peaks
#' @returns numeric, the extremes of the signal's domain
#' @export
signalDomain <- function(signal, n=3){
  if ("NMRSignal1D" %in% is(signal)){
    offset <- signal@shape$params$fwhm * n * c(-1,1)
    peaksRange <- range(sapply(signal@peaks, function(aPeak) aPeak@x))
    return(peaksRange + offset)
  }
  cat(crayon::red("signalDomain >>"
                  ,"Argument is not a S4 fusion::NMRSignal1D object"))
  stop()
}