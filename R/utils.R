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
#' @details This function intends to provide a pragmatic answer to the question
#' of where does a signal "start" and "end".
#' @returns numeric, the chemical shift coordinates of the boundaries of the signal's domain
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

#' Change the chemical shift of a \code{\linkS4class{NMRSignal1D}}
#' 
#' @param signal, \code{\linkS4class{NMRSignal1D}}
#' @param to, numeric, optional, new chemical shift coordinate
#' @param by, numeric, optional, shift to apply
#' @param hertz, logic, is \code{by} given in Hz (TRUE), or in chemical shift units
#' (FALSE, defauilt)?
#' @param SF, spectrometer field, in MHz. Only used if \code{hertz=TRUE}. Default 600 
#' @details If \code{to} is provided, the signal is shifted to the given chemical shift.
#' Otherwise, if \code{by} is provided, the signal is shifted by the given ppm.
#' Else, the signal is left unshifted
#' @return \code{\linkS4class{NMRSignal1D}}, the shifted signal
#' @export
shiftSignal <- function(signal, to, by=0, hertz=FALSE, SF=600){
  if (!missing(to)) by <- to[1] - signal@chemicalShift
  #Note how careles I am about combining "to" with "hertz=TRUE", because I am
  #not sure what would be the appropriate output anyway
  if (hertz) by <- by[1] / SF
  signal@peaks <- lapply(signal@peaks, function(p){
    p@x <- p@x + by[1]
    return(p)
  })
  signal@chemicalShift <- signal@chemicalShift + by
  return(signal)
}

#' Scale the intensity of a \code{\linkS4class{NMRSignal1D}}
#' 
#' @param signal, \code{\linkS4class{NMRSignal1D}}
#' @param to, numeric, optional, new maximum intensity
#' @param by, optional, height scaling factor
#' @details If \code{to} is provided, the signal is scaled so that its maximum
#' intensity matches the parameter. Otherwise, if \code{by} is provided, the 
#' signal is scaled by the given factor. Else, the signal is left unscaled
#' @return \code{\linkS4class{NMRSignal1D}}, the scaled signal
#' @export
scaleSignal <- function(signal, to, by=1){
  if (!missing(to)) by <- to / max(sapply(signal@peaks, function(p) p@y))
  signal@peaks <- lapply(signal@peaks, function(p){
    p@y <- p@y * by
    return(p)
  })
  return(signal)
}