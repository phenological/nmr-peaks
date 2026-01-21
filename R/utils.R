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

#' Simulate a signal from resonance parameters (c.shift and j coupling). No second order effects.
#' 
#' @param cshift, numeric, chemical shift. Compulsory argument.
#' @param j, numeric, coupling constants (Hz). The default NULL produces a singlet
#' @param multiplicity, character. Use letters s,d,t,q,p,x,h,o,n for
#' singlet, doublet, triplet, etc. If omitted will be estimated from j, see details
#' @param frequency, numeric, default 600 (MHz)
#' @param linewidth, numeric, specified as full-width at half-maximum. Default 1 (Hz)
#' @param mu, numeric, gaussian-to-lorentzian ratio in the pseudo-voigt
#' combination. Default 0.85
#' @param intensity, numeric, intensity of the signal (summit). Default 1
#' @param nbAtoms, numeric, number of resonating atoms. Default 1
#' @param fwhm, numeric, optional. Quality of life alias for \code{linewidth}.
#' The corresponding NMRSignalModel1D parameter in this library is called fwhm,
#'  which may confuse veterans. This parameter allows them to auto-pilot when
#'  calling this function.
#' @details This is rather a class constructor on a different signature; I could not
#' figure out the S4 way of doing it.
#' 
#' The signal can be specified following customary NMR practice i.e. passing
#'  a coupling constant for each group of magnetically equivalent nuclei, 
#'  along with the corresponding multiplicities. Then, the function also accepts 
#'  more unusual method of giving a j for each individual proton. For instance, 
#' \code{j=c(3,2),multiplicity="td"} produces the same outcome as \code{j=c(3,3,2)}.
#'  
#' Be warned: internally, the function runs very few checks on the arguments 
#'  and does its best to produce a reasonable output.
#' 
#' 
#' @return a \code{\linkS4class{NMRSignal1D}}
#' @export
simulateSignal <- function(cshift, j=NULL, multiplicity = NULL,frequency=600
                           ,linewidth=1, mu=0.85, intensity=1, nbAtoms=1
                           ,fwhm){
  
  #Arg checks
  ##If no j, only "s" is an acceptable multiplicity
  if (is.null(j) & !is.null(multiplicity)){
    if (multiplicity != "s"){
      cat(crayon::red("nmrPeaks::simulateSignal >> missing coupling constant(s) j"))
      stop("Abort")   
    }
  }
  #multiplicity string dictionary
  multDict <- c("s","d","t","q","p","x","h","o","n")
  
  #Write cshift and height to output peak coordinates vector
  #If there is coupling it will be overwritten further down
  px <- cshift
  py <- intensity
  
  #compute linewidth in ppm
  if(missing(fwhm)){
    fwhm <- linewidth / frequency
  } else{
    fwhm <- fwhm / frequency
  }
  
  #parse multiplicity string to numbers
  ##If not provided, compute from j and make multiplicity string
  if (is.null(multiplicity)){
    if(length(j) > 0){
      m <- unname(table(j) + 1)
    } else multiplicity <- "s"
  } else{
    #If provided, parse string
    m <- unname(sapply(strsplit(multiplicity,"")[[1]], function(s) which(multDict == s)))
  }
  
  #Split multiplet
  if(length(j) > 0){
    #canonize
    ##unique j. This can be inefficient but ensures consistency
    j <- unique(j)
    ##decreasing j order. Cosmetic only?
    o <- order(j,decreasing = TRUE)
    j <- j[o]
    m <- m[o]
    
    #Remake multiplicity string. Inefficient but ensures consistency
    multiplicity <- paste0(multDict[m],collapse="")
    
    #Expand multiplicities
    j <- unlist(mapply(function(j,m) rep(j,m-1), j, m,SIMPLIFY=F))
    
    #Compute peak chemical shifts
    copulate <- function(x,y){
      c(sapply(x,function(x) c(x+y/2,x-y/2)))
    }
    px <- Reduce(copulate,j/frequency,init=cshift)
    
    #Group and compute peak height
    px <-table(px)
    py <- unname(px)
    px <- as.numeric(names(px))
    
    #scale height
    py <- py / max(py)
    py <- py * intensity
  }
  
  
  
  #Create NMRSignal1D
  new("NMRSignal1D"
      ,nbAtoms = nbAtoms
      ,multiplicity = multiplicity
      ,chemicalShift = cshift
      ,shape = list(name="pseudoVoigt", params=list(mu=mu,fwhm=fwhm))
      ,peaks = lapply(1:length(px), function(j){
        new("NMRPeak1D", x=px[j], y=py[j])
      })
  )
}

