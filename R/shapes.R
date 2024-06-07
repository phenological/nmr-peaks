#' Calculates a Gaussian bell centered at mean, with max height and a given  FWHM
#' 
#' @param xx The range of x values to compute the function (domain)
#' @param mean The center of the gaussian
#' @param max The max height of the curve
#' @param fwhm The full width at half maximum of the shape. It is not the same as the SD
#' @returns A numeric vector of length equal to length(xx)
#' @export
fgaussian <- function(xx, mean, max, fwhm) {
  max * exp(-4*log(2)*(((xx - mean) / fwhm)^2))
}

#' Calculates a Lorentzian bell centered at mean, with max height and a given  FWHM
#' 
#' @param xx The range of x values to compute the function (domain)
#' @param mean The center of the Lorentzian
#' @param max The max height of the curve
#' @param fwhm The full width at half maximum of the shape. It is not the same as the SD
#' @returns A numeric vector of length equal to length(xx)
#' @export
florentzian <- function(xx, mean, max, fwhm) {
  gamma2 <- fwhm/2
  gamma2 <- gamma2^ 2
  max * gamma2/((xx - mean)^2 + gamma2)
}

#' Calculates a pseudoVoigt bell centered at mean, with max height and a given  FWHM. 
#' A pseudoVoigt is an approximation of the Voigt function(Convolution of a Gaussian and a Lorentzian)
#' as a linear combination of the 2 aforementioned functions.
#' 
#' @param xx The range of x values to compute the function (domain)
#' @param mu The linear combination factor:  `mu*lorentzian + (1-mu)*gaussian`
#' @param mean The center of the Voigt
#' @param max The max height of the curve
#' @param fwhm The full width at half maximum of the shape. It is not the same as the SD
#' @returns A numeric vector of length equal to length(xx)
#' @export
fpseudoVoigt <- function(xx, mu, mean, max, fwhm) {
  return ((1 - mu) * fgaussian(xx, mean, max, fwhm) + mu * florentzian(xx, mean, max, fwhm))
}
