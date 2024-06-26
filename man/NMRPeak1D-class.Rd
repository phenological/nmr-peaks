% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMRClassDefinitions.R
\docType{class}
\name{NMRPeak1D-class}
\alias{NMRPeak1D-class}
\title{An S4 class for a 1D NMR peak}
\description{
An S4 class for a 1D NMR peak
}
\details{
A peak is a NMR feature with a unique maximum. This library models
such features as pseudoVoigt profiles, but this class makes no assumptions
about the peak profile other than it can be expressed in terms of ordinate,
height and width.
}
\section{Slots}{

\describe{
\item{\code{x}}{chemical shift of the peak}

\item{\code{y}}{height of the peak}

\item{\code{fwhm}}{Full width at half maximum of peak}

\item{\code{shape}}{(optional) A peak shape for the different peaks conforming the signal. It uses the parent shape if void}

\item{\code{type}}{The Class Name. Used for moving between S4 and JSON}
}}

