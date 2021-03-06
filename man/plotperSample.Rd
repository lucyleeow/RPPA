% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barplots.R
\name{plotperSample}
\alias{plotperSample}
\alias{plotperAB}
\title{Barplots per sample/antibody}
\usage{
plotperSample(tidydf, RFIcol = "RFI", pdfoutput = TRUE)

plotperAB(tidydf, RFIcol = "RFI", pdfoutput = TRUE)
}
\arguments{
\item{tidydf}{Tidy dataframe of RFI values to plot.}

\item{RFIcol}{Column name of the RFI column as string.}

\item{pdfoutput}{Single logical indicating whether the plots should be
output to pdf file or the current graphics device.}
}
\description{
Creates separate barplots of the RFI values for EACH sample or antibody and 
writes to either pdf or the current graphics device.
}
\section{Functions}{
\itemize{
\item \code{plotperSample}: Create a barplot of the RFI values of each 
antibody. One barplot is created for each sample. Write to a file called 
'plotPerSample.pdf' or the current graphics device.

\item \code{plotperAB}: Creates a barplot of the RFI values of each 
sample. One barplot is created for each antibody. Writes to either a pdf called 
'plotPerAB.pdf' or the current graphics device.
}}

