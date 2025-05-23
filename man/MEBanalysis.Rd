% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MEBanalysis.R
\name{MEBanalysis}
\alias{MEBanalysis}
\title{Analysis of ME-designs in Block Set-up}
\usage{
MEBanalysis(data)
}
\arguments{
\item{data}{Columns of dataset should be in order of block, line1,line2, cross number and response.}
}
\value{
Returns the ANOVA table of gca and sca effect analysis.
}
\description{
For a given field data it will provide analysis result through ANOVA table including gca and sca effect analysis.
}
\examples{
library(MEDesigns)
MEBanalysis(MEdata)
}
