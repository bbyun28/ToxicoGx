% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloadTSet.R
\name{availableTSets}
\alias{availableTSets}
\title{Return a table of ToxicoSets available for download}
\usage{
availableTSets(
  saveDir = tempdir(),
  myfn = "availableToxicoSets.csv",
  verbose = TRUE
)
}
\arguments{
\item{saveDir}{\code{character} Directory to save the table of tSets}

\item{myfn}{\code{character} The filename for the table of tSets}

\item{verbose}{\code{bool} Should status messages be printed during download.}
}
\value{
A data.frame with details about the available ToxicoSet objects
}
\description{
The function fetches a table of all ToxicoSets available for download from
the PharmacoGx server. The table includes the names of the PharamcoSet, the
types of data available in the object, and the date of last update.
}
\examples{
if (interactive()){
availabletSets()
}

}
