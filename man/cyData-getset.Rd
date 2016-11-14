\name{cyData-getset}

\alias{nmarkers}
\alias{nmarkers,cyData-method}
\alias{ncells}
\alias{ncells,cyData-method}

\alias{markerData}
\alias{markerData,cyData-method}
\alias{markerData<-}
\alias{markerData<-,cyData-method}
\alias{medianIntensities}
\alias{medianIntensities,cyData-method}
\alias{medianIntensities<-}
\alias{medianIntensities<-,cyData-method}
\alias{cellIntensities}
\alias{cellIntensities,cyData-method}
\alias{cellIntensities<-}
\alias{cellIntensities<-,cyData-method}

\title{cyData getters and setters}
\description{Methods to get and set data slots in the cyData object.}

\usage{
\S4method{markerData}{cyData}(x)
\S4method{markerData}{cyData}(x) <- value

\S4method{medianIntensities}{cyData}(x)
\S4method{medianIntensities}{cyData}(x) <- value

\S4method{cellIntensities}{cyData}(x)
\S4method{cellIntensities}{cyData}(x) <- value

\S4method{nmarkers}{cyData}(x)
\S4method{ncells}{cyData}(x)
}

\arguments{
\item{x}{A cyData object.}
\item{value}{An object to use to replace the appropriate field in \code{x}}
}

\details{
These are methods to get or set the slots in a cyData object. 
Some caution is advised when using \code{cellIntensities}, as several methods in this package rely on a particular ordering of cells in this slot.
}

\value{
\code{nmarkers} will return the number of markers, while \code{ncells} will return the number of cells.
Each other method will get or set the named slot.
}

\author{
Aaron Lun
}

\examples{
example(cyData) # Mocking up an object.

nmarkers(cyd)
markerData(cyd)
markerData(cyd)$stuff <- runif(nmarkers(cyd))
markerData(cyd)

head(medianIntensities(cyd))
cellIntensities(cyd)[,1:10]
}
