\name{cyData-getset}

\alias{nmarkers}
\alias{nmarkers,cyData-method}
\alias{ncells}
\alias{ncells,cyData-method}

\alias{markernames}
\alias{markernames,cyData-method}
\alias{markernames<-}
\alias{markernames<-,cyData-method}
\alias{markernames<-,cyData,ANY-method}
\alias{sampleNames}
\alias{sampleNames,cyData-method}
\alias{sampleNames<-}
\alias{sampleNames<-,cyData-method}
\alias{sampleNames<-,cyData,ANY-method}

\alias{markerData}
\alias{markerData,cyData-method}
\alias{markerData<-}
\alias{markerData<-,cyData-method}
\alias{intensities}
\alias{intensities,cyData-method}
\alias{intensities<-}
\alias{intensities<-,cyData-method}
\alias{cellAssignments}
\alias{cellAssignments,cyData-method}
\alias{cellAssignments<-}
\alias{cellAssignments<-,cyData-method}

\alias{cellIntensities}
\alias{cellIntensities,cyData-method}
\alias{cellIntensities<-}
\alias{cellIntensities<-,cyData-method}
\alias{cellData}
\alias{cellData,cyData-method}
\alias{cellData<-}
\alias{cellData<-,cyData-method}

\title{cyData getters and setters}
\description{Methods to get and set data slots in the cyData object.}

\usage{
\S4method{markerData}{cyData}(x)
\S4method{markerData}{cyData}(x) <- value

\S4method{intensities}{cyData}(x)
\S4method{intensities}{cyData}(x) <- value

\S4method{cellAssignments}{cyData}(x)
\S4method{cellAssignments}{cyData}(x) <- value

\S4method{cellIntensities}{cyData}(x)
\S4method{cellIntensities}{cyData}(x) <- value

\S4method{cellData}{cyData}(x)
\S4method{cellData}{cyData}(x) <- value

\S4method{nmarkers}{cyData}(x)
\S4method{ncells}{cyData}(x)

\S4method{markernames}{cyData}(object)
\S4method{markernames}{cyData}(object) <- value

\S4method{sampleNames}{cyData}(object)
\S4method{sampleNames}{cyData}(object) <- value
}

\arguments{
\item{x, object}{A cyData object.}
\item{value}{
For \code{markernames} and \code{sampleNames}, a character vector specifying the names of the markers and samples, respectively.
Otherwise, an appropriate object to use to replace the corresponding slot of \code{x}.
}
}

\details{
\code{markerData}, \code{intensities}, \code{cellAssignments}, \code{cellIntensities} and \code{cellData} are methods to get or set the slots in a cyData object.
general, setting should not be performed unless adding annotation for particular markers or cells (or if you really know what you're doing). 
Caution is especially advised when using \code{cellIntensities<-} and/or \code{cellData<-}, as several methods in this package rely on a particular ordering of cells.

Several convenience functions are also provided.
\code{nmarkers} will return the number of markers, while \code{ncells} will return the number of cells.
\code{markernames} will return the marker names (i.e., the row names of \code{markerData}), while \code{sampleNames} will return the sample names (i.e., the column names of \code{object}).
The respective replacement functions will modify these names.
}

\value{
For the setters, a new cyData object will be generated with modified entries for the appropriate slots.
For the getters, an object of the class corresponding to each slot will be returned.
\code{nmarkers} and \code{ncells} will return integer scalars, while \code{markerNames} and \code{sampleNames} will return character vectors.
}

\author{
Aaron Lun
}

\examples{
example(cyData) # Mocking up an object.

nmarkers(cyd)
markernames(cyd)
sampleNames(cyd)

markerData(cyd)
markerData(cyd)$stuff <- runif(nmarkers(cyd))
markerData(cyd)

head(intensities(cyd))
cellIntensities(cyd)[,1:10]
}
