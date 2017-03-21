\name{CyData-getset}

\alias{nmarkers}
\alias{nmarkers,CyData-method}
\alias{ncells}
\alias{ncells,CyData-method}

\alias{markernames}
\alias{markernames,CyData-method}
\alias{markernames<-}
\alias{markernames<-,CyData-method}
\alias{markernames<-,CyData,ANY-method}
\alias{sampleNames}
\alias{sampleNames,CyData-method}
\alias{sampleNames<-}
\alias{sampleNames<-,CyData-method}
\alias{sampleNames<-,CyData,ANY-method}

\alias{markerData}
\alias{markerData,CyData-method}
\alias{markerData<-}
\alias{markerData<-,CyData-method}
\alias{intensities}
\alias{intensities,CyData-method}
\alias{intensities<-}
\alias{intensities<-,CyData-method}
\alias{cellAssignments}
\alias{cellAssignments,CyData-method}
\alias{cellAssignments<-}
\alias{cellAssignments<-,CyData-method}

\alias{cellIntensities}
\alias{cellIntensities,CyData-method}
\alias{cellIntensities<-}
\alias{cellIntensities<-,CyData-method}
\alias{cellData}
\alias{cellData,CyData-method}
\alias{cellData<-}
\alias{cellData<-,CyData-method}

\title{CyData getters and setters}
\description{Methods to get and set data slots in the CyData object.}

\usage{
\S4method{markerData}{CyData}(x)
\S4method{markerData}{CyData}(x) <- value

\S4method{intensities}{CyData}(x)
\S4method{intensities}{CyData}(x) <- value

\S4method{cellAssignments}{CyData}(x)
\S4method{cellAssignments}{CyData}(x) <- value

\S4method{cellIntensities}{CyData}(x)
\S4method{cellIntensities}{CyData}(x) <- value

\S4method{cellData}{CyData}(x)
\S4method{cellData}{CyData}(x) <- value

\S4method{nmarkers}{CyData}(x)
\S4method{ncells}{CyData}(x)

\S4method{markernames}{CyData}(object)
\S4method{markernames}{CyData}(object) <- value

\S4method{sampleNames}{CyData}(object)
\S4method{sampleNames}{CyData}(object) <- value
}

\arguments{
\item{x, object}{A CyData object.}
\item{value}{
For \code{markernames} and \code{sampleNames}, a character vector specifying the names of the markers and samples, respectively.
Otherwise, an appropriate object to use to replace the corresponding slot of \code{x}.
}
}

\details{
\code{markerData}, \code{intensities}, \code{cellAssignments}, \code{cellIntensities} and \code{cellData} are methods to get or set the slots in a CyData object.
general, setting should not be performed unless adding annotation for particular markers or cells (or if you really know what you're doing). 
Caution is especially advised when using \code{cellIntensities<-} and/or \code{cellData<-}, as several methods in this package rely on a particular ordering of cells.

Several convenience functions are also provided.
\code{nmarkers} will return the number of markers, while \code{ncells} will return the number of cells.
\code{markernames} will return the marker names (i.e., the row names of \code{markerData}), while \code{sampleNames} will return the sample names (i.e., the column names of \code{object}).
The respective replacement functions will modify these names.
}

\value{
For the setters, a new CyData object will be generated with modified entries for the appropriate slots.
For the getters, an object of the class corresponding to each slot will be returned.
\code{nmarkers} and \code{ncells} will return integer scalars, while \code{markerNames} and \code{sampleNames} will return character vectors.
}

\author{
Aaron Lun
}

\examples{
example(CyData) # Mocking up an object.

nmarkers(cyd)
markernames(cyd)
sampleNames(cyd)

markerData(cyd)
markerData(cyd)$stuff <- runif(nmarkers(cyd))
markerData(cyd)

head(intensities(cyd))
cellIntensities(cyd)[,1:10]
}
