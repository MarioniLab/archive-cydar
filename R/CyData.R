# This sets up the CyData reference class.

setClass("CyData", contains="SummarizedExperiment", 
         slots=c(markerData="DataFrame",
                 cellIntensities="matrix",
                 cellData="DataFrame",
                 cellAssignments="list",
                 intensities="matrix"))

setValidity2("CyData", function(object) {
    txt <- character(0)
    if (storage.mode(.raw_cellIntensities(object))!="double") {
        txt <- c(txt, "cellIntensities should be stored in double-precision")
    }
    if (storage.mode(.raw_intensities(object))!="double") {
        txt <- c(txt, "intensities should be stored in double-precision")
    }
    if (nrow(cellData(object))!=ncol(.raw_cellIntensities(object))) {
        txt <- c(txt, "number of rows in 'cellData' and columns in 'cellIntensities' should be equal")
    }
    if (nrow(markerData(object))!=nrow(.raw_cellIntensities(object))) {
        txt <- c(txt, "'markerData' and 'cellIntensities' must have the same number of rows")
    }
    if (ncol(.raw_intensities(object))!=nrow(markerData(object))) { 
        txt <- c(txt, "number of rows in 'markerData' and columns in 'intensities' must be equal")
    }
    if (nrow(.raw_intensities(object))!=nrow(object)) { 
        txt <- c(txt, "'intensities' and 'object' must have the same number of rows")
    }
    if (nrow(.raw_intensities(object))!=length(.raw_cellAssignments(object))) {
        txt <- c(txt, "number of rows in 'intensities' and length of 'cellAssignments' should be equal")
    }
    if (length(txt)) return(txt)
    return(TRUE)
})

scat <- function(fmt, vals=character(), exdent=2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}

setMethod("show", signature("CyData"), function(object) {
    callNextMethod()
    scat("markers(%d): %s\n", rownames(markerData(object)))
    scat("markerData names(%d): %s\n", colnames(markerData(object)))
    cat(sprintf("cells: %i\n", ncol(.raw_cellIntensities(object))))
    scat("cellData names(%d): %s\n", colnames(cellData(object)))
})

#############################################
# Defining a reasonably helpful constructor.

CyData <- function(markerData, intensities=NULL, cellAssignments=NULL, cellIntensities=NULL, cellData=NULL, assays=NULL, ...) {
    marker.names <- rownames(markerData)
    if (is.null(marker.names)) {
        stop("rownames of 'markerData' must contain marker names")
    }

    if (is.null(cellIntensities)) { 
        cellIntensities <- matrix(0, nrow(markerData), 0)
    } else {
        if (!is.double(cellIntensities)) storage.mode(cellIntensities) <- "double"
        dimnames(cellIntensities) <- NULL
    }

    if (is.null(cellData)) {
        cellData <- DataFrame(row.names=seq_len(ncol(cellIntensities)))
    }

    if (is.null(assays)) {
        assays <- matrix(0L, 0, 0)
    }
    se <- SummarizedExperiment(assays, ...)

    if (is.null(intensities)) {
        intensities <- matrix(0, nrow(assays), nrow(markerData))
    } else {
        if (!is.double(cellIntensities)) storage.mode(intensities) <- "double"
        dimnames(intensities) <- NULL
    }

    if (is.null(cellAssignments)) {
        cellAssignments <- rep(list(integer(0)), nrow(assays))  
    } else {
        names(cellAssignments) <- NULL
    }

    new("CyData", markerData=markerData, cellIntensities=cellIntensities, cellData=cellData, 
        intensities=intensities, cellAssignments=cellAssignments, se)
}

#############################################
# Subsetting, as the SE class doesn't use extractROWS or replaceROWs.

setMethod("[", c("CyData", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) { 
        if (is.character(i)) { 
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i, rownames(x), fmt)
        }
        x@intensities <- .raw_intensities(x)[i,,drop=FALSE] # direct slot assignment, otherwise validity checks cause errors.
        x@cellAssignments <- .raw_cellAssignments(x)[i]
    }
    callNextMethod()
})

.check_identical_refdata <- function(x, y, check.group.values=TRUE) {
    if (!identical(markerData(x), markerData(y))) {
        stop("'markerData' are not identical")
    }
    if (!identical(.raw_cellIntensities(x), .raw_cellIntensities(y))) {
        stop("'cellIntensities' are not identical")
    }
    if (check.group.values) {
        if (!identical(.raw_intensities(x), .raw_intensities(y))) {
            stop("'intensities' are not identical")
        } else if (!identical(.raw_cellAssignments(x), .raw_cellAssignments(y))) { 
            stop("'cellAssignments' are not identical")
        }
    }
    return(TRUE)
}

setMethod("[<-", c("CyData", "ANY", "ANY", "CyData"), function(x, i, j, ..., value) {
    if (!missing(i)) {
        if (is.character(i)) { 
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i, rownames(x), fmt)
        }
        if (!missing(j)) { 
            stop("simultaneous row/column replacement is not supported")
        }
        .check_identical_refdata(x, value, check.group.values=FALSE)
        intensities(x)[i,] <- .raw_intensities(value) 
        cellAssignments(x)[i] <- .raw_cellAssignments(value)
    }
    if (!missing(j)) { 
        .check_identical_refdata(x, value, check.group.values=TRUE)
    }
    callNextMethod(x=x, i=i, j=j, ..., value=value)
})

setMethod("subset", "CyData", function(x, i, j) {
    x[i, j]
})

#############################################
# Defining some getters and setters.

setGeneric("markerData", function(x) standardGeneric("markerData"))
setMethod("markerData", "CyData", function(x) x@markerData)

setGeneric("markerData<-", function(x, value) standardGeneric("markerData<-"))
setReplaceMethod("markerData", "CyData", function(x, value){
    x@markerData <- value
    validObject(x)
    return(x)
})

.raw_cellIntensities <- function(x) x@cellIntensities

setGeneric("cellIntensities", function(x) standardGeneric("cellIntensities"))
setMethod("cellIntensities", "CyData", function(x) {
    out <- .raw_cellIntensities(x)
    rownames(out) <- rownames(markerData(x))
    colnames(out) <- rownames(cellData(x))
    return(out)
})

setGeneric("cellIntensities<-", function(x, value) standardGeneric("cellIntensities<-"))
setReplaceMethod("cellIntensities", "CyData", function(x, value){
    dimnames(value) <- NULL
    x@cellIntensities <- value
    validObject(x)
    return(x)
})

.raw_cellAssignments <- function(x) x@cellAssignments

setGeneric("cellAssignments", function(x) standardGeneric("cellAssignments"))
setMethod("cellAssignments", "CyData", function(x) {
    out <- .raw_cellAssignments(x)
    names(out) <- rownames(x)
    return(out)
})

setGeneric("cellAssignments<-", function(x, value) standardGeneric("cellAssignments<-"))
setReplaceMethod("cellAssignments", "CyData", function(x, value){
    names(value) <- NULL
    x@cellAssignments <- value
    validObject(x)
    return(x)
})

setGeneric("cellData", function(x) standardGeneric("cellData"))
setMethod("cellData", "CyData", function(x) x@cellData)

setGeneric("cellData<-", function(x, value) standardGeneric("cellData<-"))
setReplaceMethod("cellData", "CyData", function(x, value){
    x@cellData <- value
    validObject(x)
    return(x)
})

.raw_intensities <- function(x) x@intensities

setGeneric("intensities", function(x) standardGeneric("intensities"))
setMethod("intensities", "CyData", function(x) {
    out <- .raw_intensities(x)
    rownames(out) <- rownames(x)
    colnames(out) <- rownames(markerData(x))
    return(out)
})

setGeneric("intensities<-", function(x, value) standardGeneric("intensities<-"))
setReplaceMethod("intensities", "CyData", function(x, value){
    dimnames(value) <- NULL
    x@intensities <- value
    validObject(x)
    return(x)
})

setGeneric("nmarkers", function(x) standardGeneric("nmarkers"))
setMethod("nmarkers", "CyData", function(x) nrow(markerData(x)))

setGeneric("ncells", function(x) standardGeneric("ncells"))
setMethod("ncells", "CyData", function(x) ncol(.raw_cellIntensities(x)))

setMethod("markernames", "CyData", function(object) {
    return(rownames(markerData(object)))
})

setReplaceMethod("markernames", c("CyData", "ANY"), function(object, value) {
    rownames(markerData(object)) <- value
    return(object) 
})

setMethod("sampleNames", "CyData", function(object) {
    return(colnames(object))          
})

setReplaceMethod("sampleNames", c("CyData", "ANY"), function(object, value) {
    colnames(object) <- value
    return(object) 
})

#############################################
# Combining objects.

setMethod("cbind", "CyData", function(..., deparse.level=1) {
    args <- unname(list(...))
    ref <- args[[1]]
    for (x in args[-1]) {
        .check_identical_refdata(x, ref)
    }
    
    base <- do.call(cbind, lapply(args, function(x) { as(x, "SummarizedExperiment") }))
    new("CyData", base, intensities=.raw_intensities(ref), cellAssignments=.raw_cellAssignments(ref),
        cellData=cellData(ref), cellIntensities=.raw_cellIntensities(ref), markerData=markerData(ref))
})

setMethod("rbind", "CyData", function(..., deparse.level=1) {
    args <- unname(list(...))
    ref <- args[[1]]
    for (x in args[-1]) {
        .check_identical_refdata(x, ref, check.group.values=FALSE)
    }

    base <- do.call(rbind, lapply(args, function(x) { as(x, "SummarizedExperiment") }))
    new.intensities <- do.call(rbind, lapply(args, .raw_intensities))
    new.assignments <- unlist(lapply(args, .raw_cellAssignments), recursive=FALSE)
    new("CyData", base, intensities=new.intensities, cellAssignments=new.assignments,
        cellData=cellData(ref), cellIntensities=.raw_cellIntensities(ref), markerData=markerData(ref))
})

setMethod("c", "CyData", function(x, ..., recursive = FALSE) {
    rbind(x, ...)
})

#############################################

