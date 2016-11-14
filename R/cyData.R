# This sets up the cyData reference class.

setClass("cyData", contains="SummarizedExperiment", 
         slots=c(markerData="DataFrame",
                 cellIntensities="matrix",
                 cellData="DataFrame",
                 intensities="matrix"))

setValidity2("cyData", function(object) {
    if (storage.mode(object@cellIntensities)!="double") {
        return("cellIntensities should be stored in double-precision")
    }
    if (storage.mode(object@intensities)!="double") {
        return("intensities should be stored in double-precision")
    }
    if (nrow(object@cellData)!=ncol(object@cellIntensities)) {
        return("number of rows in 'cellData' and columns in 'cellIntensities' should be equal")
    }
    if (nrow(object@markerData)!=nrow(object@cellIntensities)) {
        return("'markerData' and 'cellIntensities' must have same number of rows")
    }
    if (ncol(object@intensities)!=nrow(object@markerData)) { 
        return("number of rows in 'markerData' and columns in 'intensities' must be equal")
    }
    if (nrow(object@intensities)!=nrow(object)) { 
        return("'intensities' and 'object' must have same number of rows")
    }
    return(TRUE)
})

setMethod("parallelSlotNames", "cyData", function(x) {
    c("intensities", callNextMethod()) 
})

scat <- function(fmt, vals=character(), exdent=2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}

setMethod("show", signature("cyData"), function(object) {
    callNextMethod()
    scat("markers(%d): %s\n", rownames(object@markerData))
    scat("markerData names(%d): %s\n", colnames(object@markerData))
    cat(sprintf("cells: %i\n", ncol(object@cellIntensities)))
    scat("cellData names(%d): %s\n", colnames(object@cellData))
})

#############################################
# Defining a reasonably helpful constructor.

cyData <- function(markerData, intensities=NULL, cellIntensities=NULL, cellData=NULL, assays=NULL, ...) {
    marker.names <- rownames(markerData)
    if (is.null(marker.names)) {
        stop("rownames of 'markerData' must contain marker names")
    }

    if (is.null(cellIntensities)) { 
        cellIntensities <- matrix(0, nrow(markerData), 0)
    } else {
        if (!is.double(cellIntensities)) storage.mode(cellIntensities) <- "double"
        rownames(cellIntensities) <- NULL
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
        rownames(intensities) <- NULL
    }

    new("cyData", markerData=markerData, cellIntensities=cellIntensities, 
        cellData=cellData, intensities=intensities, se)
}

#############################################
# Subsetting, as the SE class doesn't use extractROWS or replaceROWs.

setMethod("[", c("cyData", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) { x@intensities <- x@intensities[i,,drop=FALSE] }
    callNextMethod()
})

.check_identical_refdata <- function(x, y, check.medians=TRUE) {
    if (!identical(x@markerData, y@markerData)) {
        stop("'markerData' are not identical")
    }
    if (!identical(x@cellIntensities, y@cellIntensities)) {
        stop("'cellIntensities' are not identical")
    }
    if (check.medians && !identical(x@intensities, y@intensities)) {
        stop("'intensities' are not identical")
    }
    return(TRUE)
}

setMethod("[<-", c("cyData", "ANY", "ANY", "cyData"), function(x, i, j, ..., value) {
    .check_identical_refdata(x, value, check.medians=FALSE)
    if (!missing(i)) { x@intensities[i,] <- value@intensities }
    callNextMethod(x=x, i=i, j=j, ..., value=value)
})

setMethod("subset", "cyData", function(x, i, j) {
    x[i, j]
})

#############################################
# Defining some getters and setters.

setGeneric("markerData", function(x) standardGeneric("markerData"))
setMethod("markerData", "cyData", function(x) x@markerData)

setGeneric("markerData<-", function(x, value) standardGeneric("markerData<-"))
setReplaceMethod("markerData", "cyData", function(x, value){
    x@markerData <- value
    validObject(x)
    return(x)
})

setGeneric("cellIntensities", function(x) standardGeneric("cellIntensities"))
setMethod("cellIntensities", "cyData", function(x) x@cellIntensities)

setGeneric("cellIntensities<-", function(x, value) standardGeneric("cellIntensities<-"))
setReplaceMethod("cellIntensities", "cyData", function(x, value){
    x@cellIntensities <- value
    validObject(x)
    return(x)
})

setGeneric("cellData", function(x) standardGeneric("cellData"))
setMethod("cellData", "cyData", function(x) x@cellData)

setGeneric("cellData<-", function(x, value) standardGeneric("cellData<-"))
setReplaceMethod("cellData", "cyData", function(x, value){
    x@cellData <- value
    validObject(x)
    return(x)
})

setGeneric("intensities", function(x) standardGeneric("intensities"))
setMethod("intensities", "cyData", function(x) x@intensities)

setGeneric("intensities<-", function(x, value) standardGeneric("intensities<-"))
setReplaceMethod("intensities", "cyData", function(x, value){
    x@intensities <- value
    validObject(x)
    return(x)
})

setGeneric("nmarkers", function(x) standardGeneric("nmarkers"))
setMethod("nmarkers", "cyData", function(x) nrow(x@markerData))

setGeneric("ncells", function(x) standardGeneric("ncells"))
setMethod("ncells", "cyData", function(x) ncol(x@cellIntensities))

#############################################
# Combining objects.

setMethod("cbind", "cyData", function(..., deparse.level=1) {
    args <- unname(list(...))
    ref <- args[[1]]
    for (x in args[-1]) {
        .check_identical_refdata(x, ref)
    }
    
    base <- do.call(cbind, lapply(args, function(x) { as(x, "SummarizedExperiment") }))
    new("cyData", base, intensities=ref@intensities, cellData=ref@cellData,
        cellIntensities=ref@cellIntensities, markerData=ref@markerData)
})

setMethod("rbind", "cyData", function(..., deparse.level=1) {
    args <- unname(list(...))
    ref <- args[[1]]
    for (x in args[-1]) {
        .check_identical_refdata(x, ref, check.medians=FALSE)
    }

    base <- do.call(rbind, lapply(args, function(x) { as(x, "SummarizedExperiment") }))
    new.meds <- do.call(rbind, lapply(args, intensities))
    new("cyData", base, intensities=new.meds, cellData=ref@cellData,
        cellIntensities=ref@cellIntensities, markerData=ref@markerData)
})

setMethod("c", "cyData", function(x, ..., recursive = FALSE) {
    rbind(x, ...)
})

#############################################

