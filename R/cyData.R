# This sets up the cyData reference class.

setClass("cyData", contains="SummarizedExperiment", 
         slots=c(markerData="DataFrame",
                 cellIntensities="matrix",
                 medianIntensities="matrix"))

setValidity2("cyData", function(object) {
    if (storage.mode(object@cellIntensities)!="double") {
        return("cellIntensities should be stored in double-precision")
    }
    if (storage.mode(object@medianIntensities)!="double") {
        return("medianIntensities should be stored in double-precision")
    }
    if (nrow(object@markerData)!=nrow(object@cellIntensities)) {
        return("'markerData' and 'cellIntensities' must have same number of rows")
    }
    if (ncol(object@medianIntensities)!=nrow(object@markerData)) { 
        return("number of rows in 'markerData' and columns in 'medianIntensities' must be equal")
    }
    if (nrow(object@medianIntensities)!=nrow(object)) { 
        return("'medianIntensities' and 'object' must have same number of rows")
    }
    return(TRUE)
})

setMethod("parallelSlotNames", "cyData", function(x) {
    c("medianIntensities", callNextMethod()) 
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
    cat(sprintf("cells: %i\n", ncol(object@cellIntensities)))
})

#############################################
# Defining a reasonably helpful constructor.

cyData <- function(markerData, medianIntensities=NULL, cellIntensities=NULL, assays=NULL, ...) {
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

    if (is.null(assays)) {
        assays <- matrix(0, 0, 0)
    }
    se <- SummarizedExperiment(assays, ...)

    if (is.null(medianIntensities)) {
        medianIntensities <- matrix(0, nrow(assays), nrow(markerData))
    } else {
        if (!is.double(cellIntensities)) storage.mode(medianIntensities) <- "double"
        rownames(medianIntensities) <- NULL
    }

    new("cyData", markerData=markerData, cellIntensities=cellIntensities,
        medianIntensities=medianIntensities, se)
}

#############################################
# Subsetting, as the SE class doesn't use extractROWS or replaceROWs.

setMethod("[", c("cyData", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) { x@medianIntensities <- x@medianIntensities[i,,drop=FALSE] }
    callNextMethod()
})

.check_identical_refdata <- function(x, y, check.medians=TRUE) {
    if (!identical(x@markerData, y@markerData)) {
        stop("'markerData' are not identical")
    }
    if (!identical(x@cellIntensities, y@cellIntensities)) {
        stop("'cellIntensities' are not identical")
    }
    if (check.medians && !identical(x@medianIntensities, y@medianIntensities)) {
        stop("'medianIntensities' are not identical")
    }
    return(TRUE)
}

setMethod("[<-", c("cyData", "ANY", "ANY", "cyData"), function(x, i, j, ..., value) {
    .check_identical_refdata(x, value, check.medians=FALSE)
    if (!missing(i)) { x@medianIntensities[i,] <- value@medianIntensities }
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

setGeneric("medianIntensities", function(x) standardGeneric("medianIntensities"))
setMethod("medianIntensities", "cyData", function(x) x@medianIntensities)

setGeneric("medianIntensities<-", function(x, value) standardGeneric("medianIntensities<-"))
setReplaceMethod("medianIntensities", "cyData", function(x, value){
    x@medianIntensities <- value
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
    new("cyData", base, medianIntensities=ref@medianIntensities,
        cellIntensities=ref@cellIntensities, markerData=ref@markerData)
})

setMethod("rbind", "cyData", function(..., deparse.level=1) {
    args <- unname(list(...))
    ref <- args[[1]]
    for (x in args[-1]) {
        .check_identical_refdata(x, ref, check.medians=FALSE)
    }

    base <- do.call(rbind, lapply(args, function(x) { as(x, "SummarizedExperiment") }))
    new.meds <- do.call(rbind, lapply(args, medianIntensities))
    new("cyData", base, medianIntensities=new.meds,
        cellIntensities=ref@cellIntensities, markerData=ref@markerData)
})

setMethod("c", "cyData", function(x, ..., recursive = FALSE) {
    rbind(x, ...)
})

#############################################

