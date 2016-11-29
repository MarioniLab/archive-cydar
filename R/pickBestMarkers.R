pickBestMarkers <- function(x, chosen, markers=NULL, downsample=10, p=0.05, naive=FALSE)
# The idea is to fit a binomial GLM to the coordinates, using the marker intensity as covariates 
# and whether it belongs in the selected set as a response. We then calculate the contamination 
# and recovery rate for each increasing set of markers, based on the number of points that fall
# into a box defined by the boundaries of the chosen points (we set 'p' to protect against outliers).
# 
# written by Aaron Lun
# created 10 June 2016
# last modified 29 November 2016
{
    .check_cell_data(x, check.clusters=FALSE)
    if (is.null(markers)) markers <- markerData(x)$used
    used <- .chosen_markers(markers, markernames(x))

    # Figuring out which cells were counted.
    cur.assignments <- cellAssignments(x)[chosen]
    was.counted <- logical(nrow(cellData(x)))
    for (g in cur.assignments) {
        expanded <- unpackIndices(list(g))[[1]]
        was.counted[expanded] <- TRUE
    }

    # Downsampling cells for speed.
    downsample <- as.integer(downsample)
    selected <- ((cellData(x)$cell.id - 1L) %% downsample)==0L
    coordinates <- t(cellIntensities(x)[used,selected,drop=FALSE])
    was.counted <- was.counted[selected]
    
    # Fitting a binomial GLM with LASSO.
    out <- glmnet(coordinates, # as the marker intensities are the covariates.
           factor(was.counted), # response is whether it is in or not.
           family="binomial", intercept=TRUE)

    # Specifying the ranges with which we will count collections.
    inside <- coordinates[was.counted,,drop=FALSE]
    inranges <- apply(inside, 2, quantile, p=c(p, 1-p))
    outside <- coordinates[!was.counted,,drop=FALSE]

    # Going through and collecting markers.
    currently.collected <- logical(ncol(coordinates))
    inside.box <- !logical(nrow(outside))
    self.inside.box <- !logical(nrow(inside))
    output <- list()
    marker.as.coef <- coef(out)[rownames(coef(out))!="(Intercept)",]

    for (m in seq_len(ncol(coef(out)))) {
        new.markers <- marker.as.coef[,m]!=0 & !currently.collected
        new.marker.names <- rownames(marker.as.coef)[new.markers]

        if (any(new.markers)) {
            currently.collected[new.markers] <- TRUE
            for (mi in which(new.markers)) { 
                inside.box <- inside.box & outside[,mi] >= inranges[1,mi] & outside[,mi] <= inranges[2,mi]
                self.inside.box <- self.inside.box & inside[,mi] >= inranges[1,mi] & inside[,mi] <= inranges[2,mi]
            }

            nfp <- sum(inside.box)
            ntp <- sum(self.inside.box)
            output[[length(output)+1]] <- data.frame(new.marker.names, nfp/(ntp+nfp), ntp/nrow(inside), 
                                                     inranges[1,new.markers], inranges[2,new.markers], m)
        }
    }

    output <- do.call(rbind, output)
    rownames(output) <- NULL
    colnames(output) <- c("Marker", "Contamination", "Recovery", "Lower", "Upper", "Iteration")
    return(output)
} 

