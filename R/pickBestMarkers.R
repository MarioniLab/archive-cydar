pickBestMarkers <- function(cell.data, chosen, tol=0.5, downsample=10, p=0.05, naive=FALSE)
# The idea is to fit a binomial GLM to the coordinates, using the marker intensity as covariates 
# and whether it belongs in the selected set as a response. We then calculate the contamination 
# and recovery rate for each increasing set of markers, based on the number of points that fall
# into a box defined by the boundaries of the chosen points (we set 'p' to protect against outliers).
# 
# written by Aaron Lun
# created 10 June 2016
# last modified 11 August 2016
{
    .check_cell_data(cell.data)
    markers <- attributes(cell.data)$markers
    distance <- as.double(tol) * sqrt(length(markers))
   
    # Downsampling cells for speed.
    downsample <- as.integer(downsample)
    cell.id <- attributes(cell.data)$cell.id
    selected <- which((cell.id %% downsample) == 0L)
    selected <- cell.data[,selected]

    # Organizing the hypersphere centers. 
    index <- as.integer(sub("^c", "", chosen))
    coords <- cell.data[,index,drop=FALSE]
    if (naive) {
        cluster.centers <- cluster.info <- NULL 
    } else {
        rownames(coords) <- markers
        coords <- prepareCellData(list(A=t(coords)))
        cluster.centers <- attributes(coords)$cluster.centers
        cluster.info <- attributes(coords)$cluster.info
    }

    # Selecting the centers.
    was.counted <- .Call(cxx_find_counted, coords, cluster.centers, cluster.info, selected, distance)
    if (is.character(was.counted)) {
        stop(was.counted)
    }
    
    # Fitting a binomial GLM with LASSO.
    coordinates <- t(selected)
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
    lose.intercept <- rownames(coef(out))!="(Intercept)"

    for (m in seq_len(ncol(coef(out)))) {
        new.markers <- coef(out)[lose.intercept,m]!=0 & !currently.collected

        if (any(new.markers)) {
            currently.collected[new.markers] <- TRUE
            for (mi in which(new.markers)) { 
                inside.box <- inside.box & outside[,mi] >= inranges[1,mi] & outside[,mi] <= inranges[2,mi]
                self.inside.box <- self.inside.box & inside[,mi] >= inranges[1,mi] & inside[,mi] <= inranges[2,mi]
            }

            nfp <- sum(inside.box)
            ntp <- sum(self.inside.box)
            output[[length(output)+1]] <- data.frame(markers[new.markers], nfp/(ntp+nfp), ntp/nrow(inside), 
                                                     inranges[1,new.markers], inranges[2,new.markers], m)
        }
    }

    output <- do.call(rbind, output)
    rownames(output) <- NULL
    colnames(output) <- c("Marker", "Contamination", "Recovery", "Lower", "Upper", "Iteration")
    return(output)
} 

