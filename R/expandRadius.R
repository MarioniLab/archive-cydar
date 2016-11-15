expandRadius <- function(x, design=NULL, tol=0.5) 
# This computes the standard deviation in the mean intensities,
# in order to compute the radius expansion required for handling
# randomly-distributed intensity shifts between samples.
#
# written by Aaron Lun
# created 27 October 2016   
# last modified 14 November 2016 
{
    .check_cell_data(x, check.clusters=FALSE)
    ci <- cellIntensities(x)
    sample.id <- cellData(x)$sample.id

    # Computing mean intensities for all markers in all samples.
    all.means <- list()
    for (s in seq_len(ncol(x))) { 
        all.means[[s]] <- rowMeans(ci[,sample.id==s])
    }
    all.means <- do.call(rbind, all.means)
    
    # Fitting a linear model to estimate variance of shift process per marker.
    if (is.null(design)) { 
        design <- matrix(1, nrow=nrow(all.means), ncol=1)
    }
    fit <- lm.fit(y=all.means, x=design)

    # Taking the mean variance across all markers, representing squared shift.
    shift2 <- mean(fit$effects[-fit$qr$pivot[seq_len(fit$qr$rank)],]^2) 
    
    # Adding to default tolerance to get per-marker radius.
    # This is done based on squared distance, hence the squaring and rooting.
    # Doubling the shift to account for distances between cells, not just distance from cell to mean.
    return(sqrt(2*shift2 + tol^2))
}
