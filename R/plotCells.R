plotCellLogFC <- function(x, y, logFC, max.logFC=NULL, zero.col=0.8, length.out=100, pch=16, ...) 
# Visualizes the cells with appropriate coloration, using a PCA plot by default.
#
# written by Aaron Lun
# created 20 April 2016
# last modified 8 June 2016
{
    if (is.null(max.logFC)) {
        max.logFC <- max(abs(logFC))
    } else {
        logFC[logFC < -max.logFC] <- -max.logFC
        logFC[logFC > max.logFC] <- max.logFC
    }
    logFC.col <- color.scale(logFC, c(0,zero.col,1), c(0,zero.col,0), c(1,zero.col,0), 
                             xrange=c(-max.logFC, max.logFC))

    if (pch %in% 21:25) {           
        plot(x, y, pch=pch, bg=logFC.col, ...)
    } else {
        plot(x, y, pch=pch, col=logFC.col, ...)
    }
                          
    values <- seq(from=-max.logFC, to=max.logFC, length.out=length.out)
    stored <- color.scale(values, c(0,zero.col,1), c(0,zero.col,0), c(1,zero.col,0), 
                          xrange=c(-max.logFC, max.logFC))
    names(stored) <- values
    return(invisible(stored))
}

plotCellIntensity <- function(x, y, intensity, irange=NULL, length.out=100, pch=16, ...) 
# Visualizes the cells with appropriate coloration, using a PCA plot by default.
#
# written by Aaron Lun
# created 20 April 2016
# last modified 8 June 2016
{
    if (is.null(irange)) { 
        irange <- range(intensity)
    } else {
        intensity[intensity > irange[2]] <- irange[2]
        intensity[intensity < irange[1]] <- irange[1]
    }

    all.cols <- viridis(n=length.out)
    mdpts <- seq(from=irange[1], to=irange[2], length.out=length.out)
    actual.threshold <- mdpts + (mdpts[2] - mdpts[1])/2
    ix <- pmin(length.out, findInterval(intensity, actual.threshold) + 1)
    cur.cols <- all.cols[ix]

    if (pch %in% 21:25) {           
        plot(x, y, pch=pch, bg=cur.cols, ...)
    } else {
        plot(x, y, pch=pch, col=cur.cols, ...)
    }
                          
    names(all.cols) <- mdpts 
    return(invisible(all.cols))
}
