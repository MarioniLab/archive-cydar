multiIntHist <- function(collected, cols=NULL,
                         xlab="Intensity", ylab="Density", ..., 
                         lwd=2, lty=1, pch=16, cex=2) 
# This plots multiple intensity distributions as histogram outlines on a single plot.
# The density at zero is also shown in points at zero.
#
# written by Aaron Lun
# created 30 January 2017   
# last modified 21 March 2017 
{   
    if (is.null(cols)) { 
        cols <- grDevices::rainbow(length(collected))
    } else {
        cols <- rep(cols, length.out=length(collected))
    }
    
    max.x <- max(sapply(collected, max))
    breaks <- seq(from=0, to=max.x, length.out=50)
    h.all <- h.zero <- vector("list", length(collected))
    for (batch in seq_along(collected)) {
        cur.vals <- collected[[batch]] 
        is.zero <- cur.vals < 1e-6 
        cur.vals <- cur.vals[!is.zero]
        hout <- graphics::hist(cur.vals, breaks=breaks, plot=FALSE)
        hout$density <- hout$density * length(is.zero)/length(cur.vals)
        h.all[[batch]] <- hout 
        h.zero[[batch]] <- sum(is.zero)/length(is.zero)
    }

    max.y <- max(sapply(h.all, FUN=function(incoming) max(incoming$density)))
    max.y <- max(max.y, unlist(h.zero))
    plot(0,0,type="n", xlab=xlab, ylab=ylab, xlim=c(0, max.x), ylim=c(0, max.y), ...)
    for (batch in seq_along(collected)) {
        lines(rep(breaks, each=2), c(0, rep(h.all[[batch]]$density, each=2), 0), 
              col=cols[batch], lwd=lwd, lty=lty)
        points(0, jitter(h.zero[[batch]]), col=cols[batch], pch=pch, cex=cex)
    }
    invisible(NULL)
}


