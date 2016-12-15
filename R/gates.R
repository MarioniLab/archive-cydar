outlierGate <- function(x, name, nmads=3, type=c("both", "upper", "lower")) 
# Constructs a gate to remove outliers in 1 dimension, 
# based on intensities that are more than 'nmads' from the median. 
#
# written by Aaron Lun
# created 15 December 2016
{
    intensity <- exprs(x)[,name]
    center <- median(intensity)
    dev <- mad(intensity, center=center)
    type <- match.arg(type)
    
    upper.threshold <- Inf
    lower.threshold <- -Inf
    if (type=="both" || type=="lower") {
        lower.threshold <- center - nmads*dev        
    }
    if (type=="both" || type=="upper") {
        upper.threshold <- center + nmads*dev
    }
    
    gate <- list(c(lower.threshold, upper.threshold))
    names(gate) <- name
    rectangleGate(filterId=paste0(name, "_outlierGate"), .gate=gate) 
}

dnaGate <- function(x, name1, name2, tol=0.5, nmads=3, type=c("both", "lower"))
# Constructs a gate to remove non-cells, doublets, and
# cells with differences in the two DNA channels.
# This assumes that most events correspond to singlets.
#
# written by Aaron Lun
# created 15 December 2016    
{
    ex1 <- exprs(x)[,name1]
    bound1 <- .get_LR_bounds(ex1, nmads=nmads)
    lower.dna1 <- bound1$left
    upper.dna1 <- bound1$right

    ex2 <- exprs(x)[,name2]
    bound2 <- .get_LR_bounds(ex2, nmads=nmads)
    lower.dna2 <- bound2$left
    upper.dna2 <- bound2$right
    
    # Calculating the lines corresponding to the boundaries.
    grad.par <- (upper.dna2 - lower.dna2) / (upper.dna1 - lower.dna1)
    intercept <- lower.dna2 - grad.par * lower.dna1
    lower.int.par <- intercept - tol
    upper.int.par <- intercept + tol

    grad.per <- -1/grad.par
    lower.int.per <- lower.dna2 - grad.per*lower.dna1
    upper.int.per <- upper.dna2 - grad.per*upper.dna1
   
    # Effectively eliminate the upper bound if necessary.
    type <- match.arg(type)
    if (type=="lower") upper.int.per <- 10000

    # Finding the vertices.
    all.vertices <- rbind(.find_vertex(grad.par, lower.int.par, grad.per, lower.int.per),
                          .find_vertex(grad.per, lower.int.per, grad.par, upper.int.par),
                          .find_vertex(grad.par, upper.int.par, grad.per, upper.int.per),
                          .find_vertex(grad.per, upper.int.per, grad.par, lower.int.par))

    colnames(all.vertices) <- c(name1, name2)
    polygonGate(filterId="dnaGate", .gate=all.vertices)
}

.get_LR_bounds <- function(x, nmads) {
    center <- median(x)
    dev <- mad(x, center=center)
    ref.lower <- center - nmads*dev
    ref.upper <- center + nmads*dev

    dens <- density(x)
    max.x <- which.max(dens$y)
    first.deriv <- diff(dens$y)
    local.min <- which(c(TRUE, first.deriv < 0) & c(first.deriv > 0, TRUE))

    lower <- dens$x[max(local.min[local.min < max.x])]
    if (lower < ref.lower) lower <- ref.lower
    upper <- dens$x[min(local.min[local.min > max.x])]
    if (upper > ref.upper) upper <- ref.upper
    return(list(left=lower, right=upper))
}

.find_vertex <- function(m1, b1, m2, b2) {
    x <- (b2-b1)/(m1-m2)
    return(c(x=x, y=m1*x+b1))
}

