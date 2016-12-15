poolCells <- function(x, equalize=TRUE, n=NULL)
# This function creates a flowFrame from cells pooled from multiple samples.
# The idea is to estimate transformation and gating parameters from the pool,
# so that you don't have sample-specific parameters that interfere with comparisons.
#
# written by Aaron Lun
# created 15 December 2016
{
    on.exit({gc()}) # Getting rid of huge memory structures that have built up.
    cell.data <- .pull_out_data(x)

    if (equalize) { 
        if (is.null(n)) n <- min(sapply(cell.data$exprs, nrow))
        n <- as.integer(n)
        if (n <= 0L) stop("'n' must be a positive integer")
        for (s in seq_along(cell.data$samples)) {
            cur.exprs <- cell.data$exprs[[s]]
            i <- round(seq(1, nrow(cur.exprs), length.out=n))
            cell.data$exprs[[s]] <- cur.exprs[i,,drop=FALSE]            
        }
    }
    new.exprs <- do.call(rbind, cell.data$exprs)

    param <- data.frame(name=cell.data$markers,
                        desc=cell.data$markers,
                        minRange=apply(new.exprs, 2, min),
                        maxRange=apply(new.exprs, 2, max),
                        row.names=paste0("$P", seq_along(cell.data$markers)))
    param$range <- param$maxRange - param$minRange
    
    flowFrame(new.exprs, AnnotatedDataFrame(param))
}
