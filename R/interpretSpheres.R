interpretSpheres <- function(x, markers=NULL, labels=NULL, metrics=NULL, num.per.row=6, plot.height=100, xlim=NULL, 
                             p=0.01, add.plot=NULL, add.plot.height=500, run=TRUE, ...)   
# This creates a Shiny app to assist interpretation of the hyperspheres.
#
# written by Aaron Lun
# created 1 November 2016    
# last modified 3 December 2017
{
    coords <- intensities(x)
    nrows <- ceiling(ncol(coords)/num.per.row)
    plot.height <- plot.height*nrows
    collim <- intensityRanges(x, p=p)

    coords <- as.matrix(coords)
    N <- nrow(coords)
    cell.int <- cellIntensities(x)
    cell.assign <- cellAssignments(x)

    # Checking marker order.
    if (!is.null(markers)) { 
        m.order <- match(markers, colnames(coords))
        stopifnot(!any(is.na(m.order)))
    } else {
        m.order <- seq_len(ncol(coords))
    }

    # Using existing labels if provided.
    if (is.null(labels)) {
        labels <- character(N)
    } else {
        stopifnot(length(labels)==N)
    } 

    # Checking metrics if provided.
    if (!is.null(metrics)) {
        stopifnot(nrow(metrics)==N)
    }

    # Setting up the internal plotting and storage.
    collected <- new.env()
    collected$labels <- labels

    collected$current <- 1
    collected$interval <- 0
    collected$history <- rep(NA_integer_, 5)
    all.dens <- .prepareDensity(x, ...)

    # Main panel arguments.
    main.args <- list(plotOutput("histograms", height = plot.height), hr())
    main.args <- append(main.args, 
                        list(
                             fluidRow(
                                      column(
                                             textInput("label", "Label:"),
                                             actionButton("addlabel", "Add label"),
                                             width=4
                                             ),
                                      column(
                                             textInput("gotonum", "Go to sphere:"),
                                             actionButton("go", "Go"),
                                             actionButton("previous", "Previous"),
                                             actionButton("continue", "Next"),
                                             width=4
                                             ),
                                      column( 
                                             actionButton("finish", "Save to R"),
                                             width=4
                                             )
                                      )
                             )
                        )
    if (!is.null(add.plot)) {
        main.args <- append(main.args, list(hr(), plotOutput("addplots", height = add.plot.height)))
    }

    # Generating the page layout.
    ui <- pageWithSidebar(
                          headerPanel("Interpreting hypersphere coordinates"),
                          sidebarPanel(
                                       h4("Metrics"),
                                       tableOutput("metrics"),
                                       hr(size=30),
                                       h4("History:"),
                                       tableOutput("history"), 
                                       hr(size=30),
                                       h4("Closest labelled:"),
                                       tableOutput("closest"), 
                                       hr(size=30),
                                       sliderInput("intbar", "Intensity interval for current hypersphere (%):",
                                                   min = 0, max = 100, value = 95),
                                       textInput("extraplot", "Add more hyperspheres:", value=""),
                                       actionButton("addtoplot", "Add to plot"),
                                       actionButton("clearplot", "Clear")
                                       ),
                          do.call(mainPanel, main.args)
                          )

    # Setting up updating functions.
    updatePlot <- function(output, extras=NULL) {
        output$histograms <- renderPlot({
            par(mfrow=c(nrows, num.per.row), mar=c(2.1, 1.1, 2.1, 1.1))
            .generateDensity(coords, all.dens, xlim=xlim, collim=collim, ordering=m.order,
                             current=collected$current, extras=extras, 
                             interval=collected$interval, cell.int=cell.int, cell.assign=cell.assign)
        })
        if (!is.null(add.plot)) {
            output$addplots <- renderPlot({ add.plot(collected$current, x) }) 
        }
    }

    updateMetrics <- function(output) {
        output$metrics <- renderTable({
            out.metrics <- c("Number", "Label")
            out.vals <- c(collected$current, collected$labels[collected$current])

            if (!is.null(metrics)) { 
                extra.metrics <- colnames(metrics)
                extra.vals <- list() 
                for (met in seq_along(extra.metrics)) {
                    my.val <- metrics[collected$current,met]
                    if (is.double(my.val)) {
                        my.val <- format(my.val, digits=5)
                    } else {
                        my.val <- as.character(my.val)
                    }
                    extra.vals[[met]] <- my.val
                }
                out.metrics <- c(out.metrics, extra.metrics)
                out.vals <- c(out.vals, unlist(extra.vals))
            }

            data.frame(paste0("<b>", out.metrics, "</b>"), out.vals)
        }, colnames=FALSE, sanitize.text.function=function(x) {x})
    }

    updateHistory <- function(output, advance=TRUE) {
        output$history <- renderTable({
            if (advance) {
                collected$history <- c(collected$current, collected$history)[1:5]
            }
            data.frame(Number=as.character(collected$history),
                Label=collected$labels[collected$history])
        })
    }

    updateClosest <- function(output) {
        output$closest <- renderTable({
            is.anno <- which(collected$labels!="")
            if (length(is.anno)) { 
                all.dist <- sqrt(colSums((t(coords[is.anno,,drop=FALSE]) - coords[collected$current,])^2))
            } else {
                all.dist <- numeric(0)
            }
            o <- order(all.dist)[1:5]
            closest <- is.anno[o]
            data.frame(Distance=all.dist[o], Number=as.character(closest), Label=collected$labels[closest])
        })
    }

    # Setting up the server actions.
    server <- function(input, output) {
        updatePlot(output)
        updateMetrics(output)
        updateHistory(output)
        updateClosest(output)

        observeEvent(input$previous, {
            if (collected$current!=1L) { 
                collected$current <- collected$current - 1L
            } else {
                warning("index cannot be negative")
            }
            updatePlot(output)
            updateMetrics(output)
            updateHistory(output)
            updateClosest(output)
        })

        observeEvent(input$continue, {
            if (collected$current!=N) {
                collected$current <- collected$current + 1L
            } else {
                warning("index greater than number of hyperspheres")
            }
            updatePlot(output)
            updateMetrics(output)
            updateHistory(output)  
            updateClosest(output)  
        })

        observeEvent(input$intbar, {
            collected$interval <- input$intbar/100
            updatePlot(output)
        })

        observeEvent(input$go, {
            candidate <- as.integer(input$gotonum)
            if (is.finite(candidate) && candidate >= 1L && candidate <= N) {
                collected$current <- candidate
            } else {
                warning("index out of range of number of hyperspheres")
            }
            updatePlot(output)
            updateMetrics(output)
            updateHistory(output)
            updateClosest(output)
        })

        observeEvent(input$addtoplot, {
            extras <- as.integer(unlist(strsplit(input$extraplot, " ") ))
            extras <- extras[!is.na(extras)]
            invalid <- !is.finite(extras) && extras < 1L && extras > N
            if (any(invalid)) {
                warning("indices out of range of number of hyperspheres")
                extras <- extras[!invalid]
            }
            updatePlot(output, extras=extras)
        })

        observeEvent(input$clearplot, {
            updatePlot(output)
        })

        observeEvent(input$addlabel, {
            collected$labels[collected$current] <- input$label
            updateHistory(output, advance=FALSE)
            updateMetrics(output)
        })

        observeEvent(input$finish, {
            stopApp(collected$labels)
        })
    }
    
    app <- shinyApp(ui, server)
    if (run) {
        return(runApp(app))
    } else {
        return(app)
    }
}

.prepareDensity <- function(x, ...) 
# Computing the density once to speed up plotting.
{ 
    ci <- cellIntensities(x)
    all.markers <- rownames(markerData(x))
    collected <- list()

    for (m in seq_along(all.markers)) {
        cur.intensities <- ci[m,]
        collected[[m]] <- density(cur.intensities, ...)
    }

    names(collected) <- all.markers
    return(collected)
}

.generateDensity <- function(coords, density.data, ordering, collim, xlim, current, extras, 
                             interval, cell.int, cell.assign, ...) 
# Plotting the densities and adding the point corresponding to the current coordinates.
{
    current.coords <- coords[current,]
    all.markers <- names(current.coords)
    all.cols <- viridis(256)

    if (interval > 0) { 
        cur.assign <- unpackIndices(cell.assign[current])
        cur.cell.int <- cell.int[,cur.assign[[1]],drop=FALSE]
        half.left <- (1-interval)/2
        prob.interval <- c(half.left, interval + half.left)
    }

    for (m in ordering) {
        curdex <- round(approx(collim[,m], c(1, 256), xout=current.coords[m], rule=2)$y)
        if (is.null(xlim)) { 
            xlim <- collim[,m]
        }
    
        curdens <- density.data[[m]]
        plot(0, 0, type="n", xlab="", ylab="", yaxt="n", bty="n", main=all.markers[m], xlim=xlim, ylim=c(0, max(curdens$y)), ...)
        my.x <- c(curdens$x[1]-10, curdens$x, curdens$x[length(curdens$x)]+10)
        my.y <- c(0, curdens$y, 0)
        polygon(my.x, my.y, col=all.cols[curdex], border=NA)
        lines(my.x, my.y)
    
        curpos <- current.coords[m]
        curpos <- pmin(xlim[2], pmax(xlim[1], curpos))
        cury <- approx(my.x, my.y, curpos, rule=2)$y
        par(xpd=TRUE)
        points(curpos, cury, pch=16, col="red", cex=1.5)
        
        if (interval > 0) {
            q.int <- quantile(cur.cell.int[m,], prob.interval)
            segments(q.int[1], cury, q.int[2], col="red")
        }

        if (length(extras)) { 
            text(curpos, cury, pos=3, current, col="red")
            for (ex in extras) {
                expos <- coords[ex,m]
                expos <- pmin(xlim[2], pmax(xlim[1], expos))
                exy <- approx(my.x, my.y, expos, rule=2)$y
                points(expos, exy, pch=16, col="red", cex=1.5)
                text(expos, exy, pos=3, ex, col="red")
            }
        }
        par(xpd=FALSE)
    }

    return(invisible())
}

