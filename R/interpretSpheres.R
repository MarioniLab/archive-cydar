interpretSpheres <- function(x, markers=NULL, labels=NULL, num.per.row=6, plot.height=100, xlim=NULL, p=0.01, run=TRUE, ...)   
# This creates a Shiny app to assist interpretation of the hyperspheres.
#
# written by Aaron Lun
# created 1 November 2016    
# last modified 14 November 2016
{
    coords <- medianIntensities(x)
    nrows <- ceiling(ncol(coords)/num.per.row)
    plot.height <- plot.height*nrows
    collim <- intensityRanges(x, p=p)

    coords <- as.matrix(coords)
    if (is.null(rownames(coords))) {
        sphere.names <- seq_len(nrow(coords))
    } else {
        sphere.names <- rownames(coords)
    }

    # Checking marker order.
    if (!is.null(markers)) { 
        m.order <- match(markers, colnames(coords))
        stopifnot(!any(is.na(m.order)))
    } else {
        m.order <- seq_len(ncol(coords))
    }

    # Using existing labels if provided.
    if (is.null(labels)) {
        labels <- character(length(sphere.names))
    } else {
        stopifnot(length(labels)==length(sphere.names))
    } 

    # Setting up the internal plotting and storage.
    collected <- new.env()
    collected$labels <- labels

    collected$current <- 1
    collected$history <- rep(NA_integer_, 5)
    all.dens <- .prepareDensity(x, ...)

    # Generating the page layout.
    ui <- pageWithSidebar(
                    headerPanel("Interpreting hypersphere coordinates"),
                    sidebarPanel(
                                 h4("History:"),
                                 tableOutput("history"),
                                 h4("Closest labelled:"),
                                 tableOutput("closest")
                                 ),
                    mainPanel(
                              plotOutput("histograms", height = plot.height), 
                              hr(),
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

    # Setting up updating functions.
    updatePlot <- function(output) {
        output$histograms <- renderPlot({
            par(mfrow=c(nrows, num.per.row), mar=c(2.1, 1.1, 2.1, 1.1))
            .generateDensity(coords[collected$current,], all.dens, xlim=xlim, collim=collim, ordering=m.order) 
        })
    }

    updateHistory <- function(output, advance=TRUE) {
        output$history <- renderTable({
            if (advance) {
                collected$history <- c(collected$current, collected$history)[1:5]
            }
            data.frame(Number=as.character(collected$history),
                Name=sphere.names[collected$history],
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
            data.frame(Distance=all.dist[o], Number=as.character(closest), 
                Name=sphere.names[closest], Label=collected$labels[closest])
        })

    }

    # Setting up the server actions.
    server <- function(input, output) {
        updatePlot(output)
        updateHistory(output)
        updateClosest(output)

        observeEvent(input$previous, {
            if (collected$current!=1L) { 
                collected$current <- collected$current - 1L
            }
            updatePlot(output)
            updateHistory(output)
            updateClosest(output)
        })

        observeEvent(input$continue, {
            if (collected$current!=length(sphere.names)) {
                collected$current <- collected$current + 1L
            }
            updatePlot(output)
            updateHistory(output)  
            updateClosest(output)  
        })

        observeEvent(input$addlabel, {
            collected$labels[collected$current] <- input$label
            updateHistory(output, advance=FALSE)
        })

        observeEvent(input$go, {
            candidate <- as.integer(input$gotonum)
            if (is.finite(candidate) && candidate >= 1L && candidate <= length(sphere.names)) {
                collected$current <- candidate
                updatePlot(output)
                updateHistory(output)
                updateClosest(output)
            }
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

.generateDensity <- function(current.coords, density.data, ordering, collim, xlim, ...) 
# Plotting the densities and adding the point corresponding to the current coordinates.
{
    all.markers <- names(current.coords)
    all.cols <- viridis(256)

    for (m in ordering) {
        curdex <- round(approx(collim[,m], c(1, 256), xout=current.coords[m], rule=2)$y)
        if (is.null(xlim)) { 
            xlim <- collim[,m]
        }
    
        curdens <- density.data[[m]]
        plot(0, 0, type="n", xlab="", ylab="", yaxt="n", bty="n", main=all.markers[m], xlim=xlim, ylim=c(0, max(curdens$y)))
        my.x <- c(curdens$x[1]-10, curdens$x, curdens$x[length(curdens$x)]+10)
        my.y <- c(0, curdens$y, 0)
        polygon(my.x, my.y, col=all.cols[curdex], border=NA)
        lines(my.x, my.y)
    
        curpos <- current.coords[m]
        curpos <- pmin(xlim[2], pmax(xlim[1], curpos))
        points(curpos, approx(my.x, my.y, curpos, rule=2)$y, pch=16, col="red", cex=1.5)
    }

    return(invisible())
}

