extract_ordisurf_data <- function(ordisurf.obj) {
    ## Extracts xyz data for an environmental variable from an ordisurf object. 
    ## Source: https://oliviarata.wordpress.com/2014/07/17/ordinations-in-ggplot2-v2-ordisurf/
    ##
    ## Arguments: ordisurf.obj - ordisurf results
    ## Output: tibble of x, y and z values for the variable, where x = NMDS1 
    ##         coordinate, y = NMDS2 coordinate, z = value of the variable. 
    ## Dependencies: tidyverse 
    
    ## extract the ordisurf grid object (a list)
    ordi.grid <- ordisurf.obj$grid
  
    ## get x and y values from the grid
    ordi.xy <- expand.grid(x = ordi.grid$x, y = ordi.grid$y)
    
    ## get z values from the matrix  
    ordi.xy$z <- as.vector(ordi.grid$z)
    
    ## get rid of the NAs
    ordi.final <- as_data_frame(na.omit(ordi.xy))

    return(ordi.final)
}  


plot_mds_ordisurf <- function(mds.obj, ordisurf.obj, 
                              col.isolines = c("blue", "red"), 
                              col.points = "grey30", col.others = "grey30") {
    ## Plots an MDS ordination with overlain ordisurf for an environmental 
    ## variable using base R graphics (the best - and most painless - label 
    ## positioning for the isolines.
    ## 
    ## Arguments: mds.obj - MDS result (from vegan)
    ##            ordisurf.obj - ordisurf result (from vegan)
    ##            col.isolines - start and end colour for the isolines colour
    ##              palette (character vector, length = 2)
    ##            col.points - colour for the points (default: grey30)
    ##            col.others - colour for the other plot elements, e.g. axes, 
    ##              axis ticks and labels, etc. (default: grey30)
    ## Dependencies: vegan
    ## Returns: plot object (base R graphics)
    
    ## make an empty plot, without axes or axis labels
    plot(mds.obj, display = "sites", 
         type = "n", axes = FALSE, ann = FALSE)
    
    ## display the stations as points
    points(mds.obj, display = "sites", 
           col = col.points, cex = 0.7, pch = 16)  
        
    ## create custom colour palette for plotting the isolines
    isoline.cols <- colorRampPalette(col.isolines)
    
    ## get the (max) number of isolines that will be plotted automatically, to 
    ## make sure the colour ramp covers the whole extent of the isolines. 
    ## Function pretty() is called under the hood by the plot.ordisurf methods 
    ## for determining the best number and level of those lines, so we're going
    ## to apply it here.
    ncols <- length(pretty(ordisurf.obj$grid$z, n = 10))
    
    ## overlay the ordisurf object on the plot 
    plot(ordisurf.obj, add = TRUE, nlevels = 10, col = isoline.cols(ncols), labcex = 0.75)
    
    ## add in all the extra plot elements
    axis(1, col.axis = col.others, col = col.others, col.ticks = NULL)
    axis(2, col.axis = col.others, col = col.others, col.ticks = NULL)
    box(col = col.others)
    title(xlab = "NMDS1", ylab = "NMDS2", col.lab = col.others)
    
}