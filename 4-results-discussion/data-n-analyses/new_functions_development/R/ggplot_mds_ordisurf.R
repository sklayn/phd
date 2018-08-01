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


ordisurf_ggplot <- function(mds.plot, ordi.df) {
  ## an ordination in ggplot, with environmental variables overlaid as 
  ## surfaces (calculated by ordisurf in package vegan)
  
  p <- mds.plot + 
    ## plot the surface representing the environmental variable's values. Change 
    ## the number of bins to the desired number.
    geom_contour(data = ordisurf.df, 
                 lwd = 0.7,
                 aes(x = x, y = y, z = z, colour = ..level..),  
                 bins = 10)
  return(p)
} 
