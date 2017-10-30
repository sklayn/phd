plot_mds_ordisurf <- function(mds.obj, 
                              ordisurf.obj, 
                              st.labels = NULL,
                              col.isolines = c("blue", "red"), 
                              col.points = "grey28",
                              col.others = "grey28") {
  ## plot the mds with overlain ordisurf for an environmental variable - with base 
  ## graphics, because the isolines with labels are better-looking than in 
  ## ggplot2 or lattice.
  ## 
  
  # create an empty plot, without axes or axis labels
  plot(mds.obj, display = "sites", type = "n", 
       axes = FALSE,
       ann = FALSE)
  
  if(is.null(st.labels)) {
    # display the stations as points, if labels are not provided as input
    points(mds.obj, display = "sites", col = col.points, cex = 0.7)  
    
  } else {
    # display the stations as text, using the labels provided 
    text(mds.obj, display = "sites", labels = st.labels, col = col.points, cex = 0.7)
  }
  
  # create custom color palette for plotting the isolines
  isoline.cols <- colorRampPalette(col.isolines)
  
  # get the (max) number of isolines that will be plotted automatically, to make
  # sure the color ramp covers the whole extent of the isolines. Function pretty()
  # is called under the hood by the plot.ordisurf methods for determining the 
  # best number and level of those lines.
  ncols <- length(pretty(ordisurf.obj$grid$z, n = 10))
  
  # overlay the ordisurf object 
  plot(ordisurf.obj, add = TRUE, nlevels = 10, col = isoline.cols(ncols))
  
  # add in all the extra plot elements
  axis(1, col.axis = col.others, col = col.others, col.ticks = NULL)
  axis(2, col.axis = col.others, col = col.others, col.ticks = NULL)
  box(col = col.others)
  title(xlab = "NMDS1", ylab = "NMDS2", col.lab = col.others)
  
}