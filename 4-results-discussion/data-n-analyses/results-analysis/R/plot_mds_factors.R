extract_envfit_scores <- function(envfit.obj, pval = 0.05, r2 = FALSE) {
  # helper function to extract NMDS vector scores from envfit objects according 
  # to the specified p-value (p < 0.05 by default). Optionally also extracts 
  # the r2 value.
  
  # extract the NMDS scores for the fitted environmental parameters 
  vector.scrs <- as.data.frame(scores(envfit.obj, display = "vectors"))
  
  # convert all the variables characterizing the vectors in the envfit object 
  # to a list 
  list.scrs <- as.list(envfit.obj$vectors)
  
  if(r2) {
    # get r2 values 
    vector.scrs$r2 <- list.scrs$r    
  }
  
  # get the p-values and subset the data frame according to them
  vector.scrs$pvals <- list.scrs$pvals
  sign.vectors <- subset(vector.scrs, pvals < pval)
  
  # clean up a little: get rid of the row names - add them as a variable in their
  # own column instead
  sign.vectors <- cbind(vars = rownames(sign.vectors), sign.vectors)
  rownames(sign.vectors) <- 1:nrow(sign.vectors)
  
  return(sign.vectors)
}


plot_mds <- function(mds.obj, stations, label.groups = TRUE) {
  ## Plots MDS result in ggplot2. 
  ## Arguments: mds.obj - mds result object.
  ##            stations - station labels to use in grouping the points on 
  ##              the plot; need to be in the same order as those used for 
  ##              the mds.
  ## Returns a ggplot object.
  
  # import libraries
  library(ggplot2)
  library(plyr)
  library(vegan)
  
  # extract the MDS scores and create data frame for plotting
  scrs <- as.data.frame(scores(mds.obj, display = "sites"))
  
  # add labels to use for grouping
  scrs <- cbind(scrs, stations)
  names(scrs) <- c("NMDS1", "NMDS2", "stations")
  
  # make a custom colour scale (to match text label colours with point colours) 
  custom.cols <- brewer.pal(n = length(levels(scrs$stations)), name = "Set2")
  names(custom.cols) <- levels(scrs$stations) 
  custom.col.scale <- scale_color_manual(name = "", values = custom.cols)
  
  p <- ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = stations)) + 
    geom_point(size = 2, position = "identity", show.legend = FALSE) +
    stat_ellipse(show.legend = FALSE) +
    custom.col.scale + 
    theme_bw() + 
    
  # annotate the groups (inside the ellipses) if desired (only if plot doesn't 
  # become too cluttered)
  if(label.groups) {
    # first, calculate the label placement (=> mean of each group)
    mds.group.labs <- ddply(scrs, 
                            .(stations), 
                            summarize, 
                            NMDS1.m = mean(NMDS1), 
                            NMDS2.m = mean(NMDS2)) 
    
    # add labels to plot, matching their colors with the group colors
    p <- p + annotate("text", 
                      x = mds.group.labs$NMDS1.m, 
                      y = mds.group.labs$NMDS2.m, 
                      label = mds.group.labs$stations, 
                      colour = custom.cols, 
                      size = 5) 
  }
  
  return(p)
}



plot_mds_factors <- function(mds.obj, envfit.obj, stations, p = 0.05) {
	## Plots MDS with fitted environmental factors as vectors (factors with 
  ## significant correlations).
  ## Arguments: mds.obj - mds result object.
  ##            envfit.obj - corresponding envfit result object.
  ##            stations - station labels to use in grouping the points on 
  ##              the plot; need to be in the same order as those used for 
  ##              the mds and envfit.
  ##            p - p-value for the desired level of significance of the 
  ##              fitted variables.
  ## NB: Calls a helper function (extract_envfit_scores) to extract envfit 
  ## vector scores from an envfit object. 
  ## 
  ## Returns a ggplot object.

	# import libraries
	library(ggplot2)
	library(grid)
	library(vegan)

	# extract the MDS scores and create data frame for plotting
	scrs <- as.data.frame(scores(mds.obj, display = "sites"))

	# add labels to use for grouping
	scrs <- cbind(scrs, stations) 		
	 
	vector.scrs <- extract_envfit_scores(envfit.obj, pval = p)

	p <- ggplot(scrs) + 
	        geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = stations)) +
	        scale_color_brewer(palette = "Set2") +
	        coord_fixed() +    # need aspect ratio of 1!
	        geom_segment(data = sign.vectors, 
      					       aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
      					       arrow = arrow(length = unit(0.25, "cm")), 
      					       colour = "grey") +
      	  geom_text(data = sign.vectors, aes(x = NMDS1, y = NMDS2, label = vars), size = 5) +
      	  theme_bw()
	
	return(p)
}