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
	        coord_fixed() +    # need aspect ratio of 1!
	        geom_segment(data = sign.vectors, 
      					       aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
      					       arrow = arrow(length = unit(0.25, "cm")), 
      					       colour = "grey") +
      	  geom_text(data = sign.vectors, aes(x = NMDS1, y = NMDS2, label = env.vars), size = 5) +
      	  theme_bw() + 
      	  scale_fill_discrete(name = "Stations") # change the label for the grouping in the legend
  
	return(p)
}