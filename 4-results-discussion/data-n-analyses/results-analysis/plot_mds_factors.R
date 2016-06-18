plot_mds_factors <- function(zoo.mds, zoo.envfit) {
	## Plots MDS with fitted environmental factors as vectors (factors with 
  ## significant correlations).

	# import libraries
	library(ggplot2)
	library(grid)
	library(vegan)

	# extract the MDS scores and create data frame for plotting
	scrs <- as.data.frame(scores(zoo.mds, display = "sites"))

	# add vector of labels to use for grouping - FIX THIS - CALLS "SITES" IN GLOBAL ENVIRONMENT!
	scrs <- cbind(scrs, sites) 		
	 
	vector.scrs <- extract_envfit_scores(zoo.envfit)

	ggplot(scrs) + 
	  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = sites)) +
	  coord_fixed() +    # need aspect ratio of 1!
	  geom_segment(data = sign.vectors, 
					       aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
					       arrow = arrow(length = unit(0.25, "cm")), 
					       colour = "grey") +
	  geom_text(data = sign.vectors, aes(x = NMDS1, y = NMDS2, label = env.vars), size = 5) +
	  theme_bw() + 
	  scale_fill_discrete(name = "Stations") # change the label for the grouping in the legend

}