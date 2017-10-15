plot_most_abnd_sp <- function(abnd.df, gr.factor, nb.sp = 15) {
  ## plots transformed abundance of the specified number of most abundant species
  ## and colors the points by the levels of the supplied factor.
  
  library(reshape2)
  library(ggplot2)
  
  # calculate the total abundance of each species, then sort in descending order
  tot.abnd <- colSums(abnd.df)
  tot.abnd.sorted <- sort(tot.abnd, decreasing = TRUE, index.return = TRUE) # make sure to return column indices, not names!
  
  # extract the given number of species (columns) from the abundance data frame - by their index
  abnd.sub <- abnd.df[, tot.abnd.sorted$ix[1:nb.sp]]
  
  # add the grouping factor column to the data frame, if it's not already there
  abnd.sub <- cbind(abnd.sub, group = gr.factor)
  
  # convert the df to long format -> easier for ggplot2 to handle  
  abnd.melted <- melt(abnd.sub, id.vars = "group") 
  
  # this doesn't keep row names, which here denote stations/replicates (does keep
  # row order, so easy to recover, if needed)
  
  ## use custom transformation for abundance values (log(y/min + 1)) - as in the 
  ## mvabund package
  # calculate the minimum non-0 value in the dataset 
  min.val <- min(abnd.melted$value[abnd.melted$value > 0]) 
  abnd.melted$value.tr <- log(abnd.melted$value / min.val + 1)
  
  ggplot(abnd.melted, aes_string(x = "variable", 
                                 y = "value.tr", 
                                 colour = "group")) + 
    geom_point(size = 2.5, alpha = 0.75) + 
    # make points somewhat transparent so that points plotted on top of one 
    # another are still visible
    scale_color_brewer(name = "Station", palette = "Set2") + 
    # reverse the order of x axis, so highest-contributing species are on top
    scale_x_discrete(name = "", limits = rev(levels(abnd.melted$variable))) + 
    scale_y_continuous(name = "Abundance (log(y/min + 1))") + 
    coord_flip() + # put species on y axis - easier to read
    theme_bw() + 
    theme(axis.text.x = element_text(size = rel(1.3)), 
          axis.text.y = element_text(size = rel(1.3)), 
          legend.text = element_text(size = rel(1.2)), 
          legend.title = element_text(size = rel(1.2)))
}
