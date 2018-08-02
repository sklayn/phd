#### Plot ordisurf results with ggplot
#### Source: https://oliviarata.wordpress.com/2014/07/17/ordinations-in-ggplot2-v2-ordisurf/

library(vegan)

#data - from the vegan package:
data(mite)
data(mite.env)
set.seed(505) #set seed for reproducibility
names(mite.env)
 
#NMDS
meta.nmds.mite <- metaMDS(mite) # no transformation of species data is made here prior to bray 
                                # curtis dissimilarities being calculated. (Bray Curtis is the default in R).
meta.nmds.mite #prints very basic dsecription
stressplot(meta.nmds.mite) # To gain the stress plot for stress values for your MDS
 
#ordisurf:
ordi <- ordisurf(meta.nmds.mite ~ mite.env$SubsDens) #created the ordisurf object
ordi.grid <- ordi$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
ordi.mite.na #looks ready for plotting!
 
 
#data for plotting
##NMDS points
mite.NMDS.data <- mite.env #there are other ways of doing this. But this is the way I do it for ease of plotting
mite.NMDS.data$NMDS1 <- meta.nmds.mite$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
mite.NMDS.data$NMDS2 <- meta.nmds.mite$points[ ,2]
 
#setting up a plot in base graphics
plot(mite.NMDS.data$NMDS1, mite.NMDS.data$NMDS2,
     xlab="NMDS1", ylab="NMDS2",  
     pch = c(1,16)[as.factor(mite.NMDS.data$Topo)]) #symbols coded by factor
legend(title = "Topo", "bottomleft", levels(mite.NMDS.data$Topo), pch = c(1,16))
ordisurf(meta.nmds.mite~mite.env$SubsDens, col = "grey50", add = TRUE)
with(ordi.grid, contour(x = x, y = y, z = z, add=TRUE)) # you can use this as a check to make sure your ordigrid is doing the same thing as an ordisurf (ie, you've extracted the components correctly)
 
## Plotting in ggplot2
mite.ordisurf.ggplot<-ggplot(mite.NMDS.data, aes(x = NMDS1, y = NMDS2))+
  stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, colour = rev(..level..)),
               binwidth = 2)+ #can change the binwidth depending on how many contours you want
  geom_point(size = 3, alpha = 0.8,  aes(shape = Topo)) + #plots the NMDS points, with shape by topo type
  theme_bw() + #for aesthetics
  scale_shape_manual("Topo Type",  values = c(1,16)) + #sets the name of the legend for shape, and says which symbols we want (equivalent to the pch command in plot)
  labs(colour = "Substrate Density")+ #another way to set the labels, in this case, for the colour legend
  scale_colour_gradient(high = "darkgreen", low = "darkolivegreen1")+ #here we set the high and low of the colour scale.  Can delete to go back to the standard blue, or specify others
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = "bottom", #legend at the bottom
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "centre")
mite.ordisurf.ggplot
#ggsave(plot = mite.ordisurf.ggplot, file= "mite.ordisurf.ggplot.png") #to save = see ?ggsave. does not work for base graphics



##############################################################################################################################

## My cleaned-up version (unused, because want labels on contour lines, but no time to look into making them more presentable). Also, way too complicated and unnecessary.. 

envdata_ordisurf_clean <- function(mds.obj, env.vars) {
  ## prepares environemntal data for plotting as surfaces overlaid on an 
  ## ordination (package vegan)
  ## mds.obj - ordination object; env.vars - data frame of environmental 
  ## variables (observations x variables)
  ## output: named list of data frames containing x, y and z values from
  ## ordisurf, where each element corresponds to a variable.
  
  library(vegan)
  
  # apply ordisurf sequentially to all environmental variables
  ordi.list.all <- apply(env.vars, MARGIN = 2, FUN = function(x) ordi <- ordisurf(mds.obj ~ x, plot = F))
  
  ## prepare the ordisurf data for plotting
  # extract the ordisurf grid objects
  ordi.list.grids <- lapply(ordi.list.all, function(x) x$grid)
  
  # they are lists themsleves, though - cannot be plotted directly.
  # Get x and y values
  ordi.list.xy <- lapply(ordi.list.grids, function(m) expand.grid(x = m$x, y = m$y))
  
  # get z values from the matrix  
  ordi.list.xyz <- mapply(function(m, n){
    m$z <- as.vector(n$z)
    return(m)
  }, 
  ordi.list.xy,
  ordi.list.grids,
  SIMPLIFY = F) 
  
  # get rid of the NAs
  ordi.list.final <- lapply(ordi.list.xyz, function(x) data.frame(na.omit(x)))
  
  return(ordi.list.final)
  
}  


ordisurf_ggplot <- function(mds.obj, ordisurf.df, stations) {
  ## an ordination in ggplot, with environmental variables overlaid as 
  ## surfaces (calculated by ordisurf in package vegan)
  
  # extract the NMDS coordinates of the points in new columns, but in the same 
  # data frame as the environmental variables (for ease of plotting in ggplot)  
  
  nmds.data <- as.data.frame(scores(mds.obj, display = "sites"))
  
  p <- ggplot() + 
    # plot the surface representing the environmental variable's values. Change 
    # the number of bins to the desired number.
    geom_contour(data = ordisurf.df, 
                 lwd = 0.7,
                 aes(x = x, y = y, z = z, colour = ..level..),  
                 bins = 10) + 
    # plot the points corresponding to the station/replicates
    geom_text(data = nmds.data, 
              aes(x = NMDS1, y = NMDS2), 
              label = stations, 
              size = 3, 
              alpha = 0.8, 
              check_overlap = T) + 
    
    theme_bw()
  
  return(p)
} 
