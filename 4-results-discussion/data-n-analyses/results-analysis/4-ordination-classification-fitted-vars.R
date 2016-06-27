## Ordintaion (nMDS) on the abundance data
# set seed in case we have to rerun the computation
set.seed(1)
mds.sand <- metaMDS(num.zoo.abnd.sand)

# basic summary of the MDS
mds.sand

# diagnostic plots for the MDS
pdf(file = file.path(figs.dir, "explor_diagnostic-plots-mds_sand.pdf"), 
    paper = "a4r", 
    width = 12, 
    height = 12, 
    useDingbats = FALSE)

par(mfrow = c(1, 2))

# stressplot
stressplot(mds.sand)

# goodness-of-fit plot 
# first plot the nMDS ordination with sites
plot(mds.sand, display = 'sites', type = 't', main = 'Goodness of fit') 
# then, add the points with size reflecting goodness of fit (smaller = better fit)
points(mds.sand, display = 'sites', cex = goodness(mds.sand)*200) 

dev.off()

par(mfrow = c(1, 1))

# save the MDS result 
saveRDS(mds.sand, file = file.path(save.dir, "mds_sand.rds"))

# plot and save the MDS (ggplot2)
pdf(file = file.path(figs.dir, "explor_mds_sand.pdf"), useDingbats = FALSE)
plot_mds(mds.sand, factors.zoo.sand$stations)
dev.off()


## Fit environmental variables to the ordination (datasets from environmental 
## data script)

# envfit (vegan) on the imputed datasets (env. variables subset only)
# set seed if we have to repeat the calculation
set.seed(1)
env.mds.sand <- lapply(env.imp.all.sand, function(x) {
  # subset each data frame to only the variable columns
  envfit(mds.sand, 
         subset(x, select = -c(station, month, year)), 
         permutations = 999)
})

# apply Bonferroni correction for multiple testing to the p-values, because there is a large number 
# of tested variables (calls custom function!).
env.mds.sand <- lapply(env.mds.sand, p_adjust_envfit)

## SAVE THE LIST OF ENVFITS
saveRDS(env.mds.sand, file = file.path(save.dir, "envfit-mds-clean_sand.rds"))

# extract the most significant variables and positions from the list of envfits 
# (calls custom helper function)
env.sand.sign.count <- lapply(env.mds.sand, sign_vars_pos_envfit)

# get the most significant environmental variables according to the envfits 
# performed on the imputed datasets (custom helper function)
sign.vars.sand <- sign_vars_freq(env.sand.sign.count, target.freq = 10)

## PLOTTING VARIABLES AS SURFACES OVERLAID ON THE MDS
# get the original numeric values of the environmental variables chosen for plotting
# from all the imputed datasets

# only use env. variables (no factors) from each df
vars.for.ordisurf.sand <- lapply(lapply(env.imp.all.sand, function(y) subset(y, select = -c(station, month, year))), # only env.variables from each df
                                 function(x) { 
                                   x <- subset(x, select = names(x) %in% levels(sign.vars.sand[[1]])) # subset according to vector of chosen sign. variables
                                   # add a numeric identifier to be used later for aggregating
                                   # to avoid having to reorder later
                                   x$id <- rownames(x)
                                   x$id <- as.numeric(x$id)
                                   return(x)
                                 }) 

# combine into a single data frame and average each variable by id (= replicate)
vars.for.ordisurf.sand <- do.call("rbind", vars.for.ordisurf.sand) 
vars.for.ordisurf.sand <- ddply(vars.for.ordisurf.sand, .(id), colwise(mean)) 

# apply ordisurf sequentially to all environmental variables except the first (= id), 
# which serves no purpose (get back a list of ordisurf objects where each element
# is an environmental variable)
ordisurf.list.all.sand <- apply(vars.for.ordisurf.sand[-1], 
                                MARGIN = 2, 
                                FUN = function(x) ordi <- ordisurf(mds.sand ~ x, plot = FALSE)) 

# check out the summaries of the fits
lapply(ordisurf.list.all.sand, summary)

## SAVE FITTED ORDISURFS AS WELL (in case)
saveRDS(ordisurf.list.all.sand, file = file.path(save.dir, "ordisurf-mds_sand.rds"))

# rearrange list to have the plots in the desired order.
# Here: on the first row of the plot will be the sediment parameters, on the 
# second - the water column parameters, and on the third - the heavy metals.
names(ordisurf.list.all.sand)
ordisurf.list.all.sand <- ordisurf.list.all.sand[c("sorting", "mean.grain.size", "dist.innermost", "depth",  
                                                   "Secchi.depth", "O2.bottom", "O2.average", "heavy.metals.noFe")]

var.labels.sand <- c("sorting", "mean grain size", "distance to innermost station", "depth", 
                     "Secchi depth", "O2 bottom", "O2 average", "heavy metals (no Fe)")

# set the file name and properties for the output graph
pdf(file = file.path(figs.dir, "mds_ordisurf_sand_most_sign_vars.pdf"), 
    paper = "a4",
    width = 12,
    height = 12,
    useDingbats = FALSE)

# modify par to fit all plots on one page (here, 4 plots per row)
par(mfrow = c(4, 2))

# plot all variables, using the custom plot_mds_ordisurf function, and adding the
# corresponding main title (variable name) on each subplot
mapply(function(m, n) {
  plot_mds_ordisurf(mds.sand, m)
  title(main = n, col.main = "grey28")
}, 
ordisurf.list.all.sand, 
var.labels.sand)

dev.off()

# return the graphics device to the original settings
par(mfrow = c(1, 1))

# clean up workspace 
rm(env.sand.sign.count, 
   ordisurf.list.all.sand, 
   sign.vars.sand, 
   var.labels.sand, 
   vars.for.ordisurf.sand)


## Classification of the communities - FINISH THIS!!
dendr.sand <- hclust(vegdist(sqrt(num.zoo.abnd.sand)), "average")

# plot and examine the dendrogram 
plot(dendr.sand, hang = -1)
rect.hclust(dendr.sand, k = 4) # here, it appears to have about 4 meaningful groups
gr.dendr.sand <- cutree(dendr.sand, k = 4)


# reorder dendrogram by some variable
dendr.by.O2 <- with(sign.vars.mean, reorder(dendr.sand, O2.bottom))
plot(dendr.by.O2, label = round(sign.vars.mean$O2.bottom, 2), hang = -1, main = "O2 bottom")

# numerical analysis of the grouping
anova(lm(O2.bottom ~ gr.dendr.sand, data = sign.vars.mean))


## ANOSIM on the groups identified by the MDS and the classification
## SIMPER to id the most significant species


