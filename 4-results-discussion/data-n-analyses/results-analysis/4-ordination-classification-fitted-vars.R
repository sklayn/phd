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

###############
# repeat each row 3 times to match the rows of the mds object and fix row names
# (very ugly, but working, hack)
env.imp.all.sand <- env.imp.all.sand[rep(seq_len(nrow(env.imp.all.sand)), each = 3), ]
rownames(env.imp.all.sand) <- 1:nrow(env.imp.all.sand)

# SAVE THIS IF NECESSARY FOR FUTURE USE 
saveRDS(env.imp.all.sand, file.path(save.dir, "env-data-imputed-clean_sand.rds"))


set.seed(1)
env.mds.sand <- dlply(env.imp.all.sand[, !names(env.imp.all.sand) %in% c("station", "month", "year")], 
                      .(.imp), 
                      function(x) {
                        x$.imp <- NULL
                        envfit(mds.sand, x, permutations = 999)
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

# only use env. variables (no factors) from each the imputed data df
sign.vars.mean.sand <- ddply(env.imp.all.sand,
                             .(.imp),
                             # subset according to vector of chosen significant variables
                             function(x) {
                               x <- subset(x, select = names(x) %in% levels(sign.vars.sand[[1]]))
                               x$id <- as.numeric(rownames(x))
                               return(x)
                               }
                            ) 

# average each variable by id (= replicate)
sign.vars.mean.sand <- ddply(sign.vars.mean.sand, .(id), colwise(mean, .cols = is.numeric)) 

# apply ordisurf sequentially to all environmental variables except the first (= id), 
# which serves no purpose (get back a list of ordisurf objects where each element
# is an environmental variable)
ordisurf.list.all.sand <- apply(sign.vars.mean.sand[-1], 
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
   sign.vars.mean.sand)


## Classification of the communities 

# dendrogram of dissimilarities between samples
set.seed(1)
dendr.sand <- hclust(vegdist(sqrt(num.zoo.abnd.sand), method = "bray"), 
                     "average")

# add station names as labels 
dendr.sand$labels <- factors.zoo.sand$stations

# plot and examine the dendrogram
pdf(file.path(figs.dir, "explor_dendrogram-sand.pdf"), useDingbats = FALSE)
plot(dendr.sand, hang = -1, 
     main = "", ylab = "Distance (Bray-Curtis)", xlab = "")
rect.hclust(dendr.sand, k = 4) # here, it appears there are about 4 distinct groups

dev.off()

gr.dendr.sand <- cutree(dendr.sand, k = 4) # or the tree can be cut at any given height, too

# reorder dendrogram by some variable (variables - same as aggregated data frame 
# used for ordisurf & plotting over ordination)
## => plots can be redone - leaves colored by cluster; values of the variable used
##    for rearranging the tree plotted below as a colored bar..   
dendr.sand.by.O2 <- with(sign.vars.mean, reorder(dendr.sand, O2.bottom))

pdf(file.path(figs.dir, "explor_dendrogram-by-O2_sand.pdf"), useDingbats = FALSE)
plot(dendr.sand.by.O2, hang = -1, 
     main = "O2 bottom", xlab = "", ylab = "Distance (Bray-Curtis)")
rect.hclust(dendr.sand.by.O2, k = 4)
dev.off()

# numerical analysis of the grouping
anova(lm(O2.bottom ~ gr.dendr.sand, data = sign.vars.mean))



## ANOSIM - groups = stations.  
## This is a non-parametric permutation procedure applied to the rank (similarity)
## matrix underlying the ordination or classification of the samples. R statistic:
##  -1 to 1; 1 = all replicates within sites are more similar to each other than
## to any other replicate from a different site; 0 = H0 is true (the same average
## similarities between and within sites). Usually 0 < R < 1 => some degree of
## difference observed between sites.
## Best course of analysis: 1) global ANOSIM - overall difference between groups; 
## if significant - 2) where does the main between-group difference come from? 
## => examine R values for each pairwise comparison: large = complete separation, 
## small - little or no difference.
anosim.sand <- anosim(vegdist(sqrt(num.zoo.abnd.sand), method = "bray"), 
                      grouping = factors.zoo.sand$stations)

anosim.sand


## import and reclassify environmental variables to use for grouping - can be used 
## to colour/order classification or ordination, too 
env.qualit <- read.csv(file.path(data.dir, "env-qualit-vars_sand.csv"), header = TRUE)
str(env.qualit)
names(env.qualit)

# copy each row (to match replicates in the zoo abundancee date frame) 
env.qualit <- env.qualit[rep(seq_len(nrow(env.qualit)), each=3),]
rownames(env.qualit) <- 1:nrow(env.qualit)


## repeat the ANOSIM using these new groups 
anosim.sand2 <- apply(env.qualit, 2, function(x) 
                                      anos <- anosim(vegdist(sqrt(num.zoo.abnd.sand), method = "bray"), 
                                                             grouping = x))
anosim.sand2

## SIMPER to id the species with the highest contribution to the differences
## between groups. 
## Good discriminating species - high contribution + small sd. 
simper.sand <- simper(sqrt(num.zoo.abnd.sand), 
                      group = factors.zoo.sand$stations)
summary(simper.sand, ordered = TRUE)

simper.sand2 <- apply(env.qualit, 2, function(x) simp <- simper(sqrt(num.zoo.abnd.sand), 
                                                                group = x))
