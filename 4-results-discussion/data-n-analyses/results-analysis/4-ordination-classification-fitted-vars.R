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

## envfit (vegan) on the (imputed) environmental data. 
# use the long-term water column data - summarized by station
water.sand.by.st <- ddply(water.sand.imp.df[, !names(water.sand.imp.df) %in% c("year", "month")], 
                          .(station), 
                          colwise(mean, .cols = is.numeric))

# repeat each row - to match the number of replicate zoobenthic samples 
# (unfortunately, hardcoded here and below for expediency)
water.sand.envfit <- water.sand.by.st[rep(seq_len(nrow(water.sand.by.st)), each = 9), ]
rownames(water.sand.envfit) <- 1:nrow(water.sand.envfit)

# for the other environmental parameters, summarize the imputed datasets by month
# and year (these were only measured in 2013-2014 anyway) 
other.env.sand.by.st <- ddply(other.env.sand.imp.df, 
                              .(station, year, month), 
                              colwise(mean, .cols = is.numeric))

# repeat each row 3 times to match the number of replicate zoobenthic samples
other.env.sand.envfit <- other.env.sand.st[rep(seq_len(nrow(other.env.sand.st)), each = 3), ]
rownames(other.env.sand.envfit) <- 1:nrow(other.env.sand.envfit)

# quick workspace cleanup
rm(water.sand.by.st, other.env.sand.by.st)

# add the heavy metals 
# import heavy metal data, if not already imported
# heavy.metals.sand <- read.csv(file.path(data.dir, "heavy-metals-sand.csv"), header = TRUE) 

# repeat each row - to match the number of replicate zoobenthic samples
heavy.metals.sand.envfit <- heavy.metals.sand[rep(seq_len(nrow(heavy.metals.sand)), each = 9), ]
rownames(heavy.metals.sand.envfit) <- 1:nrow(heavy.metals.sand.envfit)

# remove the month and year from the heavy metals dataset (only measured once, so 
# not relevant)
heavy.metals.sand.envfit <- heavy.metals.sand.envfit[, !names(heavy.metals.sand.envfit) %in% c("month", "year", "station")]

# join all 3 (water column, sediments and heavy metals) together
env.all.sand.envfit <- cbind(other.env.sand.envfit,
                             heavy.metals.sand.envfit, 
                             water.sand.envfit[, !names(water.sand.envfit) == "station"])

## SAVE IN CASE WE NEED IT 
saveRDS(env.all.sand.envfit, file = file.path(save.dir, "env-all-for-envfit_clean_sand.rds"))

# set seed if we have to repeat the calculation
set.seed(1)
envfit.mds.sand <- envfit(mds.sand,
                       env.all.sand.envfit[, !names(env.all.sand.envfit) %in% c("station", "month", "year")], 
                       permutations = 999)
                      
# apply Bonferroni correction for multiple testing to the p-values, because there is a large number 
# of tested variables (calls custom function!).
envfit.mds.sand <- p_adjust_envfit(envfit.mds.sand)

## SAVE THE ENVFIT
saveRDS(envfit.mds.sand, file = file.path(save.dir, "envfit-mds_sand.rds"))

# extract significant variables & order by r2
envfit.sand.sign.vars <- extract_envfit_scores(envfit.mds.sand, p = 0.05, r2 = TRUE)
envfit.sand.sign.vars <- arrange(envfit.sand.sign.vars, pvals, desc(r2))


## PLOT VARIABLES AS SURFACES OVERLAID ON THE MDS (ORDISURF)

# ordisurf needs the original numeric values of the environmental variables chosen
# for plotting (= the significant variables from envfit, p < 0.05)
sign.vars.sand.ordisurf <- subset(env.all.sand.envfit, 
                                  select = names(env.all.sand.envfit) %in% envfit.sand.sign.vars$vars)

dim(sign.vars.sand.ordisurf) # just to check that everything we need is there

# apply ordisurf sequentially to all significant environmental variables, 
# which serves no purpose (get back a list of ordisurf objects where each element
# is an environmental variable)
ordisurf.list.all.sand <- apply(sign.vars.sand.ordisurf, 
                                MARGIN = 2, 
                                FUN = function(x) ordi <- ordisurf(mds.sand ~ x, plot = FALSE)) 

# check out the summaries of the fits
lapply(ordisurf.list.all.sand, summary)

## SAVE FITTED ORDISURFS AS WELL (in case)
saveRDS(ordisurf.list.all.sand, file = file.path(save.dir, "ordisurf-mds_sand.rds"))

# rearrange list to have the plots in the desired order.
# Here: first the sediment parameters, then the water column parameters, and 
# finally the heavy metals.
names(ordisurf.list.all.sand)
ordisurf.list.all.sand <- ordisurf.list.all.sand[c("Ninorg", "Ntot", "NH4", "NO2", "PO4", 
                                                   "O2.bottom", "O2.average", "Secchi.depth",
                                                   "seston", "dist.innermost", "depth", "sand",
                                                   "gravel", "sorting", "mean.grain.size",
                                                   "org.matter", "heavy.metals.noFe", "heavy.metals.all",
                                                   "Mn", "Pb", "Ni", "Cu")]

var.labels.sand <- c("N-inorganic", "N-total", "NH4", "NO2", "PO4", 
                     "O2 bottom", "O2 average", "Secchi depth",
                     "seston", "distance to innermost station", "depth", "% sand",
                     "% gravel", "sorting", "mean grain size",
                     "organic matter", "heavy metals (no Fe)", "total heavy metals",
                     "Mn", "Pb", "Ni", "Cu")

# set the file name and properties for the output graph --> REDO TO FIT MULTIPLE PLOTS/PAGE!
pdf(file = file.path(figs.dir, "mds_ordisurf_sand_most_sign_vars.pdf"), 
    paper = "a4",
    useDingbats = FALSE)

# plot all variables, using the custom plot_mds_ordisurf function, and adding the
# corresponding main title (variable name) on each subplot
mapply(function(m, n) {
        plot_mds_ordisurf(mds.sand, m)
        title(main = n, col.main = "grey28")
      }, 
      ordisurf.list.all.sand, 
      var.labels.sand)

dev.off()


# clean up workspace 
rm(envfit.sand.sign.vars, 
   ordisurf.list.all.sand, 
   sign.vars.sand.ordisurf, 
   var.labels.sand,
   heavy.metals.sand.envfit, 
   water.sand.envfit, 
   other.env.sand.envfit)


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
