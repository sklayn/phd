## Graphical community structure analyses - diversity profiles, dominance curves
## (continued from community analysis script - references objects defined there)

## Diversity profiles - graphical representation of the shape of the community;
# show how the perceived diversity changes as the emphasis shifts from common to rare
# species - we can judge their respective contributions in the community composition.
diversity.profiles.sand <- diversity_profiles(num.zoo.abnd.sand, q = 50)
# diversity.profiles.zostera <- diversity_profiles(num.zoo.abnd.zostera, q = 50)

# plot the diversity profiles (all samples - in panels by station); save to file.
# !NB set useDingbats = FALSE, otherwise points might get transformed to letters 
# when the pdf file is opened in another program!
pdf(file = file.path(figs.dir, "classical-diversity-profiles_all_sand.pdf"), useDingbats = FALSE)
plot_div_profiles(diversity.profiles.sand, stations.sand, one.panel = FALSE)
dev.off()

# pdf(file = file.path(figs.dir, "classical-diversity-profiles_all_zostera.pdf"), useDingbats = FALSE)
# plot_div_profiles(diversity.profiles.zostera, stations.zostera, one.panel = FALSE)
# dev.off()


# save diversity profiles to file 
write.csv(diversity.profiles.sand, 
          file.path(save.dir, "div-profiles-sand_classical.csv"))

# write.csv(diversity.profiles.zostera, 
#           file.path(save.dir, "div-profiles-zostera_classical.csv"))



# add a measure of similarity between species/taxa. 
# If not already in the workspace, import taxonomic data from which distances will
# be derived. 
# zoo.taxa <- read.csv(file.path(data.dir, "zoo-taxonomy.csv"), header = T, row.names = 1)

# calculate diversity profiles including a measure of similarity between species -
# allows for a lot more meaningful ecological and biodiversity comparisons; then
# calculate and plot diversity profiles again (Leinster & Cobbold, 2012)
weighted.profiles.sand <- weighted_div_profiles(num.zoo.abnd.sand, zoo.taxa)
# weighted.profiles.zostera <- weighted_div_profiles(num.zoo.abnd.zostera, zoo.taxa)

# plot the profiles in panels by station and save
pdf(file = file.path(figs.dir, "weighted-diversity-profiles_all_sand.pdf"), useDingbats = FALSE)
plot_div_profiles(weighted.profiles.sand, stations.sand, one.panel = FALSE)
dev.off()

# pdf(file = file.path(figs.dir, "weighted-diversity-profiles_all_zostera.pdf"), useDingbats = FALSE)
# plot_div_profiles(weighted.profiles.zostera, stations.zostera, one.panel = FALSE)
# dev.off()

# same procedure, but using the data averaged by station - first column is the station names!
aver.profiles.sand <- weighted_div_profiles(summary.abnd.sand[-1], zoo.taxa)
# aver.profiles.zostera <- weighted_div_profiles(summary.abnd.zostera[-1], zoo.taxa)

# plot the profiles; save as pdf 
pdf(file = file.path(figs.dir, "weighted-diversity-profiles_aver_sand.pdf"), useDingbats = FALSE)
plot_div_profiles(aver.profiles.sand, stations.sand, one.panel = TRUE)
dev.off()

# pdf(file = file.path(figs.dir, "weighted-diversity-profiles_aver_zostera.pdf"), useDingbats = FALSE)
# plot_div_profiles(aver.profiles.zostera, stations.zostera, one.panel = TRUE)
# dev.off()


# plot all the weighted diversity profiles by station, and add the average profiles
# on the same graph
pdf(file = file.path(figs.dir, "weighted-div-profiles_allaver_sand.pdf"), useDingbats = FALSE)
plot_div_profiles_w_aver(weighted.profiles.sand, 
                         aver.profiles.sand, 
                         stations.sand, 
                         col.profiles = c("skyblue", "royalblue"))

dev.off()


# pdf(file = file.path(figs.dir, "weighted-div-profiles_allaver_zostera.pdf"), useDingbats = FALSE)
# plot_div_profiles_w_aver(weighted.profiles.zostera, 
#                          aver.profiles.zostera, 
#                          stations.zostera, 
#                          col.profiles = c("skyblue", "royalblue"))
# 
# dev.off()



# Abundance-biomass comparison (ABC) curves. These are cumulative k-dominance 
# curves plotted from the species abundance distribution and from the species 
# biomass distribution -> y-axis scale: 0 to 100. Used to determine visually 
# the level of disturbance of the communities. Under stable conditions, the 
# competitive dominants in the communities are K-strategists - large-bodied, 
# conservative species with long life-span, usually dominant in terms of biomass;
# when the communities suffer a perturbation, the smaller, opportunistic 
# species are favored - small body size and short life-span, can be numerically
# significant, but rarely represent a large proportion of the biomass. Thus, under
# undisturbed conditions, the biomass curve is elevated over the abundance 
# curve; under moderate disturbance, the biomass and the abundance curves are 
# coincident and may cross each other one or several times; under severe disturbance,
# the abundance curve lies above the biomass curve throughout its length.

# import the biomass data (if not already imported by then)
zoo.biomass.sand <- import_zoo_data(data.dir = data.dir, 
                                    zoo.data = "zoo-biomass-sand.csv", 
                                    station.names = stations.sand, 
                                    repl = 3)

# quick check for import errors, etc.
str(zoo.biomass.sand)
names(zoo.biomass.sand)

## SAVE THE CLEANED BIOMASS DATA
write.csv(zoo.biomass.sand, file.path(save.dir, "biomass_sand_clean.csv"), row.names = FALSE)

# make a subset of only the numeric biomass data for the species present in the 
# current dataset, and another of the mean biomass per station  
num.zoo.biomass.sand <- zoo.biomass.sand[sapply(zoo.biomass.sand, is.numeric)]
num.zoo.biomass.sand <- num.zoo.biomass.sand[colSums(num.zoo.biomass.sand) > 0]

summary.biomass <- ddply(zoo.biomass.sand, .(stations), colwise(mean, .cols = is.numeric))

# calculate the ABC curves on the transposed abundance and biomass data (the 
# function accepts input as species x samples tables).
abc.sand <- mapply(abc, 
                   as.data.frame(t(num.zoo.abnd.sand)), 
                   as.data.frame(t(num.zoo.biomass.sand)), 
                   SIMPLIFY = FALSE)

# calculate and plot the average ABC curves - by station
abc.aver.sand <- mapply(abc, 
                        as.data.frame(t(summary.abnd[-1])), 
                        as.data.frame(t(summary.biomass[-1])),
                        SIMPLIFY = FALSE)

# each sublist (= station/replicate) consists of 2 elements: a data frame with 
# the cumulative abundance and biomass for the ranked species, and a W value. 
# For the plots, we are interested in the first element of each of those lists.
abc.plots <- lapply(sapply(abc.aver.sand, function(x) x[1]), plot_dom_curves)
do.call(grid.arrange, abc.plots)

pdf(file = file.path(figs.dir, "abc_aver_sand.pdf"), useDingbats = FALSE)
plot_dom_curves_facets(abc.aver.sand, stations = stations.sand)
dev.off()

# you can do a transformation of the y-axis for your k-dominance curves: because
# very often, they approach a cumulative frequency of 100 % for a large part of 
# their length, in highly dominated communities - often after the first 2-3 top-ranked
# species, making the forms of the abundance and biomass curves difficult to 
# distinguish. So, you can transform the y-axis so that the cumulative values are
# closer to linearity - e.g. using a modified logistic transformation.

# define a custom transformation and its inverse for the y axis (modified logistic 
# transformation) using the scales package (for ggplot graphics).
modif_logistic_trans <- function() trans_new("modif_logistic", 
                                             function(x) log((1 + x)/(101 - x)), 
                                             function(x) exp((1 + x)/(101 - x)))

# plot the curves again
pdf(file = file.path(figs.dir, "abc-transformed-y_aver_sand.pdf"), useDingbats = FALSE)
plot_dom_curves_facets(abc.aver.sand, stations = stations.sand, trasf.y = TRUE)
dev.off()


# analyze the W statistic. The standardized sum of the difference curves 
# (resulting from the subtraction of the abundance Ai from the biomass Bi for 
# any given rank i) is called the W statistic.
# W is strongly positive when the biomass curve is above the abundance curve 
# throughout their lengths; strongly negative when the curves are transposed; 
# close to 0 in intermediate cases, when the curves are intertwined.

# get the w values from the ABC results (stored in the second element of each
# sample/replicate's sublist) 
w <- sapply(abc.sand, function(x) x[2])
w <- as.data.frame(t(as.data.frame(w)))
row.names(w) <- c(1:nrow(w))

# check if ANOVA assumptions are met

# perform ANOVA on W; factors - stations & years 



## Partial dominance curves. Used in order to avoid problems of over-dominance of
# species affecting the curves and their interpretation. Partial dominance curves
# compute the dominance of the second ranked species over the remainder, ignoring
# the first-ranked species, and so on with the remaining species. Thus, earlier 
# values do not affect later points on the curve. Partial dominance ABC curves 
# for undisturbed macrobenthic communities have the biomass curve over the abundance
# curve (mostly) throughout its length. The abundance curve is smoother; the biomass
# curve shows a slight and steady decline before its inevitable final rise.
# Under polluted conditions, there is still a change in position of partial dominance
# curves, and the abundance curve is much more "ragged": the implication is that the
# disturbance is affecting the whole suite of species in the community.

# calculate the partial dominance curves, using the transposed numeric data (same
# as the ABC curves above)
part.dominance.sand <- mapply(partial_dominance_curves, 
                              as.data.frame(t(num.zoo.abnd.sand)), 
                              as.data.frame(t(num.zoo.biomass.sand)), 
                              SIMPLIFY = FALSE)


# calculate the average curves/station, then plot
pd.aver.sand <- mapply(partial_dominance_curves, 
                       as.data.frame(t(summary.abnd[-1])), 
                       as.data.frame(t(summary.biomass[-1])),
                       SIMPLIFY = FALSE)

l1 <- lapply(pd.aver.sand, plot_dom_curves)
do.call(grid.arrange, l1)

pdf(file = file.path(figs.dir, "part-dom-curves_aver_sand.pdf"), useDingbats = FALSE)
plot_dom_curves_facets(pd.aver.sand, stations = stations.sand)
dev.off()
