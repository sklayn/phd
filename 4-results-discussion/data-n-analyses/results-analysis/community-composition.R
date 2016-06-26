### Community composition analysis 

# define the working subdirectories
data.dir <- "data"
functions.dir <- "R"
save.dir <- "output"
figs.dir <- "figs"

# import libraries
library(gdata)
library(ggplot2)
library(gridExtra)
library(lattice)
library(latticeExtra)
library(plyr)
library(reshape2)
library(R.utils)
library(vegan)
library(viridis)

# source the files containing the necessary funcitons (function from package R.utils)-> MAKE ONE FILE W/ ALL CUSTOM FUNCTIONS??
sourceDirectory(path = file.path(functions.dir))


# import the datasets

# specify the stations, in the order they appear in the input dataset (separately, because vector 
# will be used by other functions later, so better have it typed out only once at the beginning).

# stations.2012 <- c("Konski1", "Konski2", "Ribka1", "Ribka2", "Gradina1", "Gradina2")
stations.zostera <- c("Poda", "Otmanli", "Gradina", "Ropotamo", "Vromos") 
stations.sand <- c("Kraimorie", "Chukalya", "Akin", "Sozopol", "Agalina", "Paraskeva") 


# TEST TO SEE IF WORKS WITH CURRENT IMPORT FUNCTION! (prob. not because of number at the end of station names)
# zoo.abnd.2012 <- import_zoo_data(data.dir = data.dir, 
#                             zoo.data = "zoo-abnd-2012.csv", 
#                             stations = stations.2012, 
#                             repl = 4)

# add a column for habitat type (2 adjacent habitats sampled in 2012 - seagrass and bare sand);
# place it in the beginning with the other factors
# zoo.abnd.2012 <- data.frame(append(zoo.abnd.2012, 
#                                    list(habitat = rep(c("zostera", "sand"), each = 24)), 
#                                    after = match("stations", names(zoo.abnd.2012))))


zoo.abnd.sand <- import_zoo_data(data.dir = data.dir, 
                                 zoo.data = "zoo-abnd-sand.csv", 
                                 station.names = stations.sand, 
                                 repl = 3)

zoo.abnd.zostera <- import_zoo_data(data.dir = data.dir, 
                                    zoo.data = "zoo-abnd-zostera.csv", 
                                    station.names = stations.zostera, 
                                    repl = 4)



# after-import sanity checks (better safe than sorry!)
str(zoo.abnd.sand)  # structure & variable classes: factors are factors, numeric variables
                    # are numeric, etc.

str(zoo.abnd.zostera)
# str(zoo.abnd.2012)

names(zoo.abnd.sand)  # variable (species) names
names(zoo.abnd.zostera)
# names(zoo.abnd.2012)


## SAVE THE CLEANED AND REARRANGED DATA
write.csv(zoo.abnd.sand, 
          file = file.path(save.dir, "zoo-sand-clean.csv"), 
          row.names = FALSE)

write.csv(zoo.abnd.zostera, 
          file = file.path(save.dir, "zoo-zostera-clean.csv"), 
          row.names = FALSE)

# make subsets of the numeric columns (the species data), and the factor columns -
# to avoid retyping the command, because these subsets will be 
# used a lot in all the subsequent analyses. 
num.zoo.abnd.sand <- zoo.abnd.sand[sapply(zoo.abnd.sand, is.numeric)]
factors.zoo.sand <- zoo.abnd.sand[sapply(zoo.abnd.sand, is.factor)]

num.zoo.abnd.zostera <- zoo.abnd.zostera[sapply(zoo.abnd.zostera, is.numeric)]
factors.zoo.zostera <- zoo.abnd.zostera[sapply(zoo.abnd.zostera, is.factor)]

# summary of numeric abundance data by station (mean) - for summary plots, etc.
summary.abnd.sand <- ddply(zoo.abnd.sand, .(stations), colwise(mean, .cols = is.numeric))

summary.abnd.zostera <- ddply(zoo.abnd.zostera, .(stations), colwise(mean, .cols = is.numeric))

## Basic taxonomic composition and structure
# import taxonomic data (species names are rows, successively higher taxonomic
# ranks are columns)
zoo.taxa <- read.csv(file.path(data.dir, "zoo-taxonomy.csv"), header = T, row.names = 1)

# get the present taxa only (abundance > 0 in all samples)
community.sand <- num.zoo.abnd.sand[colSums(num.zoo.abnd.sand) > 0] 
current.taxa.sand <- subset(zoo.taxa, row.names(zoo.taxa) %in% names(community.sand))

community.zostera <- num.zoo.abnd.zostera[colSums(num.zoo.abnd.zostera) > 0] 
current.taxa.zostera <- subset(zoo.taxa, row.names(zoo.taxa) %in% names(community.zostera))

# explore the taxonomic composition of the community: number of taxa per phylum/class, etc.
table(current.taxa.sand$class)
table(current.taxa.sand$phylum)

table(current.taxa.zostera$class)
table(current.taxa.zostera$phylum)

# plot comparison of nb taxa/class and nb taxa/phylum side by side, and save for reference
pdf(file = file.path(figs.dir, "explor_nb-taxa_sand.pdf"), useDingbats = FALSE)
p.class <- barchart(sort(table(current.taxa.sand$class)), # sort for easier comparison
                                         main = "Per class", 
                                         xlab = "Number of taxa", 
                                         col = "skyblue")

p.phyl <- barchart(sort(table(current.taxa.sand$phylum)), # sort for easier comparison
                   main = "Per phylum", 
                   xlab = "Number of taxa", 
                   col = "skyblue")

# have to use grid.arrange() from gridExtra to put both plots on the same graphics device, 
# because barchart() is a lattice function.
grid.arrange(p.class, p.phyl, nrow = 2)

dev.off()
rm(p.class, p.phyl)

# idem for zostera communities
pdf(file = file.path(figs.dir, "explor_nb-taxa_zostera.pdf"), useDingbats = FALSE)
p.class <- barchart(sort(table(current.taxa.zostera$class)), # sort for easier comparison
                    main = "Per class", 
                    xlab = "Number of taxa", 
                    col = "lightgreen")

p.phyl <- barchart(sort(table(current.taxa.zostera$phylum)), # sort for easier comparison
                   main = "Per phylum", 
                   xlab = "Number of taxa", 
                   col = "lightgreen")

# have to use grid.arrange() from gridExtra to put both plots on the same graphics device, 
# because barchart() is a lattice function.
grid.arrange(p.class, p.phyl, nrow = 2)

dev.off()
rm(p.class, p.phyl)


# add new column with the most commonly used larger taxonomic groups from the literature
current.taxa.sand$group <- with(current.taxa.sand, 
                                ifelse(class == "Polychaeta", "Polychaeta", 
                                ifelse(class == "Bivalvia" | 
                                       class == "Gastropoda" | 
                                       class == "Polyplacophora", "Mollusca", 
                                ifelse(class == "Malacostraca", "Crustacea", "Varia"))))

table(current.taxa.sand$group)


current.taxa.zostera$group <- with(current.taxa.zostera, 
                                   ifelse(class == "Polychaeta", "Polychaeta", 
                                   ifelse(class == "Bivalvia" | 
                                          class == "Gastropoda" | 
                                          class == "Polyplacophora", "Mollusca", 
                                   ifelse(class == "Malacostraca", "Crustacea", "Varia"))))

table(current.taxa.zostera$group)

# plot the number of taxa in each of these taxonomic groups
pdf(file = file.path(figs.dir, "nb-taxa_sand.pdf"), useDingbats = FALSE)
barchart(sort(table(current.taxa.sand$group)), 
         main = "Number of taxa",
         xlab = "", 
         col = "skyblue")
dev.off()


pdf(file = file.path(figs.dir, "nb-taxa_zostera.pdf"), useDingbats = FALSE)
barchart(sort(table(current.taxa.zostera$group)), 
         main = "Number of taxa",
         xlab = "", 
         col = "lightgreen")
dev.off()


## calculate the contribution of each taxonomic group to the numeric community 
## composition by station and by year
# first calculate proportions by taxonomic group in each station/replicate
tax.group.props.sand <- tax_group_contribution(community.sand, current.taxa.sand$group)
tax.group.props.sand <- cbind(factors.zoo.sand, tax.group.props.sand)

tax.group.props.zostera <- tax_group_contribution(community.zostera, current.taxa.zostera$group)
tax.group.props.zostera <- cbind(factors.zoo.zostera, tax.group.props.zostera)

# save to file (in case) 
write.csv(tax.group.props.sand, 
          file = file.path(save.dir, "tax-gr-proportions_sand-clean.csv"), 
          row.names = FALSE)

write.csv(tax.group.props.zostera, 
          file = file.path(save.dir, "tax-gr-proportions_zostera-clean.csv"), 
          row.names = FALSE)



# plot the contribution of taxonomic groups: 
# by station
pdf(file = file.path(figs.dir, "tax-gr-contrib_stations_sand.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.group.props.sand)
dev.off()

pdf(file = file.path(figs.dir, "tax-gr-contrib_stations_zostera.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.group.props.zostera)
dev.off()

# by station and year
pdf(file = file.path(figs.dir, "tax-gr-contrib_st-yrs_sand.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.group.props.sand, by.years = TRUE)
dev.off()

pdf(file = file.path(figs.dir, "tax-gr-contrib_st-yrs_zostera.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.group.props.zostera, by.years = TRUE)
dev.off() 


## Diversity indices

## Alpha diversity (Whittaker, 1960) - description of the diversity in one spot
diversity.sand <- alpha_diversity(zoo.abnd.sand)

# save the data frame of diversity indices and measures in a file
write.csv(diversity.sand, 
          file = file.path(save.dir, "diversity.sand.csv"), 
          row.names = FALSE)


# calculate and plot diversity profiles - graphical representation of the shape of the community;
# show how the perceived diversity changes as the emphasis shifts from common to rare
# species - we can judge their respective contributions in the community composition.
diversity.profiles.sand <- diversity_profiles(num.zoo.abnd.sand, q = 50)

# plot the diversity profiles (all samples - in panels by station); save to file.
# !NB set useDingbats = FALSE, otherwise points might get transformed to letters 
# when the pdf file is opened in another program!
pdf(file = file.path(figs.dir, "classical-diversity-profiles_all_sand.pdf"), useDingbats = FALSE)
plot_div_profiles(diversity.profiles.sand, stations.sand)
dev.off()

# save diversity profiles to file 
write.csv(diversity.profiles.sand, 
          file.path(save.dir, "div-profiles-sand_classical.csv"))

# add a measure of similarity between species/taxa. 
# If not already in the workspace, import taxonomic data from which distances will
# be derived. 
# zoo.taxa <- read.csv(file.path(data.dir, "zoo-taxonomy.csv"), header = T, row.names = 1)

# calculate diversity profiles including a measure of similarity between species -
# allows for a lot more meaningful ecological and biodiversity comparisons; then
# calculate and plot diversity profiles again (Leinster & Cobbold, 2012)
weighted.profiles.sand <- weighted_div_profiles(num.zoo.abnd.sand, zoo.taxa)

# plot the profiles in panels by station and save
pdf(file = file.path(figs.dir, "weighted-diversity-profiles_all_sand.pdf"), useDingbats = FALSE)
plot_div_profiles(weighted.profiles.sand, stations.sand, one.panel = FALSE)
dev.off()

# same procedure, but using the data averaged by station - first column is the station names!
aver.profiles.sand <- weighted_div_profiles(summary.abnd[-1], zoo.taxa)

# plot the profiles; save as pdf 
pdf(file = file.path(figs.dir, "weighted-diversity-profiles_aver_sand.pdf"), useDingbats = FALSE)
plot_div_profiles(aver.profiles.sand, stations.sand, one.panel = TRUE)
dev.off()

# plot all the weighted diversity profiles by station, and add the average profiles
# on the same graph
pdf(file = file.path(figs.dir, "weighted-div-profiles_allaver_sand.pdf"), useDingbats = FALSE)
plot_div_profiles_w_aver(weighted.profiles.sand, 
                         aver.profiles.sand, 
                         stations.sand, 
                         col.profiles = c("skyblue", "royalblue"))

dev.off()


## Taxonomic diversity

# NEEDS UPDATED MASTER SPECIES LIST FOR BIOGEOGRAPHIC AREA

# re-import abundances, but including the full species list for the area (Black Sea 
# macroinvertebrate fauna) 
tax.zoo.abnd <- import_zoo_data(data.dir = data.dir, 
                                 zoo.data = "td-zoo-abnd-sand.csv", 
                                 stations = stations.sand, 
                                 repl = 3)

# import taxonomy data for the full species list (genus, family, order..)
taxa.list <- read.csv(file = file.path(data.dir, "td-zoo-taxonomy.csv"), 
                      header = TRUE, 
                      row.names = 1)

# calcualate teh taxonomic distance for the species list, using variable step 
taxa.dist <- taxa2dist(taxa.list, varstep = TRUE)

# calculate the taxonomic distinctness indices - on the numeric columns (species)
tax.distinctness <- taxondive(tax.zoo.abnd[sapply(tax.zoo.abnd, is.numeric)], 
                              taxa.dist)

# convert the taxonomic distinctness indices list to a data frame
tax.dist.table <- as.data.frame(lapply(tax.distinctness, 
                                       function(x) {vals <- x; return(vals)}))  

# add identifying factors to the data frame - for easier plotting and subsetting
tax.dist.table <- cbind(factors.zoo.sand, tax.dist.table)

# save tax.distinctness indices to file 
write.csv(tax.dist.table, 
          file = file.path(save.dir, "tax.distinctness.sand.csv"), 
          row.names = FALSE)

# aggregate taxonomic distinctness indices by station and year
tax.dist.table$stations <- reorder.factor(tax.dist.table$stations, 
                                            new.order = stations.sand)
tax.dist.summary <- ddply(tax.dist.table, .(stations), 
                          colwise(mean, .cols = is.numeric))


# plot selected taxonomic diversity indices:
# average tax.diversity based on presence/absence - Delta + 
# (reflects mean tax.breadth of the species lists)
pdf(file = file.path(figs.dir, "Dplus_sand.pdf"), useDingbats = FALSE)
ggplot(tax.dist.summary, mapping = aes(x = Species, y = Dplus)) +
  geom_point(size = 2, show.legend = FALSE) +
  geom_text(aes(label = stations), nudge_x = 0.1, nudge_y = -0.2) + 
  geom_hline(aes(yintercept = EDplus), colour = "gray") +
  geom_text(aes(40, EDplus, label = "EDplus"), nudge_y = 0.2, colour = "gray") + 
  theme_bw() 
dev.off() 

# variation in tax.distinctness Lambda + (also based on presence/
# absence, reflects unevenness in the tax.hierarchy)
pdf(file = file.path(figs.dir, "Lambda_sand.pdf"), useDingbats = FALSE)
ggplot(tax.dist.summary, mapping = aes(x = Species, y = Lambda)) +
  geom_point(size = 2, show.legend = FALSE) +
  geom_text(aes(label = stations), nudge_y = -4) +
  theme_bw()
dev.off()

# plot the (quantitative) average tax.distinctness Delta * 
# NEEDS STH MEANINGFUL AS X - DISTANCE FROM BURGAS, LUSI, other index/
# environmental variable - e.g. those identified as most significant in the 
# analyses of the environmental data...

# example - with distance from the innermost (= most impacted) station (log scale). 
# NB dist.inn won't be in workspace! 
pdf(file = file.path(figs.dir, "explor_AvTD-dist.innermost_sand.pdf"), useDingbats = FALSE)
plot(tax.dist.table$Dstar ~ dist.inn, 
     type = "n", 
     xlim = c(0.5, 100), 
     log = "x", 
     xlab = "log(dist.innermost.st)", ylab = "AvTD*")
text(tax.dist.table$Dstar, labels = tax.dist.table$stations)
dev.off()



## Species abundance models 
## Species accumulation and beta diversity

# species accumulation curve
sp.accumulation <- specaccum(num.zoo.abnd.sand)

# plot sac with 95% confidence intervals
pdf(file = file.path(figs.dir, "sp-accumulaton_sand.pdf"), useDingbats = FALSE)
plot(sp.accumulation, ylab = "Cumulative species richness", ci.type = "polygon", ci.col = "skyblue")
dev.off()

# beta diversity - from pairwise comparison of sites (= Sorensen index of 
# dissimilarity)
beta.div <- vegdist(num.zoo.abnd.sand, binary = TRUE)
beta.div


## Species pool
# number of unseen species - various estimates
specpool(x = num.zoo.abnd.sand)



## Graphical community structure analyses - ABC curves, k-dominance & partial 
## dominance curves

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

# import the biomass data
zoo.biomass.sand <- import_zoo_data(data.dir, "zoo-biomass-sand.csv", stations.sand, 3)

# quick check for import errors, etc.
str(zoo.biomass.sand)
names(zoo.biomass.sand)

# make a subset of only the numeric biomass data, and another of the mean biomass
# per station  
num.zoo.biomass.sand <- zoo.biomass.sand[sapply(zoo.biomass.sand, is.numeric)]
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



# Partial dominance curves. Used in order to avoid problems of over-dominance of
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



## nMDS on the abundance data
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
# first plot the NMDS ordination with sites
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