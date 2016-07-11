### PCA on environmental parameters 
### (continued from environmental data script; references objects from it)

# source: http://www.sthda.com/english/wiki/principal-component-analysis-how-to-reveal-the-most-important-variables-in-your-data-r-software-and-data-mining#at_pco=smlwn-1.0&at_si=56ab2ae7c9bd86b4&at_ab=per-2&at_pos=0&at_tot=1
library(plyr)
library(vegan)
library(FactoMineR)
library(factoextra)

# if not installed:
# devtools::install_github("kassambara/factoextra")
# install.packages("FactoMineR", "factoextra")

################################################################################################################
### explore PCA results - basic procedure

# 1. Check the variances of the principal components. Amount of variation retained
# by each PC = eigenvalue. First PC = direction with maximum amount of variation 
# in the dataset.
pca.result$eig

# 2. Visualize the importance of the PCs (-> scree plot)
fviz_screeplot(pca.result)

# 3. Plot correlations/loadings of the variables with the PCs = variable loadings.
# Variables can be plotted as points in the component space using their loadings
# as coordinates.

# look at variable coordinates
pca.result$var$coord

# visualize the variables on the factor map. Correlation circle can help visualize
# the most correlated variables (variables that group together). 
fviz_pca_var(pca.result) + theme_bw()

# 4. Explore the quality of the representation for variables on the factor map 
# (cos2 = squared loadings for variables = cor * cor = coord * coord). The closer
# a variable to the circle of correlations, the better its representation on the 
# factor map, and the more important it is to interpret these components. If a 
# variable is perfectly represented by only 2 components, the sum of the cos2 = 1,
# and the variables will be positioned on the circle of correlations. For some 
# variables - more than 2 components required to perfectly represent the data; 
# then - variables positioned inside circle of correlations. Variables close to 
# the center of the plot - less important for the first components.

# rearrange the data frame successively by PC axes so that it's easier to read
# and find the variables most correlated with each axis (custom function)
pca_repres_quality(pca.res = pca.result, choice = "var", param = "cos2")

# plot the variables factor map, and color the vectors according to the amount of 
# correlation to the PCs. Can select the minimum value of cos2 and the PC axes 
# to plot.
# PC1-PC2
fviz_pca_var(pca.result, axes = c(1, 2), select.var = list(cos2 = 0.5)) + 
  scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.7) +
  theme_minimal()

# PC1-PC3    
fviz_pca_var(pca.result, axes = c(1, 3), select.var = list(cos2 = 0.5)) + 
  scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.7) +
  theme_minimal()

# 5. Check the contributions of the variables to the PCs. Variables that are correlated
# with PC1 and PC2 - the most important in explaining the variability in the dataset.
# Variables not correlated with any PC or correlated with the last dimensions - 
# low contribution; could be removed to simplify overall analysis.
# % contribution of a variable in accounting for the variability in a given PC = 
# (variable.cos2 * 100) / (total cos2 of component)
pca.result$var$contrib

# rearrange the data frame successively by PC axes so that it's easier to read
# and find the variables with the highest contribution to each axis
pca_repres_quality(pca.result, choice = "var", param = "contrib")

# 6. Visualize the most important variables associated with a given PC. The red 
# line on the graph represents the expected average contribution (if all 
# variable contributions were uniform). For a given component, any variable with
# contribution larger than that cutoff could be considered as important.

# variable contributions on PC1, PC2 and PC3 (+ all 3 axes together)
plot_pca_var_contrib(pca.result, 
                     choice = "var", 
                     axes = c(1, 3))


# if there are many variables in the dataset, you can show only the top n
# contributing variables:
fviz_contrib(pca.result, choice = "var", axes = 1, top = 10)

# 7. Color variables on the variable factor map according to their contributions - 
# highlights the most important variables in explaining the variations retained 
# by the PCs.
fviz_pca_var(pca.result, col.var = "contrib") + theme_bw()

# variable contributions to PC1-2
fviz_pca_var(pca.result, axes = c(1, 2), col.var = "contrib") + theme_bw()
# # PC1-3
fviz_pca_var(pca.result, axes = c(1, 3), col.var = "contrib") + theme_bw()


## 8. Dimension description - identify the most correlated variables with a given PC.
# description of PC 1 to 5 -> WILL ONLY WORK WITH MORE CATEGORIES (STATION REPLICATES, YEARS,..)
dimdesc(pca.result, axes = c(1:5))

## 9. Explore individuals (= stations/species/..)
pca.result$ind$coord

# plot & color by station (or year, or other group..)
fviz_pca_ind(pca.result,
             geom = "text",
             habillage = "station", 
             jitter = list(what = "label", width = 0.2, height = 0.3)) + 
  theme_bw() + 
  theme(legend.position = "none")

# quality of representation for individuals on the PCs (cos2)
pca.result$ind$cos2
pca_repres_quality(pca.result, choice = "ind", param = "cos2")

# PC1-2
fviz_pca_ind(pca.result, axes = c(1, 2), col.ind = "cos2") + 
  scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.5) +
  theme_minimal()

# PC1-3
fviz_pca_ind(pca.result, axes = c(1, 3), col.ind = "cos2") + 
  scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.5) +
  theme_minimal()

# contribution of individuals to the PCs
pca.result$ind$contrib
pca_repres_quality(pca.result, choice = "ind", param = "contrib")

# visualize the most contributing individuals (stations) associated with a given PC
# PC1-3 + all 3 together (custom function)
plot_pca_var_contrib(pca.result, 
                     choice = "ind", 
                     axes = c(1, 3))

# top 10 individuals (stations) contributing to PC1 (for ex.)
fviz_contrib(pca.result, choice = "ind", axes = 1, top = 10)

# individuals (stations) map colored according to their contribution
fviz_pca_ind(pca.result, col.ind = "contrib") + 
  scale_color_gradient2(low = "skyblue", mid = "steelblue", high = "dark blue", midpoint = 50) +
  theme_bw()

##############################################################################################################################


## 1. PCA only on 2013-2014 data (with previously imputed missing values) - 
## 18 unique values per station (3 samplings x 6 stations) - NO GOOD; DROPPED - 
## DATASET TOO VARIABLE W/ TOO FEW MEASURES TO ACCOUNT FOR IT 

## 2. PCA only on 2013-2014 data (with previously imputed missing values) - 
## averaged by STATION - 6 unique values - NO GOOD EITHER BECAUSE MIXED WATER & 
## SEDIMENT + WATER VERY VARIABLE BECAUSE TOO FEW VALUES MEASURED (ONLY SEVERAL 
## POINT SAMPLINGS FOR THESE YEARS) 

## 3. PCA only on 2013-2014 data (with previously imputed missing values) - 
## averaged by STATION & YEAR - 12 unique values - ALSO DROPPED - SEE PCA 2

## 4. PCA on long-term data - averaged by STATION + 2 SEPARATE DATSETS FOR WATER
## AND SEDIMENTS (MAYBE 3, IF HEAVY METALS ON THEIR OWN) - 6 unique values/dataset

# summarize water column parameters by station (excluding the month and year) 
# -> for use in PCA
water.sand.by.st <- ddply(water.sand.imp.df[, !names(water.sand.imp.df) %in% c("year", "month")], 
                          .(station), 
                          colwise(mean, .cols = is.numeric))

# summarize other environmental parameters by station (mean values) -> for use in 
# PCA
other.env.sand.by.st <- ddply(other.env.sand.imp.df[, !names(other.env.sand.imp.df) %in% c("year", "month")], 
                              .(station), 
                              colwise(mean, .cols = is.numeric))

# make one dataset of only water column parameters, and another of sediment 
# parameters - for 2 different PCAs.
# Move temperature, salinity and O2 (average and bottom) + distance to innermost 
# station to water column parameters
water.sand.pca <- join(water.sand.by.st,
                       other.env.sand.by.st[, names(other.env.sand.by.st) %in% c("station",
                                                                                 "salinity", 
                                                                                 "temperature",
                                                                                 "temp.bottom", 
                                                                                 "O2.bottom", 
                                                                                 "O2.average", 
                                                                                 "dist.innermost")], 
                       by = "station")

other.env.sand.pca <- other.env.sand.by.st[, !names(other.env.sand.by.st) %in% c("salinity", 
                                                                                 "temperature",
                                                                                 "temp.bottom", 
                                                                                 "O2.bottom", 
                                                                                 "O2.average")]

# import heavy metals data - for separate PCA 
heavy.metals.sand <- read.csv(file.path(data.dir, "heavy-metals-sand.csv"), header = TRUE) 

# reorder factor station according to sand station order
heavy.metals.sand$station <- reorder.factor(heavy.metals.sand$station, 
                                            new.order = stations.sand)


# remove the no longer necessary data frames and subsets
rm(other.env.sand.by.st, water.sand.by.st)


### perform variable pruning to reduce variable redundancy in the dataset
## prune water column parameters - NO PRUNING DONE HERE, BECAUSE NOT SURE THE 
## VARIABLE SELECTION PROCEDURE IS CORRECT (OUTCOME IS NOT BINARY!). KEPT THE 
## WRAPPER FUNCTION IN CASE, AND ONLY ELIMINATED HIGHLY INTERCORRELATED VARIABLES. 

# find and eliminate highly correlated variables in the datasets
filter_correlated_vars <- function(env.df, cor.threshold = 0.75) {
  ## wrapper for several functions finding highly correlated variables 
  ## in the dataset and eliminating them
  
  library(caret)
  
  # calculate the correlation matrix (on the numeric variables only)
  env.vars.cors <- cor(env.df[, sapply(env.df, is.numeric)])
  
  # find and eliminate the highly correlated variables
  highly.cor.vars <- findCorrelation(env.vars.cors, cutoff = cor.threshold, names = FALSE)
  env.vars.fin <- env.df[, -highly.cor.vars]
  
  return(env.vars.fin)
}

water.vars.fin <- filter_correlated_vars(water.sand.pca, cor.threshold = 0.75)

### standardize the data and perform the PCA
water.vars.std <- as.data.frame(scale(water.vars.fin[, sapply(water.vars.fin, is.numeric)], 
                                      center = TRUE, scale = TRUE))

water.vars.std <- cbind(station = water.sand.pca$station, water.vars.std)

pca.water.sand.pruned <- PCA(water.vars.std, scale.unit = FALSE, quali.sup = 1, graph = FALSE)
summary(pca.water.sand.pruned)

## save for reference
saveRDS(pca.water.sand.pruned, file.path(save.dir, "pca_water-column_sand_pruned.rds"))

## clean up 
rm(water.vars.fin, water.vars.std)

## plot & check
# PCA biplots - save as pdf!
pdf(file.path(figs.dir, "explor_pca_water-column_sand_pruned_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca.water.sand.pruned, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca.water.sand.pruned, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca.water.sand.pruned, file.name = "explor_pca_water-column_sand_pruned.pdf")



### same for other environmental parameters
# prune & filter
other.env.pruned <- filter_correlated_vars(other.env.sand.pca, cor.threshold = 0.75)

# standardize and perform PCA
other.env.pruned.std <- as.data.frame(scale(other.env.pruned[, sapply(other.env.pruned, is.numeric)], 
                                            center = TRUE, scale = TRUE))

other.env.pruned.std <- cbind(station = other.env.sand.pca$station, other.env.pruned.std)

pca.other.env.sand.pruned <- PCA(other.env.pruned.std, scale.unit = FALSE, 
                                 quali.sup = 1, graph = FALSE)
summary(pca.other.env.sand.pruned)

## save PCA result
saveRDS(pca.other.env.sand.pruned, file.path(save.dir, "pca_other-env_sand_pruned.rds"))

## clean up workspace
rm(other.env.pruned.std, other.env.pruned)

## plot & check
# PCA biplots
pdf(file.path(figs.dir, "explor_pca_sediment-params_sand_pruned_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca.other.env.sand.pruned, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca.other.env.sand.pruned, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca.other.env.sand.pruned, file.name = "explor_pca_sediment-params_sand_pruned.pdf")


### same for the heavy metals
names(heavy.metals.sand)
# drop the month and year column
heavy.metals.pca <- heavy.metals.sand[, setdiff(names(heavy.metals.sand), c("month", "year"))]

# prune and filter correlated variables
heavy.metals.pruned <- filter_correlated_vars(heavy.metals.pca, cor.threshold = 0.75)

# standardize and perform PCA
heavy.metals.pruned.std <- as.data.frame(scale(heavy.metals.pruned[, sapply(heavy.metals.pruned, is.numeric)], 
                                               center = TRUE, scale = TRUE))

heavy.metals.pruned.std <- cbind(station = heavy.metals.pca$station, heavy.metals.pruned.std)

pca.heavy.metals.sand.pruned <- PCA(heavy.metals.pruned.std, scale.unit = FALSE, 
                                    quali.sup = 1, graph = FALSE)
summary(pca.heavy.metals.sand.pruned)

# save
saveRDS(pca.heavy.metals.sand.pruned, file.path(save.dir, "pca_heavy-metals_sand_pruned.rds"))

# clean workspace
rm(heavy.metals.pruned, heavy.metals.pruned.std, heavy.metals.pca)


# plot PCA biplots - save as pdf!
pdf(file.path(figs.dir, "explor_pca_heavy-metals_sand_pruned_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca.heavy.metals.sand.pruned, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca.heavy.metals.sand.pruned, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()


dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca.heavy.metals.sand.pruned, 
                    file.name = "explor_pca_heavy-metals_sand_pruned.pdf")



######################################################################################################################
### perform PCAs on the full datasets 
## PCA on water column parameters
# standardize the data (only the numeric variables)
water.sand.pca.std <- as.data.frame(scale(water.sand.pca[, sapply(water.sand.pca, is.numeric)], 
                                   center = TRUE, scale = TRUE))


# add back the factors (useful for labeling later)
water.sand.pca.std <- cbind(station = water.sand.pca$station, water.sand.pca.std)

# run the PCA
pca.water.sand <- PCA(water.sand.pca.std, scale.unit = FALSE, quali.sup = 1, graph = FALSE)

# SAVE FOR FUTURE REFERENCE
saveRDS(pca.water.sand, file.path(save.dir, "pca_water-column_sand.rds"))

# clean up workspace
rm(water.sand.pca.std)

# PCA biplots - save as pdf!
pdf(file.path(figs.dir, "explor_pca_water-column_sand_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca.water.sand, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca.water.sand, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca.water.sand, file.name = "explor_pca_water-column_sand.pdf")


## PCA on sediment parameters
# standardize the data (only the numeric variables)
other.env.sand.pca.std <- as.data.frame(scale(other.env.sand.pca[, sapply(other.env.sand.pca, is.numeric)], 
                                          center = TRUE, scale = TRUE))

# add back the factors (useful for labeling later)
other.env.sand.pca.std <- cbind(station = other.env.sand.pca$station, other.env.sand.pca.std)

# run the PCA
pca.other.env.sand <- PCA(other.env.sand.pca.std, scale.unit = FALSE, quali.sup = 1, graph = FALSE)

# SAVE FOR FUTURE REFERENCE
saveRDS(pca.other.env.sand, file.path(save.dir, "pca_sediment-params_sand.rds"))

# clean up workspace
rm(other.env.sand.pca.std)

# PCA biplots - save as pdf!
pdf(file.path(figs.dir, "explor_pca_sediment-params_sand_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca.other.env.sand, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca.other.env.sand, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca.other.env.sand, file.name = "explor_pca_sediment-params_sand.pdf")


## PCA on heavy metals
# standardize the data (only the numeric variables)
heavy.metals.sand.pca.std <- as.data.frame(scale(heavy.metals.sand[, !names(heavy.metals.sand) %in% c("station", "month", "year")], 
                                           center = TRUE, scale = TRUE))

# add back the factors (useful for labeling later)
heavy.metals.sand.pca.std <- cbind(station = heavy.metals.sand$station, heavy.metals.sand.pca.std)

# run the PCA
pca.heavy.metals.sand <- PCA(heavy.metals.sand.pca.std, scale.unit = FALSE, quali.sup = 1, graph = FALSE)

# SAVE FOR FUTURE REFERENCE
saveRDS(pca.heavy.metals.sand, file.path(save.dir, "pca_heavy-metals_sand.rds"))

# clean up workspace
rm(heavy.metals.sand.pca.std)

# PCA biplots - save as pdf!
pdf(file.path(figs.dir, "explor_pca_heavy-metals_sand_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca.heavy.metals.sand, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca.heavy.metals.sand, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey50", 
                repel = TRUE) + 
  theme_bw()

dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca.heavy.metals.sand, file.name = "explor_pca_heavy-metals_sand.pdf")


