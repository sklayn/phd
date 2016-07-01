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


prepare_data_pca <- function(env.data, summarize.by) {
  ## helper for PCA. Prepares a data frame of environmental data for performing a
  ## PCA: summarizes according to the desired factors; removes the factors; 
  ## standardizes the data.
  
  library(plyr)
  
  # summarize (average) the data by the levels of the specified factor 
  env.pca <- ddply(env.data, 
                   .variables = summarize.by, 
                   colwise(mean, .cols = is.numeric))
  
  # standardize the data (only the variables)
  env.pca.std <- as.data.frame(scale(env.pca[, !names(env.pca) %in% c("station", "year", "month")], 
                                     center = TRUE, scale = TRUE))
  
  # add back the factors (useful for labeling later)
  env.pca.std <- cbind(env.pca[, names(env.pca) %in% c("station", "year", "month")], 
                       env.pca.std)
  
  return(env.pca.std)
}


## 1. only on 2013-2014 data (with previously imputed missing values) - 
## 18 unique values per station (3 samplings x 6 stations)
env.sand.pca1 <- prepare_data_pca(env.imp.all.sand, summarize.by = c("station", "year", "month"))
pca1.sand <- PCA(env.sand.pca1, quali.sup = c(1:3), graph = FALSE)

# SAVE FOR FUTURE REFERENCE
saveRDS(pca1.sand, file.path(save.dir, "pca1_sand.rds"))

# PCA biplot - save as pdf!
pdf(file.path(figs.dir, "explor_pca1_sand_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca1.sand, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey20", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca1.sand, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey20", 
                repel = TRUE) + 
  theme_bw()

dev.off()


# explore PCA results
# Variances of the principal components. 
# Amount of variation retained by each PC = eigenvalue. First PC = direction with
# maximum amount of variation in the dataset.
pca1.sand$eig

## For pca1.sand - pretty bad representation (many PCs needed to explain the 
## observed variablility...)

# visualize the importance of the PCs (scree plot)
fviz_screeplot(pca1.sand)


# Plot correlations/loadings of the variables with the PCs = variable loadings.
# Variables can be plotted as points in the component space using their loadings
# as coordinates.

# look at variable coordinates
pca1.sand$var$coord

# visualize the variables on the factor map. Correlation circle can help visualize
# the most correlated variables (variables that group together). 
fviz_pca_var(pca1.sand) + theme_bw()

# Explore the quality of the representation for variables on the factor map 
# (cos2 = squared loadings for variables = cor * cor = coord * coord). The closer
# a variable to the circle of correlations, the better its representation on the 
# factor map, and the more important it is to interpret these components. If a 
# variable is perfectly represented by only 2 components, the sum of the cos2 = 1,
# and the variables will be positioned on the circle of correlations. For some 
# variables - more than 2 components required to perfectly represent the data; 
# then - variables positioned inside circle of correlations. Variables close to 
# the center of the plot - less important for the first components.

pca_repres_quality <- function(pca.res, choice = c("var", "ind"), param = c("cos2", "contrib")) {
  ## Extracts data frame of PCA parameters (cos2 and contribution) and their 
  ## relation to the PCs, then arranges it successively for each PC in 
  ## descending order (of contribution, correlation, etc.) for easier reading.
  
  library(plyr)
  
  # select the variables
  if(choice == "var") {
    # get their correlation with the axes
    vars.cos2 <- as.data.frame(pca.res$var$cos2)
    vars.contrib <- as.data.frame(pca.res$var$contrib)
  
  # select the individuals
  } else {
    vars.cos2 <- as.data.frame(pca.res$ind$cos2)
    vars.contrib <- as.data.frame(pca.res$ind$contrib)
  }
  
  # sort each dimension in descending order
  vars.cos2 <- cbind(x = row.names(vars.cos2), vars.cos2)
  vars.contrib <- cbind(x = row.names(vars.contrib), vars.contrib)
  
  # return only the selected (ordered) parameter 
  if (param == "cos2"){
    apply(vars.cos2[-1], 2, function(m) arrange(vars.cos2, desc(m)))
  } else {
    apply(vars.contrib[-1], 2, function(m) arrange(vars.contrib, desc(m)))
  }
  
}


# rearrange the data frame successively by PC axes so that it's easier to read
# and find the variables most correlated with each axis
pca_repres_quality(pca.res = pca1.sand, choice = "var", param = "cos2")

## For pca1.sand: PC1 - heavy metals; PC2 - Secchi depth, N compounds, chl-a; 
## PC3 - sorting, organic matter, mean grain size. 


# plot the variables factor map, and color the vectors according to the amount of 
# correlation to the PCs. Can select the minimum value of cos2 and the PC axes 
# to plot.
# PC1-PC2
fviz_pca_var(pca1.sand, axes = c(1, 2), select.var = list(cos2 = 0.5)) + 
  scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.7) +
  theme_minimal()

# PC1-PC3    
fviz_pca_var(pca1.sand, axes = c(1, 3), select.var = list(cos2 = 0.5)) + 
  scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.7) +
  theme_minimal()


# Contributions of the variables to the PCs. Variables that are correlated with
# PC1 and PC2 - the most important in explaining the variability in the dataset.
# Variables not correlated with any PC or correlated with the last dimensions - 
# low contribution; could be removed to simplify overall analysis.
# % contribution of a variable in accounting for the variability in a given PC = 
# (variable.cos2 * 100) / (total cos2 of component)
pca1.sand$var$contrib

# rearrange the data frame successively by PC axes so that it's easier to read
# and find the variables with the highest contribution to each axis
pca_repres_quality(pca.res = pca1.sand, choice = "var", param = "contrib")

## For pca1.sand: PC1 - heavy metals; PC2 - Secchi depth; PC3 - sediment 
## parameters. 

# visualize the most important variables associated with a given PC. The red 
# line on the graph represents the expected average contribution (if all 
# variable contributions were uniform). For a given component, any variable with
# contribution larger than that cutoff could be considered as important.

# variable contributions on PC1, PC2 and PC3 (+ all 3 axes together)
plot_pca_var_contrib <- function(pca.res, choice = c("var", "ind"), axes = c(start : end)) {
  ## helper plotting variable/individuals contributions (as bars) to the specified
  ## axes. 

  # plot contributions for each axis separately
  for(i in axes){
    print(fviz_contrib(pca.res, choice = choice, axes = i))
  }
  
  # plot the overall contributions to the specified axes
  print(fviz_contrib(pca.res, choice = choice, axes = axes))
}


plot_pca_var_contrib(pca1.sand, 
                     choice = "var", 
                     axes = c(1, 3))


# if there are many variables in the dataset, you can show only the top n
# contributing variables:
# fviz_contrib(pca1.sand, choice = "var", axes = 1, top = 10)

# color variables on the variable factor map according to their contributions - 
# highlights the most important variables in explaining the variations retained 
# by the PCs.
fviz_pca_var(pca1.sand, col.var = "contrib") + theme_bw()

# variable contributions to PC1-2
fviz_pca_var(pca1.sand, axes = c(1, 2), col.var = "contrib") + theme_bw()
# PC1-3
fviz_pca_var(pca1.sand, axes = c(1, 3), col.var = "contrib") + theme_bw()


## Dimension description - identify the most correlated variables with a given PC.
# description of PC 1 to 5
dimdesc(pca1.sand, axes = c(1:5))


## EXplore individuals (= stations/species/..)
pca1.sand$ind$coord

# plot & color by station (or year, or other group..)
fviz_pca_ind(pca1.sand,
             geom = "text",
             habillage = "station", 
             jitter = list(what = "label", width = 0.2, height = 0.3)) + 
  theme_bw() + 
  theme(legend.position = "none")

# quality of representation for individuals on the PCs (cos2)
pca1.sand$ind$cos2
pca_repres_quality(pca1.sand, choice = "ind", param = "cos2")

# PC1-2
fviz_pca_ind(pca1.sand, axes = c(1, 2), col.ind = "cos2") + 
  scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.5) +
  theme_minimal()

# PC1-3
fviz_pca_ind(pca1.sand, axes = c(1, 3), col.ind = "cos2") + 
  scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.5) +
  theme_minimal()

# contribution of individuals to the PCs
pca1.sand$ind$contrib
pca_repres_quality(pca1.sand, choice = "ind", param = "contrib")

# visualize the most contributing individuals (stations) associated with a given PC
# PC1-3 + all 3 together
plot_pca_var_contrib(pca1.sand, 
                     choice = "ind", 
                     axes = c(1, 3))

# top 10 individuals (stations) contributing to PC1 (for ex.)
# fviz_contrib(pca1.sand, choice = "ind", axes = 1, top = 10)

# individuals (stations) map colored according to their contribution
fviz_pca_ind(pca1.sand, col.ind = "contrib") + 
  scale_color_gradient2(low = "skyblue", mid = "steelblue", high = "dark blue", midpoint = 50) +
  theme_bw()



## put all diagnostic plots for the PCA in one file 
plot_pca_diagnostic <- function(pca.res, file.name) {
  ## quick PCA diagnostic plots in one pdf file. Crap with many hard-coded 
  ## choices and values + calls other custom functions. Only acceptable to avoid 
  ## retyping, has to be modified if at all useful for later work. 
  
  pdf(file.path(figs.dir, file.name), 
      paper = "a4r", width = 12, height = 12,
      useDingbats = FALSE)
  
  # scree plot
  print(fviz_screeplot(pca.res) + theme_bw())
  
  ## VARIABLES
  # variables factor map, with vectors colored according to the amount of 
  # correlation to the PCs
  # PC1-PC2
  print(fviz_pca_var(pca.res, axes = c(1, 2), select.var = list(cos2 = 0.5), col.var = "cos2") + 
          scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.7) +
          theme_minimal()
  )
  
  # PC1-PC3    
  print(fviz_pca_var(pca.res, axes = c(1, 3), select.var = list(cos2 = 0.5), col.var = "cos2") + 
          scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.7) +
          theme_minimal()
  )
  # variable contributions to the PC axes (1-3 + all 3 together)
  plot_pca_var_contrib(pca.res, 
                       choice = "var", 
                       axes = c(1:3))
  
  # variable map with contribution to PCs
  # variable contributions to PC1-2
  print(fviz_pca_var(pca.res, axes = c(1, 2), 
                     select.var = list(contrib = 50), 
                     col.var = "contrib") 
        + theme_bw())
  # PC1-3
  print(fviz_pca_var(pca.res, axes = c(1, 3), 
                     select.var = list(contrib = 50), 
                     col.var = "contrib") 
        + theme_bw())
  
  
  ## INDIVIDUALS
  # individuals map - PC1-2
  print(fviz_pca_ind(pca.res,
                     axes = c(1, 2),
                     geom = "text",
                     habillage = "station", 
                     jitter = list(what = "label", width = 0.2, height = 0.3)) + 
          theme_bw() + 
          theme(legend.position = "none")
  )
  # individuals map - PC1-3
  print(fviz_pca_ind(pca.res,
                     axes = c(1, 3),
                     geom = "text",
                     habillage = "station", 
                     jitter = list(what = "label", width = 0.2, height = 0.3)) + 
          theme_bw() + 
          theme(legend.position = "none")
  )
  
  # individuals map, colored by representation quality (cos2)
  # PC1-2
  print(fviz_pca_ind(pca.res, axes = c(1, 2), col.ind = "cos2") + 
          scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.5) +
          theme_minimal()
  )
  # PC1-3
  print(fviz_pca_ind(pca.res, axes = c(1, 3), col.ind = "cos2") + 
          scale_color_gradient2(low = "skyblue", mid = "navyblue", high = "red", midpoint = 0.5) +
          theme_minimal()
  )
  
  # contribution of individuals to the PCs
  # PC1-3 + all 3 together
  plot_pca_var_contrib(pca.res, 
                       choice = "ind", 
                       axes = c(1:3))
  
  
  dev.off()
}


plot_pca_diagnostic(pca.res = pca1.sand, file.name = "explor_pca1_sand.pdf")


## 2. only on 2013-2014 data (with previously imputed missing values) - 
## averaged by STATION - 6 unique values 
env.sand.pca2 <- prepare_data_pca(env.imp.all.sand, summarize.by = "station")
pca2.sand <- PCA(env.sand.pca2, quali.sup = c(1:3), graph = FALSE)

# SAVE FOR FUTURE REFERENCE
saveRDS(pca2.sand, file.path(save.dir, "pca2_sand.rds"))

# PCA biplots - save as pdf!
pdf(file.path(figs.dir, "explor_pca2_sand_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca2.sand, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey20", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca2.sand, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey20", 
                repel = TRUE) + 
  theme_bw()

dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca2.sand, file.name = "explor_pca2_sand.pdf")



## 3. only on 2013-2014 data (with previously imputed missing values) - 
## averaged by STATION & YEAR - 12 unique values
env.sand.pca3 <- prepare_data_pca(env.imp.all.sand, summarize.by = c("station", "year"))
pca3.sand <- PCA(env.sand.pca3, quali.sup = c(1:3), graph = FALSE)

# SAVE FOR FUTURE REFERENCE
saveRDS(pca3.sand, file.path(save.dir, "pca3_sand.rds"))

# PCA biplots - save as pdf!
pdf(file.path(figs.dir, "explor_pca3_sand_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca3.sand, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey20", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca3.sand, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey20", 
                repel = TRUE) + 
  theme_bw()

dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca3.sand, file.name = "explor_pca3_sand.pdf")


## 4. on long-term data - averaged by STATION - 6 unique values

# remove duplicates (there is technically only 1 valid row per station; all the
# others are obtained by copying it - so that the dataset can be used with the 
# replicates ,etc.)
env.LT.sand.pca <- env.imp.all.LT.sand[!duplicated(env.imp.all.LT.sand), ] 
row.names(env.LT.sand.pca) <- 1:nrow(env.LT.sand.pca)

# standardize the data (only the variables)
env.LT.sand.pca.std <- as.data.frame(scale(env.LT.sand.pca[, !names(env.LT.sand.pca) == "station"], 
                                   center = TRUE, scale = TRUE))

# add back the factors (useful for labeling later)
env.LT.sand.pca.std <- cbind(station = env.LT.sand.pca$station, env.LT.sand.pca.std)

# run the PCA
pca.LT.sand <- PCA(env.LT.sand.pca.std, quali.sup = 1, graph = FALSE)


# SAVE FOR FUTURE REFERENCE
saveRDS(pca.LT.sand, file.path(save.dir, "pca_LT_sand.rds"))

# PCA biplots - save as pdf!
pdf(file.path(figs.dir, "explor_pca_LT_sand_biplot.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)

# PC1-2
fviz_pca_biplot(pca.LT.sand, 
                axes = c(1, 2),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey20", 
                repel = TRUE) + 
  theme_bw()

# PC1-3
fviz_pca_biplot(pca.LT.sand, 
                axes = c(1, 3),
                label = "var", 
                habillage = "station", 
                select.var = list(cos2 = 0.5), 
                col.var = "grey20", 
                repel = TRUE) + 
  theme_bw()

dev.off()

# various PCA diagnostic & exploratory plots
plot_pca_diagnostic(pca.LT.sand, file.name = "explor_pca_LT_sand.pdf")
