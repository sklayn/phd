### Community composition & structure analysis 

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

stations.sand <- c("Kraimorie", "Chukalya", "Akin", "Sozopol", "Agalina", "Paraskeva") 
# stations.zostera <- c("Poda", "Otmanli", "Gradina", "Ropotamo", "Vromos") 
# stations.2012 <- c("Konski1", "Konski2", "Ribka1", "Ribka2", "Gradina1", "Gradina2")

zoo.abnd.sand <- import_zoo_data(data.dir = data.dir, 
                                 zoo.data = "zoo-abnd-sand.csv", 
                                 station.names = stations.sand, 
                                 repl = 3)

# zoo.abnd.zostera <- import_zoo_data(data.dir = data.dir, 
#                                     zoo.data = "zoo-abnd-zostera.csv", 
#                                     station.names = stations.zostera, 
#                                     repl = 4)

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


# after-import sanity checks (better safe than sorry!)
str(zoo.abnd.sand)  # structure & variable classes: factors are factors, numeric variables
                    # are numeric, etc.

# str(zoo.abnd.zostera)
# str(zoo.abnd.2012)

names(zoo.abnd.sand)  # variable (species) names
# names(zoo.abnd.zostera)
# names(zoo.abnd.2012)


## SAVE THE CLEANED AND REARRANGED DATA
write.csv(zoo.abnd.sand, 
          file = file.path(save.dir, "zoo-sand-clean.csv"), 
          row.names = FALSE)

# write.csv(zoo.abnd.zostera, 
#           file = file.path(save.dir, "zoo-zostera-clean.csv"), 
#           row.names = FALSE)

# make subsets of the numeric columns (the species data), and the factor columns -
# to avoid retyping the command, because these subsets will be 
# used a lot in all the subsequent analyses. 
factors.zoo.sand <- zoo.abnd.sand[sapply(zoo.abnd.sand, is.factor)]
num.zoo.abnd.sand <- zoo.abnd.sand[sapply(zoo.abnd.sand, is.numeric)]

# leave only the species actually present in the current dataset (to simplify later
# analyses)
num.zoo.abnd.sand <- num.zoo.abnd.sand[colSums(num.zoo.abnd.sand) > 0] 


# factors.zoo.zostera <- zoo.abnd.zostera[sapply(zoo.abnd.zostera, is.factor)]
# num.zoo.abnd.zostera <- zoo.abnd.zostera[sapply(zoo.abnd.zostera, is.numeric)]
# num.zoo.abnd.zostera <- num.zoo.abnd.zostera[colSums(num.zoo.abnd.zostera) > 0] 

# summary of numeric abundance data by station (mean) - for summary plots, etc.
summary.abnd.sand <- ddply(zoo.abnd.sand, .(stations), colwise(mean, .cols = is.numeric))
# summary.abnd.zostera <- ddply(zoo.abnd.zostera, .(stations), colwise(mean, .cols = is.numeric))

## Basic taxonomic composition and structure
# import taxonomic data (species names are rows, successively higher taxonomic
# ranks are columns)
zoo.taxa <- read.csv(file.path(data.dir, "zoo-taxonomy.csv"), header = T, row.names = 1)

# get the present taxa only (abundance > 0 in all samples)
current.taxa.sand <- subset(zoo.taxa, row.names(zoo.taxa) %in% names(num.zoo.abnd.sand))

# current.taxa.zostera <- subset(zoo.taxa, row.names(zoo.taxa) %in% names(num.zoo.abnd.zostera))

# explore the taxonomic composition of the community: number of taxa per phylum/class, etc.
table(current.taxa.sand$class)
table(current.taxa.sand$phylum)

# table(current.taxa.zostera$class)
# table(current.taxa.zostera$phylum)

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
# pdf(file = file.path(figs.dir, "explor_nb-taxa_zostera.pdf"), useDingbats = FALSE)
# p.class <- barchart(sort(table(current.taxa.zostera$class)), # sort for easier comparison
#                     main = "Per class", 
#                     xlab = "Number of taxa", 
#                     col = "lightgreen")
# 
# p.phyl <- barchart(sort(table(current.taxa.zostera$phylum)), # sort for easier comparison
#                    main = "Per phylum", 
#                    xlab = "Number of taxa", 
#                    col = "lightgreen")
# 
# grid.arrange(p.class, p.phyl, nrow = 2)
# 
# dev.off()
# rm(p.class, p.phyl)


# add new column with the most commonly used larger taxonomic groups from the literature
current.taxa.sand$group <- with(current.taxa.sand, 
                                ifelse(class == "Polychaeta", "Polychaeta", 
                                ifelse(class == "Bivalvia" | 
                                       class == "Gastropoda" | 
                                       class == "Polyplacophora", "Mollusca", 
                                ifelse(class == "Malacostraca", "Crustacea", "Varia"))))

table(current.taxa.sand$group)


# current.taxa.zostera$group <- with(current.taxa.zostera, 
#                                    ifelse(class == "Polychaeta", "Polychaeta", 
#                                    ifelse(class == "Bivalvia" | 
#                                           class == "Gastropoda" | 
#                                           class == "Polyplacophora", "Mollusca", 
#                                    ifelse(class == "Malacostraca", "Crustacea", "Varia"))))
# 
# table(current.taxa.zostera$group)

# plot the number of taxa in each of these taxonomic groups
pdf(file = file.path(figs.dir, "nb-taxa_sand.pdf"), useDingbats = FALSE)
barchart(sort(table(current.taxa.sand$group)), 
         main = "Number of taxa",
         xlab = "", 
         col = "skyblue")
dev.off()


# pdf(file = file.path(figs.dir, "nb-taxa_zostera.pdf"), useDingbats = FALSE)
# barchart(sort(table(current.taxa.zostera$group)), 
#          main = "Number of taxa",
#          xlab = "", 
#          col = "lightgreen")
# dev.off()


## calculate the contribution of each taxonomic group to the numeric community 
## composition by station and by year
# first calculate proportions by taxonomic group in each station/replicate
tax.group.props.sand <- tax_group_contribution(num.zoo.abnd.sand, current.taxa.sand$group)
# tax.group.props.sand <- cbind(factors.zoo.sand, tax.group.props.sand)

tax.group.props.zostera <- tax_group_contribution(num.zoo.abnd.zostera, current.taxa.zostera$group)
# tax.group.props.zostera <- cbind(factors.zoo.zostera, tax.group.props.zostera)

# save to file (in case) 
write.csv(tax.group.props.sand, 
          file = file.path(save.dir, "tax-gr-proportions_sand-clean.csv"), 
          row.names = FALSE)

# write.csv(tax.group.props.zostera, 
#           file = file.path(save.dir, "tax-gr-proportions_zostera-clean.csv"), 
#           row.names = FALSE)



# plot the contribution of taxonomic groups: 
# by station
pdf(file = file.path(figs.dir, "tax-gr-contrib_stations_sand.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.group.props.sand)
dev.off()

# pdf(file = file.path(figs.dir, "tax-gr-contrib_stations_zostera.pdf"), useDingbats = FALSE)
# plot_tax_group_contribution(tax.group.props.zostera)
# dev.off()

# by station and year
pdf(file = file.path(figs.dir, "tax-gr-contrib_st-yrs_sand.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.group.props.sand, by.years = TRUE)
dev.off()

# pdf(file = file.path(figs.dir, "tax-gr-contrib_st-yrs_zostera.pdf"), useDingbats = FALSE)
# plot_tax_group_contribution(tax.group.props.zostera, by.years = TRUE)
# dev.off() 


## frequency of occurrence of the species - by station
zoo.freq <- cbind(factors.zoo.sand$stations, num.zoo.abnd.sand)
str(zoo.freq)
names(zoo.freq) <- c("stations", names(num.zoo.abnd.sand))

# calculate relative frequencies by station



### Diversity indices ###

## Alpha diversity (Whittaker, 1960) - description of the diversity in one spot
diversity.sand <- alpha_diversity(zoo.abnd.sand)
# diversity.zostera <- alpha_diversity(zoo.abnd.zostera)
  
# save the data frame of diversity indices and measures in a file
write.csv(diversity.sand, 
          file = file.path(save.dir, "diversity.sand.csv"), 
          row.names = FALSE)

# write.csv(diversity.zostera, 
#           file = file.path(save.dir, "diversity.zostera.csv"), 
#           row.names = FALSE)


# compare the differences in diversity indices by station (ANOVA)



## Taxonomic diversity

# NB! NEEDS UPDATED MASTER SPECIES LIST FOR BIOGEOGRAPHIC AREA

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

# species accumulation curves
sp.accumulation <- specaccum(num.zoo.abnd.sand)

# plot sac with 95% confidence intervals
pdf(file = file.path(figs.dir, "sp-accumulaton_sand.pdf"), useDingbats = FALSE)
plot(sp.accumulation, ylab = "Cumulative species richness", ci.type = "polygon", ci.col = "skyblue")
dev.off()

# ranked abundance diagrams



## Species pool
# number of unseen species - various estimates
specpool(x = num.zoo.abnd.sand)


