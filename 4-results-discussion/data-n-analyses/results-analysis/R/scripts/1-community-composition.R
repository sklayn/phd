##### Community composition & structure analysis #####

# define the working subdirectories
data.dir <- "data"
functions.dir <- "R"
save.dir <- "output"
figures.dir <- "figs"

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

# TEST TO SEE IF WORKS WITH CURRENT IMPORT FUNCTION! (it should - splits labels on . now)
# zoo.abnd.2012 <- import_zoo_data(data.dir = data.dir, 
#                             zoo.data = "zoo-abnd-2012.csv", 
#                             stations = stations.2012, 
#                             repl = 4)
# add a column for habitat type (2 adjacent habitats sampled in 2012 - seagrass and bare sand);
# place it in the beginning with the other factors
# zoo.abnd.2012 <- data.frame(append(zoo.abnd.2012, 
#                                    list(habitat = rep(c("zostera", "sand"), each = 24)), 
#                                    after = match("stations", names(zoo.abnd.2012))))


# also import the biomass data (if not already imported by then)
zoo.biomass.sand <- import_zoo_data(data.dir = data.dir, 
                                    zoo.data = "zoo-biomass-sand.csv", 
                                    station.names = stations.sand, 
                                    repl = 3)



# after-import sanity checks (better safe than sorry!)
str(zoo.abnd.sand)  # structure & variable classes: factors are factors, numeric variables
                    # are numeric, etc.

# str(zoo.abnd.zostera)
# str(zoo.abnd.2012)

names(zoo.abnd.sand)  # variable (species) names
# names(zoo.abnd.zostera)
# names(zoo.abnd.2012)


# quick check for import errors, etc. for the biomass, too
str(zoo.biomass.sand)
names(zoo.biomass.sand)


## SAVE THE CLEANED AND REARRANGED DATA
write.csv(zoo.abnd.sand, 
          file = file.path(save.dir, "zoo-sand-clean.csv"), 
          row.names = FALSE)

# write.csv(zoo.abnd.zostera, 
#           file = file.path(save.dir, "zoo-zostera-clean.csv"), 
#           row.names = FALSE)


## SAVE THE CLEANED BIOMASS DATA
write.csv(zoo.biomass.sand, 
          file.path(save.dir, "biomass_sand_clean.csv"), 
          row.names = FALSE)


# make subsets of the numeric columns (the species data), and the factor columns -
# to avoid retyping the command, because these subsets will be 
# used a lot in all the subsequent analyses. 
factors.zoo.sand <- zoo.abnd.sand[sapply(zoo.abnd.sand, is.factor)]
num.zoo.abnd.sand <- zoo.abnd.sand[sapply(zoo.abnd.sand, is.numeric)]

# biomass
num.biomass.sand <- zoo.biomass.sand[sapply(zoo.biomass.sand, is.numeric)]

# leave only the species actually present in the current dataset (to simplify later
# analyses)
num.zoo.abnd.sand <- num.zoo.abnd.sand[colSums(num.zoo.abnd.sand) > 0] 
num.biomass.sand <- num.biomass.sand[colSums(num.biomass.sand) > 0]

# factors.zoo.zostera <- zoo.abnd.zostera[sapply(zoo.abnd.zostera, is.factor)]
# num.zoo.abnd.zostera <- zoo.abnd.zostera[sapply(zoo.abnd.zostera, is.numeric)]
# num.zoo.abnd.zostera <- num.zoo.abnd.zostera[colSums(num.zoo.abnd.zostera) > 0] 

# average nb species, abundance & biomass for the whole dataset, min-max values of 
# same + where/when they were found (CUSTOM FUNCTION)
summary_zoo_params(num.zoo.abnd.sand, factors.zoo.sand, "abnd") # abundance
summary_zoo_params(num.biomass.sand, factors.zoo.sand, "biomass") # biomass
summary_zoo_params(num.zoo.abnd.sand, factors.zoo.sand, "nb.sp") # number of species


# summary of numeric abundance data by station (mean) - for summary plots, etc.
summary.abnd.sand <- ddply(zoo.abnd.sand, .(stations), colwise(mean, .cols = is.numeric))
# summary.abnd.zostera <- ddply(zoo.abnd.zostera, .(stations), colwise(mean, .cols = is.numeric))

## average nb species, abundance & biomass per station/year
# NB exclude (mostly invalid) taxa - defined several rows below 
sp.summaries <- num.zoo.abnd.sand[, -which(colnames(num.zoo.abnd.sand) %in% sp2exclude)]
# nb species 
aver.nb.sp <- data.frame(station = factors.zoo.sand$stations,
                         year = factors.zoo.sand$years,
                         nb.sp = specnumber(sp.summaries))

# aver nb species by station and year
ddply(aver.nb.sp, 
      .(station, year), 
      summarize, aver.nb.sp = mean(nb.sp), 
      nb.sp.sd = sd(nb.sp))

# abundance - DO NOT exclude those taxa from here, just from the number of taxa
aver.abnd.sand <- data.frame(station = factors.zoo.sand$stations,
                             year = factors.zoo.sand$years,
                             tot.abnd = rowSums(zoo.abnd.sand[sapply(zoo.abnd.sand, is.numeric)]))

# average abundance by station and year
ddply(aver.abnd.sand, 
      .(station, year), 
      summarize, aver.abnd = mean(tot.abnd), 
      abnd.sd = sd(tot.abnd))

# same crap, different biomass
aver.biomass.sand <- data.frame(station = factors.zoo.sand$stations,
                             year = factors.zoo.sand$years,
                             tot.biomass = rowSums(zoo.biomass.sand[sapply(zoo.biomass.sand, is.numeric)]))

# average biomass by station and year
ddply(aver.biomass.sand, .(station, year), 
      summarize, 
      aver.biomass = mean(tot.biomass), 
      biomass.sd = sd(tot.biomass))

rm(aver.biomass.sand, aver.abnd.sand, aver.nb.sp, sp.summaries) # crap cleanup

## Basic taxonomic composition and structure
# import taxonomic data (species names are rows, successively higher taxonomic
# ranks are columns)
zoo.taxa <- read.csv(file.path(data.dir, "zoo-taxonomy.csv"), 
                     header = TRUE, 
                     row.names = 1)

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
                                ifelse(subclass == "Eumalacostraca" | #subclasses because otherwise Cirripedia will get stuck in Varia
                                       subclass == "Thecostraca", "Crustacea", "Varia"))))

table(current.taxa.sand$group)


# current.taxa.zostera$group <- with(current.taxa.zostera, 
#                                    ifelse(class == "Polychaeta", "Polychaeta", 
#                                    ifelse(class == "Bivalvia" | 
#                                           class == "Gastropoda" | 
#                                           class == "Polyplacophora", "Mollusca", 
#                                    ifelse(subclass == "Eumalacostraca" |
#                                           subclass == "Thecostraca", "Crustacea", "Varia"))))
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

## REAL number of taxa (without the species identified to genus/family level 
## because of bad preservation, the larvae,..)
sp2exclude <- c("Abra sp.", "Eteone sp.", "Glycera sp.", "Hydrobia sp.", "Iphinoe sp.", 
                "Micrpodeutopus sp.", "Microphthalmus sp.", "Gammaridae", "Amphipoda", 
                "Cardiidae", "Decapoda larvae", "Protodrilus sp.", "Maldanidae", 
                "Nepthtyidae", "Nereididae", "Pisces larvae", "Polychaeta larvae", "Spionidae")

sort(table(current.taxa.sand[-which(row.names(current.taxa.sand) %in% sp2exclude), ]$group))

# save this plot
pdf(file = file.path(figs.dir, "nb-taxa_sand_actual.pdf"), useDingbats = FALSE)
barchart(sort(table(current.taxa.sand[-which(row.names(current.taxa.sand) %in% sp2exclude), ]$group)),
         xlab = "", 
         scales = list(tck = c(1,0), x = list(cex = 1.2), y = list(cex = 1.2)),
         col = "royalblue")
dev.off()

## calculate the contribution of each taxonomic group to the numeric community 
## composition by station and by year
# first calculate proportions by taxonomic group in each station/replicate
tax.group.props.sand <- tax_group_contribution(num.zoo.abnd.sand, 
                                               current.taxa.sand$group)
tax.group.props.sand <- cbind(factors.zoo.sand, tax.group.props.sand)

# tax.group.props.zostera <- tax_group_contribution(num.zoo.abnd.zostera, 
#                                                   current.taxa.zostera)  
# tax.group.props.zostera <- cbind(factors.zoo.zostera, tax.group.props.zostera)

# idem for biomass
tax.gr.props.biomass.sand <- tax_group_contribution(num.zoo.biomass.sand, 
                                                    current.taxa.sand$group)
tax.gr.props.biomass.sand <- cbind(factors.zoo.sand, tax.gr.props.biomass.sand)


# save to file (in case) 
write.csv(tax.group.props.sand, 
          file = file.path(save.dir, "tax-gr-proportions_sand-clean.csv"), 
          row.names = FALSE)

# write.csv(tax.group.props.zostera, 
#           file = file.path(save.dir, "tax-gr-proportions_zostera-clean.csv"), 
#           row.names = FALSE)

# idem biomass
write.csv(tax.gr.props.biomass.sand, 
          file = file.path(save.dir, "tax-gr-proportions-biomass_sand_clean.csv"), 
          row.names = FALSE)

# plot the contribution of taxonomic groups: 
# add factors station and year (expected by plotting function) if not already in df

# by station
pdf(file = file.path(figs.dir, "tax-gr-contrib_stations_sand.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.group.props.sand)
dev.off()

# pdf(file = file.path(figs.dir, "tax-gr-contrib_stations_zostera.pdf"), useDingbats = FALSE)
# plot_tax_group_contribution(tax.group.props.zostera)
# dev.off()

# for biomass, too
pdf(file = file.path(figs.dir, "tax-gr-contrib-biomass_st_sand.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.gr.props.biomass.sand)
dev.off()


# by station and year
pdf(file = file.path(figs.dir, "tax-gr-contrib_st-yrs_sand.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.group.props.sand, by.years = TRUE)
dev.off()

# pdf(file = file.path(figs.dir, "tax-gr-contrib_st-yrs_zostera.pdf"), useDingbats = FALSE)
# plot_tax_group_contribution(tax.group.props.zostera, by.years = TRUE)
# dev.off() 

#... and for biomass
pdf(file = file.path(figs.dir, "tax-gr-contrib-biomass_st-yrs_sand.pdf"), useDingbats = FALSE)
plot_tax_group_contribution(tax.gr.props.biomass.sand, by.years = TRUE)
dev.off()

## calculate the overall proportions of the taxonomic groups (all stations/years pooled)
tax.gr.props.overall <- apply(tax.group.props.sand[-c(1:3)], 2, mean)
# convert to data frame fr easier plotting 
tax.gr.props.overall <- data.frame(tax.gr = names(tax.gr.props.overall), 
                                   prop = tax.gr.props.overall, 
                                   row.names = 1:4)
sum(tax.gr.props.overall$prop) # check calculation
# rearrange in descending order; also the factor, because otherwise new order 
# won't be recognized by ggplot 
tax.gr.props.overall <- arrange(tax.gr.props.overall, desc(prop))
tax.gr.props.overall$tax.gr <- factor(tax.gr.props.overall$tax.gr, 
                                      levels = tax.gr.props.overall$tax.gr)

# plot
pdf(file.path(figs.dir, "tax-gr-contribution_overall_sand.pdf"), useDingbats = FALSE)
ggplot(tax.gr.props.overall, aes(x = tax.gr, y = prop)) + 
  geom_bar(stat = "identity", fill = "royal blue") + 
  theme_bw() + 
  labs(x = "", y = "%") +
  theme(axis.text.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.3)))
dev.off()

rm(tax.gr.props.overall)  # only temporary anyway

##... rinse and repeat for biomass
tax.gr.props.biomass.overall <- apply(tax.gr.props.biomass.sand[-c(1:2)], 2, mean)
# convert to data frame fr easier plotting 
tax.gr.props.biomass.overall <- data.frame(tax.gr = names(tax.gr.props.biomass.overall), 
                                           prop = tax.gr.props.biomass.overall, 
                                           row.names = 1:4)
sum(tax.gr.props.biomass.overall$prop) # check calculation... all good!

# rearrange in descending order; also the factor, because otherwise new order 
# won't be recognized by ggplot 
tax.gr.props.biomass.overall <- arrange(tax.gr.props.biomass.overall, desc(prop))
tax.gr.props.biomass.overall$tax.gr <- factor(tax.gr.props.biomass.overall$tax.gr, 
                                              levels = tax.gr.props.biomass.overall$tax.gr)

# plot
pdf(file.path(figs.dir, "tax-gr-contribution_biomass_overall_sand.pdf"), useDingbats = FALSE)
ggplot(tax.gr.props.biomass.overall, aes(x = tax.gr, y = prop*100)) + 
  geom_bar(stat = "identity", fill = "dark orange") + 
  theme_bw() + 
  labs(x = "", y = "%") +
  theme(axis.text.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.3)))
dev.off()

rm(tax.gr.props.biomass.overall)  # only temporary anyway


## number of families/genera/sp per taxonomic group - fugly!
## Mollusca
# nb families
unique(subset(current.taxa.sand, group == "Mollusca")$family)
# number of species of class Bivalvia & Gastropoda
length(row.names(subset(current.taxa.sand, group == "Mollusca" & class == "Bivalvia")))
length(row.names(subset(current.taxa.sand, group == "Mollusca" & class == "Gastropoda")))
# which species, to check & remove invalid ones (unfortunately, manually)
row.names(subset(current.taxa.sand, group == "Mollusca" & class == "Bivalvia"))
row.names(subset(current.taxa.sand, group == "Mollusca" & class == "Gastropoda"))

# number of species in each family (descending order)
## NB still has to be fixed - for ex. remove sp considered invalid beforehand.. Or check 
## manually post factum
sort(table(subset(current.taxa.sand, group == "Mollusca" & class == "Bivalvia")$family), 
     decreasing = T)

## Polychaeta
# nb families
unique(subset(current.taxa.sand, group == "Polychaeta")$family)

# which species, to check & remove invalid ones (unfortunately, manually)
row.names(subset(current.taxa.sand, group == "Polychaeta"))

# number of species in each family (descending order)
## NB still has to be fixed - for ex. remove sp considered invalid beforehand.. Or check 
## manually post factum
sort(table(subset(current.taxa.sand, group == "Polychaeta")$family), 
     decreasing = T)
row.names(subset(current.taxa.sand, group == "Polychaeta" & family == "Nereididae"))
row.names(subset(current.taxa.sand, group == "Polychaeta" & family == "Phyllodocidae")) #etc...

## Crustacea
# nb families
unique(subset(current.taxa.sand, group == "Crustacea")$family)

# which species, to check & remove invalid ones (unfortunately, manually)
row.names(subset(current.taxa.sand, group == "Crustacea"))

# number of sp in each order
sort(table(subset(current.taxa.sand, group == "Crustacea")$order), decreasing = TRUE)
row.names(subset(current.taxa.sand, group == "Crustacea" & order == "Amphipoda"))
row.names(subset(current.taxa.sand, group == "Crustacea" & order == "Decapoda"))
row.names(subset(current.taxa.sand, group == "Crustacea" & order == "Cumacea"))
row.names(subset(current.taxa.sand, group == "Crustacea" & order == "Isopoda"))

# number of species in each family (descending order)
## NB still has to be fixed - for ex. remove sp considered invalid beforehand.. Or check 
## manually post factum
sort(table(subset(current.taxa.sand, group == "Crustacea")$family), 
     decreasing = T)
row.names(subset(current.taxa.sand, group == "Crustacea" & family == "Bodotriidae"))
row.names(subset(current.taxa.sand, group == "Crustacea" & family == "Corophiidae")) #etc...


## Varia
# nb families
unique(subset(current.taxa.sand, group == "Varia")$class)

# which species, to check & remove invalid ones (unfortunately, manually)
row.names(subset(current.taxa.sand, group == "Varia"))

# number of sp in each order
sort(table(subset(current.taxa.sand, group == "Varia")$class), decreasing = TRUE)


rm(sp2exclude)


## frequency of occurrence of the species - all sites combined
zoo.freq <- cbind(num.zoo.abnd.sand, factors.zoo.sand$stations)
names(zoo.freq) <- c(names(num.zoo.abnd.sand), "station")
str(zoo.freq)

# aggregate by station
zoo.freq <- ddply(zoo.freq, .(station), colwise(mean, .cols = is.numeric))

# calculate frequency of occurrence of each species
fr <- apply(zoo.freq[sapply(zoo.freq, is.numeric)], 2, function(x) sum(x > 0) / length(x))

# sort in descending order of frequency of occurrence
sort(fr, decreasing = TRUE)

rm(fr)

### see which species are the most abundant at each station
## split abundance data frame into separate dfs for each station
abnd.dfs.st <- split(zoo.abnd.sand, zoo.abnd.sand$stations)
str(abnd.dfs.st)
# get rid of factor columns
abnd.dfs.st <- lapply(abnd.dfs.st, "[", -c(1:3))

## arrange each df's columns (species) in descending order of abundance
# calculate the total abundance of each species, then sort in descending order
tot.abnd.dfs <- lapply(abnd.dfs.st, colSums)

# make sure to return column indices, not names!
tot.abnd.sorted <- lapply(tot.abnd.dfs, sort, decreasing = TRUE, index.return = TRUE)

# extract the given number of species (columns) from the abundance data frame by their index
abnd.subs <- mapply(function(m, n) m[, n$ix[1:10]], 
                    abnd.dfs.st, 
                    tot.abnd.sorted, 
                    SIMPLIFY = FALSE)


# convert the df to long format -> easier for ggplot2 to handle  
abnd.melted <- lapply(abnd.subs, melt)

## use custom transformation for abundance values (log(y/min + 1)) - as in the 
## mvabund package
# calculate the minimum non-0 value in the dataset & transform y
min.vals <- lapply(abnd.melted, function(x) min(x$value[x$value > 0]))
abnd.fin <- mapply(function(m, n) {
                      value.tr <- log(m$value / n + 1)
                      m$value.tr <- value.tr
                      return(m)
                    }, 
                    abnd.melted,
                    min.vals,
                    SIMPLIFY = FALSE)  

# plot all subsets, and add the corresponding name (the reason we're using mapply)
plots <- mapply(function(n, m) {
                  ggplot(n, aes_string(x = "variable", y = "value.tr")) + 
                    geom_point(size = 2.5) + 
                    # reverse the order of x axis, so highest-contributing species are on top
                    scale_x_discrete(name = "", limits = rev(levels(n$variable))) + 
                    scale_y_continuous(name = "Abundance (log(y/min + 1))") + 
                    coord_flip() + # put species on y axis - easier to read
                    labs(title = m) + 
                    theme_bw() + 
                    theme(axis.text.x = element_text(size = rel(1.3)), 
                          axis.text.y = element_text(size = rel(1.3)))
                }, 
                abnd.fin, 
                names(abnd.fin), 
                SIMPLIFY = FALSE)

# arrange them all on a single plotting device (gridExtra)  
n <- length(plots)
nCol <- floor(sqrt(n))
p <- do.call("grid.arrange", c(plots, ncol = nCol))

# keep for exploratory purposes; otherwise not very useful (y scales not matching,
# overlap, ..)
ggsave(file.path(figs.dir, "explor_most-abnd-sp-stations_sand.pdf"), p, 
       height = 15, width = 15,
       dpi = 300)

# get rid of the clutter
rm(p, n, nCol, plots, abnd.fin, min.val, abnd.melted, abnd.subs, tot.abnd.sorted, 
   tot.abnd.dfs, abnd.dfs.st)

### see which species contribute most to the biomass at each station
## split biomass data frame into separate dfs for each station
biomass.dfs.st <- split(zoo.biomass.sand, zoo.biomass.sand$stations)
str(biomass.dfs.st)
# get rid of factor columns
biomass.dfs.st <- lapply(biomass.dfs.st, "[", -c(1:3))

## arrange each df's columns (species) in descending order of biomass
# calculate the total biomass of each species, then sort in descending order
tot.biomass.dfs <- lapply(biomass.dfs.st, colSums)

# make sure to return column indices, not names!
tot.biomass.sorted <- lapply(tot.biomass.dfs, sort, decreasing = TRUE, index.return = TRUE)

# extract the given number of species (columns) from the biomass data frame by their index
biomass.subs <- mapply(function(m, n) m[, n$ix[1:10]], 
                      biomass.dfs.st, 
                      tot.biomass.sorted, 
                      SIMPLIFY = FALSE)


# convert the df to long format -> easier for ggplot2 to handle  
biomass.melted <- lapply(biomass.subs, melt)

## use custom transformation for abundance values (log(y/min + 1)) - as in the 
## mvabund package
# calculate the minimum non-0 value in the dataset & transform y
min.vals <- lapply(biomass.melted, function(x) min(x$value[x$value > 0]))
biomass.fin <- mapply(function(m, n) {
                        value.tr <- log(m$value / n + 1)
                        m$value.tr <- value.tr
                        return(m)
                      }, 
                      biomass.melted,
                      min.vals,
                      SIMPLIFY = FALSE)  

# plot all subsets, and add the corresponding name (the reason we're using mapply)
plots <- mapply(function(n, m) {
                    ggplot(n, aes_string(x = "variable", y = "value.tr")) + 
                      geom_point(size = 2.5) + 
                      # reverse the order of x axis, so highest-contributing species are on top
                      scale_x_discrete(name = "", limits = rev(levels(n$variable))) + 
                      scale_y_continuous(name = "Abundance (log(y/min + 1))") + 
                      coord_flip() + # put species on y axis - easier to read
                      labs(title = m) + 
                      theme_bw() + 
                      theme(axis.text.x = element_text(size = rel(1.3)), 
                            axis.text.y = element_text(size = rel(1.3)))
                  }, 
                  biomass.fin, 
                  names(biomass.fin), 
                  SIMPLIFY = FALSE)

# arrange them all on a single plotting device (gridExtra)  
n <- length(plots)
nCol <- floor(sqrt(n))
p <- do.call("grid.arrange", c(plots, ncol = nCol))

# keep for exploratory purposes; otherwise not very useful (y scales not matching,
# overlap, ..)
ggsave(file.path(figs.dir, "explor_most-biomass-sp-stations_sand.pdf"), p, 
       height = 15, width = 15,
       dpi = 300)

# get rid of the clutter
rm(p, n, nCol, plots, biomass.fin, min.vals, biomass.melted, biomass.subs, 
   tot.biomass.sorted, tot.biomass.dfs, biomass.dfs.st)


## more useful: plot most abundant species OVERALL; shape/colour by station
p <- plot_most_abnd_sp(num.zoo.abnd.sand, as.factor(as.numeric(factors.zoo.sand$stations)), 
                       nb.sp = 20)
ggsave(file.path(figs.dir, "most-abnd-sp_stations_sand.pdf"),
       p,
       width = 15, height = 15)

rm(p)

##... idem biomass: most contributing species OVERALL; shape/colour by station
### NB have to fix labels on plot - hardcoded for abundance!!
p <- plot_most_abnd_sp(num.zoo.biomass.sand, as.factor(as.numeric(factors.zoo.sand$stations)), 
                       nb.sp = 20)
ggsave(file.path(figs.dir, "most-biomass-sp_stations_sand.pdf"),
       p,
       width = 15, height = 15)

rm(p)




##### Diversity indices #####

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



#### Taxonomic diversity ####

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

# calculate the taxonomic distance for the species list, using variable step 
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


