####  Dissimilarities and environment ####

### continues from the previous scripts. tries to relate the environmental parameters
### to the observed groups of stations based on zoobenthic abundances.


## PERMANOVA - multivariate ANOVA based on dissimilarities (=> adonis in vegan)

# make a distance matrix - to be used in all computations
sand.dist <- vegdist(log(num.zoo.abnd.sand + 1))

# check the multivariate homogeneity of groups dispersions (variances)
# => assumption of PERMANOVA
sand.betadisper <- betadisper(sand.dist, group = gr.dendr.sand)  

anova(sand.betadisper)
TukeyHSD(sand.betadisper)

### there are significant differences in dispersion between groups => might
### contribute to a significant effect detected by PERMANOVA..

plot(meta.tst <- metaMDS(sand.dist, zerodist = ignore), type = "n")
points(meta.tst, select = which(gr.dendr.sand == 1), col = "red")
points(meta.tst, select = which(gr.dendr.sand == 2), col = "blue")
points(meta.tst, select = which(gr.dendr.sand == 3), col = "black")
points(meta.tst, select = which(gr.dendr.sand == 4), col = "green")
ordispider(meta.tst, group = gr.dendr.sand)

# Either proceed as is and report the difference in dispersion, or try averaging
# the data by samplings (e.g. combine the 3 replicates / sampling to obtain 18 rows
# in total)


#############################################################################################
### stratification of permutations - because repeated samplings, etc. - 
### have to account for that in the tests (rows are NOT independent)

# make a data frame with the factors & abundances - all replicates (54 rows)!
tst.permanova <- cbind(num.zoo.abnd.sand, factors.zoo.sand)

# convert replicates to numeric factor representing repeated samplings, then average 
# by replicate and station 
tst.permanova$replicates <- rep(1:3, each = 3, times = 6)
tst.permanova <- ddply(tst.permanova, .(replicates, stations), colwise(mean, .cols = is.numeric))
tst.permanova <- arrange(tst.permanova, stations)

# add in the grouping factor = groups identified in the MDS and dendrogram: 
# 1 = Kraimorie-Chukalya; 2 = Akin; 3 = Sozopol-Paraskeva; 4 = Agalina.
tst.groups <- c(rep(1, 6), rep(2, 3), rep(3, 3), rep(4, 3), rep(3, 3))
tst.permanova$groups <- tst.groups

# make a distance matrix (Bray-Curtis) and check the multivariate dispersion of 
# the groups
tst.dist <- vegdist(tst.permanova[-c(1, 2, 150)])
anova(betadisper(tst.dist, group = tst.permanova$groups))

# examine visually, through an MDS
plot(meta.tst2 <- metaMDS(tst.dist, zerodist = ignore), type = "n")
points(meta.tst2, select = which(tst.permanova$groups == 1), col = "red")
points(meta.tst2, select = which(tst.permanova$groups == 2), col = "blue")
points(meta.tst2, select = which(tst.permanova$groups == 3), col = "black")
points(meta.tst2, select = which(tst.permanova$groups == 4), col = "green")
ordispider(meta.tst2, group = tst.permanova$groups)


## looks fine this time. Go on to construct the permutation restrictions and run 
## the PERMANOVA.

### compute the true R2-value
num.tst.permanova <- tst.permanova[-which(names(tst.permanova) %in% c("stations", "replicates", "groups"))]

# remove singleton spp (rarely informative but will improve run time)
num.tst.permanova <- num.tst.permanova[ , which(specnumber(num.tst.permanova, MARGIN = 2) > 1)] 

fit <- adonis(num.tst.permanova ~ tst.permanova$groups, permutations = 1)
fit

### number of perms
B <- 1999

### setting up frame which will be populated by
### random r2 values:
pop <- rep(NA, B + 1)

### the first entry will be the true r2:
pop[1] <- fit$aov.tab[1, 5]

### set up a "permControl" object which will control how the samples are shuffled. 
### We want to test the effect of site, so we shuffle sites within years (each site
### was revisited 3 times, so there are 3 blocks of samples)
ctrl <- how(blocks = tst.permanova$replicates, within = Within(type = "series", mirror = FALSE)) 

### Number of observations:
nobs <- nrow(tst.permanova)

### check permutation (...rows represent the sample id):
### ..they are ok!
### => within each repeated sample (= sites) timepoints are shuffled,
### (e.g., for site 1: 1,2,3 - 2,3,1 - 1,3,2).
### Repeat several times to see if it's working properly. 
shuffle(nobs, control = ctrl)

### loop:
### in adonis(...) you need to put permutations = 1, otherwise
### adonis will not run
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand <- adonis(num.tst.permanova ~ tst.permanova$groups[idx], permutations = 1)
  pop[i] <- fit.rand$aov.tab[1, 5]
}

### get the p-value:
pval <- sum(pop >= pop[1])/(B + 1)

### the sign. p-value supports the H1 (-> there is a group effect).
### ..and the fact that samples are not iid is allowed by
### the customized perms - so this p-value is trustworthy as opposed
### to tests not acknowledging dependency of data points..

#########################################################################################################
## using "strata" argument from within the function: 
adonis(num.tst.permanova ~ tst.permanova$groups, strata = c(rep(1:3, times = 6)))

## looks ok... 

#### plot the groups and color prettily according to different factors (for ex. LUSI, O2, other significant...)




## cca - but biased towards rare species, so probably not the best..
cca(num.zoo.abnd.sand ~ sand + O2.bottom + Cu + Pb + env.qualit$LUSI.3000.impact, 
    data = env.all.sand.envfit)
anova(cca(num.zoo.abnd.sand ~ sand + O2.bottom + Cu + Pb + env.qualit$LUSI.3000.impact, 
          data = env.all.sand.envfit), by = "term")

## dbrda
capscale(num.zoo.abnd.sand ~ sand + O2.bottom + Cu + Pb + env.qualit$LUSI.3000.impact, 
         data = env.all.sand.envfit)

### both give distributions that are probably distorted (first - triangle, 
### second - arch...) - either change constraints - choose more judiciously, or 
### don't use either analysis..

