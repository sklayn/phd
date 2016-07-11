## significance pruning of insignificant (mostly noise) explanatory variables
## - way too many to use for GLMs, etc. (choke the computer, parallel processing or no)
library(vtreat)
library(ggExtra)
library(corrplot)
library(caret)
library(Hmisc)

# prepare data
env.sand.pruning <- env.all.sand.envfit[-which(names(env.all.sand.envfit) %in% c("month", "year"))] 

# add LUSI impact and recode to make binary (high/low impact stations) - will 
# try pruning data variables based on that 
env.sand.pruning$LUSI.3000.impact <- env.qualit$LUSI.3000.impact
env.sand.pruning$LUSI.3000.impact.binary <- with(env.sand.pruning, ifelse(LUSI.3000.impact == "moderate" | 
                                                           LUSI.3000.impact == "high", 
                                                         "high", 
                                                         "low"))

env.sand.pruning$LUSI.3000.impact.binary <- as.factor(env.sand.pruning$LUSI.3000.impact.binary)

### in the end does not work so well, prob. too unbalanced - most of the stations 
### are under "high" impact; there are only a few under "low". Will stick with the 
### groups from the MDS or with the stations themselves..

# keep only unique observations (= no repeating of rows to match the number of 
# replicate zoobenthic samples)
env.sand.pruning <- env.sand.pruning[seq(1, nrow(env.sand.pruning), by = 3), ]

# add the groups 
env.sand.pruning$group <- tst.groups
env.sand.pruning$group <- as.factor(env.sand.pruning$group)

### only x-aware scaling - prune variables, scale & fit model

################################################################################################################
## helper functions - CHECK THESE SOMETIME! HARD-CODED Y in one place, ETC.

# get the theoretical significance glm model (wrt deviance)
get_significance <- function(model) {
  delta_deviance <- model$null.deviance - model$deviance
  df <- model$df.null - model$df.residual
  sig <- pchisq(delta_deviance, df, lower.tail=FALSE)
}

# return a frame of the deviance scores on the permuted data
permutation_test <- function(dataf, ycol, nperm) {
  nrows <- dim(dataf)[1]
  y <- dataf[[ycol]]
  X <- dataf[, setdiff(colnames(dataf), ycol), drop=FALSE]
  varnames = colnames(X)
  fmla <- paste("y", paste(varnames, collapse = "+"), sep = "~")
  deviances <- numeric(nperm)
  for(i in seq_len(nperm)) {
    # random order of rows
    ord <- sample.int(nrows, size=nrows, replace=FALSE)
    model <- glm(fmla, data = cbind(y = y[ord], X),
                 family = binomial(link = "logit"))
    #print(summary(model))
    deviances[[i]] <- model$deviance
  }
  deviances
}

score_variable <- function(dframe, ycol, var, nperm,
                           title = '') {
  df <- data.frame(y = dframe[[ycol]], x = dframe[[var]])
  
  mod <- glm("y~x", data = df,
             family = binomial(link="logit"))
  vdev <- mod$deviance
  vperm <- permutation_test(df, ycol, nperm)
  
  print(summary(mod))
  
  # count how many times vdev >= deviances from perm test
  num <- sum(vperm <= vdev)
  vscore <- num / nperm
  
  print(ggplot(data.frame(nullperm = vperm), aes(x = nullperm)) +
          geom_density() + 
          geom_vline(xintercept = vdev, color = 'red') +
          ggtitle(paste(title, "left tail area ~", vscore)))
  
  print(paste("Chi-sq estimate:", get_significance(mod)))
  
}

# get the signal scores for the variables
# in a data set
# assume output is a binary variable named y (?)
get_chiscores <- function(dframe, varnames) {
  nvar <- length(varnames)
  scores <- numeric(nvar)
  for(i in seq_len(nvar)) {
    model <- glm(paste("group~", varnames[i]), 
                 dframe,
                 family = binomial(link = "logit"))
    scores[i] <- get_significance(model)
  }
  
  sframe <- data.frame(var = varnames,
                       scores = scores, 
                       stringsAsFactors = FALSE)
  sframe
}

#
# Plot the scores of each variable
# frm has columns var and scores
# (output of get_chiscores)
#
scoreplot <- function(frm, threshold, sort = 1) {
  n <- dim(frm)[1]
  frm$var <- reorder(frm$var, frm$scores * sort, FUN = sum)
  frm$goodvar <- frm$scores < threshold
  
  ggplot(frm, aes(x = var, y = scores, ymin = 0, ymax = scores, color = goodvar)) +
    geom_pointrange() +
    geom_hline(yintercept = threshold, color = "red", linetype = 2) +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "darkgray")) +
    theme(legend.position = "none")
}

################################################################################################

### variable pruning based on significance. y = groups determined by MDS/classification 
## (-> correspond more or less to level of impact)

# remove the LUSI columns (decided they're not useful here). Keep the stations, in case
# (but have to be excluded in the next lines).
env.sand.pruning <- env.sand.pruning[-which(names(env.sand.pruning) %in% 
                                    c("LUSI.3000.impact", "LUSI.3000.impact.binary"))]

sframe <- get_chiscores(dframe = env.sand.pruning, 
                        varnames = setdiff(colnames(env.sand.pruning[-1]), "group"))


threshold <- 0.05  # be generous in accepting a variable (1/20 false positive rate)
scoreplot(sframe, threshold)
vars.sel <- sframe[sframe$scores < threshold,]$var

# find correlated variables among those selected -> could be eliminated without 
# losing too much information
sand.pruned.vars <- subset(env.sand.pruning, select = vars.sel)

# use this function to calculate the correlations, because it also gives the p-values
sand.pruned.vars.cors <- rcorr(as.matrix(sand.pruned.vars))

# plot correlation matrix to inspect visually (and save as pdf in case needed later)
pdf(file = file.path(figs.dir, "explor_corrplot_env-vars-pruned.pdf"), 
    paper = "a4r", width = 12, height = 12, 
    useDingbats = FALSE)
corrplot(sand.pruned.vars.cors$r, type ="lower", order = "hclust", tl.col="black", tl.srt = 45)
dev.off()


# get rid of the most strongly correlated variables (redundant and cumbersome if fitting 
# a model)

# set the cutoff value of r2 pretty high, anyway. NB names = FALSE, to return column index
highly.cor.vars <- findCorrelation(sand.pruned.vars.cors$r, cutoff = 0.8, names = FALSE)
sand.pruned.vars.fin <- sand.pruned.vars[, -highly.cor.vars]
names(sand.pruned.vars.fin)

# remove one of the oxygen variables, too - no need for two, even if they are not correlated
sand.pruned.vars.fin$O2.average <- NULL


# PCA with pruned dataset, but only x-aware (variables scaled and centered to cope 
# with different unit variances) => this is just for checking if the 
# stations/groups get separated as expected, and with a satisfactory amount of 
# variance explained by the retained variables. 
tst.princ <- prcomp(sand.pruned.vars.fin, center = TRUE, scale. = TRUE) 

# check the results
summary(tst.princ)
biplot(tst.princ)


# plot variable magnitudes -> to check (visually) PC significance 
# df2 <- data.frame(pc = 1:length(tst.princ$sdev), 
#                   magnitude = tst.princ$sdev)
# 
# ggplot(df2, mapping = aes(x = pc, y = magnitude)) + geom_point()
# 
# 


## cleanup
rm(sframe, highly.cor.vars, threshold, sand.pruned.vars, 
   sand.pruned.vars.cors, vars.sel, tst.princ, env.sand.pruning)



############################################################################################################
#############################################################################################################
### y-aware scaling (prior to PCA,..) => some other time, esp. with well-formulated 
### binary categorical outcome / numeric outcome..

## categorical outcome (binary)
# design treatment plan
treatmentsC <- designTreatmentsC(env.sand.pruning, 
                                 varlist = names(env.sand.pruning)[-c(1, 34)], 
                                 outcomename = "LUSI.3000.impact.binary", 
                                 outcometarget = "high")

scoreFrame <- treatmentsC$scoreFrame
ggplot(scoreFrame, mapping = aes(x = varName, y = sig)) + geom_point() + coord_flip()

# prepare the treated frames, with y-aware scaling
tst.y.scaled <- prepare(treatmentsC, env.sand.pruning, pruneSig = 0.05, scale = TRUE) 
# varRestriction = paste(vars.sel, "_clean", sep = ""))

vars <- setdiff(colnames(tst.y.scaled), "LUSI.3000.impact.binary")

# prcomp defaults to scale. = FALSE, but we already scaled/centered in vtreat- which we don't want to lose.
dmTst <- as.matrix(tst.y.scaled[, vars])
tst.princ <- prcomp(dmTst, center = FALSE, scale. = FALSE)

ggplot(data = data.frame(pc = 1:length(tst.princ$sdev), 
                         magnitude = tst.princ$sdev), 
       aes(x = pc, y = magnitude)) + geom_point()


#####################################################################################################
#####################################################################################################