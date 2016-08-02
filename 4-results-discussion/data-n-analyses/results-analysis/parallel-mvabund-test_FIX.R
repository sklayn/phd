library(foreach)
library(doParallel)

### run mvabund permutations in parallel - TESTING! ###

## DOES NOT work as I want it -> blazing fast, but currently performs several 
## anovas, each with a certain number of the permutations I want, in parallel, 
## but how to combine the result? 

# set up the number of iterations
nb.boots = 250

# find out how many cores the computer has, and make the cluster
nb.cores <- detectCores()
clust <- makeForkCluster(nb.cores - 2) # in case
registerDoParallel(clust)

#start time - for info
strt<-Sys.time()

# run the analysis
ls <- foreach(i = rep(nb.boots, 4)) %dopar% {
  to.ls <- anova.manyglm(sand.glm.lusi, 
                         test = "LR",
                         p.uni = "adjusted",
                         nBoot = nb.boots,
                         show.time = "all")
}

print(Sys.time()-strt)
stopCluster(clust)


## try with manylm / anova.manylm
