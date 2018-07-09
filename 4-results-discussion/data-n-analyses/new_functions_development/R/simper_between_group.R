simper_between_group <- function (comm, group, permutations = 0, trace = FALSE, parallel = getOption("mc.cores"), 
    ...) 
{
    EPS <- sqrt(.Machine$double.eps) ## this is the square root of the smallest positive floating-point number x such that 1 + x != 1 on the current machine R is running on. Used later to calculate the p-value from the permutations.
    
    ## check for empty sites (should not be included in the analysis) 
    if (any(rowSums(comm, na.rm = TRUE) == 0)) 
        warning("you have empty rows: results may be meaningless")
    
    pfun <- function(x, comm, comp, i, contrp) {
        ## helper calculating contribution of species to Bray-Curtis distance between groups for one permutation from the matrix of permutations - will be used for more efficient calculation of a large number of permutations.
        
        ## get the groups in the current permutation 
        groupp <- group[perm[x, ]]
        
        ## calculate the Bray-Curtis distance and the contribution of each species to it for each pair of groups 
        ga <- comm[groupp == comp[i, 1], , drop = FALSE]
        gb <- comm[groupp == comp[i, 2], , drop = FALSE]
        n.a <- nrow(ga)
        n.b <- nrow(gb)
        for (j in seq_len(n.b)) {
            for (k in seq_len(n.a)) {
                mdp <- abs(ga[k, , drop = FALSE] - gb[j, , drop = FALSE])
                mep <- ga[k, , drop = FALSE] + gb[j, , drop = FALSE]
                contrp[(j - 1) * n.a + k, ] <- mdp/sum(mep)
            }
        }
        ## calculate the average contribution of each species 
        colMeans(contrp)
    }
    
    ## transform input abundances to matrix 
    comm <- as.matrix(comm)
    
    ## get matrix of all unique pairwise combinations of groups
    comp <- t(combn(unique(as.character(group)), 2))
    
    ## make list to hold SIMPER results
    outlist <- NULL
    
    ## get number of columns (=species) and number of rows (=sites) in input abundance matrix 
    P <- ncol(comm)
    nobs <- nrow(comm)
    
    ## make a permutation matrix using the input number of permutations and the number of sites, where each row gives indices of observations for a permutation.  
    perm <- getPermuteMatrix(permutations, nobs, ...)
    
    ## check if the generated permutation matrix is correct (each permutation should contain all observations, only scrambled, i.e. the nb rows/sites in the input community matrix should equal the nb columns in the permutation matrix)   
    if (ncol(perm) != nobs) 
        stop(gettextf("'permutations' have %d columns, but data have %d rows", 
            ncol(perm), nobs))

    nperm <- nrow(perm)    
    if (nperm > 0) 
        ## make a matrix to hold the results for each species' contribution to the group dissimilarity in each permutation    
        perm.contr <- matrix(nrow = P, ncol = nperm)
    
    ## check if option for parallel computation of permutations is selected; if yes, make a cluster and run the permutations in parallel 
    if (is.null(parallel)) 
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    isParal <- hasClus || parallel > 1
    isMulticore <- .Platform$OS.type == "unix" && !hasClus
    
    if (isParal && !isMulticore && !hasClus) {
        parallel <- makeCluster(parallel)
    }
    
    ### 1. Calculate the observed contribution of species to the Bray-Curtis dissimilarity between groups 
    ## for each pairwise combination of groups,  
    for (i in seq_len(nrow(comp))) {
        ## get the abundance matrix in the pair of groups 
        group.a <- comm[group == comp[i, 1], , drop = FALSE]
        group.b <- comm[group == comp[i, 2], , drop = FALSE]
        
        ## get the number of sites/observations in each group
        n.a <- nrow(group.a)
        n.b <- nrow(group.b)
        
        ## make a matrix to hold the species contributions to each of the combinations of observations in the groups
        contr <- matrix(ncol = P, nrow = n.a * n.b)
        
        ## for each observation in the two groups,  
        for (j in seq_len(n.b)) {
            for (k in seq_len(n.a)) {
                ## calculate the contribution of each species to the overall Bray-Curtis dissimilarity between groups, and put the result in the matrix of species contributions 
                md <- abs(group.a[k, , drop = FALSE] - group.b[j, 
                  , drop = FALSE])
                me <- group.a[k, , drop = FALSE] + group.b[j, 
                  , drop = FALSE]
                contr[(j - 1) * n.a + k, ] <- md/sum(me)  ## this complicated row index is so that the result is a triangular matrix
            }
        }
        
        ## calculate the observed average contribution of each species to each combination of observations
        average <- colMeans(contr)

        ### 2. Do permutations if such are required at input (-> for significance testing)
        if (nperm > 0) {
            if (trace) 
                cat("Permuting", paste(comp[i, 1], comp[i, 2], 
                  sep = "_"), "\n")
            ## make a matrix to hold the results of each (same way and shape as in the calculation of the observed contributions
            contrp <- matrix(ncol = P, nrow = n.a * n.b)
            
            ## if the parallel option is selected, do the calculations in parallel
            if (isParal) {
                if (isMulticore) {
                  perm.contr <- mclapply(seq_len(nperm), function(d) pfun(d, 
                    comm, comp, i, contrp), mc.cores = parallel)
                  perm.contr <- do.call(cbind, perm.contr)
                }
                else {
                  perm.contr <- parSapply(parallel, seq_len(nperm), 
                    function(d) pfun(d, comm, comp, i, contrp))
                }
            }
            else {
                ## for each permutation, perform the calculation of contributions using the pre-defined helper function 
                perm.contr <- sapply(1:nperm, function(d) pfun(d, 
                  comm, comp, i, contrp))
            }
            
            ## calculate the p-value from the permutations
            p <- (rowSums(apply(perm.contr, 2, function(x) x >= 
                average - EPS)) + 1)/(nperm + 1)
        }
        else {
            p <- NULL ## a p-value cannot be estimated if no permutations are performed
        }
        
        
        ### 3. Prepare the results and output
        
        ## overall dissimilarity = sum of the observed average dissimilarities between groups
        overall <- sum(average)
        
        ## standard deviation and ratio of each species' contribution  
        sdi <- apply(contr, 2, sd)
        ratio <- average/sdi
        
        ## abundance of each species in each of the groups 
        ava <- colMeans(group.a)
        avb <- colMeans(group.b)
        
        ## sort the average contribution in decreasing order and calculate the cumulative contribution
        ord <- order(average, decreasing = TRUE)
        cusum <- cumsum(average[ord]/overall)
        
        ## put these values into an output list with item names = names of the pair of groups being compared
        out <- list(species = colnames(comm), average = average, 
            overall = overall, sd = sdi, ratio = ratio, ava = ava, 
            avb = avb, ord = ord, cusum = cusum, p = p)
        
        outlist[[paste(comp[i, 1], "_", comp[i, 2], sep = "")]] <- out
    }
    
    ## stop the parallel cluster, if one exists 
    if (isParal && !isMulticore && !hasClus) 
        stopCluster(parallel)
    
    ## set the attributes and class of the output, and done! 
    attr(outlist, "permutations") <- nperm
    attr(outlist, "control") <- attr(perm, "control")
    class(outlist) <- "simper"
    outlist

}

