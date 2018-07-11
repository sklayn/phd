simper_within_group_2 <- function (comm, group, ...) {

    ## check for empty sites - all-0 species abundances (should not be included in the analysis!) 
    if (any(rowSums(comm, na.rm = TRUE) == 0)) 
        warning("you have empty rows: results may be meaningless")
    
    ## transform input abundances to matrix & get the double square root of the abundances. FROM NOW ON THE CALCULATIONS WILL BE PERFORMED ON THE TRANSFORMED ABUNDANCES! 
    comm <- as.matrix(comm)
    comm.tr <- sqrt(sqrt(comm))
    
    ## make list to hold SIMPER results
    outlist <- NULL
    
    
    ### 1. Calculate the observed contribution of species to the Bray-Curtis similarity within groups 
    ## make named lists to hold the results for each group
    observed.contr <- vector("list", length = length(unique(as.character(group))))
    names(observed.contr) <- unique(as.character(group))
    
    overall.group.similarity <- vector("list", length = length(unique(as.character(group))))
    names(overall.group.similarity) <- unique(as.character(group))
    
    for (i in unique(as.character(group))) {
        
        ## subset the (transformed) community matrix to get only the current group
        current.comm <- comm.tr[group == i, , drop = FALSE]
        
        ## make a matrix of all pairwise combinations of samples within the group
        pair.comb <- t(combn(nrow(current.comm), 2))
        
        ## make the row and column names of the contribution matrices 
        contr.rows <- paste(pair.comb[, 1], pair.comb[, 2], sep = "_")
        contr.cols <- colnames(current.comm)
        
        ## matrix to hold the numerator of the species contributions to each of the combinations of observations in the groups
        contr.md <- matrix(ncol = ncol(current.comm), nrow = nrow(pair.comb))
        rownames(contr.md) <- contr.rows
        colnames(contr.md) <- contr.cols
        
        ## matrix to hold the denominator of the species contributions to each of the combinations of observations in the groups
        contr.me <- matrix(ncol = ncol(current.comm), nrow = nrow(pair.comb))
        rownames(contr.me) <- contr.rows
        colnames(contr.me) <- contr.cols
        
        ## matrix for the species contribution to the similarity between pairs of samples within a group  
        contr.sp <- matrix(ncol = ncol(current.comm), nrow = nrow(pair.comb))
        rownames(contr.sp) <- contr.rows
        colnames(contr.sp) <- contr.cols
        
        for (j in 1:nrow(pair.comb)){
            ## get the abundance matrix for the current pair of samples 
            group.a <- current.comm[pair.comb[j, 1], , drop = FALSE]
            group.b <- current.comm[pair.comb[j, 2], , drop = FALSE]
            
            for (k in 1:ncol(contr.md)) {
                ## calculate the numerator and denominator of the contribution of each species to the overall Bray-Curtis dissimilarity between groups, and put the result in the matrix of species contributions 
                md <- min(group.a[, k, drop = FALSE], group.b[, k, drop = FALSE])
                me <- group.a[, k, drop = FALSE] + group.b[, k, drop = FALSE]
                
                ## fill in the corresponding matrices 
                contr.md[j, k] <- md
                contr.me[j, k] <- me
                
                ## contribution of each species to the similarity between pairs of samples within a group
                contr.sp[j, k] <- sum(md)/sum(me) 
                
            }
        }
        
        ## fill in lists
        observed.contr[[i]] <- contr.sp ## THIS IS NOT IT 
        overall.group.similarity[[i]] <- 100 * (2 * rowSums(contr.md, na.rm = TRUE)) / (rowSums(contr.me, na.rm = TRUE)) ## THIS, HOWEVER, WORKS, INCREDIBLE AS IT IS...
    }
    
    
    ## calculate the observed average contribution of each species to each combination of observations within a group
    Si.average <- lapply(observed.contr, colMeans) ## this is now a list of vectors - 1 vector per group, with 140 values in each - 1 for each species
    
    
    ### 2. Prepare the results and output
    
    ## overall similarity of a group = sum of the observed average similarities within groups
    Si.group.overall <- lapply(observed.contr, function(x) 200 * rowSums(x, na.rm = TRUE))
    
    ## standard deviation and ratio of each species' contribution  
    sdi <- lapply(observed.contr, function(x) apply(x, 2, sd))
    ratio <- mapply(function(x, y) x/y, Si.average, sdi, SIMPLIFY = FALSE) # MAYBE SIMPLIFY = TRUE? DEPENDS.. 
    
    ## UNTRANSFORMED abundance of each species in each of the groups 
    av.abnd <- lapply(unique(as.character(group)), function(x) colMeans(comm[group == x, , drop = FALSE]))
    names(av.abnd) <- unique(as.character(group))
    
    ## sort the average contribution in decreasing order and calculate the % contribution
    ord <- lapply(Si.average, order, decreasing = TRUE) 
    Si.average.ordered <- mapply(function(x, y) x[y], Si.average, ord, SIMPLIFY = FALSE)
    aver.prop <- mapply(function(x, y) x/y, Si.average.ordered, Si.group.overall, SIMPLIFY = FALSE)
    
    cumul.percent <- lapply(aver.prop, cumsum)
    
    
    ## put these values into an output list with item names = names of the pair of groups being compared
    
    
    
    outlist
    
}



