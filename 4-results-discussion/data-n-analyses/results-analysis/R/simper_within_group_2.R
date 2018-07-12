simper_within_group <- function (comm, group, ...) {
    ### Homebrewed rough companion to the vegan function simper. Calculates the 
    ### contribution of each species to the overall similarity within a group. 
    ### The average contribution of the ith species could be defined by taking 
    ### the average (Si.av), over all pairs of samples within a group, of the ith 
    ### term in the Bray-Curtis similarity equation. A species is considered 
    ### typical of a group if it is found at a consistent abundance throughout 
    ### this group, i.e. the standard deviation of its contribution is low, and 
    ### therefore the ratio Si.av/SD(Si) is high.
    ### This function gives the same result as the function SIMPER in PRIMER-E.
    ###
    ### Arguments: comm - community abundance matrix (untransformed), of the form 
    ###                   sites x species (vegan style).
    ###            group - character vector of the groups.
    ### Output: list of tibbles containing the similarity breakdown for each 
    ###         group. Each tibble contains: species in descending order by their 
    ###         contribution to the group similarity, their average 
    ###         (untransformed) abundance within the group, Si.av, Si.av/SD(Si), 
    ###         the % contribution of each species to the overall group 
    ###         similarity, and the cumulative %. The list also contains the 
    ###         group names and the overall group similarity (as attributes).
    ### Dependencies: tidyverse, plyr

    library(tidyverse)
    
    ## check for empty sites - all-0 species abundances (should not be included 
    ## in the analysis!) 
    if (any(rowSums(comm, na.rm = TRUE) == 0)) 
        warning("you have empty rows: results may be meaningless")
    
    ## transform input abundances to matrix just in case and get the double 
    ## square root of the abundances. FROM NOW ON THE CALCULATIONS WILL BE 
    ## PERFORMED ON THE TRANSFORMED ABUNDANCES! 
    comm <- as.matrix(comm)
    comm.tr <- sqrt(sqrt(comm))
    
    
    ### 1. Calculate the observed contribution of species to the Bray-Curtis 
    ###    similarity within groups 
    
    ## make named list to hold the results for each group
    observed.contr <- vector("list", length = length(unique(as.character(group))))
    names(observed.contr) <- unique(as.character(group))
    
    ## for each group,
    for (i in unique(as.character(group))) {
        
        ## subset the (transformed) community matrix to get only the current group
        current.comm <- comm.tr[group == i, , drop = FALSE]
        
        ## make a matrix of all pairwise combinations of samples within the group
        pair.comb <- t(combn(nrow(current.comm), 2))

        ## matrix for the species contribution to the similarity between pairs of
        ## samples within a group  
        contr.sp <- matrix(ncol = ncol(current.comm), nrow = nrow(pair.comb))
        
        ## row names of contribution matrix = sample combination being compared; 
        ## column names = species
        rownames(contr.sp) <- paste(pair.comb[, 1], pair.comb[, 2], sep = "_")
        colnames(contr.sp) <- colnames(current.comm) 
        
        for (j in 1:nrow(pair.comb)){
            ## get the abundance matrix for the current pair of samples 
            group.a <- current.comm[pair.comb[j, 1], , drop = FALSE]
            group.b <- current.comm[pair.comb[j, 2], , drop = FALSE]
            
            ## calculate the numerator and denominator of the Bray-Curtis 
            ## similarity equation
            md <- pmin(group.a, group.b)
            me <- group.a + group.b

            ## calculate the contribution of each species to the similarity 
            ## between pairs of samples
            contr.sp[j, ] <- 200 * md/sum(me) 

        }
        
        ## fill in results list
        observed.contr[[i]] <- contr.sp 
    }
    
    
    ### 2. Calculate all derived results for output 
    
    ## observed average contribution of each species to the group similarity
    Si.average <- lapply(observed.contr, colMeans) ## this is now a list of vectors - 1 vector per group, with 1 value per species
    
    ## overall group similarity = the mean of all pairwise simiarities, which are 
    ## themselves the sums of the species contributions (= row sums)
    S.group.overall <- lapply(observed.contr, function(x) mean((rowSums(x))))
    
    ## standard deviation and ratio of each species' contribution  
    sdi <- lapply(observed.contr, function(x) apply(x, 2, sd))
    ratio <- mapply(function(x, y) x/y, Si.average, sdi, SIMPLIFY = FALSE)
    
    ## UNTRANSFORMED abundance of each species in each of the groups 
    av.abnd <- lapply(unique(as.character(group)), 
                      function(x) colMeans(comm[group == x, , drop = FALSE]))
    names(av.abnd) <- unique(as.character(group))
    
    ## sort the average contribution in descending order - this will be the order
    ## for all the output variables!
    ## NB IT'S VERY IMPORTANT TO MAKE SURE THAT ALL VARIABLES ARE IN THE SAME 
    ## ORDER; OTHERWISE THE OUTPUT MIGHT BE SCRAMBLED... 
    ord <- lapply(Si.average, order, decreasing = TRUE) 
    Si.average.ordered <- mapply(function(x, y) x[y], 
                                 Si.average, 
                                 ord, 
                                 SIMPLIFY = FALSE)
    
    ## apply the SAME order to the ratio of the average similarity and its 
    ## standard deviation and the average (untransformed) abundance
    ratio.ordered <- mapply(function(x, y) x[y], ratio, ord, SIMPLIFY = FALSE)
    av.abnd.ordered <- mapply(function(x, y) x[y], av.abnd, ord, SIMPLIFY = FALSE)

    ## % contribution and cumulative % contribution of each species    
    aver.contrib.percent <- mapply(function(x, y) 100 * x/y, 
                                   Si.average.ordered, 
                                   S.group.overall, 
                                   SIMPLIFY = FALSE)
    cumul.contrib.percent <- lapply(aver.contrib.percent, cumsum)
    
    
    ### 3. Format output 
    
    ## make a large list of everything we want to output 
    all.list <- list("av.abnd.ordered" = av.abnd.ordered, 
                     "Si.average.ordered" = Si.average.ordered, 
                     "ratio.ordered" = ratio.ordered, 
                     "aver.contrib.percent" = aver.contrib.percent, 
                     "cumul.contrib.percent" = cumul.contrib.percent)
    
    ## make output list to hold final well-formatted output
    outlist <- vector("list", length = length(unique(as.character(group))))
    names(outlist) <- unique(as.character(group))
    
    for (l in unique(as.character(group))) {
        ## subset each group from the incredibly clumsy list of lists, 
        ## and collapse this group sublist to a tibble
        group.out.df <- plyr::ldply(lapply(all.list, function(x) x[[l]])) %>% 
            ## this results in a tibble with rows = variables (= each of the 
            ## sublists) and cols = species. Make it a long tibble.
            gather(key = species, value = var, -.id) %>% 
            ## spread so that each variable is in its own column
            spread(key = .id, value = var) %>% 
            ## sort on descending average contribution to the similarity
            arrange(desc(Si.average.ordered)) %>% 
            ## add the group name to the tibble %>% 
            mutate(group = l) %>%
            ## rename and rearrange columns
            rename(av.abnd = av.abnd.ordered, 
                   Si.av = Si.average.ordered, 
                   `Si.av/SD` = ratio.ordered, 
                   `%contrib` = aver.contrib.percent, 
                   `%cum` = cumul.contrib.percent) %>% 
            select(group, species, av.abnd, Si.av, `Si.av/SD`, `%contrib`, `%cum`)
        
        ## add formatted results tibble to output list
        outlist[[l]] <- group.out.df
    }
    
    ## add the overall average group similarity as an attribute to the list
    attr(outlist, "overall.group.similarity") <- unlist(S.group.overall)
    
    ## voila, return result! 
    return(outlist)
}

