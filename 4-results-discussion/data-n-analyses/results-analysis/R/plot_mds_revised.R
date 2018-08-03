plot_mds <- function(mds.obj, groups = NULL) {
    ## Plots MDS result in ggplot2. 
    ## Arguments: mds.obj - MDS result object.
    ##            groups - optional factor for coloring & labeling the points 
    ##              on the plot; needs to match in length and order the rows of 
    ##              the original MDS community dataset.
    ## Returns: ggplot object.
    ## Dependencies: vegan, tidyverse 

    ## extract the MDS scores and create tibble for plotting. NB in this case 
    ## the rownames are not preserved as a variable - no need since if desired, 
    ## labels/groups are to be supplied in the function call.  
    scrs <- as_data_frame(scores(mds.obj, display = "sites"))
    
    if(!is.null(groups)) {
        ## if there are groups, add them to the tibble
        scrs <- scrs %>% 
            mutate(group = groups) 
        
        ## plot the ordination in ggplot with points coloured by group
        p <- ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = group)) + 
            geom_point(size = 2, position = "identity")
        
    } else {
        ## plot the ordination in ggplot
        p <- ggplot(scrs, aes(x = NMDS1, y = NMDS2)) + 
            geom_point(size = 2, position = "identity")
    }

    return(p)
}



extract_envfit_scores <- function(envfit.obj, pval = 0.05, r2 = FALSE) {
    ## Helper function to extract NMDS vector scores from envfit objects 
    ## according to the specified p-value (p < 0.05 by default). Optionally 
    ## also extracts the r2 value.
    ## Arguments: envfit.obj - object containing the results of an envfit call.
    ##            pval - significance (p-value) by which to limit the returned 
    ##             factors.
    ##            r2 - correlation of the factors to the ordination.
    ## Returns: tibble of NMDS coordinates of vectors for the fitted factors 
    ##          that are below the specified significance level, their p-value 
    ##          and (optionally) r2 value.
    ## Dependencies: vegan, tidyverse 
    
    ## extract the NMDS scores for the fitted environmental parameters, keeping 
    ## the row names (= environmental parameters)
    vector.scrs <- rownames_to_column(as.data.frame(envfit.obj$vectors$arrows), 
                                      var = "param")

    ## if requested, extract the correlation coefficient (r2) of the environmental
    ## parameters with the ordination
    if(r2) {
        vector.scrs <- vector.scrs %>% 
            mutate(r2 = envfit.obj$vectors$r)
    }
    
    ## get the p-values and subset the tibble according to the required threshold
    vector.scrs <- vector.scrs %>% 
        mutate(p = envfit.obj$vectors$pvals) %>% 
        filter(p < pval)

    return(vector.scrs)
}



plot_envfit <- function(mds.plot, envfit.scrs, param.labels, label.col) {
    ## Plots fitted environmental factors of envfit as vectors over an 
    ## existing MDS plot (made with ggplot).
    ## 
    ## Arguments: mds.plot - MDS plot made with ggplot.
    ##            envfit.scrs - tibble with corresponding envfit results. 
    ##            param.labels - name of the column containing the parameter 
    ##             names in the envfit scores tibble.
    ##            label.col - colour to use for the parameter vectors and labels. 
    ## Dependencies: tidyverse, ggrepel 
    ## Returns: ggplot object.

    var.labels <- enquo(param.labels)
    
    p <- mds.plot + 
        
        ## need aspect ratio of 1 - to make the length of the vectors 
        ## proportional to their correlation with the ordination!
        coord_fixed() +    
        
        ## plot the vectors over the ordination. Each starts at (0, 0) and ends 
        ## at (NMDS1, NMDS2)
        geom_segment(data = envfit.scrs,
                     aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                     arrow = arrow(length = unit(0.25, "cm")), 
                     colour = label.col, 
                     inherit.aes = FALSE) +
        
        ## label the vectors 
        ggrepel::geom_text_repel(data = envfit.scrs, 
                                 aes(x = NMDS1, y = NMDS2, label = !!var.labels), 
                                 colour = label.col,
                                 inherit.aes = FALSE) 
    
    return(p)
}