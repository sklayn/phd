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



plot_pca_diagnostic <- function(pca.res, file.name) {
  ## quick PCA diagnostic plots in one pdf file. Crap with many hard-coded 
  ## choices and values + calls other custom functions. Only acceptable to avoid 
  ## retyping, has to be modified if at all useful for later work. 
  
  pdf(file.path(figs.dir, file.name), 
      paper = "a4r", width = 12, height = 12,
      useDingbats = FALSE)
  
  # NB - all the print statements around the plots - because ggplot objects; 
  # otherwise not written in the file
  
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
