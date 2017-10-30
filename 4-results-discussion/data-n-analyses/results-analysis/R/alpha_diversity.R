alpha_diversity <- function(zoo.data) {
  ## uses vegan package functions to calculate species richness S, Shannon H', 
  ## Pielou J', and effective Shannon diversity (the Hill number), and returns 
  ## them in a tibble.
  ## Arguments: zoo.data - numeric abundance data; samples x species (vegan format)
  ## Dependencies: vegan, dplyr
  
  # number of species per sample/replicate
  sp_richness = specnumber(zoo.data)
  
  # Shannon-Wiener biodiversity index
  shannon_h = diversity(zoo.data) 
 
  # make a tibble to return
  diversity_sand <- tibble(sp_richness = sp_richness, 
                           shannon_h = shannon_h, 
                           evenness_j = shannon_h / log(sp_richness), # Pielou evenness
                           h_effective = exp(shannon_h)) # effective diversity
  
  return(diversity_sand)
}