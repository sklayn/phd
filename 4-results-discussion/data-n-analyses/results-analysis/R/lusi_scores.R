lusi_scores <- function(land_use) {
  ## calculates Land Use Simplified Index (LUSI) (Flo et al. 2011; Romero, 2011) 
  ## scores based on land use areas in the watersheds.
  ## Arguments: land_use - tibble of areas (in ha?) for each type of land use in the 
  ##              watersheds (within the chosen buffer; here - 3000). For the LUSI 
  ##              score only the urban, (irrigated) agricultural land, and industrial
  ##              land are considered, as representative of the main types of anthropogenic 
  ##              pressure on coastal ecosystems. However, the rest of the land
  ##              uses are also needed in order to calculate the correct % use in the 
  ##              watershed. 
  ## Output: tibble of urban, agricultural, and industrial scores for each watershed.
  ## Dependencies: tidyverse
  
  library(tidyverse)
  
  ## calculate % land use in the watersheds
  land_use_props <- land.use %>%
    mutate(total_area = rowSums(select(., -watershed))) %>% 
    mutate_if(is.numeric, funs((./total_area) * 100))

  ## calculate the scores corresponding to each pressure considered in the index: 
  ## score assignment as in Flo et al. 2011
  land_use_scores <- land_use_props %>% 
    # urban
    mutate(urban_score = case_when(urban <= 33 ~ 1, 
                                   urban > 33 & urban <= 66 ~ 2, 
                                   TRUE ~ 3)) %>% 
    # agricultural
    mutate(agric_score = case_when(agricultural <= 10 ~ 0, 
                                   agricultural > 10 & agricultural <= 40 ~ 1,
                                   TRUE ~ 2)) %>% 
    # industrial
    mutate(ind_score = case_when(industrial <= 10 ~ 0,
                                 TRUE ~ 1))
  
  ## return the LUSI scores by watershed
  return(land_use_scores %>% select(watershed, contains("score")))
}
