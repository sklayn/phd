## Land Use Simplified Index (LUSI) - with corrections for additional sources of pressure

lusi <- function(land.use.areas, 
                 other.pressures.scores, 
                 correction.coast) {
  ## calculates a modified version of the Land Use Simplified Index (LUSI) (Flo 
  ## et al. 2011; Romero, 2011) with corrections for more types of coastal pressures. 
  ## Arguments: land.use.area - data frame of areas (in ?) corresponding to each 
  ##              type of land use for the considered water body. For the LUSI 
  ##              score only the urban, (irrigated) agricultural land, and industrial
  ##              land are considered, as representing the main types of anthropogenic 
  ##              pressure on coastal ecosystems. However, the rest of the land
  ##              uses are also needed in order to calculate the correct % use. 
  ##            other.pressures.scores - data frame of correction scores for other 
  ##              significant pressures (rivers/point sources, harbors and tourism, 
  ##              influence of adjacent water bodies) on the considered water body 
  ##              to be included as well.   
  ##            correction.coast - correction factor for the shape of the coastline, 
  ##              which influences water residency times and thus the amount of 
  ##              pressure exerted on coastal ecosystems.
  ##            NB Row order should match in all 3 data frames!  
  ## Output: LUSI score (numeric).
  
  # convert the land use areas to %
  land.use.props <- (land.use.areas / rowSums(land.use.areas)) * 100
  
  ## calculate the scores corresponding to each pressure considered in the index: 
  ## score assignment as in Flo et al. 2011
  # urban
  urban.score <- apply(land.use.props, 1, 
                       function(x) {
                          if(x["urban"] < 33) {
                            urban.score = 1
                          } else if(x["urban"] > 33 & x["urban"] < 66) {
                            urban.score = 2
                          } else {
                            urban.score = 3
                          }
                         
                          return(urban.score)
                       })
  
  # agricultural
  agric.score <- apply(land.use.props, 1, 
                       function(x) {
                         if(x["agriculture"] < 10) {
                           agric.score = 0
                         } else if(x["agriculture"] > 10 & x["agriculture"] < 40) {
                           agric.score = 1
                         } else {
                           agric.score = 2
                         }
                         
                         return(agric.score)
                       })
  
  
  # industrial
  indust.score <- apply(land.use.props, 1, 
                       function(x) {
                         if(x["industrial"] < 10) {
                           indust.score = 0
                         } else {
                           indust.score = 1
                         }
                         
                         return(indust.score)
                       })
  
  
  # calculate LUSI using these scores and the additional scores for the other pressures
  lusi <- (urban.score + agric.score + indust.score + rowSums(other.pressures.scores)) * correction.coast
  names(lusi) <- "LUSI"
  
  return(lusi)

}