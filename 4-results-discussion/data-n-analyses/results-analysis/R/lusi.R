## Land Use Simplified Index (LUSI) - with corrections for more sources of pressure


### HAVE TO INCLUDE ALL LAND USES! OTHERWISE WRONG % HENCE WRONG SCORE!

lusi <- function(land.use.areas = c(urban, agric, indust), 
                 other.pressures.scores = c(rivers, harbors, adj.wb), 
                 correction.coast) {
  ## calculates a modified version of the Land Use Simplified Index (LUSI) (Flo 
  ## et al. 2011, Romero, 2011) with corrections for more types of coastal pressures. 
  ## Arguments: land.use.area - vector of areas (in ?) corresponding to each 
  ##              type of land use for the considered water body. For LUSI only the 
  ##              urban, (irrigated) agricultural land, and industrial land are 
  ##              considered, as corresponding to the main types of anthropogenic 
  ##              pressure on coastal ecosystems. Values should be in this order!
  ##            other.pressures.scores - vector of correction scores for other significant 
  ##              pressures (rivers, harbors and tourism, influence of adjacent water 
  ##              bodies) on the considered water body to be included as well.  
  ##            correction.coast - correction factor for the shape of the coastline, 
  ##              which influences water residency times and thus the amount of 
  ##              pressure exerted on coastal ecosystems.
  ## Output: LUSI score (numeric).
  
  # convert the land use areas to %
  land.use.props <- (land.use.areas / sum(land.use.areas)) * 100
  
  ## calculate the scores corresponding to each one, as in Flo et al. 2011
  # urban
  if(land.use.props$urban < 33) {
    urban.score = 0
  } else if(land.use.props$urban > 33 & land.use.props$urban < 66) {
    urban.score = 1
  } else {
    urban.score = 2
  }
  
  # agricultural
  if(land.use.props$agricultural < 10) {
    agric.score = 0
  } else if(land.use.props$agricultural > 10 & land.use.props$agricultural < 40) {
    agric.score = 1
  } else {
    agric.score = 2
  }
  
  # industrial
  if(land.use.props$industrial < 10) {
    indust.score = 0
  } else {
    indust.score = 1
  }
  
  # calculate LUSI using these scores and the additional scores for the other pressures
  lusi <- (urban.score + agric.score + indust.score + sum(other.pressures.scores)) * correction.coast

  return(lusi)

}