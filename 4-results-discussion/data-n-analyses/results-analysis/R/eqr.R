eqr <- function (index.values, ref.high, ref.bad) {
    # Calculates Ecological Quality Ratio (EQR) values of a biotic index, given 
    # reference high and bad values. EQR is defined as the ratio of the observed 
    # value versus the value of the same metric under the reference conditions.
    
    # make empty vector to hold the EQRs
    eqr <- c(rep(NA, length(index.values))) 
    
    # calculate EQR
    for (i in 1:length(index.values)) {
      eqr[i] <- round((index.values[i] - ref.bad) / (ref.high - ref.bad), 3)
    }

    return(eqr)
}
