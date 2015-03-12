# return_est.R

# Function: return_est
# Input: Risk of outcome in exposed (EY1)
#        Risk of outcome in unexposed (EY0)
#        Desired metric (metric)
# Output: Estimate (est)

return_est <- function(EY1, EY0, metric) {
  
  # Risk difference
  if (metric == "rd") {
    est  <- EY1 - EY0
    return(est)
    
    # Relative risk  
  } else if (metric == "rr") {
    est <- EY1/EY0
    return(est)
    
    # Odds ratio  
  } else if (metric == "or") {
    est <- (EY1/(1-EY1)) / (EY0/(1-EY0))
    return(est)
  }
}


