# an_unadj.R

# Function: an_unadj
# Input: an_data (format as pop)
#        match_weights (numeric vector, same # obs as an_data)
#        samp_weights (numeric vector, same # obs as an_data)
#        metric (string, "rr", "rd"  or "rr")
# Output: est (numeric scalar)

an_unadj <- function(an_data, match_weights, samp_weights, metric) {
  
  # Combine weights. This will be a vector of 1's if no weighting is necessary.
  weights <- match_weights * samp_weights
  
  EY1 <- weighted.mean(an_data$outcome[an_data$A_1==1], weights[an_data$A_1==1])
  EY0 <- weighted.mean(an_data$outcome[an_data$A_1==0], weights[an_data$A_1==0])

  return(return_est(EY1, EY0, metric))
}