# q_models.R

# Function: q_models
# Inputs: an_data (format same as pop),
#         subclass (numeric vector, same # obs as an_data),
#         match_weights (numeric vector, same # obs as an_data),
#         samp_weights (numeric vector, same # obs as an_data),
#         analyses (data frame of all analysis specifications )
# Output: q (data fame of all formulated estimates of the outcome model-- n*3 rows per method-- 1st n obs are predictions for observed data, 2nd n obs are predictions when A_1=1 for everyone and 3rd n obs are predictions when A_1=0 for everyone)
#          
# This function estimates the probability of the outcome given treatment and covariates

q_models <- function(an_data, subclass, match_weights, samp_weights, analyses) {
  
  # Collapse analyses to all unique combinations of q estimates, for analyses that require them
  runs <- analyses[,c("method","covs","q_para","q_dist","q_link","q_covform","q_slfull","q_sldefault","q_sllibs")]
  runs <- runs[runs$method %in% c("gcomp","tmle")]
  runs <- runs[!duplicated(runs),]
  
  # If there are weights, and the method is gcomp or tmle, then the q model estimation should include those weights.

  # For datasets that were not pre-processed with subclassification:
  if (is.na(subclass[1])) {
    for (i in 1:nrow(runs)) {
      temp_run <- run[i,]
  
      temp_q <- est_q(an_data = an_data, match_weights = match_weights, samp_weights = samp_weights, 
                      covs = temp_run$covs[[1]], 
                      para = temp_run$q_para, dist = temp_run$q_dist, link = temp_run$q_link, 
                      covform = temp_run$q_covform, slfull = temp_run$q_slfull, 
                      sldefault = temp_run$q_sldefault, sllibs = temp_run$q_sllibs[[1]])
      
      q_df <- data.frame(temp_q, temp_run, row.names=NULL) 
      if (i ==1) { q <- q_df
      } else { q <- rbind(q, q_df)  }
    }
    q$subclass <- NA
  
  # For datasets pre-processed with subclassification-- need to estimate q within each subclass
  } else {
    nsubclass <- length(unique(subclass))
    
    for (i in 1:nrow(runs)) {
      temp_run <- run[i,]
    
      # Loop through subclasses and generate estimates
      for (s in 1:nsubclass) {
        sub_data <- an_data[subclass == s,]
        sub_match_weights <- match_weights[subclass == s]
        sub_samp_weights <- samp_weights[subclass == s]
        temp_est <- est_q(an_data = sub_data, match_weights = sub_match_weights, samp_weights = sub_samp_weights,
                          covs = temp_run$covs[[1]], 
                          para = temp_run$q_para, dist = temp_run$q_dist, link = temp_run$q_link, 
                          covform = temp_run$q_covform, slfull = temp_run$q_slfull, 
                          sldefault = temp_run$q_sldefault, sllibs = temp_run$q_sllibs[[1]])
        temp_est$subclass <- s
        if (s ==1 ) { temp_q <- temp_est } else { temp_q <- c(temp_q, temp_est) }
      }
      
      q_df <- data.frame(temp_q, temp_run, row.names=NULL) 
      if (i ==1) { q <- q_df
      } else { q <- rbind(q, q_df)  } 
    }
  }
  
  names(q)[1] <- "q"
  return(q)
}

