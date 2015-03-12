# analyze_sub.R

# Runs the particular analyses, depending on the type of analysis, within the analyze function
# To be called separately for cases with and without subclassification

# Function: analyze_sub
# Input: an_data (format as pop),
#        subclass (numeric vector, same # obs as an_data),
#        match_weights (numeric vector, same # obs as an_data),
#        samp_weights (numeric vector, same # obs as an_data),
#        match_estimand ("ate" or "att"),
#        a single row of analyses (as a list)
# Output: est (ATE or ATT estimate, a numeric scalar)

analyze_sub <- function(an_data, subclass, match_weights, samp_weights, match_estimand, analyses, g, q) {
  
  # Unadjusted----------------------------------------
  if (analyses$method == "unadj")   { 
    est <- an_unadj(an_data = an_data, match_weights = match_weights, samp_weights = samp_weights, metric = analyses$metric)
    
  # G computation-------------------------------------
  } else if (analyses$method == "gcomp")   { 
    q.covs <- as.character(lapply(q$covs, paste, collapse=","))
    q.sllibs <- as.character(lapply(q$sllibs, paste, collapse=","))
    q <- q$q[q$method == "gcomp" & q.covs == paste(analyses$covs[[1]], collapse=",") & q$q_para == analyses$q_para 
             & q$q_dist == analyses$q_dist
             & q$q_link == analyses$q_link & q$q_covform == analyses$q_covform
             & q$q_slfull == analyses$q_slfull & q$q_sldefault == analyses$q_sldefault
             & q.sllibs == paste(analyses$q_sllibs[[1]], collapse=",")]
    est <- an_gcomp(A = an_data$A_1, q = q, estimand = analyses$estimand, metric = analyses$metric)
  
  # Propensity weighting -----------------------------
  } else if (analyses$method == "pweight") { 
    g.covs <- as.character(lapply(g$covs, paste, collapse=","))
    g.sllibs <- as.character(lapply(g$sllibs, paste, collapse=","))
    g <- g$g[g$method == "pweight" & g.covs == paste(analyses$covs[[1]], collapse = ",") 
             & g$pscore_para == analyses$pscore_para 
             & g$pscore_link == analyses$pscore_link & g$pscore_covform == analyses$pscore_covform
             & g$pscore_slfull == analyses$pscore_slfull & g$pscore_sldefault == analyses$pscore_sldefault
             & g.sllibs == paste(analyses$pscore_sllibs[[1]], collapse = ",")]
    est <- an_pweight(an_data = an_data, match_weights = match_weights, samp_weights = samp_weights, estimand = analyses$estimand,
                      match_estimand = match_estimand, 
                      metric = analyses$metric, pweight_type = analyses$pweight_type, g = g)
  
  # TMLE ---------------------------------------------
  } else if (analyses$method == "tmle")    { 
    g.covs <- as.character(lapply(g$covs, paste, collapse=","))
    g.sllibs <- as.character(lapply(g$sllibs, paste, collapse=","))
    q.covs <- as.character(lapply(q$covs, paste, collapse=","))
    q.sllibs <- as.character(lapply(q$sllibs, paste, collapse=","))
    g <- g$g[g$method == "tmle" & g.covs == paste(analyses$covs[[1]], collapse = ",") 
             & g$pscore_para == analyses$pscore_para 
             & g$pscore_link == analyses$pscore_link & g$pscore_covform == analyses$pscore_covform
             & g$pscore_slfull == analyses$pscore_slfull & g$pscore_sldefault == analyses$pscore_sldefault
             & g.sllibs == paste(analyses$pscore_sllibs[[1]], collapse = ",")]
    q <- q$q[q$method == "tmle" & q.covs == paste(analyses$covs[[1]], collapse=",") & q$q_para == analyses$q_para 
             & q$q_dist == analyses$q_dist
             & q$q_link == analyses$q_link & q$q_covform == analyses$q_covform
             & q$q_slfull == analyses$q_slfull & q$q_sldefault == analyses$q_sldefault
             & q.sllibs == paste(analyses$q_sllibs[[1]], collapse=",")]
    est <- an_tmle(an_data = an_data, match_weights = match_weights, samp_weights = samp_weights, covs = analyses$covs[[1]], 
                   estimand = analyses$estimand,
                   match_estimand = match_estimand, 
                   metric = analyses$metric, g = g, q = q) 
  } 
  return(est)
}



