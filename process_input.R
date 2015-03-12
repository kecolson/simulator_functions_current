# Process user inputs

# Test input
# input <- list(sample_size=500, all_exposed=F, all_controls=F, ate=T, att=T, effect_metric=1, outcome_diff=0,
#               match=1, match_mode="basic", match_basic_covs=c("W1","W2"),
#               an_gcomp=T, an_pweight=F, an_aiptw=F, an_tmle=T, an_covs=c("W1","W2"),an_pweight_type=NA,
#               analysis_pscore_mode="basic",
#               analysis_q_mode="basic")

process_input <- function() {
  
  dfSLlibs <- c("SL.glm", "SL.step", "SL.rpartPrune", "SL.earth", "SL.nnet", "SL.gam")
  
  input <- list(sample_size = input$sample_size,
                all_exposed = input$all_exposed,
                all_controls = input$all_controls,
                ate = input$ate,
                att = input$att,
                effect_metric = input$effect_metric,
                outcome_diff = input$outcome_diff,
                match = input$match,
                match_mode = input$match_mode,
                match_basic_more = input$match_basic_more,
                match_basic_covs = input$match_basic_covs,
                match_exact = input$match_exact,
                match_nonexact = input$match_nonexact,
                match_exact_covs = input$match_exact_covs,
                match_nonexact_covs = input$match_nonexact_covs,
                match_distance_pscore = input$match_distance_pscore,
                match_distance_mahal = input$match_distance_mahal,
                match_distance_class = input$match_distance_class,
                match_distance_nnet = input$match_distance_nnet,
                match_method_nn = input$match_method_nn,
                match_method_opt = input$match_method_opt,
                match_method_gen = input$match_method_gen,
                match_method_sub = input$match_method_sub,
                match_method_full = input$match_method_full,
                match_nn_ratio = input$match_nn_ratio,
                match_nn_replace = input$match_nn_replace,
                match_nn_exact_covs = input$match_nn_exact_covs,
                match_gen_ratio = input$match_gen_ratio,
                match_sub_pscore_subclasses = input$match_sub_pscore_subclasses,
                match_discard = input$match_discard,
                match_discard_reest = input$match_discard_reest,
                match_pscore_para = input$match_pscore_para,
                match_pscore_sl = input$match_pscore_sl,
                match_pscore_link = input$match_pscore_link,
                match_pscore_backtrans = input$match_pscore_backtrans,
                match_pscore_covform = input$match_pscore_covform,
                match_pscore_slfull = input$match_pscore_slfull,
                match_pscore_sldefault = input$match_pscore_sldefault,
                match_pscore_sllibs = input$match_pscore_sllibs,
                an_gcomp  = input$an_gcomp ,
                an_pweight = input$an_pweight,
                an_aiptw = input$an_aiptw,
                an_tmle = input$an_tmle,
                an_bcm = input$an_bcm,
                an_covs = input$an_covs,
                an_pweight_type = input$an_pweight_type,
                analysis_pscore_mode = input$analysis_pscore_mode,
                analysis_pscore_para = input$analysis_pscore_para,
                analysis_pscore_sl = input$analysis_pscore_sl,
                analysis_pscore_link = input$analysis_pscore_link,
                analysis_pscore_covform = input$analysis_pscore_covform,
                analysis_pscore_slfull = input$analysis_pscore_slfull,
                analysis_pscore_sldefault = input$analysis_pscore_sldefault,
                analysis_pscore_sllibs = input$analysis_pscore_sllibs,
                analysis_q_mode = input$analysis_q_mode,
                analysis_q_para = input$analysis_q_para,
                analysis_q_sl = input$analysis_q_sl,
                analysis_q_dist_gaus = input$analysis_q_dist_gaus,
                analysis_q_dist_bin = input$analysis_q_dist_bin,
                analysis_q_dist_gam = input$analysis_q_dist_gam,
                analysis_q_dist_pois = input$analysis_q_dist_pois,
                analysis_q_gaus_link = input$analysis_q_gaus_link,
                analysis_q_bin_link = input$analysis_q_bin_link,
                analysis_q_gam_link = input$analysis_q_gam_link,
                analysis_q_pois_link = input$analysis_q_pois_link,
                analysis_q_covform = input$analysis_q_covform,
                analysis_q_slfull = input$analysis_q_slfull,
                analysis_q_sldefault = input$analysis_q_sldefault,
                analysis_q_sllibs = input$analysis_q_sllibs)
  
  #   CONVERT ANY NULL'S TO NA'S HERE?
  
  ######################
  # Overall design 
  ######################
  
  designs <- expand.grid(sample_size  = input$sample_size,   
                         all_exposed  = input$all_exposed,                 
                         all_controls = input$all_controls               
  )
  
  ######################
  # Matching
  ######################
  
  estimand <- c("ate","att")[c(input$ate,input$att)]
  
  # 1. No matching. Set all matching variables to NA.
  # Always have a version with no matching
  
  m0 <- expand.grid(estimand = NA, match = 0, match_covs = list(NA), match_exact = NA,
                    match_distance = NA, match_method = NA, match_nonexact_ratio = NA, match_replace = NA, match_firstexact_covs = list(NA), 
                    match_sub_pscore_subclasses = NA,
                    match_discard = NA, match_discard_reest = NA, 
                    match_pscore_para = NA, match_pscore_link = NA, match_pscore_backtrans = NA, match_pscore_covform = NA, 
                    match_pscore_slfull = NA, match_pscore_sldefault = NA, match_pscore_sllibs = list(NA))
  
  if (input$match == 0) {
      matches <- m0
  }
  
  # 2.A. Basic matching. Set default matching methods accordingly.
  
  if (input$match == 1) {
      if (input$match_mode == "basic") {
        
        # Default 1: nearest neighbor 1:1 matching with replacement and parametric propensity score with main terms
        m1 <- expand.grid(estimand = estimand, match = 1, match_covs = list(input$match_basic_covs), match_exact = F, 
                         match_distance = "pscore", match_method = "nn", match_nonexact_ratio = 1, match_replace = T, match_firstexact_covs = list(NA), 
                         match_sub_pscore_subclasses = NA,
                         match_discard = 0, match_discard_reest = NA, 
                         match_pscore_para = T, match_pscore_link = "logit", match_pscore_backtrans = F, match_pscore_covform = "mainterm", 
                         match_pscore_slfull = NA, match_pscore_sldefault = NA, match_pscore_sllibs = list(NA))
        
        # Default 2: nearest neighbor 1:1 matching with replacement and propensity score estimated using default SuperLearner settings
        m2 <- expand.grid(estimand = estimand, match = 1, match_covs = list(input$match_basic_covs), match_exact = F, 
                         match_distance = "pscore", match_method = "nn", match_nonexact_ratio = 1, match_replace = T, match_firstexact_covs = list(NA), 
                         match_sub_pscore_subclasses = NA,
                         match_discard = 0, match_discard_reest = NA, 
                         match_pscore_para = F, match_pscore_link = NA, match_pscore_backtrans = NA, match_pscore_covform = NA, 
                         match_pscore_slfull = 1, match_pscore_sldefault = 1, match_pscore_sllibs = list(dfSLlibs))
        
        # Default 3: full matching with parametric propensity score with main terms
        m3 <- expand.grid(estimand = estimand, match = 1, match_covs = list(input$match_basic_covs), match_exact = F,  
                         match_distance = "pscore", match_method = "full", match_nonexact_ratio = NA, match_replace = NA, match_firstexact_covs = list(NA), 
                         match_sub_pscore_subclasses = NA,
                         match_discard = 0, match_discard_reest = NA, 
                         match_pscore_para = T, match_pscore_link = "logit", match_pscore_backtrans = F, match_pscore_covform = "mainterm", 
                         match_pscore_slfull = NA, match_pscore_sldefault = NA, match_pscore_sllibs = list(NA))
        
        # Default 4: full matching with default SuperLearner propensity score
        m4 <- expand.grid(estimand = estimand, match = 1, match_covs = list(input$match_basic_covs), match_exact = F, 
                         match_distance = "pscore", match_method = "full", match_nonexact_ratio = NA, match_replace = NA, match_firstexact_covs = list(NA), 
                         match_sub_pscore_subclasses = NA,
                         match_discard = 0, match_discard_reest = NA, 
                         match_pscore_para = F, match_pscore_link = NA, match_pscore_backtrans = NA, match_pscore_covform = NA, 
                         match_pscore_slfull = 1, match_pscore_sldefault = 1, match_pscore_sllibs = list(dfSLlibs))
        
        matches <- rbind(m0,m1,m2,m3,m4)
      }
    
  # 2.B. Advanced matching. 
    
      if (input$match_mode == "adv") {
          m1 <- m2 <- NA
        
  # 2.B.i. Advanced matching. Exact. Follow by cases and fill in.
        
          if (input$match_exact == T) {
              m1 <- expand.grid(estimand = estimand, match = 1, match_covs = list(input$match_exact_covs), match_exact = T, 
                                match_distance = NA, match_method = "exact", match_nonexact_ratio = NA, match_replace = NA, match_firstexact_covs = list(NA), 
                                match_sub_pscore_subclasses = NA,
                                match_discard = NA, match_discard_reest = NA, 
                                match_pscore_para = NA, match_pscore_link = NA, match_pscore_backtrans = NA, match_pscore_covform = NA, 
                                match_pscore_slfull = NA, match_pscore_sldefault = NA, match_pscore_sllibs = list(NA))
          }
          
  # 2.B.ii. Advanced matching. Non-exact. Follow by cases and fill in.
          
          if (input$match_nonexact == T) {
            
              match_distance <- c("pscore","mahal","class","nnet")[c(input$match_distance_pscore, input$match_distance_mahal, input$match_distance_class, input$match_distance_nnet)]
              match_method   <- c("nn","opt","gen","sub","full")[c(input$match_method_nn, input$match_method_opt, input$match_method_gen, input$match_method_sub, input$match_method_full)]
              match_pscore_para <- c(T,F)[c(input$match_pscore_para, input$match_pscore_sl)]
              
              m2 <- expand.grid(estimand = estimand, match = 1, match_covs = list(input$match_nonexact_covs), match_exact = F, 
                               match_distance = match_distance, match_method = match_method, match_nonexact_ratio = NA, match_replace = NA, 
                               match_firstexact_covs = list(NA), 
                               match_sub_pscore_subclasses = NA,
                               match_discard = input$match_discard, match_discard_reest = input$match_discard_reest, 
                               match_pscore_para = input$match_pscore_para, match_pscore_link = NA, match_pscore_backtrans = NA, 
                               match_pscore_covform = NA, 
                               match_pscore_slfull = NA, match_pscore_sldefault = NA, match_pscore_sllibs = list(NA))
              
              # Add details of settings for specific non-exact matching methods
              
              if (input$match_method_nn == T) {
                m2$match_nonexact_ratio [m2$match_method == "nn"]  <- input$match_nn_ratio
                m2$match_replace        [m2$match_method == "nn"]  <- input$match_nn_replace
                m2$match_firstexact_covs[m2$match_method == "nn"]  <- list(input$match_nn_exact_covs)
              }
              
              if (input$match_method_gen == T) {
                m2$match_nonexact_ratio[m2$match_method == "gen"] <- input$match_gen_ratio
              }
              
              if (input$match_method_sub == T) {
                m2$match_sub_pscore_subclasses[m2$match_method == "sub"] <- input$match_sub_pscore_subclasses
              }
              
              if (input$match_distance_pscore == T & input$match_pscore_para == T) {
                  m2$match_pscore_link     [m2$match_distance == "pscore" & m2$match_pscore_para == T] <- input$match_pscore_link
                  m2$match_pscore_backtrans[m2$match_distance == "pscore" & m2$match_pscore_para == T] <- input$match_pscore_backtrans
                  m2$match_pscore_covform  [m2$match_distance == "pscore" & m2$match_pscore_para == T] <- input$match_pscore_covform
              }
              
              if (input$match_distance_pscore == T & input$match_pscore_para == F) {
                m2$match_pscore_slfull   [m2$match_distance == "pscore" & m2$match_pscore_para == F] <- input$match_pscore_slfull
                m2$match_pscore_sldefault[m2$match_distance == "pscore" & m2$match_pscore_para == F] <- input$match_pscore_sldefault
                m2$match_discard        [m2$match_distance == "pscore" & m2$match_pscore_para == F] <- NA
                m2$match_discard_reest  [m2$match_distance == "pscore" & m2$match_pscore_para == F] <- NA 
                m2$match_firstexact_covs[m2$match_distance == "pscore" & m2$match_pscore_para == F] <- NA 
                if (input$match_pscore_sldefault == 0) m2$match_pscore_sllibs[m2$match_distance == "pscore" & m2$match_pscore_para == F] <- list(input$match_pscore_sllibs)
                if (input$match_pscore_sldefault == 1) m2$match_pscore_sllibs[m2$match_distance == "pscore" & m2$match_pscore_para == F] <- list(dfSLlibs)
                
              }
          }
        
          if (!is.na(m1) & !is.na(m2)) { matches <- rbind(m0,m1,m2)
          } else if (!is.na(m1) &  is.na(m2)) { matches <- rbind(m0,m1)
          } else if ( is.na(m1) & !is.na(m2)) { matches <- rbind(m0,m2) }
      }
  }
  
  # For now, dropping matching methods aimed at the ATE, since I haven't figured out how to do those yet
  matches <- matches[matches$estimand!="ate",]
  
  
  ######################
  # Analyses
  ######################
  
  # 1.A. Basic mode for treatment mechanism ----------
  
  if (input$analysis_pscore_mode == "basic") {
    
    # (1) Parametric GLM with Binomial distribution and logit link, and covariates as main linear terms with first order interactions
    g1 <- expand.grid(pscore_para = T, pscore_link = "logit", pscore_covform = "mainterm_inter", 
                     pscore_slfull = NA, pscore_sldefault = NA, pscore_sllibs = list(NA))
    
    # (2) Full SuperLearner with default libraries
    g2 <- expand.grid(pscore_para = F, pscore_link = NA, pscore_covform = NA, 
                     pscore_slfull = 1, pscore_sldefault = 1, pscore_sllibs = list(dfSLlibs))
    
    gmodels <- rbind(g1,g2)
  }
  
  # 1.B. Advanced mode for treatment mechanism -----------
  
  if (input$analysis_pscore_mode == "adv") {
    g1 <- g2 <- NA
  
    if (input$analysis_pscore_para == T) {
        g1 <- expand.grid(pscore_para = T, pscore_link = input$analysis_pscore_link,
                          pscore_covform = input$analysis_pscore_covform, 
                          pscore_slfull = NA, pscore_sldefault = NA, pscore_sllibs = list(NA))
    }
    
    if (input$analysis_pscore_sl == F) {
        g2 <- expand.grid(pscore_para = F, pscore_link = NA, pscore_covform = NA, 
                          pscore_slfull = input$analysis_pscore_slfull, pscore_sldefault = input$analysis_pscore_sldefault, pscore_sllibs = list(NA))
        if (input$analysis_pscore_sldefault == 1) g2$pscore_sllibs <- list(dfSLlibs)
        if (input$analysis_pscore_sldefault == 0) g2$pscore_sllibs <- list(input$analysis_pscore_sllibs)
    }
    
    if (!is.na(g1) & !is.na(g2)) { gmodels <- rbind(g1,g2)
    } else if (!is.na(g1) &  is.na(g2)) { gmodels <- g1
    } else if ( is.na(g1) & !is.na(g2)) { gmodels <- g2 }
    
  }
  
  # 2.A. Basic mode for outcome model ----------------
  
  if (input$analysis_q_mode == "basic") {
    
    # (1) Parametric GLM with Gaussian distribution and identity link for continuous outcomes and Binomial distribution with logit link for binomial outcomes, and covariates as main linear terms with first order interactions
    q1 <- expand.grid(q_para = T, q_dist = "gaus", q_link = "identity", q_covform = "mainterm_inter", 
                     q_slfull = NA, q_sldefault = NA, q_sllibs = list(NA))
    
    # (2) Full SuperLearner outcome model with default libraries
    q2 <- expand.grid(q_para = F, q_dist = NA, q_link = NA, q_covform = NA, 
                     q_slfull = 1, q_sldefault = 1, q_sllibs = list(dfSLlibs))
    
    qmodels <- rbind(q1,q2)
  }
  
  # 2.B. Advanced mode for outcome model ----------------
  
  if (input$analysis_q_mode == "adv") {
    q1 <- q2 <- NA
      
      if (input$analysis_q_para == T) {
          analysis_q_dist <- c("gaus","bin","gam","pois")[c(input$analysis_q_dist_gaus,input$analysis_q_dist_bin,input$analysis_q_dist_gam,input$analysis_q_dist_pois)]
          q1 <- expand.grid(q_para = T, q_dist = analysis_q_dist, q_link = NA, q_covform = input$analysis_q_covform, 
                            q_slfull = NA, q_sldefault = NA, q_sllibs = list(NA))
          
          q1$qlink[q1$q_dist=="gaus"] <- input$analysis_q_gaus_link
          q1$qlink[q1$q_dist=="bin" ] <- input$analysis_q_bin_link
          q1$qlink[q1$q_dist=="gam" ] <- input$analysis_q_gam_link
          q1$qlink[q1$q_dist=="pois"] <- input$analysis_q_pois_link
      }
  
      if (input$analysis_q_sl == T) {
          q2 <- expand.grid(q_para = F, q_dist = NA, q_link = NA, q_covform = NA, 
                            q_slfull = input$analysis_q_slfull, q_sldefault = input$analysis_q_sldefault, q_sllibs = list(NA))
          
          if (input$analysis_q_sldefault == 1) q2$q_sllibs <- list(dfSLlibs)
          if (input$analysis_q_sldefault == 0) q2$q_sllibs <- list(input$analysis_q_sllibs)
      }
      
      if (!is.na(q1) & !is.na(q2)) { qmodels <- rbind(q1,q2)
      } else if (!is.na(q1) & is.na(q2)) { qmodels <- q1
      } else if (is.na(q1) & !is.na(q2)) { qmodels <- q2 }
  }
  
  
  # 3. Main analysis components ----------------------
  
  method <- c("gcomp","pweight","tmle")[c(input$an_gcomp,input$an_pweight,input$an_tmle)]
  estimand <- c("ate","att")[c(input$ate,input$att)]
  metric <- c("rr","rd","or")[input$effect_metric]
  
  an <- expand.grid(estimand = estimand, metric = metric, method = c("unadj", method), diff = input$outcome_diff, covs = list(input$an_covs), pweight_type = NA)
  
  # Don't allow modified iptw for now:
  an$pweight_type[an$method == "pweight"] <- input$an_pweight_type[input$an_pweight_type != "pweight_modht"] 
  
  # 4. Bring all analyses components together -------------------------
  
  models   <- cbind( gmodels[rep(row.names(gmodels), nrow(qmodels)),], qmodels[rep(row.names(qmodels), each=nrow(gmodels)),] )
  
  analyses <- cbind( an[rep(row.names(an), nrow(models)),], models[rep(row.names(models), each=nrow(an)),] ) 
  
  
  list(designs, matches, analyses)

}



## For reference:

## Variables from ui script (10 vars): 
# Sampling: sample_size (numeric), all_exposed (0/1), all_controls (0/1)
# Target parameter: ate (T/F), att (T/F), effect metric (1/2/3)
# Timing and outcome: outcome_diff (0/1)
# Matching: match (0/1), match_mode ('basic','adv'), match_basic_more(T/F)

## Vars from matching script (28 vars):
# Basic matching: match_basic_covs (string vector)
# Advanced matching: match_exact (T/F), match_nonexact (T/F)
# Exact matching: match_exact_covs (string vector), 
# Non-exact matching: match_nonexact_covs (string vector), match_distance_{pscore, mahal, class, nnet} (T/F), match_method_{nn, opt, gen, sub, full} (T/F)
# NN matching:      match_nn_ratio (numeric),  match_nn_replace (T/F),  match_nn_exact_covs (string vector)
# Genetic matching: match_gen_ratio (numeric)
# Subclassficiation: match_sub_pscore_classes (numeric)
# Additional options: match_discard (0/1/2/3/4/5/6), match_discard_reest (T/F)

## Vars from pscore_options_match (10 vars): 
# Pscore options: match_pscore_para (T/F), match_pscore_sl (T/F)
# Parametric: match_pscore_link ('logit','probit','complementary log-log','log','Cauchy CDF'), match_pscore_backtrans (T/F), match_pscore_covform ('mainterm','mainterm_inter','quintile'), match_pscore_formula (text)
# Semi-parametric: match_pscore_slfull (0/1), match_pscore_sldefault (0/1), match_pscore_sllibs (string vector)

## Vars from analysis
# Analysis methods: an_gcomp (T/F), an_pweight (T/F), an_aiptw (T/F), an_tmle (T/F)
# Analysis covariates: an_covs (string vector)
# Additional pweight options: an_pweight_type ('pweight_ht', 'pweight_modht')

## vars from pscore_options_analysis
# Pscore options: analysis_pscore_mode ('basic','adv'), analysis_pscore_basic_more (T/F), analysis_pscore_para (T/F), analysis_pscore_sl (T/F)
# Parametric: analysis_pscore_link ('logit','probit','complementary log-log','log','Cauchy CDF'), analysis_pscore_covform ('mainterm','mainterm_inter','quintile'), analysis_pscore_formula (text)
# Semi-parametric: analysis_pscore_slfull (0/1), analysis_pscore_sldefault (0/1), analysis_pscore_sllibs (string vector)

## Vars from analysis_q_options
# Main options: analysis_q_mode('basic','adv'), analysis_q_basic_more (T/F), analysis_q_para (T/F), analysis_q_sl (T/F)
# Parametric model distribution: analysis_q_dist_gaus (T/F), analysis_q_dist_bin (T/F), analysis_q_dist_gam (T/F), analysis_q_dist_pois (T/F)
# Parametric links: analysis_q_gaus_link (string vector), analysis_q_bin_link (string vector), analysis_q_gam_link (string vector), analysis_q_pois_link (string vector)
# Parametric covariate forms: analysis_q_covform ('mainterm','mainterm_inter','quintile'), analysis_q_formula (text)
# Semi parametric: analysis_q_slfull (0/1), analysis_q_sldefault (0/1), analysis_q_sllibs (string vector)



