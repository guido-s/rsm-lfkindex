runsim <- function(params, nsim, func.n) {
  ##
  for (i in seq_len(nsim)) {
    ##
    n_selected <- 0 # reset counter
    dat.i <- NULL   # append to this, reset at beginning
    
    ## We want to select N_studies
    ##
    while (n_selected < params$N_studies) {
      ##
      draw1 <- generate_meta(
        k = 1, # generate 1 study at a time
        theta1 = params$theta_treat,
        theta2 = params$theta_control,
        sd1 = params$sd_treat,
        sd2 = params$sd_control,
        tau = params$tau,
        func.n = func.n
      )
      
      ## Selection probabilities
      ##
      probs <-
        calc_copas_probs(draw1,
                         theta1 = params$theta_treat - params$theta_control,
                         sd1 = params$sd_treat,
                         sd2 = params$sd_control,
                         rho = params$rho,
                         tau = params$tau,
                         p_small = params$p_small,
                         p_large = params$p_large,
                         n_small = params$n_small,
                         n_large = params$n_large
                         )
      
      ## Realize publish probability (TRUE = published)
      ##
      is.publ <- probs - runif(length(probs), 0, 1) > 0
      ##
      if (is.publ) {
        ## append data
        dat.i <-
          rbind(
            dat.i,
            with(draw1, data.frame(n.e, n.c, mean.e, mean.c, sd.e, sd.c))
          )
        ## Set counter
        ##
        n_selected = nrow(dat.i)
      }
    }
    
    ## Calculate final meta analysis
    ## - use REML estimator of tau2 with DL as fall-back
    ##
    
    m.i <-
      try(meta::metacont(n.e = dat.i$n.e, n.c = dat.i$n.c,
                         mean.e = dat.i$mean.e, mean.c = dat.i$mean.c,
                         sd.e = dat.i$sd.e, sd.c = dat.i$sd.c,
                         method.tau = "REML", method.tau.ci = ""),
          silent = TRUE)
    ##
    if ("try-error" %in% class(m.i))
      m.i <-
      meta::metacont(n.e = dat.i$n.e, n.c = dat.i$n.c,
                     mean.e = dat.i$mean.e, mean.c = dat.i$mean.c,
                     sd.e = dat.i$sd.e, sd.c = dat.i$sd.c,
                     method.tau = "DL", method.tau.ci = "")
    
    ## Conduct tests for funnel plot asymmetry
    ##
    test.i <- testbias(m.i)
    
    ## Results
    ##
    res.i <- data.frame(scenario = params$scenario,
                        N_studies = params$N_studies,
                        theta_treat = params$theta_treat,
                        sd_treat = params$sd_treat,
                        theta_control = params$theta_control,
                        sd_control = params$sd_control,
                        tau = params$tau,
                        rho = params$rho,
                        ##
                        idx = i,
                        ##
                        TE.common = m.i$TE.common,
                        seTE.common = m.i$seTE.common,
                        stat.common = m.i$statistic.common,
                        pval.common = m.i$pval.common,
                        ##
                        TE.random = m.i$TE.random,
                        seTE.random = m.i$seTE.random,
                        stat.random = m.i$statistic.random,
                        pval.random = m.i$pval.random,
                        ##
                        method.tau = m.i$method.tau,
                        tau2 = m.i$tau2,
                        I2 = m.i$I2,
                        Q = m.i$Q,
                        df.Q = m.i$df.Q,
                        pval.Q = m.i$pval.Q,
                        ##
                        lfkindex = test.i$lfkindex,
                        #
                        stat.egger = test.i$stat.egger,
                        df.egger = test.i$df.egger,
                        pval.egger = test.i$pval.egger,
                        #
                        stat.rank = test.i$stat.rank,
                        pval.rank = test.i$pval.rank,
                        #
                        stat.thompson = test.i$stat.thompson,
                        df.thompson = test.i$df.thompson,
                        pval.thompson = test.i$pval.thompson,
                        method.tau.thompson = test.i$method.tau
                        )
    ##
    if (i == 1)
      res <- res.i
    else
      res <- rbind(res, res.i)
  }
  ##
  cat(".")
  ##
  attr(res, "nsim") <- nsim
  ##
  res
}


## Function to draw N
##
draw_n_Schwarzer2002 <- function(x) {
  n <- rep(0, x)
  sel <- n < 30
  ##
  while (any(n < 30)) {
    n[sel] <- rlnorm(sum(sel), 3.798, sqrt(1.104058))
    n[sel] <- ifelse(n[sel] %% 2, n[sel] + 2 - n[sel] %% 2, n[sel])
    sel <- n < 30
  }
  ##
  list(e = n / 2, c = n / 2)
}
##
draw_n_unif <- function(x, min = 50, max = 500) {
  n1 <- n2 <- ceiling(runif(x, min, max))
  list(e = n1, c = n2)
}


## Generate meta-analysis object
##
generate_meta <- function(k, theta1, sd1, theta2, sd2, tau, func.n) {
  
  ## draw study Ns
  ##
  N <- do.call(func.n, list(k))
  
  ## 1. define observed TE:
  ##    TE_true + sqrt(true TE_variance + tau^2) * error
  ##    [random draw from standard normal]
  ##
  ## True variance of the sampling distribution of the difference
  ## between means (two samples)
  ##
  TE_var <- sd1^2 / N$e + sd2^2 / N$c
  TE_observed <-
    (theta1 - theta2) + sqrt(TE_var + tau^2) *
    rnorm(k, mean = 0, sd = 1)
  
  ## 2. Draw study means from sampling distribution of mean
  ##    [sd / sqrt(N)] - centered on true mean
  ##
  mean_treat_observed <-
    rnorm(k, theta1, sd1 / sqrt(N$e))
  ## draw study SDs from sampling distribution 
  sd1_observed <-
    rnorm(k, sd1, sd1 / 10) # constant variation of sd estimate
  
  ## 3. Calc M for control group as observed treatment group minus
  ##    observed TE
  ##
  mean_control_observed <- mean_treat_observed - TE_observed
  
  ## 4. Draw SD for control group randomly
  ##
  sd2_observed <-
    rnorm(k, sd2, sd2 / 10) # constant variation of sd estimate
  
  ## 5. Plug into metacont()
  ##
  res <-
    meta::metacont(n.e = N$e, 
                   mean.e = mean_treat_observed,
                   sd.e = sd1_observed,
                   n.c = N$c,
                   mean.c = mean_control_observed,
                   sd.c = sd2_observed,
                   method.tau = "DL", method.tau.ci = "")
  
  res
}


## Selection probabilites from Copas selection model
##
calc_copas_probs <- function(x,
                             theta1,
                             sd1, sd2 = 1,
                             rho, tau,
                             p_small, p_large, # range of publication probabilities
                             n_small, n_large  # range of study SEs (95%CI)
                             ) {
  
  ## Calibrate gammas for probability of inclusion to scenario sd's
  ## pooled sd for two groups with same N
  ## (because sd1 and sd2 can be different)
  sd_pool <- sqrt((sd1^2 + sd2^2) / 2)
  ## standard errors
  s_small <- sd_pool / sqrt(n_small)
  s_large <- sd_pool / sqrt(n_large)

  ## gammas
  gamma_1 <- (qnorm(p_large) - qnorm(p_small)) / (1 / s_large - 1 / s_small) 
  gamma_0 <- qnorm(p_small) - gamma_1 / s_small
  
  ## Formula from
  ## https://journal.r-project.org/archive/2009/RJ-2009-012/RJ-2009-012.pdf
  ##
  ## true variance of the sampling distribution of the difference between means (two samples)
  ##
  sigma2 <- sd1^2 / x$n.e + sd2^2 / x$n.c
  ##
  nom <- gamma_0 + gamma_1 / x$seTE +
    rho * sqrt(sigma2) * (x$TE - theta1) / (sigma2 + tau^2)
  denom  <- sqrt(1 - rho^2 * sigma2 / (sigma2 + tau^2))
  Phi <- pnorm(nom / denom, mean = 0, sd = 1) # probability of inclusion
  
  Phi
}


## Conduct tests for funnel plot asymmetry
##
testbias <- function(x, params) {
  
  lfk <- metasens::lfkindex(x)
  egger <- meta::metabias(x, method = "Egger")
  rank <- meta::metabias(x, method = "Begg")
  #
  thompson <- try(meta::metabias(x, method = "Thompson"), silent = TRUE)
  #
  if ("try-error" %in% class(thompson)) {
    thompson <-
      meta::metabias(update(x, method.tau = "DL"), method = "Thompson")
  }
  
  res <- data.frame(
    lfkindex = lfk$lfkindex, # LFK-index
    stat.egger = egger$statistic, # Egger's test
    df.egger = egger$df,
    pval.egger = egger$pval,
    stat.rank = rank$statistic, # Begg and Mazumdar test
    pval.rank = rank$pval,
    #
    stat.thompson = thompson$statistic, # Thompson-Sharp test
    df.thompson = thompson$df,
    pval.thompson = thompson$pval,
    method.tau = thompson$x$method.tau
  )
  
  res
}
