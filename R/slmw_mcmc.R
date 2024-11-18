#' @title Spatial Linear Regression Modeling Weibull Distributions

#' @description
#' Markov chain Monte Carlo (MCMC) algorithm for linear spatial model.

#' @usage
#' slmw_mcmc = function(nsim, W, S, Sstar, startingPoint = NULL,
#'                       details = FALSE, forecasting = TRUE, burn = 0,
#'                       thin = 1, verbose = +Inf)

#' @details
#' MCMC sampler for linear spatial model without the binomial random
#'  variables. Details are the same of \code{slmbw_mcmc} except the part for
#'  the binomial distributions.

#' @references
#' ...

#' @param nsim The number of simulations.
#' @param W List of sublists containing Weibull realizations.
#' @param S Matrix of covariates for observed point.
#' @param Sstar Matrix of covariates for interpolating points.
#' @param startingPoint List for starting point.
#' @param details If it is true then also values of Gamma, Delta, and Pi are
#'  collected (large memory consumption).
#' @param forecasting If it true then hyper-parameters values at interpolating
#'  points are collected.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return A named list containing a posterior sample.

#' @export
slmw_mcmc = function(
  nsim, W, S, Sstar, startingPoint = NULL, details = FALSE,
  forecasting = TRUE, burn = 0, thin = 1, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste("slmw_mcmc: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn <= 0) burn = 0



  ##### prior parameters

  ### inverse compounded gamma for variances of white noises
  k_sigmaSq = 2
  alpha_sigmaSq = 0.5
  beta_sigmaSq = 0.5

  ### Student-t for coefficient
  nu_beta = 2
  xi_beta = 0
  omegaSq_beta = 1



  ##### get dimensions

  ### number of in sample spatial points
  M = nrow(S)

  ### number of out of sample spatial points
  if (forecasting) Mstar = nrow(Sstar)

  ### dimension of spatial coordinates
  p = ncol(S)

  ### number of replications for each spatial point
  T = length(W[[1]])

  ### add intercept
  S = cbind(rep(1, M), S)
  if (forecasting) Sstar = cbind(rep(1, Mstar), Sstar)




  ########## initialize MCMC

  ### get the values of binomial random variables
  N = matrix(0, nrow = M, ncol = T)
  for (i in 1:M) N[i, ] = lengths(W[[i]])

  ### get objects for vectorization of arguments of Weibull loglikelihood
  vecW = double(length = M * T)
  vecGamma = double(length = M * T)
  vecDelta = double(length = M * T)
  vecGamma_prop = double(length = M * T)
  vecDelta_prop = double(length = M * T)

  ### get vectorized arguments for Weibull random variables
  ll = 0
  for (i in 1:T) {
    for (j in 1:M) {
      vecW[(ll + 1):(ll + N[j, i])] = W[[j]][[i]]
      ll = ll + N[j, i]
    }
  }
  log_vecW = log(vecW)



  ##### initial points of the chain

  ### check if starting point is provided
  if (is.null(startingPoint)) {

    # Weibull and binomial parameters (mapped into the real axis )
    Gamma = matrix(nrow = M, ncol = T)
    Delta = matrix(nrow = M, ncol = T)
    #Pi = qlogis(N / 365)

    for (i in 1:M) {
      for (j in 1:T) {
        if (length(W[[i]][[j]]) > 1) {
          mle_wei = suppressWarnings(
            EnvStats::eweibull(W[[i]][[j]])$parameters
          )
        } else {
          mle_wei = suppressWarnings(
            EnvStats::eweibull(c(W[[i]][[j]], 1e-1))$parameters
          )
        }
        Gamma[i, j] = log(mle_wei["shape"])
        Delta[i, j] = log(mle_wei["scale"])
      }
    }

    # starting point for the coefficients is obtained through mle of means
    projMat = tryCatch(
      expr = solve(crossprod(S)) %*% t(S),
      error = function(err) {
        svd_tmp = svd(crossprod(S), nv = 0)
        return(tcrossprod(svd_tmp$u %*% diag(1 / sqrt(svd_tmp$d))) %*% t(S))
      }
    )
    beta_gamma = c(projMat %*% rowMeans(Gamma))
    beta_delta = c(projMat %*% rowMeans(Delta))
    #beta_pi = c(projMat %*% rowMeans(Pi))

    # variances of noises
    sigmaSq_gamma = mean((Gamma - c(S %*% beta_gamma)) ^ 2)
    sigmaSq_delta = mean((Delta - c(S %*% beta_delta)) ^ 2)
    #sigmaSq_pi = mean((Pi - c(S %*% beta_pi)) ^ 2)

  } else {

    # Weibull and binomial parameters (mapped into the real axis)
    Gamma = startingPoint[["Gamma"]]
    Delta = startingPoint[["Delta"]]
    #Pi = startingPoint[["Pi"]]

    # coefficients
    beta_gamma = startingPoint[["beta_gamma"]]
    beta_delta = startingPoint[["beta_delta"]]
    #beta_pi = startingPoint[["beta_pi"]]

    # variances of noises
    sigmaSq_gamma = startingPoint[["sigmaSq_gamma"]]
    sigmaSq_delta = startingPoint[["sigmaSq_delta"]]
    #sigmaSq_pi = startingPoint[["sigmaSq_pi"]]

  }



  ##### initialize dynamic and storing variables for MCMC

  ### get vectorized arguments for Weibull
  ll = 0
  for (i in 1:T) {
    for (j in 1:M) {
      vecGamma[(ll + 1):(ll + N[j, i])] = rep(Gamma[j, i], N[j, i])
      vecDelta[(ll + 1):(ll + N[j, i])] = rep(Delta[j, i], N[j, i])
      ll = ll + N[j, i]
    }
  }
  exp_vecGamma = exp(vecGamma)

  ### get probabilities for binomial
  #probs_binom = c(plogis(Pi))
  #probs_binom[probs_binom >= 1 - sqrt(.Machine$double.eps)] = 1 - sqrt(
  #  .Machine$double.eps
  #)
  #probs_binom[probs_binom <= sqrt(.Machine$double.eps)] = sqrt(
  #  .Machine$double.eps
  #)

  ### compute current log-likelihood for Weibull
  loglike_weibull = tryCatch(
    expr = {
      sum(
        vecGamma - exp_vecGamma * vecDelta +
          (exp_vecGamma - 1) * log_vecW -
          (vecW / exp(vecDelta)) ^ exp_vecGamma
      )
    }, error = function(err) -Inf
  )
  if (is.na(loglike_weibull)) loglike_weibull = -Inf

  ### compute current log-likelihood for binomial
  #loglike_binom = sum(
  #  N * log(probs_binom) + (365 - N) * log(1 - probs_binom)
  #)

  ### storing arrays
  keep_mix_beta_gamma = matrix(nrow = p + 1, ncol = nsim)
  keep_mix_beta_delta = matrix(nrow = p + 1, ncol = nsim)
  #keep_mix_beta_pi = matrix(nrow = p + 1, ncol = nsim)

  keep_mix_sigmaSq = matrix(nrow = 2, ncol = nsim)
  rownames(keep_mix_sigmaSq) = c("gamma", "delta")#, "pi")

  keep_beta_gamma = matrix(nrow = p + 1, ncol = nsim)
  keep_beta_delta = matrix(nrow = p + 1, ncol = nsim)
  #keep_beta_pi = matrix(nrow = p + 1, ncol = nsim)

  keep_sigmaSq = matrix(nrow = 2, ncol = nsim)
  rownames(keep_sigmaSq) = c("gamma", "delta")#, "pi")

  if (details) {
    keep_Gamma = array(dim = c(M, T, nsim))
    keep_Delta = array(dim = c(M, T, nsim))
    #keep_Pi = array(dim = c(M, T, nsim))
  }

  if (forecasting) {
    keep_mu_gamma_star = matrix(nrow = Mstar, ncol = nsim)
    keep_mu_delta_star = matrix(nrow = Mstar, ncol = nsim)
    #keep_mu_pi_star = matrix(nrow = Mstar, ncol = nsim)
  }




  ########## draw the chain
  for (insim in (1 - burn):nsim) {

    ithin = 0

    ##### iterate up to thinning value
    while (ithin < thin) {

      ##########################################################################
      #                                                                        #
      #                        update mixing variables                         #
      #                                                                        #
      ##########################################################################

      ### coefficients of linear functions
      mix_beta_gamma = 1 / rgamma(
        p + 1, shape = 0.5 * (nu_beta + 1),
        rate = 0.5 * (nu_beta + (beta_gamma - xi_beta)^2 / omegaSq_beta)
      )
      mix_beta_delta = 1 / rgamma(
        p + 1, shape = 0.5 * (nu_beta + 1),
        rate = 0.5 * (nu_beta + (beta_delta - xi_beta)^2 / omegaSq_beta)
      )
      #mix_beta_pi = 1 / rgamma(
      #  p + 1, shape = 0.5 * (nu_beta + 1),
      #  rate = 0.5 * (nu_beta + (beta_pi - xi_beta)^2 / omegaSq_beta)
      #)

      ### variances of white noises
      mix_sigmaSq_gamma = rgamma(
        1, shape = alpha_sigmaSq + k_sigmaSq,
        rate = beta_sigmaSq + 1 / sigmaSq_gamma
      )
      mix_sigmaSq_delta = rgamma(
        1, shape = alpha_sigmaSq + k_sigmaSq,
        rate = beta_sigmaSq + 1 / sigmaSq_delta
      )
      #mix_sigmaSq_pi = rgamma(
      #  1, shape = alpha_sigmaSq + k_sigmaSq,
      #  rate = beta_sigmaSq + 1 / sigmaSq_pi
      #)





      ##########################################################################
      #                                                                        #
      #                         update scale variables                         #
      #                                                                        #
      ##########################################################################

      ### variances of white noises
      sigmaSq_gamma = 1 / rgamma(
        1, shape = k_sigmaSq + 0.5 * (M * T),
        rate = mix_sigmaSq_gamma + 0.5 * sum((c(S %*% beta_gamma) - Gamma)^2)
      )
      sigmaSq_delta = 1 / rgamma(
        1, shape = k_sigmaSq + 0.5 * (M * T),
        rate = mix_sigmaSq_delta + 0.5 * sum((c(S %*% beta_delta) - Delta)^2)
      )
      #sigmaSq_pi = 1 / rgamma(
      #  1, shape = k_sigmaSq + 0.5 * (M * T),
      #  rate = mix_sigmaSq_pi + 0.5 * sum((c(S %*% beta_pi) - Pi)^2)
      #)





      ##########################################################################
      #                                                                        #
      #                       update location variables                        #
      #                                                                        #
      ##########################################################################

      ##### Weibull parameters

      ### get a sample from the prior
      # shape
      beta_gamma_prior = rnorm(p + 1, mean = 0, sd = sqrt(mix_beta_gamma))
      Gamma_prior = c(S %*% beta_gamma_prior) + sqrt(sigmaSq_gamma) *
        matrix(rnorm(M * T), nrow = M, ncol = T)
      # scale
      beta_delta_prior = rnorm(p + 1, mean = 0, sd = sqrt(mix_beta_delta))
      Delta_prior = c(S %*% beta_delta_prior) + sqrt(sigmaSq_delta) *
        matrix(rnorm(M * T), nrow = M, ncol = T)

      ### get threshold for acceptance and initial angle
      strike = loglike_weibull - rgamma(1, shape = 1, rate = 1)
      angle = runif(1, min = 0, max = 2 * pi)

      ### get initial brackets
      angle_min = angle - 2 * pi
      angle_max = angle

      ### try until acceptance occurs
      while (TRUE) {
        beta_gamma_prop = beta_gamma * cos(angle) +
          beta_gamma_prior * sin(angle)
        Gamma_prop = Gamma * cos(angle) + Gamma_prior * sin(angle)
        beta_delta_prop = beta_delta * cos(angle) +
          beta_delta_prior * sin(angle)
        Delta_prop = Delta * cos(angle) + Delta_prior * sin(angle)

        ll = 0
        for (i in 1:T) {
          for (j in 1:M) {
            vecGamma_prop[(ll + 1):(ll + N[j, i])] = rep(
              Gamma_prop[j, i], N[j, i]
            )
            vecDelta_prop[(ll + 1):(ll + N[j, i])] = rep(
              Delta_prop[j, i], N[j, i]
            )
            ll = ll + N[j, i]
          }
        }
        exp_vecGamma_prop = exp(vecGamma_prop)

        loglike_weibull_prop = tryCatch(
          expr = {
            sum(
              vecGamma_prop - exp_vecGamma_prop * vecDelta_prop +
                (exp_vecGamma_prop - 1) * log_vecW -
                (vecW / exp(vecDelta_prop)) ^ exp_vecGamma_prop
            )
          }, error = function(err) return(-Inf)
        )
        if (is.na(loglike_weibull_prop)) loglike_weibull_prop = -Inf

        if (loglike_weibull_prop >= strike) {
          beta_gamma = beta_gamma_prop
          Gamma = Gamma_prop
          vecGamma = vecGamma_prop
          exp_vecGamma = exp_vecGamma_prop
          beta_delta = beta_delta_prop
          Delta = Delta_prop
          vecDelta = vecDelta_prop
          loglike_weibull = loglike_weibull_prop
          break
        } else {
          if (angle < 0) angle_min = angle else angle_max = angle
          angle = runif(1, min = angle_min, max = angle_max)
        }
      }



      ##### binomial parameters

      ### get a sample from the prior
      # shape
      #beta_pi_prior = rnorm(p + 1, mean = 0, sd = sqrt(mix_beta_pi))
      #Pi_prior = c(S %*% beta_pi_prior) + sqrt(sigmaSq_pi) *
      #  matrix(rnorm(M * T), nrow = M, ncol = T)

      ### get threshold for acceptance and initial angle
      #strike = loglike_binom - rgamma(1, shape = 1, rate = 1)
      #angle = runif(1, min = 0, max = 2 * pi)

      ### get initial brackets
      #angle_min = angle - 2 * pi
      #angle_max = angle

      ### try until acceptance occurs
      #while (TRUE) {
      #  beta_pi_prop = beta_pi * cos(angle) + beta_pi_prior * sin(angle)
      #  Pi_prop = Pi * cos(angle) + Pi_prior * sin(angle)

      #  probs_binom = c(plogis(Pi_prop))
      #  probs_binom[probs_binom >= 1 - sqrt(.Machine$double.eps)] = 1 - sqrt(
      #    .Machine$double.eps
      #  )
      #  probs_binom[probs_binom <= sqrt(.Machine$double.eps)] = sqrt(
      #    .Machine$double.eps
      #  )

      #  loglike_binom_prop = sum(
      #    N * log(probs_binom) + (365 - N) * log(1 - probs_binom)
      #  )

      #  if (loglike_binom_prop >= strike) {
      #    beta_pi = beta_pi_prop
      #    Pi = Pi_prop
      #    loglike_binom = loglike_binom_prop
      #    break
      #  } else {
      #    if (angle < 0) angle_min = angle else angle_max = angle
      #    angle = runif(1, min = angle_min, max = angle_max)
      #  }
      #}





      ##########################################################################
      #                                                                        #
      #                         end of a complete scan                         #
      #                                                                        #
      ##########################################################################

      if (insim > 0) ithin = ithin + 1 else ithin = thin

    }

    ##########################################################################
    #                                                                        #
    #                          store current values                          #
    #                                                                        #
    ##########################################################################

    ### mixing variables
    keep_mix_beta_gamma[, insim] = mix_beta_gamma
    keep_mix_beta_delta[, insim] = mix_beta_delta
    #keep_mix_beta_pi[, insim] = mix_beta_pi
    keep_mix_sigmaSq[, insim] = c(
      mix_sigmaSq_gamma, mix_sigmaSq_delta#, mix_sigmaSq_pi
    )

    ### variances of white noises
    keep_sigmaSq[, insim] = c(sigmaSq_gamma, sigmaSq_delta)#, sigmaSq_pi)

    ### coefficients
    keep_beta_gamma[, insim] = beta_gamma
    keep_beta_delta[, insim] = beta_delta
    #keep_beta_pi[, insim] = beta_pi

    ### if required, collect parameters (large memory consumption)
    if (details) {
      keep_Gamma[, , insim] = Gamma
      keep_Delta[, , insim] = Delta
      #keep_Pi[, , insim] = Pi
    }

    ### if required, forecast parameters for other points
    if (forecasting) {

      keep_mu_gamma_star[, insim] = c(Sstar %*% beta_gamma)
      keep_mu_delta_star[, insim] = c(Sstar %*% beta_delta)
      #keep_mu_pi_star[, insim] = c(Sstar %*% beta_pi)
    }





    ### print status of the chain
    if (insim %% verbose == 0) {
      print(paste(
        "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
        sep = ""
      ))
    }

  }




  ### print end time if required
  if (!is.infinite(verbose)) {
    print(paste("slmw_mcmc: end time at ", Sys.time(), sep = ""))
  }




  ########## return results

  ### initialize list
  res = list()

  ### mixing variables
  res[["mix_beta_gamma"]] = keep_mix_beta_gamma
  res[["mix_beta_delta"]] = keep_mix_beta_delta
  #res[["mix_beta_pi"]] = keep_mix_beta_pi
  res[["mix_sigmaSq"]] = keep_mix_sigmaSq

  ### variances of white noises
  res[["sigmaSq"]] = keep_sigmaSq

  ### linear coefficients
  res[["beta_gamma"]] = keep_beta_gamma
  res[["beta_delta"]] = keep_beta_delta
  #res[["beta_pi"]] = keep_beta_pi

  ### values of parameters
  if (details) {
    res[["Gamma"]] = keep_Gamma
    res[["Delta"]] = keep_Delta
    #res[["Pi"]] = keep_Pi
  }

  ### hyper-parameters values for other points
  if (forecasting) {
    res[["mu_gamma_star"]] = keep_mu_gamma_star
    res[["mu_delta_star"]] = keep_mu_delta_star
    #res[["mu_pi_star"]] = keep_mu_pi_star
  }

  ### return list of results
  return(res)

}