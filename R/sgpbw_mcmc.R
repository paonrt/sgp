#' @title Spatial Gaussian Processes Modeling Binomial and Weibull Distributions

#' @description
#' Markov chain Monte Carlo (MCMC) algorithm for semi-parametric spatial model.

#' @usage
#' sgpbw_mcmc = function(nsim, W, S, Sstar, startingPoint = NULL,
#'                       details = FALSE, forecasting = TRUE, burn = 0,
#'                       thin = 1, verbose = +Inf)

#' @details
#' The parametrization for Weibull is such that \eqn{\mathrm{E}(W) = \log
#'  (\delta) \: \Gamma \big( 1 + 1 / \log(\gamma) \big)}.
#'
#' Binomial parameters are \eqn{n = 365} and \eqn{p = (1 + \exp \big( -\pi)
#'  \big)^{-1}}.
#'
#' \code{W}: list containing sublists for each observed point. Each sublist
#'  contains a vector of Weibull realizations. There is a vector for each year,
#'  thus, the length of each vector determines the binomial draws.
#'
#' \code{S}: matrix containing covariates of observed points. One column for
#'  each covariates and one row for each observed point. Note that, \eqn{i}-th
#'  row of \code{S} contains covariates for observed point with draws collected
#'  at the \eqn{i}-th sublist of \code{W}.
#'
#' \code{Sstar}: matrix containing covariates of interpolating points.
#'  One column for each covariates and one row for each observed point.
#'
#' \code{startingPoint}: the starting point is generated in the following way:
#' \itemize{
#'  \item if \code{startingPoint == NULL} then:
#'  \enumerate{
#'    \item set \code{Pi} equal to the MLE,
#'    \item set \code{mu_Pi = rowMeans(Pi)},
#'    \item set \code{psi_Pi = mean(mu_Pi)},
#'    \item set \code{tauSq_Pi = mean(mu_Pi^2) - psi_Pi^2},
#'    \item set \code{sigmaSq_Pi = mean(rowMeans(Pi^2) - mu_Pi^2)}
#'    \item set \code{lengthScale_Pi} equal to the MLE w.r.t.
#'     \code{(mu_pi - psi_pi) / sqrt(tauSq_pi)},
#'    \item do the same for \code{Gamma} and \code{Delta} and their
#'      hyper-parameters.
#'  }
#'  \item \code{else} it must be a named list with names equal to the variables
#'    to be initialized
#' }
#'
#' Only one value every \code{thin} values is kept in the chain, so the true
#'  number of complete scans will be \code{nsim * thin + burn}. By default
#'  \code{thin = 1}, that is no thinning.
#'
#' The current time and the current number of iteration are printed one every
#' \code{verbose} iterations. Furthermore:
#' \itemize{
#'  \item if \code{verbose == +-Inf} then there is no printing,
#'  \item if \code{verbose != +-Inf} then at least start and end of simulation
#'   are reported.
#' }

#' @references
#' ...

#' @param nsim The number of simulations.
#' @param W List of sublists containing Weibull realizations.
#' @param S Matrix of covariates for observed point.
#' @param Sstar Matrix of covariates for interpolating points.
#' @param startingPoint List for starting point.
#' @param details If it is true then also values of Gamma, Delta, and Pi are
#'  collected (large memory consumption).
#' @param forecasting If it true then Gaussian processes values at interpolating
#'  points are collected.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return A named list containing a posterior sample.

#' @export
sgpbw_mcmc = function(
  nsim, W, S, Sstar, startingPoint = NULL, details = FALSE,
  forecasting = TRUE, burn = 0, thin = 1, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste("sgpbw_mcmc: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn <= 0) burn = 0



  ##### prior parameters

  ### log-normal prior for lengthscales
  meanlog_lengthScales = 0
  varlog_lengthScales = 4

  ### inverse compounded gamma for variances of Gaussian processes
  k_tauSq = 2
  alpha_tauSq = 0.5
  beta_tauSq = 0.5

  ### inverse compounded gamma for variances of white noises
  k_sigmaSq = 2
  alpha_sigmaSq = 0.5
  beta_sigmaSq = 0.5

  ### Student-t for mean of Gaussian processes
  nu_psi = 2
  xi_psi = 0
  omegaSq_psi = 1



  ##### get dimensions

  ### number of in sample spatial points
  M = nrow(S)

  ### number of out of sample spatial points
  if (forecasting) Mstar = nrow(Sstar)

  ### dimension of spatial coordinates
  p = ncol(S)

  ### number of replications for each spatial point
  T = length(W[[1]])




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

  ### get distances matrix
  distMat_list = list()
  for (i in 1:p) distMat_list[[i]] = unidist(S[, i])



  ##### initial points of the chain

  ### check if starting point is provided
  if (is.null(startingPoint)) {

    # Weibull and binomial parameters (mapped into the real axis)
    Gamma = matrix(nrow = M, ncol = T)
    Delta = matrix(nrow = M, ncol = T)
    Pi = qlogis(N / 365)

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

    # Gaussian processes
    mu_gamma = rowMeans(Gamma)
    mu_delta = rowMeans(Delta)
    mu_pi = rowMeans(Pi)

    # means of Gaussian processes
    psi_gamma = mean(mu_gamma)
    psi_delta = mean(mu_delta)
    psi_pi = mean(mu_pi)

    # variances of Gaussian processes
    tauSq_gamma = mean(mu_gamma^2) - psi_gamma^2
    tauSq_delta = mean(mu_delta^2) - psi_delta^2
    tauSq_pi = mean(mu_pi^2) - psi_pi^2

    # length scales of Gaussian processes
    lengthScale_gamma = hyperTuner_corrMatern1.5(
      (mu_gamma - psi_gamma) / sqrt(tauSq_gamma), D_list = distMat_list
    )
    lengthScale_delta = hyperTuner_corrMatern1.5(
      (mu_delta - psi_delta) / sqrt(tauSq_delta), D_list = distMat_list
    )
    lengthScale_pi = hyperTuner_corrMatern1.5(
      (mu_pi - psi_pi) / sqrt(tauSq_pi), D_list = distMat_list
    )

    # variances of noises
    sigmaSq_gamma = mean(rowMeans(Gamma^2) - mu_gamma^2)
    sigmaSq_delta = mean(rowMeans(Delta^2) - mu_delta^2)
    sigmaSq_pi = mean(rowMeans(Pi^2) - mu_pi^2)

  } else {

    # Weibull and binomial parameters (mapped into the real axis )
    Gamma = startingPoint[["Gamma"]]
    Delta = startingPoint[["Delta"]]
    Pi = startingPoint[["Pi"]]

    # Gaussian processes
    mu_gamma = startingPoint[["mu_gamma"]]
    mu_delta = startingPoint[["mu_delta"]]
    mu_pi = startingPoint[["mu_pi"]]

    # means of Gaussian processes
    psi_gamma = startingPoint[["psi_gamma"]]
    psi_delta = startingPoint[["psi_delta"]]
    psi_pi = startingPoint[["psi_pi"]]

    # variances of Gaussian processes
    tauSq_gamma = startingPoint[["tauSq_gamma"]]
    tauSq_delta = startingPoint[["tauSq_delta"]]
    tauSq_pi = startingPoint[["tauSq_pi"]]

    # length scales of Gaussian processes
    lengthScale_gamma = startingPoint[["lengthScale_gamma"]]
    lengthScale_delta = startingPoint[["lengthScale_delta"]]
    lengthScale_pi = startingPoint[["lengthScale_pi"]]

    # variances of noises
    sigmaSq_gamma = startingPoint[["sigmaSq_gamma"]]
    sigmaSq_delta = startingPoint[["sigmaSq_delta"]]
    sigmaSq_pi = startingPoint[["sigmaSq_pi"]]

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
  probs_binom = c(plogis(Pi))
  probs_binom[probs_binom >= 1 - sqrt(.Machine$double.eps)] = 1 - sqrt(
    .Machine$double.eps
  )
  probs_binom[probs_binom <= sqrt(.Machine$double.eps)] = sqrt(
    .Machine$double.eps
  )

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
  loglike_binom = sum(
    N * log(probs_binom) + (365 - N) * log(1 - probs_binom)
  )

  ### storing arrays
  keep_mix_psi = matrix(nrow = 3, ncol = nsim)
  rownames(keep_mix_psi) = c("gamma", "delta", "pi")

  keep_mix_sigmaSq = matrix(nrow = 3, ncol = nsim)
  rownames(keep_mix_sigmaSq) = c("gamma", "delta", "pi")

  keep_mix_tauSq = matrix(nrow = 3, ncol = nsim)
  rownames(keep_mix_tauSq) = c("gamma", "delta", "pi")

  keep_lengthScale_gamma = matrix(nrow = p, ncol = nsim)
  keep_lengthScale_delta = matrix(nrow = p, ncol = nsim)
  keep_lengthScale_pi = matrix(nrow = p, ncol = nsim)

  keep_sigmaSq = matrix(nrow = 3, ncol = nsim)
  rownames(keep_sigmaSq) = c("gamma", "delta", "pi")

  keep_tauSq = matrix(nrow = 3, ncol = nsim)
  rownames(keep_tauSq) = c("gamma", "delta", "pi")

  keep_psi = matrix(nrow = 3, ncol = nsim)
  rownames(keep_psi) = c("gamma", "delta", "pi")

  keep_mu_gamma = matrix(nrow = M, ncol = nsim)
  keep_mu_delta = matrix(nrow = M, ncol = nsim)
  keep_mu_pi = matrix(nrow = M, ncol = nsim)

  if (details) {
    keep_Gamma = array(dim = c(M, T, nsim))
    keep_Delta = array(dim = c(M, T, nsim))
    keep_Pi = array(dim = c(M, T, nsim))
  }

  if (forecasting) {
    keep_mu_gamma_star = matrix(nrow = Mstar, ncol = nsim)
    keep_mu_delta_star = matrix(nrow = Mstar, ncol = nsim)
    keep_mu_pi_star = matrix(nrow = Mstar, ncol = nsim)

    ### get distances matrix
    starDistMat_list = list()
    crossDistMat_list = list()
    for (i in 1:p) {
      starDistMat_list[[i]] = unidist(Sstar[, i])
      crossDistMat_list[[i]] = cross_unidist(Sstar[, i], S[, i])
    }
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

      ### means of Gaussian processes
      mix_psi_gamma = 1 / rgamma(
        1, shape = 0.5 * (nu_psi + 1),
        rate = 0.5 * (nu_psi + (psi_gamma - xi_psi)^2 / omegaSq_psi)
      )
      mix_psi_delta = 1 / rgamma(
        1, shape = 0.5 * (nu_psi + 1),
        rate = 0.5 * (nu_psi + (psi_delta - xi_psi)^2 / omegaSq_psi)
      )
      mix_psi_pi = 1 / rgamma(
        1, shape = 0.5 * (nu_psi + 1),
        rate = 0.5 * (nu_psi + (psi_pi - xi_psi)^2 / omegaSq_psi)
      )

      ### variances of white noises
      mix_sigmaSq_gamma = rgamma(
        1, shape = alpha_sigmaSq + k_sigmaSq,
        rate = beta_sigmaSq + 1 / sigmaSq_gamma
      )
      mix_sigmaSq_delta = rgamma(
        1, shape = alpha_sigmaSq + k_sigmaSq,
        rate = beta_sigmaSq + 1 / sigmaSq_delta
      )
      mix_sigmaSq_pi = rgamma(
        1, shape = alpha_sigmaSq + k_sigmaSq,
        rate = beta_sigmaSq + 1 / sigmaSq_pi
      )

      ### variances of Gaussian processes
      mix_tauSq_gamma = rgamma(
        1, shape = alpha_tauSq + k_tauSq,
        rate = beta_tauSq + 1 / tauSq_gamma
      )
      mix_tauSq_delta = rgamma(
        1, shape = alpha_tauSq + k_tauSq,
        rate = beta_tauSq + 1 / tauSq_delta
      )
      mix_tauSq_pi = rgamma(
        1, shape = alpha_tauSq + k_tauSq,
        rate = beta_tauSq + 1 / tauSq_pi
      )





      ##########################################################################
      #                                                                        #
      #                          update length scales                          #
      #                                                                        #
      ##########################################################################

      ##### log-shape of Weibull random variables

      ### get Kbar list of gamma for current values
      listKbar_gamma = getKbar_list(
        (mu_gamma - psi_gamma) / sqrt(tauSq_gamma),
        lengthScale = lengthScale_gamma, D_list = distMat_list
      )

      ### get a sample from the prior
      lengthScale_gamma_prior = rlnorm(
        p, meanlog = meanlog_lengthScales, sdlog = sqrt(varlog_lengthScales)
      )

      ### get threshold for acceptance and initial angle
      strike = listKbar_gamma$logden - rgamma(1, shape = 1, rate = 1)
      angle = runif(1, min = 0, max = 2 * pi)

      ### get initial brackets
      angle_min = angle - 2 * pi
      angle_max = angle

      ### try until acceptance occurs
      while (TRUE) {
        lengthScale_gamma_prop = exp(
          log(lengthScale_gamma) * cos(angle) +
            log(lengthScale_gamma_prior) * sin(angle)
        )

        listKbar_gamma = getKbar_list(
          (mu_gamma - psi_gamma) / sqrt(tauSq_gamma),
          lengthScale = lengthScale_gamma_prop, D_list = distMat_list
        )

        if (listKbar_gamma$logden >= strike) {
          lengthScale_gamma = lengthScale_gamma_prop
          break
        } else {
          if (angle < 0) angle_min = angle else angle_max = angle
          angle = runif(1, min = angle_min, max = angle_max)
        }
      }



      ##### log-scale of Weibull random variables

      ### get Kbar list of delta for current values
      listKbar_delta = getKbar_list(
        (mu_delta - psi_delta) / sqrt(tauSq_delta),
        lengthScale = lengthScale_delta, D_list = distMat_list
      )

      ### get a sample from the prior
      lengthScale_delta_prior = rlnorm(
        p, meanlog = meanlog_lengthScales, sdlog = sqrt(varlog_lengthScales)
      )

      ### get threshold for acceptance and initial angle
      strike = listKbar_delta$logden - rgamma(1, shape = 1, rate = 1)
      angle = runif(1, min = 0, max = 2 * pi)

      ### get initial brackets
      angle_min = angle - 2 * pi
      angle_max = angle

      ### try until acceptance occurs
      while (TRUE) {
        lengthScale_delta_prop = exp(
          log(lengthScale_delta) * cos(angle) +
            log(lengthScale_delta_prior) * sin(angle)
        )

        listKbar_delta = getKbar_list(
          (mu_delta - psi_delta) / sqrt(tauSq_delta),
          lengthScale = lengthScale_delta_prop, D_list = distMat_list
        )

        if (listKbar_delta$logden >= strike) {
          lengthScale_delta = lengthScale_delta_prop
          break
        } else {
          if (angle < 0) angle_min = angle else angle_max = angle
          angle = runif(1, min = angle_min, max = angle_max)
        }
      }



      ##### inverse logit of probability of binomial random variables

      ### get Kbar list of delta for current values
      listKbar_pi = getKbar_list(
        (mu_pi - psi_pi) / sqrt(tauSq_pi),
        lengthScale = lengthScale_pi, D_list = distMat_list
      )

      ### get a sample from the prior
      lengthScale_pi_prior = rlnorm(
        p, meanlog = meanlog_lengthScales, sdlog = sqrt(varlog_lengthScales)
      )

      ### get threshold for acceptance and initial angle
      strike = listKbar_pi$logden - rgamma(1, shape = 1, rate = 1)
      angle = runif(1, min = 0, max = 2 * pi)

      ### get initial brackets
      angle_min = angle - 2 * pi
      angle_max = angle

      ### try until acceptance occurs
      while (TRUE) {
        lengthScale_pi_prop = exp(
          log(lengthScale_pi) * cos(angle) +
            log(lengthScale_pi_prior) * sin(angle)
        )

        listKbar_pi = getKbar_list(
          (mu_pi - psi_pi) / sqrt(tauSq_pi),
          lengthScale = lengthScale_pi_prop, D_list = distMat_list
        )

        if (listKbar_pi$logden >= strike) {
          lengthScale_pi = lengthScale_pi_prop
          break
        } else {
          if (angle < 0) angle_min = angle else angle_max = angle
          angle = runif(1, min = angle_min, max = angle_max)
        }
      }





      ##########################################################################
      #                                                                        #
      #                         update scale variables                         #
      #                                                                        #
      ##########################################################################

      ### variances of white noises
      sigmaSq_gamma = 1 / rgamma(
        1, shape = k_sigmaSq + 0.5 * (M * T),
        rate = mix_sigmaSq_gamma + 0.5 * sum((mu_gamma - Gamma)^2)
      )
      sigmaSq_delta = 1 / rgamma(
        1, shape = k_sigmaSq + 0.5 * (M * T),
        rate = mix_sigmaSq_delta + 0.5 * sum((mu_delta - Delta)^2)
      )
      sigmaSq_pi = 1 / rgamma(
        1, shape = k_sigmaSq + 0.5 * (M * T),
        rate = mix_sigmaSq_pi + 0.5 * sum((mu_pi - Pi)^2)
      )

      ### variances of Gaussian processes
      tauSq_gamma = 1 / rgamma(
        1, shape = k_tauSq + 0.5 * M,
        rate = mix_tauSq_gamma + 0.5 * sum((psi_gamma - mu_gamma)^2)
      )
      tauSq_delta = 1 / rgamma(
        1, shape = k_tauSq + 0.5 * M,
        rate = mix_tauSq_delta + 0.5 * sum((psi_delta - mu_delta)^2)
      )
      tauSq_pi = 1 / rgamma(
        1, shape = k_tauSq + 0.5 * M,
        rate = mix_tauSq_pi + 0.5 * sum((psi_pi - mu_pi)^2)
      )





      ##########################################################################
      #                                                                        #
      #                       update location variables                        #
      #                                                                        #
      ##########################################################################

      ##### Weibull parameters

      ### get a sample from the prior
      # shape
      psi_gamma_prior = rnorm(1, mean = 0, sd = sqrt(mix_psi_gamma))
      mu_gamma_prior = c(
        psi_gamma_prior + sqrt(tauSq_gamma) *
          listKbar_gamma$decompKbar %*% rnorm(M)
      )
      Gamma_prior = mu_gamma_prior + sqrt(sigmaSq_gamma) *
        matrix(rnorm(M * T), nrow = M, ncol = T)
      # scale
      psi_delta_prior = rnorm(1, mean = 0, sd = sqrt(mix_psi_delta))
      mu_delta_prior = c(
        psi_delta_prior + sqrt(tauSq_delta) *
          listKbar_delta$decompKbar %*% rnorm(M)
      )
      Delta_prior = mu_delta_prior + sqrt(sigmaSq_delta) *
        matrix(rnorm(M * T), nrow = M, ncol = T)

      ### get threshold for acceptance and initial angle
      strike = loglike_weibull - rgamma(1, shape = 1, rate = 1)
      angle = runif(1, min = 0, max = 2 * pi)

      ### get initial brackets
      angle_min = angle - 2 * pi
      angle_max = angle

      ### try until acceptance occurs
      while (TRUE) {
        psi_gamma_prop = psi_gamma * cos(angle) + psi_gamma_prior * sin(angle)
        mu_gamma_prop = mu_gamma * cos(angle) + mu_gamma_prior * sin(angle)
        Gamma_prop = Gamma * cos(angle) + Gamma_prior * sin(angle)
        psi_delta_prop = psi_delta * cos(angle) + psi_delta_prior * sin(angle)
        mu_delta_prop = mu_delta * cos(angle) + mu_delta_prior * sin(angle)
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
          psi_gamma = psi_gamma_prop
          mu_gamma = mu_gamma_prop
          Gamma = Gamma_prop
          vecGamma = vecGamma_prop
          exp_vecGamma = exp_vecGamma_prop
          psi_delta = psi_delta_prop
          mu_delta = mu_delta_prop
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
      psi_pi_prior = rnorm(1, mean = 0, sd = sqrt(mix_psi_pi))
      mu_pi_prior = c(
        psi_pi_prior + sqrt(tauSq_pi) *
          listKbar_pi$decompKbar %*% rnorm(M)
      )
      Pi_prior = mu_pi_prior + sqrt(sigmaSq_pi) *
        matrix(rnorm(M * T), nrow = M, ncol = T)

      ### get threshold for acceptance and initial angle
      strike = loglike_binom - rgamma(1, shape = 1, rate = 1)
      angle = runif(1, min = 0, max = 2 * pi)

      ### get initial brackets
      angle_min = angle - 2 * pi
      angle_max = angle

      ### try until acceptance occurs
      while (TRUE) {
        psi_pi_prop = psi_pi * cos(angle) + psi_pi_prior * sin(angle)
        mu_pi_prop = mu_pi * cos(angle) + mu_pi_prior * sin(angle)
        Pi_prop = Pi * cos(angle) + Pi_prior * sin(angle)

        probs_binom = c(plogis(Pi_prop))
        probs_binom[probs_binom >= 1 - sqrt(.Machine$double.eps)] = 1 - sqrt(
          .Machine$double.eps
        )
        probs_binom[probs_binom <= sqrt(.Machine$double.eps)] = sqrt(
          .Machine$double.eps
        )

        loglike_binom_prop = sum(
          N * log(probs_binom) + (365 - N) * log(1 - probs_binom)
        )

        if (loglike_binom_prop >= strike) {
          psi_pi = psi_pi_prop
          mu_pi = mu_pi_prop
          Pi = Pi_prop
          loglike_binom = loglike_binom_prop
          break
        } else {
          if (angle < 0) angle_min = angle else angle_max = angle
          angle = runif(1, min = angle_min, max = angle_max)
        }
      }





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
    keep_mix_psi[, insim] = c(mix_psi_gamma, mix_psi_delta, mix_psi_pi)
    keep_mix_sigmaSq[, insim] = c(
      mix_sigmaSq_gamma, mix_sigmaSq_delta, mix_sigmaSq_pi
    )
    keep_mix_tauSq[, insim] = c(mix_tauSq_gamma, mix_tauSq_delta, mix_tauSq_pi)

    ### length scales
    keep_lengthScale_gamma[, insim] = lengthScale_gamma
    keep_lengthScale_delta[, insim] = lengthScale_delta
    keep_lengthScale_pi[, insim] = lengthScale_pi

    ### variances of white noises, variances and means of Gaussian processes
    keep_sigmaSq[, insim] = c(sigmaSq_gamma, sigmaSq_delta, sigmaSq_pi)
    keep_tauSq[, insim] = c(tauSq_gamma, tauSq_delta, tauSq_pi)
    keep_psi[, insim] = c(psi_gamma, psi_delta, psi_pi)

    ### Gaussian processes
    keep_mu_gamma[, insim] = mu_gamma
    keep_mu_delta[, insim] = mu_delta
    keep_mu_pi[, insim] = mu_pi

    ### if required, collect parameters (large memory consumption)
    if (details) {
      keep_Gamma[, , insim] = Gamma
      keep_Delta[, , insim] = Delta
      keep_Pi[, , insim] = Pi
    }

    ### if required, forecast parameters for other points
    if (forecasting) {
      # starKbar_gamma = getKbar_only(
      #   lengthScale_gamma, D_list = starDistMat_list
      # )
      # starKbar_delta = getKbar_only(
      #   lengthScale_delta, D_list = starDistMat_list
      # )
      # starKbar_pi = getKbar_only(
      #   lengthScale_pi, D_list = starDistMat_list
      # )

      crossKbar_gamma = getCrossKbar(
        lengthScale_gamma, crossD_list = crossDistMat_list
      )
      crossKbar_delta = getCrossKbar(
        lengthScale_delta, crossD_list = crossDistMat_list
      )
      crossKbar_pi = getCrossKbar(
        lengthScale_pi, crossD_list = crossDistMat_list
      )

      tmp = crossKbar_gamma %*% listKbar_gamma$invKbar
      keep_mu_gamma_star[, insim] = psi_gamma + tmp %*% (mu_gamma - psi_gamma) +
        sqrt(tauSq_gamma * (rep(1, Mstar) - rowSums(tmp * crossKbar_gamma))) *
          rnorm(Mstar)

      tmp = crossKbar_delta %*% listKbar_delta$invKbar
      keep_mu_delta_star[, insim] = psi_delta + tmp %*% (mu_delta - psi_delta) +
        sqrt(tauSq_delta * (rep(1, Mstar) - rowSums(tmp * crossKbar_delta))) *
          rnorm(Mstar)

      tmp = crossKbar_pi %*% listKbar_pi$invKbar
      keep_mu_pi_star[, insim] = psi_pi + tmp %*% (mu_pi - psi_pi) +
        sqrt(tauSq_pi * (rep(1, Mstar) - rowSums(tmp * crossKbar_pi))) *
          rnorm(Mstar)
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
    print(paste("sgpbw_mcmc: end time at ", Sys.time(), sep = ""))
  }




  ########## return results

  ### initialize list
  res = list()

  ### mixing variables
  res[["mix_psi"]] = keep_mix_psi
  res[["mix_sigmaSq"]] = keep_mix_sigmaSq
  res[["mix_tauSq"]] = keep_mix_tauSq

  ### length scales
  res[["lengthScale_gamma"]] = keep_lengthScale_gamma
  res[["lengthScale_delta"]] = keep_lengthScale_delta
  res[["lengthScale_pi"]] = keep_lengthScale_pi

  ### variances of white noises, variances and means of Gaussian processes
  res[["psi"]] = keep_psi
  res[["sigmaSq"]] = keep_sigmaSq
  res[["tauSq"]] = keep_tauSq

  ### Gaussian processes
  res[["mu_gamma"]] = keep_mu_gamma
  res[["mu_delta"]] = keep_mu_delta
  res[["mu_pi"]] = keep_mu_pi

  ### values of parameters
  if (details) {
    res[["Gamma"]] = keep_Gamma
    res[["Delta"]] = keep_Delta
    res[["Pi"]] = keep_Pi
  }

  ### Gaussian processes values for other points
  if (forecasting) {
    res[["mu_gamma_star"]] = keep_mu_gamma_star
    res[["mu_delta_star"]] = keep_mu_delta_star
    res[["mu_pi_star"]] = keep_mu_pi_star
  }

  ### return list of results
  return(res)

}