# matern1.5_kernel = function(x, y, length_scale) {
#   tmp = sqrt(3) * abs(x - y) / length_scale
#   return(c(
#     (1 + tmp) * exp(-tmp)
#   ))
# }

# matern1.5_corrMat = function(x, length_scale) {
#   d = length(x)
#   corrMat = matrix(0, nrow = d, ncol = d)
#   for (i in 1:(d - 1)) {
#     corrMat[(i + 1):d, i] = matern1.5_kernel(
#       x = x[(i + 1):d], y = x[i], length_scale = length_scale
#     )
#   }
#   corrMat[upper.tri(corrMat)] = t(corrMat)[upper.tri(corrMat)]
#   diag(corrMat) = rep(1, d)
#   return(corrMat)
# }

# matern1.5_crossCorrMat = function(x, y, length_scale) {
#   dx = length(x)
#   dy = length(y)
#   crossCorrMat = matrix(0, nrow = dx, ncol = dy)
#   for (i in 1:dy) {
#     crossCorrMat[, i] = matern1.5_kernel(
#       x = x, y = y[i], length_scale = length_scale
#     )
#   }
#   return(crossCorrMat)
# }

chol_w_pert = function(x, n = dim(x)[1]) {
  c = 0
  while (TRUE) {
    res = tryCatch(
      chol(x + c * diag(sqrt(.Machine$double.eps), nrow = n)),
      error = function(err) NA
    )
    c = c + 1
    if (!is.na(res[1])) break
  }
  return(res)
}

unidist = function(x) {
  d = length(x)
  D = matrix(0, nrow = d, ncol = d)
  for (i in 1:(d - 1)) {
    D[(i + 1):d, i] = abs(x[(i + 1):d] - x[i])
  }
  D[upper.tri(D)] = t(D)[upper.tri(D)]
  diag(D) = rep(0, d)
  return(D)
}

cross_unidist = function(x, y) {
  dx = length(x)
  dy = length(y)
  D = matrix(0, nrow = dx, ncol = dy)
  for (i in 1:dx) {
    D[i, ] = abs(x[i] - y)
  }
  return(D)
}

matern1.5_unidist2corrMat = function(D, length_scale) {
  corrMat = D
  tmp = sqrt(3) * D[lower.tri(D)] / length_scale
  tmp = (1 + tmp) * exp(-tmp)
  corrMat[lower.tri(D)] = tmp
  diag(corrMat) = 1
  corrMat[!lower.tri(D)] = t(corrMat)[!lower.tri(D)]
  return(corrMat)
}

matern1.5_crossUnidist2crossCorrMat = function(D, length_scale) {
  tmp = sqrt(3) * D / length_scale
  crossCorrMat = (1 + tmp) * exp(-tmp)
  return(crossCorrMat)
}

matern1.5_unidist2derivFactor = function(D, length_scale) {
  3 * D^2 / (length_scale^3 + sqrt(3) * D * length_scale^2)
}

# dinvLomax = function(x, alpha, beta) {
#   x^(alpha - 1) * alpha * beta^alpha / (1 + x * beta)^(alpha + 1)
# }

hyperTuner_corrMatern1.5 = function(x, D_list) {

  p = length(D_list)
  n = length(x)
  logden_curr = -Inf
  # lengthScales_current = rep(exp(2), p)
  # for (i in 1:p) {
  #   lengthScales_current[i] = median(-sqrt(3) / log(0.5) * D_list[[i]])
  # }
  lengthScales_current = rep(1, p)
  lengthScales_gradient = -2 / (1 + lengthScales_current)
  step_size = 0
  iterations_counter = -1

  while (TRUE) {

    lengthScales_proposal = lengthScales_current +
      step_size * lengthScales_gradient
    lengthScales_proposal[
      lengthScales_proposal <= sqrt(.Machine$double.eps)
    ] = sqrt(.Machine$double.eps)

    Kbar = matrix(1, nrow = n, ncol = n)
    for (i in 1:p) {
      Kbar = Kbar * matern1.5_unidist2corrMat(
        D = D_list[[i]], length_scale = lengthScales_proposal[i]
      )
    }

    decompKbar = chol_w_pert(Kbar, n = n)
    logdetKbar = 2 * sum(log(diag(decompKbar)))
    decomp_invKbar = backsolve(decompKbar, x = diag(1, nrow = n))
    # decompKbar = tryCatch(
    #   chol(Kbar),
    #   error = function(err) svd(Kbar, nv = 0)
    # )
    # if (is.matrix(decompKbar)) {
    #   logdetKbar = 2 * sum(log(diag(decompKbar)))
    #   decomp_invKbar = backsolve(decompKbar, x = diag(1, nrow = n))
    # } else {
    #   logdetKbar = sum(log(decompKbar$d))
    #   decomp_invKbar = t(1 / sqrt(decompKbar$d) * t(decompKbar$u))
    # }

    logden_prop = -0.5 * tcrossprod(x %*% decomp_invKbar) - 0.5 * logdetKbar

    if (logden_prop >= logden_curr) {
      logden_curr = logden_prop
      lengthScales_current = lengthScales_proposal
      step_size = 1
      iterations_counter = iterations_counter + 1

      invKbar = tcrossprod(decomp_invKbar)
      for (i in 1:p) {
        lengthScales_gradient[i] = sum(
          (tcrossprod(invKbar %*% x) - invKbar) *
            t(Kbar * matern1.5_unidist2derivFactor(
              D = D_list[[i]], length_scale = lengthScales_current[i]
            ))
        )
      }
    } else {
      step_size = 0.5 * step_size
    }

    if (step_size * sum(lengthScales_gradient^2) <= .Machine$double.eps) break
    if (iterations_counter == 1000) break

  }

  return(lengthScales_current)

}

getKbar_list = function(x, lengthScale, D_list) {

  if (any(lengthScale <= 0)) return(list(logden = -Inf))

  p = length(D_list)
  n = length(x)

  Kbar = matrix(1, nrow = n, ncol = n)
  for (i in 1:p) {
    Kbar = Kbar * matern1.5_unidist2corrMat(
      D = D_list[[i]], length_scale = lengthScale[i]
    )
  }

  decompKbar = chol_w_pert(Kbar, n = n)
  logdetKbar = 2 * sum(log(diag(decompKbar)))
  decomp_invKbar = backsolve(decompKbar, x = diag(1, nrow = n))
  # decompKbar = tryCatch(chol(Kbar), error = function(err) svd(Kbar, nv = 0))
  # if (is.matrix(decompKbar)) {
  #   logdetKbar = 2 * sum(log(diag(decompKbar)))
  #   decomp_invKbar = backsolve(decompKbar, x = diag(1, nrow = n))
  # } else {
  #   logdetKbar = sum(log(decompKbar$d))
  #   decomp_invKbar = t(1 / sqrt(decompKbar$d) * t(decompKbar$u))
  #   decompKbar = t(sqrt(decompKbar$d) * t(decompKbar$u))
  # }

  invKbar = tcrossprod(decomp_invKbar)
  logden = -0.5 * tcrossprod(x %*% decomp_invKbar) - 0.5 * logdetKbar

  return(list(
    decompKbar = decompKbar,
    logdetKbar = logdetKbar,
    invKbar = invKbar,
    decomp_invKbar = decomp_invKbar,
    logden = logden
  ))

}

# getKbar_only = function(lengthScale, D_list) {
#   p = length(D_list)
#   n = nrow(D_list[[1]])

#   Kbar = matrix(1, nrow = n, ncol = n)
#   for (i in 1:p) {
#     Kbar = Kbar * matern1.5_unidist2corrMat(
#       D = D_list[[i]], length_scale = lengthScale[i]
#     )
#   }

#   return(Kbar)
# }

getCrossKbar = function(lengthScale, crossD_list) {
  p = length(crossD_list)
  dx = nrow(crossD_list[[1]])
  dy = ncol(crossD_list[[1]])

  crossKbar = matrix(1, nrow = dx, ncol = dy)
  for (i in 1:p) {
    crossKbar = crossKbar * matern1.5_crossUnidist2crossCorrMat(
      D = crossD_list[[i]], length_scale = lengthScale[i]
    )
  }

  return(crossKbar)
}

# rescale = function(
#   x, compute_coeff = TRUE, return_coeff = FALSE, loc = NULL, scale = NULL
# ) {
#   if (!is.numeric(x)) {
#     x = as.numeric(factor(x))
#     x_obs = x[!is.na(x)]
#   }
#   if (compute_coeff) {
#     loc = mean(x_obs)
#     scale = sd(x_obs)
#   }

#   if (!return_coeff) {
#     return((x - loc) / scale)
#   } else {
#     return(list(
#       x = (x - loc) / scale,
#       loc = loc,
#       scale = scale
#     ))
#   }
# }