log_post_b_EB <- function (b_i, y_i, X_i, Z_i, Sigma,
                        betas, invD, alphas, log_dens, Data = Data, assoc = assoc,
                        extraForm, XderivY_i, ZderivY_i) {

  log_post_y <- list()
  eta_y <- list()
  deriv_eta_y <- list()
  for (o in 1:length(X_i)){
    #marg_part <- paste0("prob_class", seq_len(classes), " * (X_i[[o]] %*% betas$betas", o, seq_len(classes), ")", collapse = " + ")
    eta_y[[o]] <- as.vector(X_i[[o]] %*% eval(parse(text = paste0("betas$betas", o))) +
                       Z_i[[o]] %*% b_i[eval(parse(text = paste0("Data$RE_ind", o, "_pred")))])

    if (!is.null(extraForm)){
      if (o %in% assoc_from) {
        # deriv_eta_y[[o]] <- as.vector(XderivY_i$deriv_fixed %*% eval(parse(text = paste0("betas$betas", o, "[", extraForm$indFixed, "]" ))) +
        #                               ZderivY_i$deriv_random %*% b_i[eval(parse(text = paste0("Data$RE_ind", o, "[extraForm$indFixed]")))])

        indtest <- if (length( eval(parse(text = paste0("Data$RE_ind", o, "[extraForm$indRandom]")))) > 1) {
          ind <- eval(parse(text = paste0("Data$RE_ind", o, "[extraForm$indRandom]")))
          paste0("c(", paste0(ind, collapse = ","), ")" )
        } else { eval(parse(text = paste0("Data$RE_ind", o, "[extraForm$indRandom]"))) }
        deriv_eta_y[[o]] <- as.vector(XderivY_i$deriv_fixed %*% eval(parse(text = paste0("betas$betas", o, "[", if (length(  extraForm$indFixed) > 1) {
          paste0("c(", paste0(extraForm$indFixed, collapse = ","), ")" )
        } else { extraForm$indFixed }, "]" ))) +
                                     ZderivY_i$deriv_random %*% b_i[ eval(parse(text = indtest))  ])

      }
    }


  }


  if (assoc == TRUE){
    if (!is.null(extraForm)){
      for (o in assoc_from){
        eta_y[[assoc_to]] <- eta_y[[assoc_to]] + deriv_eta_y[[o]] * unlist(alphas)
      }
    } else{
      for (o in assoc_from){
        eta_y[[assoc_to]] <- eta_y[[assoc_to]] + eta_y[[o]] * unlist(alphas)
      }
    }
  }

  for (o in 1:length(X_i)){
    log_post_y[[o]] <- log_dens[[o]](y_i[[o]], eta_y[[o]], Sigma = Sigma[[o]])

      #c(0.5 * crossprod(b_i, invD) %*% b_i)
  }

  unique(-(sum(unlist(log_post_y)) -
    c(0.5 * crossprod(b_i, invD) %*% b_i)))

}



log_post_b_MCMC <- function (b_i, y_i, X_i, Z_i, Sigma,
                        betas, invD, alphas, log_dens, Data, assoc = assoc,
                        extraForm, XderivY_i, ZderivY_i) {
  log_post_y <- list()
  eta_y <- list()
  deriv_eta_y <- list()
  for (o in 1:length(X_i)){
    #marg_part <- paste0("prob_class", seq_len(classes), " * (X_i[[o]] %*% betas$betas", o, seq_len(classes), ")", collapse = " + ")
    eta_y[[o]] <- as.vector(X_i[[o]] %*% eval(parse(text = paste0("betas$betas", o))) +
                              Z_i[[o]] %*% b_i[eval(parse(text = paste0("Data$RE_ind", o, "_pred")))])

    if (!is.null(extraForm)){
      if (o %in% assoc_from) {
        # deriv_eta_y[[o]] <- as.vector(XderivY_i$deriv_fixed %*% eval(parse(text = paste0("betas$betas", o, "[", extraForm$indFixed, "]" ))) +
        #                                 ZderivY_i$deriv_random %*% b_i[eval(parse(text = paste0("Data$RE_ind", o, "[extraForm$indFixed]")))])
        indtest <- if (length( eval(parse(text = paste0("Data$RE_ind", o, "[extraForm$indRandom]")))) > 1) {
          ind <- eval(parse(text = paste0("Data$RE_ind", o, "[extraForm$indRandom]")))
          paste0("c(", paste0(ind, collapse = ","), ")" )
        } else { eval(parse(text = paste0("Data$RE_ind", o, "[extraForm$indRandom]"))) }
        deriv_eta_y[[o]] <- as.vector(XderivY_i$deriv_fixed %*% eval(parse(text = paste0("betas$betas", o, "[", if (length(  extraForm$indFixed) > 1) {
          paste0("c(", paste0(extraForm$indFixed, collapse = ","), ")" )
        } else { extraForm$indFixed }, "]" ))) +
          ZderivY_i$deriv_random %*% b_i[ eval(parse(text = indtest))  ])

      }
    }


  }


  if (assoc == TRUE){
    if (!is.null(extraForm)){
      for (o in assoc_from){
        eta_y[[assoc_to]] <- eta_y[[assoc_to]] + deriv_eta_y[[o]] * unlist(alphas)
      }
    } else{
      for (o in assoc_from){
        eta_y[[assoc_to]] <- eta_y[[assoc_to]] + eta_y[[o]] * unlist(alphas)
      }
    }
  }

  for (o in 1:length(X_i)){
    log_post_y[[o]] <- log_dens[[o]](y_i[[o]], eta_y[[o]], Sigma = Sigma[[o]])

  }

  (sum(unlist(log_post_y)) -
      c(0.5 * crossprod(b_i, invD) %*% b_i))
}



# dmvnorm <-
#   function (x, mu, Sigma = NULL, invSigma = NULL, log = FALSE, prop = TRUE) {
#     if (!is.matrix(x))
#       x <- rbind(x)
#     if (is.matrix(mu)) {
#       if (nrow(mu) != nrow(x))
#         stop("incorrect dimensions for 'mu'.")
#       p <- ncol(mu)
#     } else {
#       p <- length(mu)
#       mu <- rep(mu, each = nrow(x))
#     }
#     if (is.null(Sigma) && is.null(invSigma))
#       stop("'Sigma' or 'invSigma' must be given.")
#     if (!is.null(Sigma)) {
#       if (is.list(Sigma)) {
#         ev <- Sigma$values
#         evec <- Sigma$vectors
#       } else {
#         ed <- eigen(Sigma, symmetric = TRUE)
#         ev <- ed$values
#         evec <- ed$vectors
#       }
#       invSigma <- evec %*% (t(evec) / ev)
#       if (!prop)
#         logdetSigma <- sum(log(ev))
#     } else {
#       if (!prop)
#         logdetSigma <- - determinant(as.matrix(invSigma))$modulus
#     }
#     ss <- x - mu
#     quad <- 0.5 * rowSums((ss %*% invSigma) * ss)
#     if (!prop)
#       fact <- - 0.5 * (p * log(2 * pi) + logdetSigma)
#     if (log) {
#       if (!prop) as.vector(fact - quad) else as.vector(- quad)
#     } else {
#       if (!prop) as.vector(exp(fact - quad)) else as.vector(exp(- quad))
#     }
#   }

#object$family$linkinv
#function (eta)
#  .Call(C_logit_linkinv, eta)
#<environment: namespace:stats>
# test_fun<- object$family$linkinv
# > object$family$linkinv(eta)
# 866       867       868
# 0.1320198 0.1424191 0.1529611
# > exp(-1.883217)/(1 + exp(-1.883217))
# [1] 0.1320198
# >


binomial_log_dens = function (y, eta, Sigma = 0) {
  mu_y <- exp(eta) / (1 + exp(eta)) #mu_fun(eta)
  y <- as.numeric(unclass(y) - 1)
  out <- if (NCOL(y) == 2) {
    dbinom(y[, 1], y[, 1] + y[, 2], mu_y, TRUE)
  } else {
    dbinom(y, 1, mu_y, TRUE)
  }
  #attr(out, "mu_y") <- mu_y
  out
}

normal_log_dens = function (y, eta, Sigma) {
  mu_y <- (eta)
  out <- dnorm(y, mu_y, Sigma, log = TRUE)
  #attr(out, "mu_y") <- mu_y
  out
}
