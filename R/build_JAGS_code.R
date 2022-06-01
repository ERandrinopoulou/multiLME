JAGSmodel <- function(families, hc, predicted, corr_RE,
                      assoc, assoc_from, assoc_to, extraForm, Data1,
                      norm_area, time_var){


    K <- length(families)

    myt <- function(times = 1) {
      tb <- "    "
      paste(rep(tb, times), collapse = "")
    }


    long_betas <- function(k){
      longbetas <- paste0("c_betas", k, "[1:ncx", k, "]", sep = "", collapse = "")
      longbetas
    }

    long_b <- function(k){
      longb <- paste0("b", "[i, RE_ind",k, "]", sep = "", collapse = "")
      longb
    }

    long_u <- function(k){
      longu <- paste0("u", k, "[i, 1:ncz", k, "]", sep = "", collapse = "")
      longu
    }

    longCont <- function(k) {
      longStart1 <- paste0(myt(1), "for (i in 1:n", k, ") {\n" )
      longStart2 <- paste0(myt(2), "for (j in offset", k, "[i]:(offset", k, "[i+1] - 1)) {\n" )
      long2 <- paste0(myt(3), "muy", k, "[j] <- inprod(", paste0(long_betas(k), collapse = " + "), ", Xc", k, "[j, 1:ncx", k, "]) +\n", myt(4),
                      " inprod(", paste0(long_b(k), sep = "", collapse = " + "), ", Z", k, "[j, 1:ncz", k, "])\n")
      if (assoc == TRUE) {

      if (k %in% assoc_from) {
        if (!is.null(extraForm)){
          paramtr <- paste0(myt(3), "f_derivY", k, "[j] <- betas", k, "[", extraForm$indFixed, "] * XderivY[j, ] + b[i, ",
                            eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indFixed]"))), "] * ZderivY[j,]\n")
        }
       }

        if (k %in% assoc_to) {

          if (!is.null(extraForm)){

            if (norm_area == TRUE){

              ind_time_var <- c(1:dim(X1)[2])[colnames(Data1$X1) == "year"]
              mu_alpha <- paste0(" * (f_derivY", assoc_from, "[j])/X", assoc_to, "[,", ind_time_var, "]")
            } else {
              mu_alpha <- paste0(" * f_derivY", assoc_from, "[j]")
            }

          } else {
            mu_alpha <- paste0(" * muy", assoc_from, "[j]")
          }

          long2 <- paste0(myt(3), "muy", k, "[j] <- inprod(", paste0(long_betas(k), collapse = " + "), ", Xc", k, "[j, 1:ncx", k, "]) +\n", myt(4),
                          " inprod(", paste0(long_b(k), sep = "", collapse = " + "), ", Z", k, "[j, 1:ncz", k, "])")
          alphast <- paste0(" + alpha", k, assoc_from, mu_alpha)
          alphast <- paste0(alphast, collapse = "")
          long2 <- paste0(long2, if(!is.null(extraForm) & (k %in% assoc_from) ){ paste0(paramtr, collapse = "") }, alphast, "\n")
        }
      }

      if (k %in% assoc_from) {
        if (!is.null(extraForm)){
          # paramtr <- paste0(myt(3), "f_derivY", k, "[j] <- betas", k, "[", extraForm$indFixed, "] * XderivY[j, ] + b[i, ",
          #                   eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indFixed]"))), "] * ZderivY[j,]\n")
          paramtr <- paste0(myt(3), "f_derivY", k, "[j] <- betas", k, "[", if (length(  extraForm$indFixed) > 1) {
            paste0("c(", paste0(extraForm$indFixed, collapse = ","), ")" )
          } else { extraForm$indFixed }, "] %*% XderivY[j, ] + b[i, ",
          if (length( eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indRandom]")))) > 1) {
            ind <- eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indRandom]")))
            paste0("c(", paste0(ind, collapse = ","), ")" )
          } else { eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indRandom]"))) }  , "] %*% ZderivY[j,]\n")

        }
      }

      long2hier <- paste0(myt(3), "muy", k, "[j] <- inprod(", paste0(long_u(k), sep = "", collapse = " + "), ", Z", k, "[j, 1:ncz", k, "])\n")

      long3 <- paste0(myt(3), "y", k, "[j] ~ dnorm(muy", k, "[j], tau", k, ")\n")
      long3_pred <- paste0(myt(3), "y", k, "_pred[j] ~ dnorm(muy", k, "[j], tau", k, ")\n")

      longEnd <- paste0(myt(1), "}\n")
      paste0(longStart1, longStart2, if (hc == "TRUE") {paste0(long2hier)} else {paste0(long2)}, if(!is.null(extraForm) & (k %in% assoc_from) ){ paste0(paramtr, collapse = "") },
             long3,
             if (predicted == TRUE) {long3_pred}, myt(2), "}\n", myt(1), "}\n")
    }

    longBin <- function(k) {
      alphast <- " + alpha * muy2[j]"
      longStart1 <- paste0(myt(1), "for (i in 1:n", k, ") {\n" )
      longStart2 <- paste0(myt(2), "for (j in offset", k, "[i]:(offset", k, "[i+1] - 1)) {\n" )
      long2 <- paste0(myt(3), "muy", k, "[j] <- inprod(", paste0(long_betas(k), collapse = " + "), ", Xc", k, "[j, 1:ncx", k, "]) +\n", myt(4),
                      " inprod(", paste0(long_b(k), sep = "", collapse = " + "), ", Z", k, "[j, 1:ncz", k, "])\n")
      if (assoc == TRUE) {

        if (k %in% assoc_from) {
          if (!is.null(extraForm)){
          #   paramtr <- paste0(myt(3), "f_derivY", k, "[j] <- betas", k, "[", extraForm$indFixed, "] * XderivY[j, ] + b[i, ",
          #                     eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indFixed]"))), "] * ZderivY[j,]\n")
            paramtr <- paste0(myt(3), "f_derivY", k, "[j] <- betas", k, "[", if (length(  extraForm$indFixed) > 1) {
              paste0("c(", paste0(extraForm$indFixed, collapse = ","), ")" )
            } else { extraForm$indFixed }, "] %*% XderivY[j, ] + b[i, ",
            if (length( eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indRandom]")))) > 1) {
              ind <- eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indRandom]")))
              paste0("c(", paste0(ind, collapse = ","), ")" )
            } else { eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indRandom]"))) }  , "] %*% ZderivY[j,]\n")

          }
        }

        if (k %in% assoc_to) {

          if (!is.null(extraForm)){

            if (norm_area == TRUE){
              ind_time_var <- c(1:dim(X1)[2])[colnames(Data1$X1) == "year"]
              mu_alpha <- paste0(" * (f_derivY", assoc_from, "[j])/X", assoc_to, "[,", ind_time_var, "]")
            } else {
              mu_alpha <- paste0(" * f_derivY", assoc_from, "[j]")
            }

          } else {
            mu_alpha <- paste0(" * muy", assoc_from, "[j]")
          }

          long2 <- paste0(myt(3), "muy", k, "[j] <- inprod(", paste0(long_betas(k), collapse = " + "), ", Xc", k, "[j, 1:ncx", k, "]) +\n", myt(4),
                          " inprod(", paste0(long_b(k), sep = "", collapse = " + "), ", Z", k, "[j, 1:ncz", k, "])")
          alphast <- paste0(" + alpha", k, assoc_from, mu_alpha)
          alphast <- paste0(alphast, collapse = "")
          long2 <- paste0(long2, if (!is.null(extraForm) & (k %in% assoc_from) ){ paste0(paramtr, collapse = "") }, alphast, collapse = " ")
        }
      }

      if (k %in% assoc_from) {
        if (!is.null(extraForm)){
          paramtr <- paste0(myt(3), "f_derivY", k, "[j] <- betas", k, "[", extraForm$indFixed, "] * XderivY[j, ] + b[i, ",
                            eval(parse(text = paste0("Data1$RE_ind", k, "[extraForm$indFixed]"))), "] * ZderivY[j,]\n")
        }
      }

      long2hier <- paste0(myt(3), "muy", k, "[j] <- inprod(", paste0(long_u(k), sep = "", collapse = "+ "), ", Z", k, "[j, 1:ncz", k, "] )\n")

      long3 <- paste0(myt(3), "Pr", k, "[j] <- max(1.00000E-05, min(0.99999, (exp(muy", k, "[j])/(1 + exp(muy", k, "[j])))))\n")
      long4 <- paste0(myt(3), "y", k, "[j] ~ dbin(Pr", k, "[j], 1)\n")
      long4_pred <- paste0(myt(3), "y", k, "_pred[j] ~ dbin(Pr", k, "[j], 1)\n")
      paste0(longStart1, longStart2, if (hc == "TRUE") {paste0(long2hier)} else {paste0(long2)}, if(!is.null(extraForm) & (k %in% assoc_from) ){ paste0(paramtr, collapse = "") },
             long3, long4,
             if (predicted == TRUE) {long4_pred}, myt(2), "}\n", myt(1), "}\n")
    }


    long <- NULL
    for (k in 1:K){
      long[k] <- paste0(if (families[k] == "gaussian") longCont(k) else if (families[k] == "binomial") longBin(k))
    }


    longNew <- paste0(long, collapse = "")

    longPar <- paste0(longNew, collapse = "")


    prior.b <- paste0(myt(2), "b", "[i, ", "1:n_RE] ~ dmnorm(mu0", "[], inv.D", "[, ])\n",  myt(), "} \n")


    prior.b <- paste0(prior.b, collapse = "")
    beg <- paste0(myt(1), "for (i in 1:n) {\n",  sep = "", collapse = "")
    b <- paste0(beg, prior.b, collapse = "")

    # #### HC
    # if (hc == TRUE) {
    #
    #   hierlong <- "string"
    #
    #   for (j in seq_along(colmns_HC)) {
    #     ii <- colmns_HC[[j]]
    #     j.incr <- j + 0
    #     longbetasN <- "string"
    #
    #     if (length(ii) > 1) {
    #       ind.cl <- paste0("c(", paste(ii, collapse = ", "), ")")
    #       for (i in 1:classes) {
    #         longbetasN[i] <- paste0("equals(v[i], ", i, ") * betas", i, "[", ind.cl, "]", sep = "", collapse = "")
    #       }
    #
    #       hierlong[j] <- paste0(myt(1), "mu.u[i, ", j.incr, "] <- inprod((",  paste0(longbetasN, collapse = " + "),
    #                             "), Xhc[i, ", ind.cl, "])\n", sep = "", collapse = "")
    #     } else {
    #       ind.cl <- ii
    #       for (i in 1:classes) {
    #         longbetasN[i] <- paste0("equals(v[i], ", i, ") * betas", i, "[", ind.cl, "]", sep = "", collapse = "")
    #       }
    #
    #       hierlong[j] <- paste0(paste0(myt(1), "mu.u[i, ", j.incr, "] <- ", paste0(longbetasN, collapse = " + "), "\n",
    #                                    sep = "", collapse = ""))
    #     }
    #   }
    #
    #
    #   uf <- function(j){
    #     u <- ""
    #     for (i in 1:classes) {
    #       u[i] <- paste0("equals(v[i], ", i, ") * u", i, "[i, ", j, "]" )
    #     }
    #     paste0(u)
    #   }
    #
    #   prior.b.hier <- ""
    #   for (i in 1:classes) {
    #     for (j in seq_along(colmns_HC)) {
    #       prior.b.hier <- paste0(prior.b.hier, myt(), "b", i, "[i, ", j, "] <- (",  paste0(uf(j), collapse = " + "), ") - mu.u[i, ", j, "]\n", sep = "", collapse = "")
    #     }
    #   }
    #
    #
    #   prior.b.hier2 <- "string"
    #   for (i in 1:classes) {
    #     prior.b.hier2[i] <- paste0(myt(), "u", i, "[i, ] ~ dmnorm(mu.u[i, ], inv.D", i, "[, ]) \n")
    #   }
    #
    # }
    #



    # prior.betas <- "string"
    # m = 1
    # for (k in 1:K) {
    #   for (i in 1:classes) {
    #     prior.betas[m] <- paste0(myt(), "betas", k, i, "[1:(ncx", k, ")] ~ dmnorm(priorMean.betas", k, i, "[], priorTau.betas", k, i, "[, ])\n")
    #     m <- m + 1
    #   }
    # }

    prior.betas <- "string"
    m = 1
    for (k in 1:K) {
        prior.betas[m] <- paste0(myt(), "for (k", k, " in 1:ncx", k, ") {\n", myt(2), "c_betas", k, "[k", k, "] ~ dnorm(0.0, tau_betas", k, ")\n", myt(), "} \n")
        m <- m + 1
    }


    if (assoc == TRUE) {
      prior.alpha <- "string"
      m = 1
      for (k in assoc_to) {
        for (u in assoc_from){
          prior.alpha[m] <- paste0(myt(), "alpha", k, u, " ~ dnorm(0.0, tau_alpha)", "\n")
          m <- m + 1
        }
      }
    }

    prior.tau <- "string"
    m = 1
    for (k in 1:K) {
      if (families[k] == "gaussian") {
        prior.tau[m] <- paste0(myt(), "tau", k, " ~ dgamma(priorA.tau", k, ", priorB.tau", k, ")\n")
        m <- m + 1
      }
    }
    if (any(families == "gaussian")) prior.tau <- paste0(prior.tau, collapse = "") else prior.tau <- NULL

    prior.sigma <- "string"
    m = 1
    for (k in 1:K) {
      if (families[k] == "gaussian") {
        prior.sigma[m] <- paste0(myt(), "sigma", k , " <- 1/sqrt(tau", k, ")\n")
        m <- m + 1
      }
    }
    if (any(families == "gaussian")) prior.sigma <- paste0(prior.sigma, collapse = "") else prior.tau <- NULL

    prior.invD <- paste0(myt(), "inv.D[1:n_RE, 1:n_RE] ~ dwish(4*priorR.D[, ], priorK.D)\n",
                         myt(), "for (l in 1:n_RE) {\n",
                         myt(2), "priorR.D[l, l] ~ dgamma(0.5, 0.01)\n",
                         myt(), "}\n")

    prior.invD_uncorr <- paste0(myt(), "for (l in 1:n_RE) {\n",
                         myt(2), "inv.D[l, l] ~ dgamma(0.01, 0.01)\n",
                         myt(), "}\n")

    if (assoc == TRUE){
      priors <- function(k) {
        prior.betas <- paste0(prior.betas, collapse = "")
        prior.alpha <- paste0(prior.alpha, collapse = "")
        prior.invD <-  if (corr_RE == TRUE) {paste0(prior.invD, collapse = "")
        } else {paste0(prior.invD_uncorr, collapse = "")}
        paste0(prior.betas, prior.alpha, prior.tau, prior.sigma, prior.invD)
      }
    } else {
      priors <- function(k) {
        prior.betas <- paste0(prior.betas, collapse = "")
        prior.invD <-  if (corr_RE == TRUE) {paste0(prior.invD, collapse = "")
        } else {paste0(prior.invD_uncorr, collapse = "")}
        paste0(prior.betas, prior.tau, prior.sigma, prior.invD)
      }
    }


    for (k in 1:K){
      priorsPar <- paste(priors(k), collapse = "")
    }

    prior.backTransform_betas <- "string"
    m = 1
    for (k in 1:K) {
        prior.backTransform_betas[m] <- paste0(myt(1), "betas", k, "[1] = c_betas", k, "[1] - inprod(c_betas", k, "[2:ncx", k, "], means_X", k, ") \n",
                                               myt(1), "betas", k, "[2:ncx", k, "] = c_betas", k, "[2:ncx", k, "] \n")
        m <- m + 1

    }
    prior.backTransform_betas <- paste0(prior.backTransform_betas, collapse = "")

    model <- function() {
      end <- paste0("}\n", sep = "", collapse = "")
      paste0("model { \n", longPar, b,
             priorsPar, prior.backTransform_betas, end, collapse = "")
    }


    model <- model()


    filename <- file.path("mixedmodel.txt")

    write.table(model, filename, append = FALSE, quote = FALSE, sep = " ",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE, qmethod = c("escape", "double"))


}


