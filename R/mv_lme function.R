mv_lme <- function(formulas, data, families, hc,
                   predicted, engine = "JAGS",
                   parameters.b = FALSE,
                   assoc = TRUE,
                   assoc_from = 2,
                   assoc_to = 1,
                   control = list(), ...){

  ########
  # Data #
  ########

  components <- lapply(unname(formulas), extractFrames, data = data)

  components <- unlist(components, recursive = FALSE)
  n_outcomes <- length(formulas)
  names(components) <- paste0(names(components),
                              rep(seq_len(n_outcomes),
                                  each = length(components) / n_outcomes))
  colmns_HC <- components[grep("colmns_HC", names(components), fixed = TRUE)]
  colmns_nHC <- components[grep("colmns_nHC", names(components), fixed = TRUE)]
  seq_outcomes <- seq_len(n_outcomes)

  nams_vars <- c("n", "offset", "Z", "X", "Xc", "means_X", "ncx", "ncz", "y")

    #
    # "N", "id", "X", "Z", "Z_", "Zinv", "Zv", "Ztinv", "X", "Xhc", "ncx", "y", "n", "offset",
    #              "ZrowsStart", "ZrowsEnd","Xc", "Xs", "Zc", "Zs", "XhcC", "XhcS",
    #              "means_X", "SDs_X", "mean_sd_X", "means_Z", "SDs_Z", "mean_sd_Z",
    #              "means_Xhc", "SDs_Xhc", "mean_sd_Xhc", "colmns_nHC")
  vars <- paste0(rep(nams_vars, each = n_outcomes), seq_outcomes)
  #if (any(ind_td <- sapply(colmns_nHC, length))) {
  #  vars <- c(vars, paste0("X", which(ind_td > 0)))
  #}

  Data_data <- c(list(n = components$n1), components[vars])

  Data_data$n_RE <- sum(unlist(components[grep("ncz", names(components), fixed = TRUE)]))
  RE_inds <- mapply(function (sq, incr) seq_len(sq) + incr,
                    sq = components[grep("ncz", names(components), fixed = TRUE)],
                    incr = cumsum(c(0, head(sapply(colmns_HC, length), -1))),
                    SIMPLIFY = FALSE)
  names(RE_inds) <- paste0("RE_ind", seq_along(RE_inds))
  Data1 <- c(Data_data, RE_inds)

  ##########
  # Priors #
  ##########

  K <- length(families)

  datamu <- paste0("mu0", " = mu", sep = "", collapse = "")

  datamu <- paste0(datamu, sep = "", collapse = ", ")

  ncZ <- Data1$n_RE



  datainvD <- paste0("priorR.D = priorR.D,
                      priorK.D = (ncZ + 1)")

  databet <- "string"
  m <- 1
  for (k in 1:K) {
      databet[m] <- paste0("priorMean.betas", k, " = betas", k, sep = "", collapse = "")
      m <- m + 1
  }
  databet <- paste0(databet, sep = "", collapse = ", ")

  datatau <- "string"
  if (any(families == "gaussian")){
    K_gaus <- which(families == "gaussian")
    m <- 1
    for (k in K_gaus) {
      datatau[m] <- paste0("priorA.tau", k, " = priorA.tau", k, ", priorB.tau", k, " = priorB.tau", k, sep = "", collapse = "")
      m <- m + 1
    }
    datatau <- paste0(datatau, sep = "", collapse = ", ")
  }

  # CGANCE
  # databetVar <- "string"
  # m <- 1
  # for (k in 1:K) {
  #   for (i in 1:classes) {
  #     databetVar[m] <- paste0("priorTau.betas", k, i, " = diag(1/var.betas", k, i, ")", sep = "", collapse = "")
  #     m <- m + 1
  #   }
  # }
  # databetVar <- paste0(databetVar, sep = "", collapse = ", ")

  databetVar <- "string"
  m <- 1
  for (k in 1:K) {
      databetVar[m] <- paste0("tau_betas", k, " = 1/var.betas", k, sep = "", collapse = "")
      m <- m + 1
  }
  databetVar <- paste0(databetVar, sep = "", collapse = ", ")


  if (assoc == TRUE){
    dataalphaVar <- paste0("tau_alpha", " = 1/var.alpha", sep = "", collapse = "")
  }


  nZ <- Data1[grep("ncz", names(Data1), fixed = TRUE)]
  mu <- rep(0, Data1$n_RE)
  mu <- paste0("c(", paste0(mu, collapse = ","), ")")
  mu_s <- paste0("mu", " <- ", mu)
  mu_s <- paste0(mu_s, sep = "; ", collapse = "")
  eval(parse(text = mu_s))

  priorD <- "NA"
  priorD_s <- paste0("priorD", " <- ", priorD)
  priorD_s <- paste0(priorD_s, sep = "; ", collapse = "")
  eval(parse(text = priorD_s))

  priorR.D = matrix(0, ncZ, ncZ)
  diag(priorR.D) = priorD

  nX <- Data1[grep("ncx", names(Data1), fixed = TRUE)]
  betas <- lapply(nX, function(x) rep(0, x))
  betas <- lapply(betas, function(x) paste0("c(", paste0(x, collapse = ","), ")"))
  args <- c(1:K)
  betas_s <- mapply(x = args, z = betas, function(x, z) paste0("betas", x, " <- ", z))
  betas_s <- paste0(betas_s, sep = "; ", collapse = "")
  eval(parse(text = betas_s))

  # CHANGE
  # var.betas <- lapply(nX, function(x) rep(100, x))
  # var.betas <- lapply(var.betas, function(x) paste0("c(", paste0(x, collapse = ","), ")"))
  # var.betas <- rep(var.betas, classes)
  # args <- expand.grid(x = c(1:K), y = c(1:classes))
  # var.betas_s <- mapply(x = args$x, y = args$y, z = var.betas, function(x, y, z) paste0("diag(var.betas", x, y, " <- ", z, ")"))
  # var.betas_s <- paste0(var.betas_s, sep = "; ", collapse = "")
  # eval(parse(text = var.betas_s))

  args <- expand.grid(x = c(1:K))
  var.betas_s <- mapply(x = args, function(x, z) paste0("var.betas", x, " <-  100"))
  var.betas_s <- paste0(var.betas_s, sep = "; ", collapse = "")
  eval(parse(text = var.betas_s))

  if (assoc == TRUE) {
    var.alpha <-  100
  }


  if (any(families == "gaussian")){
    K_gaus <- which(families == "gaussian")
    priorA.tau <- 0.01
    args <- expand.grid(x = c(K_gaus))
    priorA.tau_s <- lapply(args$x, function(x) paste0("priorA.tau", x, " <- ", priorA.tau))
    priorA.tau_s <- paste0(priorA.tau_s, sep = "; ", collapse = "")
    eval(parse(text = priorA.tau_s))
  }

  if (any(families == "gaussian")){
    K_gaus <- which(families == "gaussian")
    priorB.tau <- 0.01
    args <- expand.grid(x = c(K_gaus))
    priorB.tau_s <- mapply(x = args$x,function(x) paste0("priorB.tau", x, " <- ", priorB.tau))
    priorB.tau_s <- paste0(priorB.tau_s, sep = "; ", collapse = "")
    eval(parse(text = priorB.tau_s))
  }



  if (assoc == TRUE) {
    data_priors <- paste0("Data2 <- list(", datamu, ",\n", datainvD, ",\n",
                        databetVar, ",\n", dataalphaVar,
                        if (any(families == "gaussian")) { paste0(",\n", datatau) }, ")")
  } else {
    data_priors <- paste0("Data2 <- list(", datamu, ",\n", datainvD, ",\n",
                          databetVar,
                          if (any(families == "gaussian")) { paste0(",\n", datatau) }, ")")
  }
  eval(parse(text = data_priors))

  Data <- c(Data1, Data2)

  ########
  # MCMC #
  ########

  con <- list(n.processors = parallel::detectCores() - 1,
              working.directory = getwd(), clear.model = TRUE,
              seed = 1L, optimize_only = FALSE, verbose = FALSE,
              n.iter = 28000L, n.burnin = 3000L, n.thin = 50L,
              n.adapt = 3000L, n.chains = 2L, seed = 1L, n.cores = 1L)



  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (!any(namc == "n.thin")) {
    con$n.thin <- if (engine == "JAGS") {
      max(1, floor((con$n.iter - con$n.burnin) * con$n.chains / 1000))
    } else {
      max(1, floor((con$n.iter - con$n.warmup) * con$n.chains / 1000))
    }
  }
  if (length(noNms <- namc[!namc %in% namC]) > 0)
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  ######################
  # Parameters to save #
  ######################

  args <- c(1:n_outcomes)
  params_betas <- paste0('betas', paste0(args))
  if (any(ind_gs <- sapply(families, function (x) x == "gaussian"))) {
    params_betas <- c(params_betas, paste0("tau", which(ind_gs)), paste0("sigma", which(ind_gs)))
  }
  params_b <- paste0('b')
  params_invD <- paste0('inv.D')

  if (assoc == TRUE){
    params_alpha <- paste0('alpha', assoc_to, assoc_from)
  }

  params_pred <- paste0('y', seq_len(n_outcomes), '_pred')


  if (assoc == TRUE){
    if (parameters.b == TRUE){
      params <- c(params_betas, params_b,
                  params_alpha, params_invD,
                  if (predicted == TRUE) {params_pred})
    } else {
      params <- c(params_betas, params_alpha, params_invD,
                  if (predicted == TRUE) {params_pred})
    }
  } else {
    if (parameters.b == TRUE){
      params <- c(params_betas, params_b,
                  params_invD,
                  if (predicted == TRUE) {params_pred})
    } else {
      params <- c(params_betas, params_invD,
                  if (predicted == TRUE) {params_pred})
    }
  }
  ############
  # Initials #
  ############

  # inits <- function () {
  #   ints <- lapply(components[grep("ncx", names(components), fixed = TRUE)],
  #                  rnorm, sd = 0.1)
  #   names(ints) <- paste0('betas', seq_len(n_outcomes))
  #   ints$u <- drop(matrix(rnorm(Data$n * Data$n_RE), Data$n,
  #                         Data$n_RE))
  #   ints$inv_D <- ints$D <- if (Data$n_RE > 1) diag(Data$n_RE) else 1
  #   if (any(ind_gs)) {
  #     nms <- which(ind_gs)
  #     taus <- rep(list(1), length(nms))
  #     names(taus) <- paste0("tau", nms)
  #     ints <- c(ints, taus)
  #   }
  #   ints
  # }

  ####################
  # Build JAGS model #
  ####################

  JAGSmodel(families = families, hc, predicted,
            assoc, assoc_from, assoc_to)

  #############
  # Fit model #
  #############


  model_name <- "mixedmodel.txt"
  # Data$priorMeanbetas1 <- rep(0, Data$ncx1)
  # Data$priorMeanbetas2 <- rep(0, Data$ncx2)
  # Data$priorMeanbetas3 <- rep(0, Data$ncx3)
  # Data$priorMeanbetas4 <- rep(0, Data$ncx4)
  # Data$priorTaubetas1 <- diag(0.01, Data$ncx1)
  # Data$priorTaubetas2 <- diag(0.01, Data$ncx2)
  # Data$priorTaubetas3 <- diag(0.01, Data$ncx3)
  # Data$priorTaubetas4 <- diag(0.01, Data$ncx4)

  # # centering
  # for(o in 1:n_outcomes){
  #   code_cen <- paste0("Data$X", o, "_mean <- as.vector(apply(Data$X", o, ", 2, mean))")
  #   eval(parse(text = code_cen))
  # }

  fit <- jags(data = Data, parameters.to.save = params, #inits = inits,
                  model.file = file.path(con$working.directory, model_name),
                  n.chains = con$n.chains, parallel = TRUE, #con$n.processors > 1,
                  n.cores = con$cores,
                  n.adapt = con$n.adapt, n.iter = con$n.iter, n.burnin = con$n.burnin,
                  n.thin = con$n.thin, seed = NULL, verbose = con$verbose)

  #############################

  out <- list(mcmc = fit$sims.list, components = components, data = data,
                  families = families, control = con, mcmc.info = fit$mcmc.info[-which(names(fit$mcmc.info) %in% c("end.values"))],
                  DIC = fit$DIC, pD = fit$pD, Rhat = fit$Rhat)

  Xnams <- lapply(components[grep("^X[0-9]", names(components))], colnames)
  for (i in seq_along(Xnams)) {
        colnames(out$mcmc[[paste0("betas", i)]]) <- Xnams[[i]]
  }
      #pat <- paste0("^Z", paste(rep("[0-9]", nchar(as.character(n_outcomes))), collapse = ""))
  Znams <- lapply(components[grep("^Z[0-9]", names(components))], colnames)
  Znams <- unlist(mapply(paste0, Znams, seq_len(n_outcomes), SIMPLIFY = FALSE))
      #dimnames(out$mcmc$D) <- dimnames(out$mcmc$inv.D) <- list(NULL, Znams, Znams)
      #dimnames(out$mcmc$b) <- list(NULL, NULL, Znams)
      # calculate statistics
  summary_fun <- function (FUN, ...) {
      out <- lapply(out$mcmc, function (x) {
          if (!is.null(dim(x)) && length(dim(x)) > 1) {
            d <- if (is.matrix(x)) 2L else c(2L, 3L)
            apply(x, d, FUN, ...)
          } else {
            FUN(x, ...)
          }
      })
        out[!sapply(out, is.null)]
  }
  out$postMeans <- summary_fun(mean, na.rm = TRUE)
  out$postModes <- summary_fun(modes)
  out$EffectiveSize <- summary_fun(effectiveSize)
  out$StErr <- summary_fun(stdErr)
  out$StDev <- summary_fun(sd, na.rm = TRUE)
  out$CIs <- summary_fun(quantile, probs = c(0.025, 0.975))
  out$Pvalues <- summary_fun(computeP)
  out$call <- formulas
  class(out) <- "mvlme"
  out

  #fit


}

