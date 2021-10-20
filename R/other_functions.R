extractFrames <- function (formula, data) {
  Terms <- terms(formula)
  term_labels <- attr(Terms, "term.labels")
  which_RE <- grep("|", term_labels, fixed = TRUE)
  namesVars <- all.vars(formula)
  respVar <- as.character(formula)[2L]
  # Fixed Effects
  formYx <- paste(term_labels[-which_RE], collapse = " + ")
  formYx <- as.formula(paste(respVar, "~", formYx))
  TermsX <- terms(formYx, data = data)
  mfX <- model.frame(TermsX, data)
  TermsX <- terms(mfX)
  X <- model.matrix(TermsX, data)
  # Random Effects
  spl <- unlist(strsplit(term_labels[which_RE], " | ", fixed = TRUE))
  idVar <- spl[2L]
  data <- data[complete.cases(data[namesVars]), ]
  id <- data[[idVar]]
  id <- match(id, unique(id))
  formYz <- paste(spl[1], collapse = " + ")
  formYz <- as.formula(paste(respVar, "~", formYz))
  TermsZ <- terms(formYz, data = data)
  mfZ <- model.frame(TermsZ, data = data)
  TermsZ <- terms(mfZ)
  Z <- model.matrix(TermsZ, data)
  # response variable
  y <- model.response(mfX)
  if (is.factor(y))
    y <- as.vector(unclass(y) - 1)
  # repeated measurements offset
  offset <- as.vector(c(1, 1 + cumsum(tapply(y, id, length))))
  # hierarchical centering
  find_positions <- function (nams1, nams2) {
    nams1 <- gsub("^", "\\^", nams1, fixed = TRUE)
    vals <- c(glob2rx(nams1), glob2rx(paste0(nams1, ":*")),
              glob2rx(paste0("*:", nams1)))
    out <- sort(unique(unlist(lapply(vals, grep, x = nams2))))
    out
  }
  check_td <- function (x, id) {
    !all(sapply(split(x, id), function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
  }
  has_interceptX <- attr(TermsX, "intercept")
  has_interceptZ <- attr(TermsZ, "intercept")
  performHC <- has_interceptX && (has_interceptX == has_interceptZ)
  if (performHC) {
    terms.labs_X <- attr(TermsX, "term.labels")
    terms.labs_Z <- attr(TermsZ, "term.labels")
    # check for time-varying covariates
    timeTerms <- if (length(terms.labs_Z))
      unlist(lapply(terms.labs_Z, FUN = function(x) grep(x, colnames(X), fixed = TRUE)))
    which_td <- unname(which(apply(X, 2, check_td, id = id)))
    all_TDterms <- unique(c(timeTerms, which_td))
    baseline <- seq_len(ncol(X))[-all_TDterms]
    ind_colmns <- c(list(baseline), lapply(colnames(Z)[-1L], find_positions, 
                                           nams2 = colnames(X)))
    ind_colmns2 <- seq_len(ncol(X))
    ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
    data.id <- data[!duplicated(id), ]
    Xhc <- if (length(terms.labs_Z)) {
      mfHC <- model.frame(TermsX, data = data.id)
      which.timevar <- unique(unlist(lapply(terms.labs_Z, 
                                            FUN = function (x) grep(x, names(mfHC), fixed = TRUE))))
      mfHC[which.timevar] <- lapply(mfHC[which.timevar], 
                                    function (x) { x[] <- 1; x })
      model.matrix(formYx, mfHC)
    } else {
      model.matrix(formYx, model.frame(TermsX, data = data.id))
    }
  }
  environment(TermsX) <- environment(TermsZ) <- NULL
  
  Xc = scale(X[, -1], center = TRUE, scale = FALSE) # except intercept
  Xc = cbind(X[,1], Xc)
  Xs = scale(X[, -1], center = TRUE, scale = TRUE)
  Xs = cbind(X[,1], Xs)
  
  
  if (ncol(X) >2) {
    means_X = apply(X[, -1], 2, mean)
    SDs_X = apply(X[, -1], 2, sd)
  } else{
    means_X = mean(X[, -1])
    SDs_X = sd(X[, -1])
  }
  mean_sd_X = means_X/SDs_X
  
  # extract results
  list(N = nrow(Z), n = length(unique(id)), offset = offset, idVar = idVar, respVar = respVar,
       id = id, y = y, X = X, Z = Z, TermsX = TermsX,
       TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
       Xc = Xc, Xhc = Xhc, colmns_HC = ind_colmns, colmns_nHC = ind_colmns2,
       means_X = means_X, ncx = ncol(X), ncz = ncol(Z))
}


extractFrames_Pred <- function (formula, data) {
  Terms <- terms(formula)
  term_labels <- attr(Terms, "term.labels")
  which_RE <- grep("|", term_labels, fixed = TRUE)
  namesVars <- all.vars(formula)
  respVar <- as.character(formula)[2L]
  # Fixed Effects
  formYx <- paste(term_labels[-which_RE], collapse = " + ")
  formYx <- as.formula(paste(respVar, "~", formYx))
  TermsX <- terms(formYx, data = data)
  mfX <- model.frame(TermsX, data)
  TermsX <- terms(mfX)
  X <- model.matrix(TermsX, data)
  # Random Effects
  spl <- unlist(strsplit(term_labels[which_RE], " | ", fixed = TRUE))
  idVar <- spl[2L]
  data <- data[complete.cases(data[namesVars]), ]
  id <- data[[idVar]]
  id <- match(id, unique(id))
  formYz <- paste(spl[1], collapse = " + ")
  formYz <- as.formula(paste(respVar, "~", formYz))
  TermsZ <- terms(formYz, data = data)
  mfZ <- model.frame(TermsZ, data = data)
  TermsZ <- terms(mfZ)
  Z <- model.matrix(TermsZ, data)
  # response variable
  y <- model.response(mfX)
  if (is.factor(y))
    y <- as.vector(unclass(y) - 1)
  # repeated measurements offset
  offset <- as.vector(c(1, 1 + cumsum(tapply(y, id, length))))
  # hierarchical centering
  find_positions <- function (nams1, nams2) {
    nams1 <- gsub("^", "\\^", nams1, fixed = TRUE)
    vals <- c(glob2rx(nams1), glob2rx(paste0(nams1, ":*")),
              glob2rx(paste0("*:", nams1)))
    out <- sort(unique(unlist(lapply(vals, grep, x = nams2))))
    out
  }
  check_td <- function (x, id) {
    !all(sapply(split(x, id), function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
  }
  has_interceptX <- attr(TermsX, "intercept")
  has_interceptZ <- attr(TermsZ, "intercept")
  performHC <- has_interceptX && (has_interceptX == has_interceptZ)
  if (performHC) {
    terms.labs_X <- attr(TermsX, "term.labels")
    terms.labs_Z <- attr(TermsZ, "term.labels")
    # check for time-varying covariates
    timeTerms <- if (length(terms.labs_Z))
      unlist(lapply(terms.labs_Z, FUN = function(x) grep(x, colnames(X), fixed = TRUE)))
    which_td <- unname(which(apply(X, 2, check_td, id = id)))
    all_TDterms <- unique(c(timeTerms, which_td))
    baseline <- seq_len(ncol(X))[-all_TDterms]
    ind_colmns <- c(list(baseline), lapply(colnames(Z)[-1L], find_positions, 
                                           nams2 = colnames(X)))
    ind_colmns2 <- seq_len(ncol(X))
    ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
    data.id <- data[!duplicated(id), ]
    Xhc <- if (length(terms.labs_Z)) {
      mfHC <- model.frame(TermsX, data = data.id)
      which.timevar <- unique(unlist(lapply(terms.labs_Z, 
                                            FUN = function (x) grep(x, names(mfHC), fixed = TRUE))))
      mfHC[which.timevar] <- lapply(mfHC[which.timevar], 
                                    function (x) { x[] <- 1; x })
      model.matrix(formYx, mfHC)
    } else {
      model.matrix(formYx, model.frame(TermsX, data = data.id))
    }
  }
  environment(TermsX) <- environment(TermsZ) <- NULL
  

  
  # extract results
  list(N = nrow(Z), n = length(unique(id)), offset = offset, idVar = idVar, respVar = respVar,
       id = id, y = y, X = X, Z = Z, TermsX = TermsX,
       TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
       Xhc = Xhc, colmns_HC = ind_colmns, colmns_nHC = ind_colmns2,
       ncx = ncol(X), ncz = ncol(Z))
}


modes <-
  function (y) {
    test <- try(d <- density(y, bw = "nrd", adjust = 3, n = 1000), silent = TRUE)
    if (!inherits(test, "try-error")) d$x[which.max(d$y)] else NA
  }

stdErr <-
  function (x) {
    x <- as.matrix(x)
    vars <- apply(x, 2L, var)
    ess <- effectiveSize(x)
    sqrt(vars / ess)
  }

effectiveSize <-
  function (x) {
    # copied (and made a bit more efficient) from the coda package
    spectrum0.ar <- function (x) {
      d <- dim(x)
      nrx <- d[1L]
      ncx <- d[2L]
      v0 <- numeric(ncx)
      res <- as.matrix(lm.fit(cbind(1, seq_len(nrx)), cbind(x, x))$residuals)
      for (i in seq_len(ncx)) {
        if (identical(all.equal(sd(res[, i]), 0), TRUE)) {
          v0[i] <- 0
        } else {
          ar.out <- ar(x[, i], aic = TRUE)
          v0[i] <- ar.out$var.pred / (1 - sum(ar.out$ar))^2
        }
      }
      v0
    }
    x <- as.matrix(x)
    spec <- spectrum0.ar(x)
    ifelse(spec == 0, 0, nrow(x) * apply(x, 2L, var) / spec)
  }

computeP <-
  function (x) {
    above <- mean(x >= 0)
    below <- mean(x < 0)
    2 * min(above, below)
  }


summary.mvlme <- function (object, classes, ...) {
  families <- object$families
  n_outcomes <- length(families)
  components <- object$components
  extract_components <- function (nam) {
    components[grep(nam, names(components), fixed = TRUE)]
  }
  respVars <- unlist(extract_components("respVar"), use.names = FALSE)
  descrpt <- data.frame(" " = unlist(extract_components("N"), use.names = FALSE),
                        row.names = respVars, check.rows = FALSE, check.names = FALSE)
  out <- list(n = components$n1, descrpt = descrpt, inv_D = object$postMeans[grep("inv.D", names(object$postMeans))],
              families = families, respVars = respVars,
              control = object$control, mcmc.info = object$mcmc.info[-which(names(object$mcmc.info) %in% c("end.values"))],
              DIC = object$DIC, pD = object$pD)
  for (i in 1:n_outcomes) {
      for (j in 1:classes) {
        out[[paste0("Outcome", i, j)]] <- data.frame("PostMean" = unlist(object$postMeans[grep(paste0("betas", i, j), names(object$postMeans))]),
                                              "StDev" = unlist(object$StDev[grep(paste0("betas", i, j), names(object$postMeans))]),
                                              "StErr"= unlist(object$StErr[grep(paste0("betas", i, j), names(object$postMeans))]),
                                              "2.5%" = unlist(lapply(object$CIs[grep(paste0("betas", i, j), names(object$postMeans))], function(x) x[1, ])),
                                              "97.5%" = unlist(lapply(object$CIs[grep(paste0("betas", i, j), names(object$postMeans))], function(x) x[2, ])),
                                              "P" = unlist(object$Pvalues[grep(paste0("betas", i, j), names(object$postMeans))]),
                                              "Rhat" = unlist(object$Rhat[grep(paste0("betas", i, j), names(object$postMeans))]),
                                              row.names = names(unlist(object$postMeans[grep(paste0("betas", i, j), names(object$postMeans))])),
                                              check.names = FALSE)
        if (families[[i]] == "gaussian") {
          D <- data.frame("PostMean" = object$postMeans[[paste0("sigma", i)]],
                          "StDev" = object$StDev[[paste0("sigma", i)]],
                          "StErr"= object$StErr[[paste0("sigma", i)]],
                          "2.5%" = object$CIs[[paste0("sigma", i)]][1],
                          "97.5%" = object$CIs[[paste0("sigma", i)]][2],
                          "P" = object$Pvalues[[paste0("sigma", i)]],
                          "Rhat" = object$Rhat[[paste0("sigma", i)]],
                          row.names = "tau", check.names = FALSE)
          out[[paste0("Outcome", i, j)]] <- rbind(out[[paste0("Outcome", i, j)]], D)
        }
      }
 

  }
  class(out) <- "summary.mvlclme"
  out
}


