DynPred_mv_lme <- function(object, newdata, families, hc,
                             level = 0.95, IdVar, timeVar,
                             M = 200, times = NULL,
                             assoc = TRUE,
                             assoc_from,
                             assoc_to,
                             extraForm = NULL,
                             extraForm_tr = NULL){

    #newdata$id <- newdata[[IdVar]]
    #K <- length(families)

    #try(setwd("Predictions"))

    ########
    # Data #
    ########

    n_outcomes <- length(object$call)
    seq_outcomes <- seq_len(n_outcomes)

    data <- object$data
    formYx <- formulas <- object$call
    term_labels <- lapply(formYx, function(x) attr(terms(x), "term.labels"))
    len_ran_ef <- lapply(term_labels, function(x) which(x == x[grep("|", x, fixed = TRUE)]))
    formYx <- mapply(function(x, y) drop.terms(terms(x), dropx = y, keep.response = TRUE),
                     x = formYx,
                     y = len_ran_ef,
                     SIMPLIFY = FALSE)
    mfX <- lapply(formYx, function(x) model.frame(terms(x), data = object$data)) # model.frame(terms(formYx), data = data)
    TermsX <- lapply(mfX, function(x) attr(x, "terms"))   #attr(mfX, "terms")
    components <- object$components
    formtest <- components[grep("TermsZ", names(components), fixed = TRUE)]
    formYz <- lapply(formtest, function(x)  formula(x)) #formula(lmeObject$modelStruct$reStruct[[1]])
    mfZ <- lapply(formYz, function(x) model.frame(x, data = data)) #model.frame(terms(formYz), data = data)
    TermsZ <- lapply(mfZ, function(x) attr(x, "terms"))   #attr(mfZ, "terms")


    ###########
    # Newdata #
    ###########

    mfX_new_pred <- mapply(function(x) model.frame(x, data = newdata, na.action = NULL),
                           x = TermsX,
                           SIMPLIFY = FALSE) #model.frame(TermsX, data = newdata_pred, na.action = NULL)
    X_pred <- mapply(function(x, y) model.matrix(x, y),
                         x = TermsX,
                         y = mfX_new_pred,
                         SIMPLIFY = FALSE) #model.matrix(formYx, mfX_new_pred)
    mfZ_new_pred <- mapply(function(x) model.frame(x, data = newdata, na.action = NULL),
                           x = TermsZ,
                           SIMPLIFY = FALSE) #model.frame(TermsZ, data = newdata_pred, na.action = NULL)
    Z_pred <- mapply(function(x, y) model.matrix(as.formula(x), y),
                         x = TermsZ,
                         y = mfZ_new_pred,
                         SIMPLIFY = FALSE) #model.matrix(formYz, mfZ_new_pred)



    times_orig  <- object$data[[timeVar]]
    if (is.null(times)) {
      times <-  seq(min(times_orig, na.rm = T), max(times_orig, na.rm = T), length.out = 100)
    }
    last_time <- max(newdata[[timeVar]])
    t <- max(newdata[[timeVar]])
    times.to.pred <- times[times > t]


    newdata_pred <- tail(newdata, n = 1L)
    newdata_pred <- newdata_pred[rep(1:nrow(newdata_pred),
                             length(times.to.pred)), ]
    newdata_pred[[timeVar]] <- unlist(times.to.pred)


    mfX_new_pred2 <- mapply(function(y, z)  model.frame(y, data = newdata_pred, na.action = NULL,
                                                     xlev = .getXlevels(y, z)),
                            y = TermsX,
                            z = mfX,
                            SIMPLIFY = FALSE) #model.frame(TermsX, data = newdata_pred, na.action = NULL)
    X_pred2 <- mapply(function(x, y) model.matrix(x, y),
                     x = TermsX,
                     y = mfX_new_pred2,
                     SIMPLIFY = FALSE) #model.matrix(formYx, mfX_new_pred)
    mfZ_new_pred2 <- mapply(function(y, z) model.frame(y, data = newdata_pred, na.action = NULL,
                                                       xlev = .getXlevels(y, z)),
                            y = TermsX,
                            z = mfX,
                           SIMPLIFY = FALSE) #model.frame(TermsZ, data = newdata_pred, na.action = NULL)
    Z_pred2 <- mapply(function(x, y) model.matrix(as.formula(x), y),
                     x = TermsZ,
                     y = mfZ_new_pred2,
                     SIMPLIFY = FALSE) #model.matrix(formYz, mfZ_new_pred)


    components_newdata <- lapply(unname(formulas), extractFrames_Pred, data = newdata)

    components_newdata <- unlist(components_newdata, recursive = FALSE)
    n_outcomes <- length(formulas)
    names(components_newdata) <- paste0(names(components_newdata),
                                rep(seq_len(n_outcomes),
                                    each = length(components_newdata) / n_outcomes))
    colmns_HC <- components_newdata[grep("colmns_HC", names(components_newdata), fixed = TRUE)]
    colmns_nHC <- components_newdata[grep("colmns_nHC", names(components_newdata), fixed = TRUE)]

    nams_vars <- c("ncx", "ncz")
    vars <- paste0(rep(nams_vars, each = n_outcomes), seq_outcomes)
    # if (any(ind_td <- sapply(colmns_nHC, length))) {
    #   vars <- c(vars, paste0("X", which(ind_td > 0)))
    # }

    Data_data <- c(list(n = components_newdata$n1), components_newdata[vars])

    Data_data$n_RE <- sum(unlist(components_newdata[grep("ncz", names(components_newdata), fixed = TRUE)]))
    RE_inds <- mapply(function (sq, incr) seq_len(sq) + incr,
                      sq = components_newdata[grep("ncz", names(components_newdata), fixed = TRUE)],
                      incr = cumsum(c(0, head(sapply(colmns_HC, length), -1))),
                      SIMPLIFY = FALSE)
    names(RE_inds) <- paste0("RE_ind", seq_along(RE_inds))
    Data1new <- c(Data_data, RE_inds)

    names(X_pred) <- paste0("X", seq_outcomes)
    names(Z_pred) <- paste0("Z", seq_outcomes)

    names(X_pred2) <- paste0("newMatF", seq_outcomes)
    names(Z_pred2) <- paste0("newMatR", seq_outcomes)


    Data1new <- c(Data1new, X_pred, Z_pred, X_pred2, Z_pred2)

    names(Data1new) <- paste0(names(Data1new), "_pred")

    Data <- c(Data1new)#, Data2)

    ########################################
    # Design matrices for parameterization #
    ########################################

    if (!is.null(extraForm)){
      mfX_derivY <- model.frame(terms(extraForm$fixed), data = newdata)
      mfZ_derivY <- model.frame(terms(extraForm$random), data = newdata)
      XderivY <- model.matrix(extraForm$fixed, mfX_derivY)
      ZderivY <- model.matrix(extraForm$random, mfZ_derivY)

      Data$deriv_fixed <- XderivY
      Data$deriv_random <- ZderivY

      mfX_derivY_new <- model.frame(terms(extraForm$fixed), data = newdata_pred)
      mfZ_derivY_new <- model.frame(terms(extraForm$random), data = newdata_pred)
      XderivY_new <- model.matrix(extraForm$fixed, mfX_derivY_new)
      ZderivY_new <- model.matrix(extraForm$random, mfZ_derivY_new)

      Data$new_deriv_fixed <- XderivY_new
      Data$new_deriv_random <- ZderivY_new

    }

    if (!is.null(extraForm_tr)){

      mfX_derivY_tr <- model.frame(terms(extraForm_tr$fixed), data = newdata)
      mfZ_derivY_tr <- model.frame(terms(extraForm_tr$random), data = newdata)
      XderivY_tr <- model.matrix(extraForm_tr$fixed, mfX_derivY_tr)
      ZderivY_tr <- model.matrix(extraForm_tr$random, mfZ_derivY_tr)

      Data$deriv_fixed <- XderivY - XderivY_tr
      Data$deriv_random <- ZderivY - ZderivY_tr

      mfX_derivY_new_tr <- model.frame(terms(extraForm_tr$fixed), data = newdata_pred)
      mfZ_derivY_new_tr <- model.frame(terms(extraForm_tr$random), data = newdata_pred)
      XderivY_new_tr <- model.matrix(extraForm_tr$fixed, mfX_derivY_new_tr)
      ZderivY_new_tr <- model.matrix(extraForm_tr$random, mfZ_derivY_new_tr)

      Data$new_deriv_fixed <- XderivY_new - XderivY_new_tr
      Data$new_deriv_random <- ZderivY_new - ZderivY_new_tr

    }

    #####################
    # Obtain parameters #
    #####################

    mcmc <- object$mcmc
    mcmc <- mcmc[!(names(mcmc) %in% "b")]

    mcmc <- mcmc[names(mcmc) != "deviance"]


    samples <- sample(nrow(mcmc[[1]]), M)
    mcmc[] <- lapply(mcmc, function (x) {
      if (length(dim(x)) == 3) {
        x[samples, , ]
      } else if (length(dim(x)) == 2) {
        x[samples, , drop = FALSE]
      } else if (all(is.null(dim(x)), length(x) > 0)) {
        x[samples]
      } else {
        x[]
      }
    })



    ##########################################
    # Updating random effects from posterior #
    ##########################################

    outcomes <- lapply(formulas, function(x) all.vars(x)[1])

    log_denss <- lapply(families, function(x) switch(x,
                                                    'binomial' = binomial_log_dens,
                                                    'gaussian' = normal_log_dens))




    b_i <- rep(0, Data$n_RE_pred)
    post_modes <- matrix(,,length(b_i))
    post_hessians <- vector("list")





    # ff <- function(b_i, y_i, X_i, Z_i, Sigma,
    #                betas, invD, log_dens, prob_class) {
    #
    #   -(
    #     sum(eval(parse(text = (paste0("(-0.5 * ((y_i[[", seq_len(n_outcomes),
    #             "]] - X_i[[", seq_len(n_outcomes), "]]%*%betas$betas", seq_len(n_outcomes), which(c(v) == 1),
    #             "-Z_i[[", seq_len(n_outcomes), "]]%*%b_i[Data$RE_ind", seq_len(n_outcomes), "]) / Sigma[[",
    #             seq_len(n_outcomes), "]])^2)",
    #                                   collapse = " + ")))
    #     ))  -
    #       (0.5*(t(b_i)%*%(invD)%*%(b_i))) )
    # }
    for (o in seq_len(n_outcomes)){
         vec <- paste0("object$postMeans$sigma", o)
         if (!is.numeric(eval(parse(text = vec)))) {
           eval(parse(text = paste0(vec, " <- 0")))
         }
    }



    #gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, i = i)
    opt_i <- optim(par = b_i, fn = log_post_b_EB, method = "BFGS",
                   y_i = lapply(outcomes, function(x) newdata[[x]]),
                   X_i = Data[grep("X", names(Data), fixed = TRUE)],#[c(1,3,4,5)],
                   #X_lis = Data[grep("newMatF", names(Data), fixed = TRUE)],
                   Z_i = Data[grep("Z", names(Data), fixed = TRUE)],
                   #Z_lis = Data[grep("newMatR", names(Data), fixed = TRUE)],
                   betas = object$postMeans[grep("betas", names(object$postMeans), fixed = TRUE)],
                   invD = Matrix::bdiag(object$postMeans[grep(paste0("inv.D"), names(mcmc), fixed = TRUE)]),
                   #object$postMeans$inv.D,  #object$postMeans$inv_D,#object$postMeans$inv.D1,
                   Sigma = object$postMeans[grep("sigma", names(object$postMeans), fixed = TRUE)],
                   alphas = object$postMeans[grep("alpha", names(object$postMeans), fixed = TRUE)],
                   Data = Data,
                   log_dens = log_denss,
                   hessian = TRUE,
                   assoc = assoc,
                   assoc_from = assoc_from,
                   assoc_to = assoc_to,
                   extraForm = extraForm,
                   extraForm_tr = extraForm_tr,

                   XderivY_i = Data[grep("deriv_fixed", names(Data), fixed = TRUE)],
                   ZderivY_i = Data[grep("deriv_random", names(Data), fixed = TRUE)]
    )

    EBs <- list()
    post_modes[1,] <- opt_i$par
    post_hessians[[1]] <- opt_i$hessian
    EBs <- list(post_modes = post_modes, post_hessians = post_hessians)

    # EB mean prediction
    # for (k in 1:n_outcomes){
    #  if(v[[1]] == 1) {
    #     marg_part <- paste0("Data$newMatF", k, "_pred %*% object$postMeans$betas", k,"1")
    #   } else if (v[[2]]==1) {
    #     marg_part <- paste0("Data$newMatF", k, "_pred %*% object$postMeans$betas", k, "2")
    #   }
    #   pred_mean <- paste0("pred_y", k, " <- c(", marg_part, ") + rowSums(Data$newMatR", k, "_pred * EBs$post_modes[rep(1, dim(Data$newMatR", k, "_pred)[1]), Data$RE_ind", k, "])")
    #   eval(parse(text = pred_mean))
    # }


    for (k in 1:n_outcomes){
      marg_part <- paste0("Data$newMatF", k, "_pred %*% object$postMeans$betas", k)
      pred_mean <- paste0("pred_y", k, " <- c(", marg_part, ") + rowSums(Data$newMatR", k, "_pred * EBs$post_modes[rep(1, dim(Data$newMatR", k, "_pred)[1]), Data$RE_ind", k, "])")
      eval(parse(text = pred_mean))
    }

    if (assoc == TRUE){
     for (k in 1:n_outcomes){

           if (k %in% assoc_from){
             if (!is.null(extraForm)) {
               pred_mean <- paste0("new_pred_y", k, " <- c(Data$new_deriv_fixed %*% object$postMeans$betas", k, "[",
                                   if (length(  extraForm$indFixed) > 1) {
                                       paste0("c(", paste0(extraForm$indFixed, collapse = ","), ")" )
                                    } else { extraForm$indFixed }
               , "]",
                                   " + rowSums(Data$new_deriv_random", " * EBs$post_modes[rep(1, dim(Data$new_deriv_random", ")[1]), ",
                                   "Data$RE_ind", k, "[extraForm$indRandom]]))")
             } else {
             marg_part <- paste0("Data$newMatF", k, "_pred %*% object$postMeans$betas", k)
             pred_mean <- paste0("new_pred_y", k, " <- c(", marg_part, ") + rowSums(Data$newMatR", k,
                                 "_pred * EBs$post_modes[rep(1, dim(Data$newMatR", k, "_pred)[1]), Data$RE_ind", k, "])")
             }
           eval(parse(text = pred_mean))
           }

      }


     for (u in assoc_to){
       for (o in assoc_from){

          pred_mean2 <- paste0("pred_y", u, " <- pred_y", u, " + new_pred_y", o, " * unlist(object$postMeans$alpha", u, o, ")")


        }
      }

      eval(parse(text = pred_mean2))


    }


    row_split_ind <- row(EBs$post_modes)
    mu <- split(EBs$post_modes, row_split_ind)


    SigmaEBs <- lapply(EBs$post_hessians, solve)

    scale = 1.6
    scale <- rep(scale, length.out = length(SigmaEBs))
    Sigmaa <- mapply("*", scale, SigmaEBs, SIMPLIFY = FALSE)

    df = 4

    EBs_proposed <- mapply(rmvt, mu = mu, Sigma = Sigmaa, SIMPLIFY = FALSE,
                           MoreArgs = list(n = M, df = df))

    ######################
    # Comments/questions #
    ######################

    # I am not using it for now, I need to check it...
    dmvt_proposed <- mapply(dmvt, x = EBs_proposed, mu = mu, Sigma = Sigmaa,
                            MoreArgs = list(df = df, log = TRUE, prop = FALSE),
                            SIMPLIFY = FALSE)

    b_current <- mu

    ## Shouldn't mu = mu be mu = EBs_proposed?
    dmvt_current <- mapply(dmvt, x = mu, mu = mu, Sigma = Sigmaa, SIMPLIFY = FALSE,
                           MoreArgs = list(df = df, log = TRUE, prop = FALSE))
    # dmvt_current <- mapply(dmvt, x = mu, mu = EBs_proposed, Sigma = Sigmaa, SIMPLIFY = FALSE,
    #                        MoreArgs = list(df = df, log = TRUE, prop = FALSE))
    # #Error in (function (x, mu, Sigma = NULL, invSigma = NULL, df, log = TRUE,  :
    # #                      dims [product 20] do not match the length of object [8000]

    b <- vector("list", M)
    success_rate <- matrix(FALSE, M, 1)#length(y_lis))

    for (k in 1:n_outcomes){
      pred_list <- paste0("eta", k, " <- list()")
      eval(parse(text = pred_list))
    }

    for (k in 1:n_outcomes){
      pred_list <- paste0("der_eta", k, " <- list()")
      eval(parse(text = pred_list))
    }



    # if(v[[1]] == 1) {
    #   log_post_b <- function(b_i, y_i, X_i, Z_i, Sigma,
    #                  betas, invD, log_dens, prob_class) {
    #
    #     sum(eval(parse(text = (paste0("-0.5 * ((y_i[[", seq_len(n_outcomes), "]] - X_i[[", seq_len(n_outcomes), "]]%*%betas$betas", seq_len(n_outcomes), "1-Z_i[[", seq_len(n_outcomes), "]]%*%b_i[Data$RE_ind", seq_len(n_outcomes), "]) / Sigma[[", seq_len(n_outcomes), "]])^2", collapse = " + ")))
    #     ))  -
    #       (0.5*(t(b_i)%*%(invD)%*%(b_i)))
    #   }
    # } else if (v[[2]]==1) {
    #   log_post_b <- function(b_i, y_i, X_i, Z_i, Sigma,
    #                  betas, invD, log_dens, prob_class) {
    #
    #     sum(eval(parse(text = (paste0("-0.5 * ((y_i[[", seq_len(n_outcomes), "]] - X_i[[", seq_len(n_outcomes), "]]%*%betas$betas", seq_len(n_outcomes), "2-Z_i[[", seq_len(n_outcomes), "]]%*%b_i[Data$RE_ind", seq_len(n_outcomes), "]) / Sigma[[", seq_len(n_outcomes), "]])^2", collapse = " + ")))
    #     ))  -
    #       (0.5*(t(b_i)%*%(invD)%*%(b_i)))
    #   }
    # }


    # log_post_b <- function(b_i, y_i, X_i, Z_i, Sigma,
    #                        betas, invD, log_dens, prob_class) {
    #
    #   sum(eval(parse(text = (paste0("-0.5 * ((y_i[[", seq_len(n_outcomes), "]] - X_i[[",
    #                                 seq_len(n_outcomes), "]]%*%betas$betas",
    #                                 seq_len(n_outcomes), which(c(v) == 1), "-Z_i[[",
    #                                 seq_len(n_outcomes), "]]%*%b_i[Data$RE_ind",
    #                                 seq_len(n_outcomes), "]) / Sigma[[", seq_len(n_outcomes), "]])^2",
    #                                 collapse = " + ")))
    #   ))  -
    #     (0.5*(t(b_i)%*%(invD)%*%(b_i)))
    # }

    for (m in seq_len(M)){
      #mcmc_sel <- lapply(mcmc, function(x) if (is.matrix(x) == TRUE) x[m, ] else x[m])

      log_post_b_current <- log_post_b_MCMC(b_i = unlist(b_current),
                                            y_i = lapply(outcomes, function(x) newdata[[x]]),
                                            X_i = Data[grep("X", names(Data), fixed = TRUE)],
                                            Z_i = Data[grep("Z", names(Data), fixed = TRUE)],
                                            betas = lapply(mcmc[grep("betas", names(mcmc), fixed = TRUE)], function(x) x[m,]),
                                            invD = matrix(as.numeric(as.data.frame(mcmc[grep(paste0("inv.D"), names(mcmc), fixed = TRUE)])[m, ,]),
                                                                Data$n_RE_pred, Data$n_RE_pred),
                                            Sigma = lapply(mcmc[grep("sigma", names(mcmc), fixed = TRUE)], function(x) x[m]),
                                            alphas = lapply(mcmc[grep("alpha", names(mcmc), fixed = TRUE)], function(x) x[m]),
                                            log_dens = log_denss, Data = Data, assoc = assoc,
                                            extraForm = extraForm,
                                            XderivY_i = Data[grep("deriv_fixed", names(Data), fixed = TRUE)],
                                            ZderivY_i = Data[grep("deriv_random", names(Data), fixed = TRUE)]
                                            )


      b_new <- lapply(EBs_proposed, function (x, m) x[m, ], m = m)
      log_post_b_new <- log_post_b_MCMC(b_i = unlist(b_new), y_i = lapply(outcomes, function(x) newdata[[x]]),
                                       X_i = Data[grep("X", names(Data), fixed = TRUE)],
                                       Z_i = Data[grep("Z", names(Data), fixed = TRUE)],
                                       betas = lapply(mcmc[grep("betas", names(mcmc), fixed = TRUE)], function(x) x[m,]),
                                       invD = matrix(as.numeric(as.data.frame(mcmc[grep(paste0("inv.D"), names(mcmc), fixed = TRUE)])[m, ,]),
                                                     Data$n_RE_pred, Data$n_RE_pred),
                                       Sigma = lapply(mcmc[grep("sigma", names(mcmc), fixed = TRUE)], function(x) x[m]),
                                       alphas = lapply(mcmc[grep("alpha", names(mcmc), fixed = TRUE)], function(x) x[m]),
                                       log_dens = log_denss, Data = Data,
                                       assoc = assoc,
                                       assoc_from = assoc_from,
                                       assoc_to = assoc_to,
                                       extraForm = extraForm,
                                       extraForm_tr = extraForm_tr,
                                       XderivY_i = Data[grep("deriv_fixed", names(Data), fixed = TRUE)],
                                       ZderivY_i = Data[grep("deriv_random", names(Data), fixed = TRUE)]
                                   )
      ############
      # New part #
      ############
      # dmvt_proposed <- mapply(dmvt, x = b_new, mu = b_current, Sigma = Sigmaa,
      #                         MoreArgs = list(df = df, log = TRUE, prop = FALSE),
      #                         SIMPLIFY = FALSE)
      #
      #
      # ## Shouldn't mu = mu be mu = EBs_proposed?
      # dmvt_current <- mapply(dmvt, x = b_current, mu = b_new, Sigma = Sigmaa, SIMPLIFY = FALSE,
      #                        MoreArgs = list(df = df, log = TRUE, prop = FALSE))

      # calc_alpha <- function (log_post_new, log_post_old, log_prop_new,
      #                         log_prop_old) {
      #   min(exp(log_post_new + log_prop_old - log_post_old - log_prop_new), 1)
      # }
      #alphas <- mapply(calc_alpha, log_post_b_new, log_post_b_current, dmvt_current,
      #                 lapply(dmvt_proposed, "[", m))

      #alphas <- min(exp(log_post_b_new + unlist(lapply(dmvt_proposed, "[", m)) - log_post_b_current - unlist(dmvt_current)), 1) WRONG??
      accept_ratio <- min(exp(log_post_b_new + unlist(dmvt_current) - log_post_b_current - unlist(lapply(dmvt_proposed, "[", m))), 1)

      accept_ratio <- min(exp(log_post_b_new - log_post_b_current), 1)

      #min(exp(log_post_b_new + dmvt_current - log_post_b_current - dmvt_proposed), 1)
      # calc_alpha <- function (log_post_new, log_post_old) {
      #   min(exp(log_post_new - log_post_old), 1)
      # }
      # alphas <- mapply(calc_alpha, log_post_b_new, log_post_b_current)

      keep_ind <- runif(length(accept_ratio)) <= accept_ratio
      if (any(keep_ind)) {
        b_current[keep_ind] <- b_new[keep_ind]
        dmvt_current[keep_ind] <- lapply(dmvt_proposed, "[", m)[keep_ind]
      }
      success_rate[m, ] <- keep_ind
      b[[m]] <- do.call("rbind", b_current)


      for (k in 1:n_outcomes){
        # marg_part <- paste0("prob_class", seq_len(classes), " * (Data$newMatF", k, "_pred %*% mcmc$betas", k, seq_len(classes), "[m,])", collapse = " + ")
          if (families[[k]] == "gaussian"){
            marg_part <- paste0("(Data$newMatF", k, "_pred %*% mcmc$betas", k,"[m,])", collapse = " + ")
            pred <- paste0("eta", k, "[[m]] <- ", marg_part, " + rowSums(Data$newMatR",k ,
                           "_pred * b[[m]][rep(1, dim(Data$newMatR",k , "_pred)[1]), Data$RE_ind", k, ", drop = FALSE])")
            eval(parse(text = pred))

            if (k %in% assoc_from){
              if (!is.null(extraForm)){
                marg_part <- paste0("(Data$new_deriv_fixed", " %*% mcmc$betas", k,"[m,",
                                          if (length(  extraForm$indFixed) > 1) {
                                            paste0("c(", paste0(extraForm$indFixed, collapse = ","), ")" )
                                          } else { extraForm$indFixed }
                                  , "])", collapse = " + ")
                pred <- paste0("der_eta", k, "[[m]] <- ", marg_part, " + rowSums(Data$new_deriv_random", " * b[[m]][rep(1, dim(Data$new_deriv_random", ")[1]), Data$RE_ind", k, "[extraForm$indRandom]",", drop = FALSE])")
              } else {
                marg_part <- paste0("(Data$newMatF", k, "_pred %*% mcmc$betas", k,"[m,])", collapse = " + ")
                pred <- paste0("der_eta", k, "[[m]] <- ", marg_part, " + rowSums(Data$newMatR",k ,
                               "_pred * b[[m]][rep(1, dim(Data$newMatR",k , "_pred)[1]), Data$RE_ind", k, ", drop = FALSE])")

              }
              eval(parse(text = pred))
            }

          }

          if (families[[k]] == "binomial"){
            marg_part <- paste0("(Data$newMatF", k, "_pred %*% mcmc$betas", k, "[m,])", collapse = " + ")
            pred2 <- paste0("eta_binom <- ", marg_part, " + rowSums(Data$newMatR",k , "_pred * b[[m]][rep(1, dim(Data$newMatR",k , "_pred)[1]), Data$RE_ind", k, ", drop = FALSE])")
            eval(parse(text = pred))

            if (!is.null(extraForm)){
              marg_part <- paste0("(Data$new_deriv_fixed", " %*% mcmc$betas", k,"[m,",
                                        if (length(  extraForm$indFixed) > 1) {
                                          paste0("c(", paste0(extraForm$indFixed, collapse = ","), ")" )
                                        } else { extraForm$indFixed }
                                  , "])", collapse = " + ")
              pred2 <- paste0("der_eta", k, "[[m]] <- ", marg_part, " + rowSums(Data$new_deriv_random", " * b[[m]][rep(1, dim(Data$new_deriv_random", ")[1]), Data$RE_ind", k, "[extraForm$indRandom]",", drop = FALSE])")
            } else {
              marg_part <- paste0("(Data$newMatF", k, "_pred %*% mcmc$betas", k,"[m,])", collapse = " + ")
              pred2 <- paste0("der_eta", k, "[[m]] <- ", marg_part, " + rowSums(Data$newMatR",k , "_pred * b[[m]][rep(1, dim(Data$newMatR",k , "_pred)[1]), Data$RE_ind", k, ", drop = FALSE])")

            }

            eval(parse(text = pred2))
            pred <- paste0("eta", k, "[[m]] <- exp(eta_binom)/(1 + exp(eta_binom))")
            }
      }

      if (assoc == TRUE){
        for (u in assoc_to){
          for (o in assoc_from){
            pred22 <- paste0("eta", u, "[[m]] <- eta", u, "[[m]] + der_eta", o, "[[m]] * mcmc$alpha", u, o, "[[m]]")
            eval(parse(text = pred22))
          }
        }
      }

    }

    #   for (k in 1:n_outcomes){
    #     marg_part <- paste0("(Data$newMatF", k, "_pred %*% mcmc$betas", k, which(c(v) == 1), "[m,])", collapse = " + ")
    #     pred <- paste0("eta", k, " [[m]] <- ", marg_part, " +
    #                 rowSums(Data$newMatR",k , "_pred * b[[m]][rep(1, dim(Data$newMatR",k , "_pred)[1]), Data$RE_ind", k, ", drop = FALSE])")
    #     eval(parse(text = pred))
    #  }
    #
    # }


    ###########
    # Results #
    ###########
    pred_eta <- paste0("pred_eta", seq_len(n_outcomes), " <- do.call(cbind, eta", seq_len(n_outcomes), ")")
    eval(parse(text = pred_eta))

    eval(parse(text = paste0("Mean", seq_len(n_outcomes), " = rowMeans(pred_eta", seq_len(n_outcomes), ", na.rm = TRUE)\n", collapse = "")))
    eval(parse(text = paste0("SD", seq_len(n_outcomes), " = apply(pred_eta", seq_len(n_outcomes), ", 1, FUN = sd, na.rm = TRUE)\n", collapse = "")))
    eval(parse(text = paste0("Median", seq_len(n_outcomes), " = apply(pred_eta", seq_len(n_outcomes), ", 1, FUN = median, na.rm = TRUE)\n", collapse = "")))
    eval(parse(text = paste0("Lower", seq_len(n_outcomes), " = apply(pred_eta", seq_len(n_outcomes), ", 1, quantile, probs = c((1 - level) / 2, (1 + level) / 2)[1], na.rm = TRUE)\n", collapse = "")))
    eval(parse(text = paste0("Upper", seq_len(n_outcomes), " = apply(pred_eta", seq_len(n_outcomes), ", 1, quantile, probs = c((1 - level) / 2, (1 + level) / 2)[2], na.rm = TRUE)\n", collapse = "")))

    ress <- paste0("Mean", seq_len(n_outcomes), ", SD", seq_len(n_outcomes), ", Median", seq_len(n_outcomes), ",
                             Lower", seq_len(n_outcomes), ", Upper", seq_len(n_outcomes), collapse = ",")
    ress2 <- paste0("data.frame(", ress, ")")

    ress3 <- eval(parse(text = ress2))

    out <- data.frame(ress3,  "Time" = times.to.pred)
    out

}


