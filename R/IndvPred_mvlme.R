
IndvPred_mvlme <- function(lmeObject, newdata, timeVar, times = NULL, M = 200L, idVar, 
                           interval = c("confidence", "prediction"),
                           all_times = FALSE,
                           level = 0.95, return_data = FALSE, seed = 1L) {
  if (!inherits(lmeObject, "mvglmer")) ##  && !inherits(lmeObject, "lmeComponents")
    stop("Use only with 'lme' or 'lmeComponents' objects.\n")
  interval <- interval #match.arg(interval)
  if (inherits(lmeObject, "mvglmer")) {
    
    #####################################################
    #### Obtain design matrices, posterior estimates ####
    #####################################################
    
    data <- lmeObject$data
    formYx <- formula(lmeObject)
    term_labels <- lapply(formYx, function(x) attr(terms(x), "term.labels"))
    len_ran_ef <- lapply(term_labels, function(x) which(x == x[grep("|", x, fixed = TRUE)]))
    formYx <- mapply(function(x, y) drop.terms(terms(x), dropx = y, keep.response = TRUE),
                     x = formYx, 
                     y = len_ran_ef, 
                     SIMPLIFY = FALSE)
    mfX <- lapply(formYx, function(x) model.frame(terms(x), data = data)) # model.frame(terms(formYx), data = data)
    TermsX <- lapply(mfX, function(x) attr(x, "terms"))   #attr(mfX, "terms")
    components <- lmeObject$components
    formtest <- components[grep("TermsZ", names(components), fixed = TRUE)]
    formYz <- lapply(formtest, function(x)  formula(x)) #formula(lmeObject$modelStruct$reStruct[[1]])
    mfZ <- lapply(formYz, function(x) model.frame(x, data = data)) #model.frame(terms(formYz), data = data)
    TermsZ <- lapply(mfZ, function(x) attr(x, "terms"))   #attr(mfZ, "terms")
    mcmcMeans <- lmeObject$postMeans
    betas <- mcmcMeans[grep("beta", names(mcmcMeans), fixed = TRUE)] #fixef(lmeObject)
    sigma <- mcmcMeans[grep("sigma", names(mcmcMeans), fixed = TRUE)] #sigma <- lmeObject$sigma
    D <- mcmcMeans$D #lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", sigma^2)[[1]]
    mcmc <- lmeObject$mcmc
    betas_mcmc <- mcmc[grep("betas", names(mcmc), fixed = TRUE)] 
    V <- var(do.call(cbind, betas_mcmc)) #V <- vcov(lmeObject)
    times_orig <- lapply(timeVar, function(x) data[[x]]) #data[[timeVar]]
    times_orig <- lapply(times_orig, function(x) x[!is.na(x)]) #times_orig[!is.na(times_orig)]
    
  } else {
    formYx <- lmeObject$formYx
    TermsX <- lmeObject$TermsX
    formYz <- lmeObject$formYz
    TermsZ <- lmeObject$TermsZ 
    idVar <- lmeObject$idVar
    betas <- lmeObject$betas
    sigma <- lmeObject$sigma
    D <- lmeObject$D
    V <- lmeObject$V
    times_orig <- lmeObject$times_orig
  }
  
  ##########################################
  #### drop missing values from newdata ####
  ##########################################
  
  all_vars <- mapply(function(x, y)  unique(c(all.vars(x), all.vars(y))),
                     x = TermsX,
                     y = TermsZ, 
                     SIMPLIFY = FALSE) #unique(c(all.vars(TermsX), all.vars(TermsZ)))
  newdata_nomiss <- lapply(all_vars, function(x) newdata[complete.cases(newdata[unlist(x)]), ])
  mfX_new <- mapply(function(x, y) model.frame(x, data = y),
                    x = TermsX, 
                    y = newdata_nomiss,
                    SIMPLIFY = FALSE) #model.frame(TermsX, data = newdata_nomiss)
  X_new <- mapply(function(x, y) model.matrix(x, y),  
                  x = formYx,
                  y = mfX_new, 
                  SIMPLIFY = FALSE) #model.matrix(formYx, mfX_new)
  mfZ_new <- mapply(function(x, y) model.frame(x, data = y),
                    x = TermsZ, 
                    y = newdata_nomiss,
                    SIMPLIFY = FALSE)  #model.frame(TermsZ, data = newdata_nomiss)
  Z_new <- mapply(function(x, y) model.matrix(as.formula(x), y),    
                  x = formYz,
                  y = mfZ_new, 
                  SIMPLIFY = FALSE) #model.matrix(formYz, mfZ_new)
  na_ind <- lapply(mfX_new, function(x) attr(mfX_new, "na.action")) #attr(mfX_new, "na.action")
  y_new <- lapply(mfX_new, function(x) model.response(x, "numeric")) #model.response(mfX_new, "numeric")
  if (length(idVar) > 1)
    stop("the current version of the function only works with a single grouping variable.\n")
  if (is.null(newdata[[idVar]]))
    stop("subject id variable not in newdata.")
  id_nomiss <- lapply(newdata_nomiss, function(x) match(x[[idVar]], unique(x[[idVar]]))) #match(newdata_nomiss[[idVar]], unique(newdata_nomiss[[idVar]]))
  n <- length(unique(id_nomiss))
  
  ############################################
  ##### Random Effects - Empirical Bayes ##### 
  ############################################
  
  modes <- matrix(0.0, n, dim(D)[1])  #matrix(0.0, n, ncol(Z_new))
  
  post_vars <- DZtVinv <- vector("list", n)
  for (i in seq_len(n)) {
    id_i <- lapply(id_nomiss, function(x) x == i) #id_i <- id_nomiss == i
    X_new_id <- mapply(function(x, y) x[y, , drop = FALSE],
                       x = X_new, 
                       y = id_i, 
                       SIMPLIFY = FALSE)  #X_new[id_i, , drop = FALSE]
    Z_new_id <- mapply(function(x, y) x[y, , drop = FALSE],
                       x = Z_new, 
                       y = id_i, 
                       SIMPLIFY = FALSE)  #Z_new[id_i, , drop = FALSE]
    Z_new_id <- JMbayes:::bdiag(Z_new_id)
    Sigmas <- mapply(function(x, y) x^2*diag(sum(id_i[[1]])),
                     x = sigma, 
                     y = id_i, 
                     SIMPLIFY = FALSE) 
    Sigmas <- JMbayes:::bdiag(Sigmas)
    Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) + 
                      Sigmas) #solve(Z_new_id %*% tcrossprod(D, Z_new_id) + sigma^2 * diag(sum(id_i)))
    DZtVinv[[i]] <- tcrossprod(D, Z_new_id) %*% Vi_inv
    Xbetas <- mapply(function(x, y, z) x - y %*% z,  
                     x = y_new, 
                     y = X_new_id, 
                     z = betas, 
                     SIMPLIFY = FALSE)
    Xbetas <- do.call(rbind, Xbetas)
    modes[i, ] <- as.vector(DZtVinv[[i]] %*% Xbetas) #c(DZtVinv[[i]] %*% (y_new[id_i] - X_new_id %*% betas))
    t1 <- DZtVinv[[i]] %*% Z_new_id %*% D
    X_new_id <- JMbayes:::bdiag(X_new_id)
    t2 <- DZtVinv[[i]] %*% X_new_id %*% V %*% 
      crossprod(X_new_id, Vi_inv) %*% Z_new_id %*% D
    post_vars[[i]] <- D - t1 + t2
  }
  
  ##############################
  #### Obtain fitted values ####
  ##############################
  
  XB <- mapply(function(x, y) x %*% y,
               x = X_new,
               y = betas,
               SIMPLIFY = FALSE)
  
  dims <- lapply(Z_new, function(x) dim(x)[2])
  end <- cumsum(dims)
  Z_new_AllCol <- do.call(cbind, Z_new)
  start <- which( colnames(Z_new_AllCol)=="(Intercept)" )
  modes_list <- mapply(function(x, y, z) modes[z, x:y],
                       x = start,
                       y = end,
                       z = id_nomiss,
                       SIMPLIFY = FALSE)
  Zb <- mapply(function(x, y) rowSums(x * y),
               x = Z_new,
               y = modes_list,
               SIMPLIFY = FALSE)
  fitted_y <- mapply(function(x, y) x + y,
                     x = XB,
                     y = Zb,
                     SIMPLIFY = FALSE)  #c(X_new %*% betas) + rowSums(Z_new * modes[id_nomiss, , drop = FALSE])
  
  ############################
  #### Obtain predictions ####
  ############################
  
  if ( any(is.null(times)) || any(sapply(times, function(x) !is.numeric(x))) ) {
    times <- lapply(times_orig, function(x) seq(min(x), max(x), length.out = 100))   #seq(min(times_orig), max(times_orig), length.out = 100)
  }
  id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
  last_time <- lapply(timeVar, function(x) tapply(newdata[[x]], id, max)) #tapply(newdata[[timeVar]], id, max)
  times_to_pred <- mapply(function(x, y) if (all_times) x else x[x > c(y)], 
                          x = times, 
                          y = last_time, 
                          SIMPLIFY = FALSE)
  
  #lapply(last_time, function (t) if (all_times) times else times[times > t])
  id_pred <- lapply(times_to_pred, function(x) rep(seq_len(n), each = sum(sapply(x, length)*n))) #rep(seq_len(n), sapply(times_to_pred, length))
  
  ##newdata_pred <- newdata_pred[id_pred, ]
  #right_rows(newdata, newdata[[timeVar]], id, times_to_pred)
  #newdata_pred[[timeVar]] <- unlist(times_to_pred)
  newdata_pred <- list()
  for (i in 1:length(TermsX)){
    #f <- lmeObject$call$formulas[[i+1]]
    #f <- drop.terms(terms(as.formula(f)), dropx = length(terms(as.formula(f))), keep.response = TRUE)
    #newdata_new <- model.frame(f, data = newdata)
    
    #newdata_new <- tail(newdata_new, n = 1)
    #timeVar_group <- timeVar[[i]]
    #newdata_new <- newdata_new[rep(1, length(times_to_pred[[i]])), , drop = FALSE]
    #newdata_new[[timeVar_group]] <- unlist(times_to_pred[[i]])
    #newdata_pred[[i]] <- newdata_new
    
    newdata_new <- newdata[unique(unlist(all_vars))]
    newdata_new <- tail(newdata_new, n = 1)
    timeVar_group <- timeVar[[i]]
    n_pred <- length(unlist(times_to_pred[[i]]))
    newdata_new <- newdata_new[rep(1, times = n_pred), ]
    newdata_new[[timeVar_group]] <- unlist(times_to_pred[[i]])
    
    #f <- lmeObject$call$formulas[[i+1]]
    #f <- drop.terms(terms(as.formula(f)), dropx = length(terms(as.formula(f))), keep.response = TRUE)
    #newdata_pred[[i]] <- model.frame(f, data = newdata_new)
    newdata_pred[[i]] <- newdata_new
  }
  
  mfX_new_pred <- mapply(function(x, y) model.frame(x, data = y, na.action = NULL),
                         x = TermsX,
                         y = newdata_pred,
                         SIMPLIFY = FALSE) #model.frame(TermsX, data = newdata_pred, na.action = NULL)
  X_new_pred <- mapply(function(x, y) model.matrix(x, y),
                       x = formYx,
                       y = mfX_new_pred,
                       SIMPLIFY = FALSE) #model.matrix(formYx, mfX_new_pred)
  mfZ_new_pred <- mapply(function(x, y) model.frame(x, data = y, na.action = NULL),
                         x = TermsZ,
                         y = newdata_pred,
                         SIMPLIFY = FALSE) #model.frame(TermsZ, data = newdata_pred, na.action = NULL)
  Z_new_pred <- mapply(function(x, y) model.matrix(as.formula(x), y),
                       x = formYz,
                       y = mfZ_new_pred,
                       SIMPLIFY = FALSE) #model.matrix(formYz, mfZ_new_pred)
  XB_pred <- mapply(function(x, y) x %*% y,
                    x = X_new_pred,
                    y = betas,
                    SIMPLIFY = FALSE)
  modes_list_pred <- mapply(function(x, y, z) modes[z, x:y],
                            x = start,
                            y = end,
                            z = id_pred,
                            SIMPLIFY = FALSE)
  Zb_pred <- mapply(function(x, y) rowSums(x * y),
                    x = Z_new_pred,
                    y = modes_list_pred,
                    SIMPLIFY = FALSE)
  predicted_y <- mapply(function(x, y) x + y,
                        x = XB_pred,
                        y = Zb_pred,
                        SIMPLIFY = FALSE)   #c(X_new_pred %*% betas) + rowSums(Z_new_pred * modes[id_pred, , drop = FALSE])
  
  
  ####################################################
  #### Obtain confidence intervals of predictions #### 
  ####################################################
  
  set.seed(seed)
  betas_M <- MASS::mvrnorm(M, unlist(betas), V) #MASS::mvrnorm(M, betas, V)
  #modes_fun <- function (betas) {
  #   t(mapply("%*%", DZtVinv, split(y_new - X_new %*% betas, id_nomiss)))
  #}
  #modes_M <- lapply(split(betas_M_list[[i]], row(betas_M_list[[i]])), modes_fun)
  modes_fun <- function(betas) {
    Xbetas <- mapply(function(x, y, z) x - y %*% z,
                     x = y_new,
                     y = X_new,
                     z = betas,
                     SIMPLIFY = FALSE)
    Xbetas <- do.call(rbind, Xbetas)
    (DZtVinv[[1]] %*% Xbetas)
  }
  dims <- lapply(X_new, function(x) dim(x)[2])
  end <- cumsum(dims)
  Z_new_AllCol <- do.call(cbind, X_new)
  start <- which(colnames(Z_new_AllCol)=="(Intercept)")
  betas_M_list <- mapply(function(x, y) betas_M[, x:y],
                         x = start,
                         y = end,
                         SIMPLIFY = FALSE)
  modes_M <- list()
  for (m in 1:M){
    betas_Mi <- lapply(betas_M_list, "[", m , )
    modes_M[[m]] <- modes_fun(betas_Mi)
  }
  modes_M <- lapply(modes_M, function(x) as.vector(x))
  #matrix_row <- function (m, i) m[i, , drop = FALSE]
  #modes_M <- lapply(seq_len(n), function (i) t(sapply(modes_M, matrix_row, i = i)))
  modes_M <- do.call(rbind, modes_M)
  b_M <- modes_M <- list(modes_M)
  for (i in seq_len(n)) {
    b_M[[i]] <- t(apply(modes_M[[i]], 1, MASS::mvrnorm, n = 1, Sigma = post_vars[[i]]))
  }
  n_pred <- lapply(predicted_y, function(x) length(x)) #length(predicted_y)
  sampled_y <- lapply(n_pred, function(x) matrix(0.0, x, M)) #matrix(0.0, n_pred, M)
  for (m in seq_len(M)) {
    betas_m <- lapply(betas_M_list, "[", m , ) #betas_M[m, ]
    b_m <- t(sapply(b_M, function (x) x[m, ]))
    #mean_m <- c(X_new_pred %*% betas_m) + 
    #  rowSums(Z_new_pred * b_m[id_pred, , drop = FALSE])
    XB_m <- mapply(function(x, y) x %*% y,
                   x = X_new_pred,
                   y = betas_m,
                   SIMPLIFY = FALSE)
    start <- grep("(Intercept)", colnames(b_m), fixed = TRUE)
    dims <- lapply(Z_new, function(x) dim(x)[2])
    end <- cumsum(dims)
    b_m_list <- mapply(function(x, y, z) b_m[z, x:y],
                       x = start,
                       y = end,
                       z = id_pred,
                       SIMPLIFY = FALSE)
    Zb_m <- mapply(function(x, y) rowSums(x * y),
                   x = Z_new_pred,
                   y = b_m_list,
                   SIMPLIFY = FALSE)
    mean_m <- mapply(function(x, y) x + y,
                     x = XB_m,
                     y = Zb_m,
                     SIMPLIFY = FALSE)
    len_pred <- lapply(mean_m, function(x) length(x))
    outcomes <- lapply(formYx, function(x) all.vars(x)[1])
    #outcomes <- names(lmeObject$components)[grep("y", names(lmeObject$components), fixed = TRUE)]
    Outcomes <- rep(unlist(outcomes), times = unlist(len_pred))
    sampled_y[[m]] <- if (interval == "confidence") data.frame(outcome = Outcomes, prediction = do.call(rbind, mean_m))#mean_m 
    else data.frame(outcome = Outcomes, prediction = do.call(rbind, mapply(function(x, y, z) as.matrix(rnorm(x, y, z)),
                                                                           x = n_pred, y = mean_m, z = sigma, SIMPLIFY = FALSE))) #rnorm(n_pred, mean_m, lmeObject$sigma)
  }
  low <- lapply(seq_len(length(Outcomes)), function(x) {
    pred_m <- lapply(sampled_y, "[", x , )
    pred_m <- do.call(rbind, pred_m)
    quantile(pred_m[,2],  probs = (1 - level) / 2)
  })
  low <- do.call(rbind, low)
  colnames(low) <- "low"
  upp <- lapply(seq_len(length(Outcomes)), function(x) {
    pred_m <- lapply(sampled_y, "[", x , )
    pred_m <- do.call(rbind, pred_m)
    quantile(pred_m[,2],  probs = 1 - (1 - level) / 2)
  })
  upp <- do.call(rbind, upp)
  colnames(upp) <- "upp"
  #low <- apply(sampled_y, 1, quantile, probs = (1 - level) / 2)
  #upp <- apply(sampled_y, 1, quantile, probs = 1 - (1 - level) / 2)
  
  ###################################
  #### Obtain output of function ####
  ################################### 
  
  rm(list = ".Random.seed", envir = globalenv())
  if (!return_data) {
    #list(times_to_pred = times_to_pred, predicted_y = predicted_y, low = low, upp = upp)
    res <- data.frame(id_pred = do.call(rbind, lapply(id_pred, function(x) as.matrix(x))), Outcomes = Outcomes, 
                      times_to_pred = do.call(rbind, lapply(times_to_pred, function(x) as.matrix(x))),
                      predicted_y = do.call(rbind, predicted_y), low = low, upp = upp)
    res
  } else {
    res <- data.frame(id_pred = do.call(rbind, lapply(id_pred, function(x) as.matrix(x))), Outcomes = Outcomes, 
                      times_to_pred = do.call(rbind, lapply(times_to_pred, function(x) as.matrix(x))),
                      predicted_y = do.call(rbind, predicted_y), low = low, upp = upp)
    
    var_names <- lapply(formula(lmeObject), function(x) all.vars(x))
    newdat <- lapply(var_names, function(x) newdata[, x])
    newdat <- mapply(cbind, newdat, "Outcome"= unique(Outcomes), SIMPLIFY = F)
    colnames <- c("y","group", "year", "id", "Outcomes")
    newdat <- lapply(newdat, setNames, colnames)
    newdat <- do.call(rbind, newdat)
    
    res <- list(predictions = res, dataset = newdat)
    res
    #out_data <- rbind(newdata, newdata_pred)
    #out_data$pred <- c(fitted_y, predicted_y)
    #out_data$low <- c(rep(NA, length(fitted_y)), low)
    #out_data$upp <- c(rep(NA, length(fitted_y)), upp)
    #out_data[order(out_data[[idVar]], out_data[[timeVar]]), ]
  }
}
