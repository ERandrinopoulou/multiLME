\name{DynPred_mv_lme}
\alias{DynPred_mv_lme}
\title{Dynamic predictions from a multivariate mixed-effects model}
\description{
  This function obtains dynamic predictions from a multivariate mixed-effects model, assuming different association structures between the longitudinal outcomes.

}
\usage{
  DynPred_mv_lme(object, newdata, families, hc, level = 0.95, IdVar,
    timeVar, M = 200, times = NULL, assoc = TRUE, assoc_from,
    assoc_to, extraForm = NULL, extraForm_tr = NULL)
}
\arguments{
  \item{object}{An mv_lme object.}
  \item{newdata}{The dataset of the new subject to predict.}
  \item{families}{Distributions for the outcomes (currently only “gaussian” is supported).}
  \item{hc}{Logical; whether hierarchical centering is applied. Default is FALSE. (Note: TRUE option under development)}
  \item{level}{The confidence interval level. Default is 0.95.}
  \item{IdVar}{The variable name of the Id column.}
  \item{time_var}{The time variable used in the longitudinal models.}
  \item{M}{The number of posterior samples. Default is 200.}
  \item{times}{The times to obtain future predictions. Default is NULL. In that case the maximum of a sequence of times until the maximum time in the whole data set will be assumed.}
  \item{assoc}{Logical; whether to include an association parameter between longitudinal outcomes. Default is TRUE.}
  \item{assoc_from}{Index of the longitudinal outcome used as a covariate in the association structure.}
  \item{assoc_to}{Index of the longitudinal outcome model to which other outcomes are added as covariates.}
  \item{extraForm}{ Specifies the functional form of the association. Default is NULL. If NULL, the underlying value association structure is assumed (see examples below for other structures).}
  \item{extraForm_tr}{Specifies the functional form for the truncated period of the association. Default is NULL.

}
}
\value{
 Returns a data frame with the mean predictions together with the prediction intervals.
}
\examples{
  # Simulate data
  n = 100

  K <- 10  # number of planned repeated measurements per subject, per outcome
  t.max <- 10 # maximum follow-up time

  # parameters for the linear mixed effects model 1
  betaa1 <- c("(Intercept)" = 3.01, "Time1" = 0.47, "Group" = -2)
  sigma1.y <- 1.017671  # measurement error standard deviation

  # parameters for the linear mixed effects model 2
  betaa2 <- c("(Intercept)" = 2.15, "Time1" = 0.91)
  sigma2.y <- 1  # measurement error standard deviation

  # association parameter
  alphaa <- 0.9 # association parameter - value

  # Variance-covariance matrix for random effects
  invDvec <- c(1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1)

  invDmat <- matrix(invDvec, 4, 4)
  D <- solve(invDmat)

  # design matrices for the longitudinal measurement model
  times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max)))))
  group <- sample(c(0,1), n, replace = TRUE)

  DF <- data.frame(year = times, drug = factor(rep(group, each = K)))

  X1 <- model.matrix(~ year + drug, data = DF)
  Z1 <- model.matrix(~ year, data = DF)

  X2 <- model.matrix(~ year, data = DF)
  Z2 <- model.matrix(~ year, data = DF)

  #simulate random effects
  library(MASS)

  b <- mvrnorm(n, rep(0, nrow(D)), D)

  # simulate longitudinal responses
  id <- rep(1:n, each = K)

  eta.y2 <- as.vector(X2 %*% betaa2 + rowSums(Z2 * b[id, 3:4]))
  y2 <- rnorm(n * K, eta.y2, sigma2.y)

  eta.y1 <- as.vector(X1 %*% betaa1 + rowSums(Z1 * b[id, 1:2])) + alphaa * eta.y2
  y1 <- rnorm(n * K, eta.y1, sigma1.y)

  dat <- DF[, ]
  dat$id <- id
  dat$y1 <- y1
  dat$y2 <- y2

  names(dat) <- c("year", "drug", "id", "y1", "y2")
  data <- dat
  # Fit models
  formulas <- list(y1 ~ year + drug + (year | id),
                   y2 ~ year + (year | id))
  families <- c("gaussian", "gaussian")


  fit <- mv_lme(formulas = formulas, data = data,
                families = families,
                corr_RE = FALSE,
                assoc = TRUE,
                assoc_from = 2,
                assoc_to = 1,
                extraForm = NULL,
                time_var = "year")
                object = fit

  newdata_all = data[data$id == 1, ]
  newdata <- newdata_all[1:6,]


  res_pred <- DynPred_mv_lme(object = fit, newdata = newdata,
                             families,
                             level = 0.95, IdVar = "id", timeVar = "year",
                             M = 200, times = NULL,
                             assoc = TRUE,
                             assoc_from = 2,
                             assoc_to = 1,
                             extraForm = NULL)



  p1 <- xyplot(res_pred$Mean1 + res_pred$Lower1 + res_pred$Upper1 ~ res_pred$Time,
         xlim = c(0, 10),
         ylim = c(min(y1), max(y1)),
         strip = FALSE, xlab = "Year", ylab = "Y1",
         lwd = c(3, 2, 2), lty = c(1, 3, 3), type = "l", col = 1,
         panel = function(x, y,  ...) {
           panel.xyplot(x, y, ...)
           panel.abline(v = max(newdata$year), col = 2, lty = 2)
           panel.points(newdata_all$year, newdata_all$y1, pch = 16, cex = 0.8, col = "blue")
           })

  print(p1)
}
