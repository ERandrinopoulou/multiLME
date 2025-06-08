\name{mv_lme}
\alias{mv_lme}
\title{Multivariate mixed-effects model}
\description{
  This function fits a multivariate mixed-effects model, assuming different association structures between the longitudinal outcomes.

}
\usage{
  mv_lme(formulas = formulas, data = data, 
                     families = families, 
                     corr_RE = FALSE,
                     assoc = TRUE,
                     assoc_from = 2,
                     assoc_to = 1,
                     extraForm = NULL,
                     time_var = "year")
}
\arguments{
  \item{formulas}{A list of formulas specifying the longitudinal submodels (excluding the association structure).}
  \item{data}{The dataset containing all variables used in the models.}
  \item{families}{Distributions for the outcomes (currently only “gaussian” is supported).}
  \item{hc}{Logical; whether hierarchical centering is applied. Default is FALSE. (Note: TRUE option under development)}
  \item{predicted}{Logical; whether to save predicted values. Default is FALSE.}
  \item{corr_RE}{Logical; whether to assume correlated random effects via a shared variance-covariance matrix. Default is TRUE.}
  \item{engine}{Currently, only "JAGS" is supported.}
  \item{parameters.b}{Logical; whether to save the random effects. Default is FALSE.}
  \item{assoc}{Logical; whether to include an association parameter between longitudinal outcomes. Default is TRUE.}
  \item{assoc_from}{Index of the longitudinal outcome used as a covariate in the association structure.}
  \item{assoc_to}{Index of the longitudinal outcome model to which other outcomes are added as covariates.}
  \item{extraForm}{Specifies the functional form of the association. Default is NULL. If NULL, the underlying value association structure is assumed (see examples below for other structures).}
  \item{extraForm_tr}{Specifies the functional form for the truncated period of the association. Default is NULL.}
  \item{norm_area}{Logical; whether to use the normalized area under the curve. Default is FALSE. (Note: feature under development for predictions)}
  \item{time_var}{The time variable used in the longitudinal models.}
  \item{control}{Control Parameters (via jagsUI::jags):}
  \item{n.processors}{parallel::detectCores() - 1}
  \item{working.directory}{getwd()}
  \item{clear.model}{TRUE}
  \item{seed}{1L}
  \item{optimize_only}{FALSE}
  \item{verbose}{FALSE}
  \item{n.iter}{28000L}
  \item{n.burnin}{3000L}
  \item{n.thin}{50L}
  \item{n.adapt}{3000L}
  \item{n.chains}{2L}
  \item{n.cores}{1L}
}
\value{
  Results can be accessed via the fit object. 
  
  Posterior means are available through fit$postMeans.
  
  Corresponding 95\% credible intervals, computed from the 2.5th and 97.5th percentiles of the posterior distributions, can be retrieved using fit$CIs.
  
  Additional output elements, such as mcmc.info, Deviance Information Criterion (DIC), effective number of parameters (pD), and convergence diagnostics (Rhat), are returned by jagsUI::jags() function.
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
}