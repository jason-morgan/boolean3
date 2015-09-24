## ----------------------------------------------------------------------------
## Test file for boolean package. VERY ROUGH. UNIT TESTS NEEDED.
## ----------------------------------------------------------------------------

library(devtools)
setwd("~/lib/R")
load_all("boolean3")

## Load the library.
## library(boolean3)
library(optimx)
library(lattice)
library(snow)
library(rlecuyer)
library(mvtnorm)
library(numDeriv)
library(rgenoud)
## library(MCMCpack)
## library(glogis)

## Set up test dataset.
set.seed(12345)
N  <- 2000
Df <- cbind(1, rmvnorm(N, mean=rep(0, 5)))

## Set coefficients
beta.a <- c(-2.00, 0.33, 0.66, 1.00)
beta.b <- c(0.00, 1.50, -0.25)

## Generate path probabilities following a normal model.
y.a <- as.vector(pnorm(tcrossprod(beta.a, Df[, 1:4])))
y.b <- as.vector(pnorm(tcrossprod(beta.b, Df[, c(1, 5, 6)])))

## AND and OR-model
or <- function(x, y) { x + y - x * y }
and <- function(x, y) { x * y }
y.star.OR  <- or(y.a, y.b)
y.star.AND <- and(y.a, y.b)

## Observed responses
y.OR <- rbinom(N, 1, y.star.OR)
y.AND <- rbinom(N, 1, y.star.AND)

## Set up data.frame for estimation
Df <- cbind(1, Df)
Df <- as.data.frame(Df)
Df[,1] <- y.OR
Df[,2] <- y.AND
names(Df) <- c("y.OR", "y.AND", "x1", "x2", "x3", "x4", "x5")


## ----------------------------------------------------------------------------
## Tests of boolean3
## ----------------------------------------------------------------------------

## Set up model, note the new syntax.

## OR model
mod.OR <- boolprep(y.OR ~ (a | b), a ~ x1 + x2 + x3, b ~ x4 + x5,
                   data = Df, family=list(binomial("probit")))

## AND model
mod.AND <- boolprep(y.AND ~ (a & b), a ~ x1 + x2 + x3, b ~ x4 + x5,
                    data = Df, family=list(binomial("probit")))


## The link function used in the model (family = binomial(logit) by default) is
## stored in the boolean object returned by this function. The link function
## now calls built-in C code, which is much faster. A pre-parsed equation for
## the boolean portion of the likelihood is also stored in the boolean
## object. The purpose of storing these in the object is that it facilitates
## easy post-estimation procedures that should always be consistent with the
## original model.
str(mod.OR$link)
str(mod.AND$link)
mod.OR$model$boolmod

## It's also possible to specify a different link function for each path. These
## should be in the same order as specified in the model formula.
mod.OR2 <- boolprep(y.OR ~ (a | b), a ~ x1 + x2 + x3, b ~ x4 + x5,
                    data = Df,
                    family=list(binomial("probit"), binomial("logit")))
str(mod.OR2$link)

## Fit the model using the nlminb optimizer.
(fit.OR <- boolean(mod.OR, method="nlminb", control=list(trace=1)))


cbind(est=coef(fit.OR), truth=c(beta.a, beta.b))

## Fit with the nlminb and nlm optimizer.
(fit1.OR <- boolean(mod.OR, method="nlm"))

## Re-fit, with BFGS and a higher maximum iterations. All of the options that
## go along with nlm(), optim(), and genoud() should be transparently passable.
(fit2.OR <- boolean(mod.OR, method="BFGS", control = list(maxit = 500)))

## Induce a convergence warning message.
(fit3.OR <- boolean(mod.OR, method="BFGS", control = list(maxit = 5)))

## Using MCMC. Then take a more detailed look at the results using tools from
## coda package. BROKEN SPECIFICATION
## (fit4.OR <- boolean(mod.OR, method="mcmc", mcmc=30000))

## Using all methods provided by optimx, increase verbosity. MCMC, SANN, and
## genoud are not used when all.methods=TRUE.
## (fit5.OR <- boolean(mod.OR, control=list(trace=0, all.methods=TRUE)))

## Using genoud.
(fit6.OR <- boolean(mod.OR, method="genoud", cluster=c("localhost", "localhost"),
                    print.level=2))
(fit6.OR <- boolean(mod.OR, method="genoud", print.level=2))

## Using SANN.
(fit7.OR <- boolean(mod.OR, method="SANN"))

## The fit is stored as "model.fit", within the boolean object. The Hessian is
## always saved and must be inverted to get the covariance matrix (the summary
## function below does this for you). *Important note*: all optimizers are set
## to *minimize* the negative of the log-likelihood. This is the preferred
## method for many optimizers, so it was easiest simply to calc -1*logLik at
## presentation. This does, however, affect the sign of the Hessian saved in the
## model.fit object; i.e., you don't want to multiply by -1 when calculating
## standard errors. If this becomes an issue, I can simply store the Hessian
## returned by the optimizers as -1 * Hessian.
str(fit.OR$model.fit)

## Create a summary object, saving and printing it. Then take a look at the
## objects stored in the summary object.
(smry <- summary(fit.OR))
str(smry)

## Extract log-likelihood and coefficient vector.
logLik(fit.OR)
coef(fit.OR)
smry$methods$nlminb$coef

## Now extract profile likelihoods (actually, estimated likelihoods, if I
## understand the likelihood theory -- need to check this) for all of the
## covariates. Since the plot is a lattice object, it can be stored as well. The
## default is to simply plot it, but storing it allows alternative plots to be
## made without recalculating the profiles.
(prof <- boolprof(fit.OR))

## Now, extract the plots for x1_a and x4_b.
plot(prof, y = c("x1_a", "x4_b"))
plot(prof, y = c(1, 3), scales = list(y = list(relation = "free")))

## You can use variable or index matching with boolprof to select particular
## covariates of interest.
boolprof(fit.OR, vars = c(1, 3))
boolprof(fit.OR, vars = c("x1_a", "x4_b"))

## Plots of predicted probabilities are available with boolprob. With boolprob,
## either vars or newdata *must* be specified.
boolprob(fit.OR, vars = c("x1_a", "x4_b"))
boolprob(fit.OR, vars = c(2, 3, 4, 6))
boolprob(fit.OR, vars = c(2:4, 6, 7))

## Specifying conf.int = TRUE gets you simulated confidence intervals. The
## number of samples to pull from the distribution of the estimated
## coefficients is controlled by n, which n=100 is default.
(prob <- boolprob(fit.OR, vars = c(2, 3, 4, 6), n = 100, conf.int = TRUE))

## Choose a different method estimate upon which to base the estimates.
(prob <- boolprob(fit.OR, method="spg", vars=c(2, 3, 4, 6), n=100,
                  conf.int=TRUE))

## As with the other components of the model, you can extract the predicted
## probabilities.
str(prob)
prob$est

## Inspect MCMC mixing.
plot(fit4$model.fit$mcmc$mcmc)

## Bootstrapping models.
(bs <- boolboot(fit.OR, n=10))

b <- coef(fit.OR)
se <- summary(fit.OR)$methods$nlminb$coef[,2]
lower <- b - 4*se
upper <- b + 4*se
cbind(lower, upper)

(bs <- boolboot(fit.OR, n=10))
cbind(lower, upper)

summary(fit.OR)

(bs <- boolboot(fit.OR, n=20, cluster=c("localhost", "localhost")))
