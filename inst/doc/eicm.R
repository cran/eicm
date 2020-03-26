## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(timeit = local({
  now = NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res = difftime(Sys.time(), now, units="secs")
      now <<- NULL
      # use options$label if you want the chunk label as well
      if(res > 60)
	      sprintf("<p style=\"text-align:right\">-> running time of the above code chunk: %.1f minutes</p>", res / 60)
      else
	      sprintf("<p style=\"text-align:right\">-> running time of the above chunk: %.1f seconds</p>", res)
    }
  }})
)

library(eicm)
set.seed(2)		# for vignette reproducibility

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # fit & select an EICM with 2 latent variables
#  m <- eicm(occurrences, n.latent=2)
#  
#  # display estimated coefficients (note they are organized in matrices)
#  coef(m$selected.model)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  m <- eicm(occurrences, n.latent=2, regularization=c(6, 0.5), penalty=4, theta.threshold=0.5)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # excluding interactions with the formula syntax
#  m <- eicm(occurrences, forbidden=list(
#      sp3 ~ sp4 + sp5,		# sp3 must not be affected by sp4 nor sp5
#      sp4 ~ .,				# sp4 must not be affected by any other
#      sp1 ~ . - sp8			# sp1 must not be affected by any other except sp8
#      ))
#  
#  # display species interaction coefficients
#  # note the zeroed coefficients are those that were excluded
#  coef(m$fitted.model)$sp

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # Excluding interactions with the matrix syntax
#  
#  # create a square matrix species x species, all zeroes
#  mask <- matrix(0, nrow=ncol(occurrences), ncol=ncol(occurrences))
#  
#  # set to 1 those interactions we want to include
#  mask[4, 2] <- 1		# species #2 may affect species #4
#  mask[6, 1] <- 1		# species #1 may affect species #6
#  mask[, 7] <- 1		# species #7 may affect any other
#  mask[2, ] <- 1		# species #2 may be affected by any other
#  
#  m <- eicm(occurrences, mask.sp=mask)
#  
#  # display species interaction coefficients
#  # note the zeroed coefficients are those that were excluded
#  coef(m$fitted.model)$sp

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # exclude interactions departing from species with 20 or less occurrences
#  m <- eicm(occurrences, mask.sp=mask, exclude.prevalence=20)

## ---- echo = TRUE, message = FALSE, eval = TRUE, results="hide", timeit=TRUE, cache=TRUE----
# Load the included true model (32 species, 30 interactions, 2 environmental predictors)
data(truemodel)

# realize the model
occurrences <- predict(truemodel, nrepetitions=1)

# Pre-define a network topology
# For illustrative purposes, let's assume we know the true network topology so we build
#  the mask from the true model coefficients:
# 0: don't estimate interaction
# 1: estimate interaction
mask <- ifelse(truemodel$model$sp == 0, 0, 1)

# And now we estimate their values
# here we discard all the predictors - note that we don't provide any predictor matrix,
# only the presence/absence data.
m <- eicm(occurrences, n.latent=2, mask.sp=mask, do.selection=FALSE,
	regularization=c(1, 0.1))

## ---- echo = TRUE, eval = TRUE, message = FALSE, results="hide", fig.cap="Estimated interaction network with a given topology. Blue: positive interactions; red: negative interactions. Line thickness is proportional to interaction strength."----
# Plot estimated species interaction network
plot(m, type="network")

## ---- echo = TRUE, eval = TRUE, message = FALSE, results="hide", fig.width=10, fig.asp=0.4, fig.fullwidth=TRUE, fig.keep="last", fig.show="hold", fig.cap="Estimated versus true coefficients of a fitted EICM model over simulated data, with the true environmental predictors dropped, and estimated as latent variables."----
# Plot the network topology comparison between true and estimated models
# (plot omitted from vignette)
plot(m, true.model=truemodel, type="network")

# Plot estimated versus true coefficients
# Note that the environmental coefficients here relate to estimated latent variables,
# not the true predictors.
plot(m, true.model=truemodel, nenv.to.plot=2, legend=FALSE)

## ---- echo = TRUE, message = FALSE, eval = TRUE-------------------------------
# display estimated species interaction coefficients (head only)
# the zeroed coefficients are those that were excluded a priori
head(round(coef(m$fitted.model)$sp, 3))

## ---- echo = TRUE, message = FALSE, eval = TRUE, results="hide", timeit=TRUE, cache=TRUE, fig.width=10, fig.asp=0.4, fig.fullwidth=TRUE, fig.cap="Estimated versus true coefficients of a fitted EICM model over simulated data in the presence of a large amount of NAs in the response matrix. Note that the true environmental predictors were dropped, and estimated as latent variables, and that their estimation is accurate even in the presence of a large amount of NAs in the occurrence data."----
# Randomly set 1/2 of the occurrence data to NA
# (8000 records out of 16000)
occurrences[sample(length(occurrences), 8000)] <- NA

# Fit the model
m <- eicm(occurrences, n.latent=2, mask.sp=mask, do.selection=FALSE,
	regularization=c(1, 0.1))

plot(m, true.model=truemodel, nenv.to.plot=2, legend=FALSE)

## ---- echo = TRUE, message = FALSE, eval = FALSE, timed=TRUE, fig.asp=3-------
#  # Confidence intervals are appended to the model object.
#  m <- confint(m$fitted.model)
#  plot(m)		# by default, plots CIs

## ---- echo = TRUE, message = FALSE, eval = TRUE-------------------------------
# load the included parameterized model
data(truemodel)

# for reference, plot the species interaction network
plot(truemodel, type="network")

## ---- echo = TRUE, message = FALSE, eval = TRUE-------------------------------
# Unconditional predictions
# let's fix environmental predictors at 0, for simplicity.
predict(truemodel, newdata=cbind(env01=0, env02=0))

## ---- echo = TRUE, message = FALSE, eval = TRUE-------------------------------
# Conditional predictions
# predict probabilities for all species conditional on the
# known presence of sp011 (compare sp014 and sp004 with the above)
predict(truemodel, newdata=cbind(env01=0, env02=0, sp011=1))

## ---- echo = TRUE, message = FALSE, eval = TRUE-------------------------------
# Propagation of indirect species effects
# predict probabilities for all species conditional on the known 
# absence (first line) and known presence (second line) of sp005 and sp023
newdata <- cbind(env01=0, env02=0, sp012=c(0, 1), sp018=c(0, 1))
newdata
predict(truemodel, newdata=newdata, nrep=100000)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  nenv <- 2			# 2 environmental predictors
#  nsp <- 20			# 20 species
#  nsamples <- 400		# 400 samples
#  
#  # random responses to environment
#  env.coefs <- matrix(runif((nenv + 1) * nsp, -4, 4), nrow=nsp)
#  
#  # let's define some true species interactions
#  sp.coefs <- matrix(0, nrow=nsp, ncol=nsp)
#  sp.coefs[3, 5] <- 2	 # this reads: species #5 affects species #3 with a "strength" of 2
#  sp.coefs[4, 1] <- 3
#  sp.coefs[1, 2] <- -2
#  
#  # random environmental predictors
#  env <- matrix(rnorm(nenv * nsamples), ncol=nenv, nrow=nsamples)
#  
#  # define the true model
#  true <- as.eicm(env=env, env.coefs=env.coefs, sp.coefs=sp.coefs)
#  
#  # realize the model
#  occurrences <- predict(true, nrepetitions=1)
#  
#  # we now have a plain presence/absence data matrix that follows from the true model

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # Generate a model that, when realized, shows a frequency distribution that
#  # follows a Beta distribution, yet accounting for all the desired effects
#  # (including species interactions)
#  
#  truemodel <- generateEICM(nspecies=32, nsamples=500, nenv=2,
#      ninteractions=30, shape1=1.5, shape2=3)
#  

## ---- echo = TRUE, eval = TRUE, fig.width=7, fig.height=7, fig.align="center", fig.cap="Frequency disrtibutions of four realizations of the same model, generated with generateEICM(...)"----
# plot frequency histogram of 4 realizations of truemodel:
# should approx. follow a Beta distribution.
par(mfrow=c(2, 2))
for(i in 1:4) {
    occurrences <- predict(truemodel, nrepetitions=1)
    spcounts <- apply(occurrences, 2, sum)

    hist(spcounts / nrow(occurrences), breaks=seq(0, 1, by=0.1), xlim=c(0, 1),
        main="Frequency distribution of one realization", xlab=
        "Frequency in samples",ylab="Number of species")
}

