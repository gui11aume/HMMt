BaumWelch.NB <- function (x, y = NULL, m, Q = NULL, alpha = NULL, theta.x = NULL,
   initialProb = NULL, maxiter = 500, tol = 1e-05, dig = 3)
# Author: Guillaume Filion.
# Date: October 8, 2009.
# Returns the estimate of the HMM with emissions distributed as multi dimensional
# correlated T variables.
# x is the vector (or a list of vectors) of observations.
# y the vector (or list of vectors) of observations of the control experiment.
# Q is the transition matrix. mu and sd and nu are the parameters of the model.
# p is the dimension of the T variables and m is the number of states.
# filter is a circular filter applied to x, to account for positional dependence.
# Its coefficients should be positive and add up to 1.

{

###############################################
#              OPTION PROCESSING              #
###############################################

   if (is.vector(x))
      x <- list(x)
   if (is.list(x)) {
      if (!all(sapply(X = x, FUN = is.vector))) {
         stop ("non vector element found in list x")
      }
      else if (!all(sapply(X = x, FUN = is.numeric))) {
         stop ("x must be numeric")
      }
      else {
         concat.x <- unlist(x)
      }
   }
   else {
      stop ("x must be a vector or a list")
   }

   if (!is.null(y)) {
      if (is.vector(y))
         y <- list(y)
      if (is.list(y)) {
         if (!all(sapply(X = y, FUN = is.vector))) {
            stop ("non vector element found in list y")
         }
         else if (!all(sapply(X = y, FUN = is.numeric))) {
            stop ("y must be numeric")
         }
         else if (!identical(sapply(X = x, FUN = length),
            sapply(X = y, FUN = length))) {
            stop ("lists x and y must have same structure")
         }
         else {
            concat.y <- unlist(y)
         }
      }
      else {
         stop ("y must be a vector or a list")
      }
   }
   else {
      concat.y <- rep(0, length(concat.x))
   }

   if (!is.null(Q)) {
      if (!all.equal(dim(Q), c(m,m)))
         stop("Q must be an m by m matrix")
   }
   else {
      Q <- matrix(0.1/(m-1), ncol = m, nrow = m)
      diag(Q) <- 0.9
   }

   if (is.null(alpha))
      alpha <- 1

   if (is.null(theta.x)) {
      init.x <- 2*(1:m)*mean(concat.x, na.rm = TRUE) / (m+1)
      init.y <- mean(concat.y, na.rm = TRUE)
      tilde.theta.x <- init.x / (init.x + alpha)
      tilde.theta.y <- init.y / (init.y + alpha)
      theta.x <- tilde.theta.x * (1 - tilde.theta.y) /
         (1 - tilde.theta.x*tilde.theta.y)
      theta.y <- tilde.theta.y * (1 - tilde.theta.x) /
         (1 - tilde.theta.x*tilde.theta.y)
      theta.alpha <- 1 - theta.x - theta.y
   }

   adjustInitialProb <- ifelse(is.null(initialProb), TRUE, FALSE)




###############################################
#            VARIABLE DEFINITIONS             #
###############################################


   n <- length(concat.x)
   blockSizes <- sapply(X = x, FUN = length)

   oldLogLikelihood <- -Inf

   phi <- matrix(NA_real_, nrow = n, ncol = m)

   # Tabulation saves a lot of time at the M step.
   counts <- c(sum(concat.x + concat.y == 0), tabulate(bin = concat.x
      + concat.y))
   levels <- (0:(length(counts)-1))[counts > 0]
   counts <- counts[counts > 0]
   



###############################################
#           FUNCTION DEFINITIONS              #
###############################################


   initial.steady.state.probabilities <- function () {
   # Compute the steady-state initial probabilities.

      spectralDecomposition <- eigen(t(Q))

      if (is.complex(spectralDecomposition$values[1]))
         return(rep(1/m,m))
      if (spectralDecomposition$values[1] > 0.99)
         return( spectralDecomposition$vectors[,1] /
            sum(spectralDecomposition$vectors[,1]))

      return(rep(1/m,m))

   }

   E.step <- function () {
   # Performs the E-step of the modified Baum-Welch algorithm

      # Emission probabilities are computed up to constant terms.
      for (i in 1:m)
         emissionProb[,i] <<- (theta.alpha[i]**alpha) * (theta.x[i]**concat.x) *
            (theta.y[i]**concat.y)

      # Emission probabilities of NAs are set to 1 for every states.
      emissionProb[is.na(emissionProb)] <<- 1

      if (adjustInitialProb)
         initialProb <- initial.steady.state.probabilities()

      logLikelihood <<- 0
      cumulative.n <- 0
      transitions <- matrix(double(m * m), nrow = m)

      counter = 0

      for (n.i in blockSizes) {

         counter = counter+1

         forwardback <- .Fortran("fwdb", as.integer(m), as.integer(n.i), initialProb, 
            emissionProb[(cumulative.n + 1):(cumulative.n + n.i),], Q, double(n.i),
            matrix(double(n.i * m), nrow = n.i), matrix(double(m^2), nrow = m), double(1),
            PACKAGE = "BaumWelchT")

         logLikelihood <<- logLikelihood + forwardback[[9]]
         phi[(cumulative.n + 1):(cumulative.n + n.i),] <<- forwardback[[7]]
         transitions <- transitions + forwardback[[8]]

         cumulative.n <- cumulative.n + n.i

      }

      Q <<- transitions / rowSums(transitions)

   }




###############################################
#                 MAIN LOOP                   #
###############################################


   for (iter in 1:maxiter) {

      cat(paste("iteration:", iter, "\n"))

      E.step()

      # Modified M-step.

      sumPhi <- colSums(phi, na.rm = TRUE)
      sumPhi.x <- colSums(phi*concat.x, na.rm = TRUE)
      sumPhi.y <- colSums(phi*concat.y, na.rm = TRUE)
      mean.x <- sumPhi.x / sumPhi
      mean.y <- sumPhi.y / sumPhi

      # Find an upper bound.
      no.upper.bound <- TRUE
      new.alpha <- alpha
      lower.bound <- 0
      while (no.upper.bound) {
         zeta <- sumPhi.x / (sumPhi.y + new.alpha*sumPhi)
         ksi <- (1+zeta) * sum((sumPhi.x + sumPhi.y + new.alpha*sumPhi) / (n*(1 + zeta)))
         if (sum(counts*digamma(new.alpha + levels)) / n - digamma(new.alpha) +
            log(new.alpha) - sum(sumPhi*log(ksi)) / n > 0) {
# Version without constraint:
#        if (sum(counts*digamma(new.alpha + levels)) / n - digamma(new.alpha) +
#           log(new.alpha) - sum(sumPhi*log(mean.x + mean.y + new.alpha)) / n
            lower.bound <- new.alpha
            new.alpha <- 2 * new.alpha
         }
         else {
            no.upper.bound <- FALSE
            upper.bound <- new.alpha
         }
      }
      # alpha is estimated with a precision of 0.01.
      while (upper.bound - lower.bound > 0.0001) {
         new.alpha <- (upper.bound + lower.bound) / 2
         zeta <- sumPhi.x / (sumPhi.y + new.alpha*sumPhi)
         ksi <- (1+zeta) * sum((sumPhi.x + sumPhi.y + new.alpha*sumPhi) / (1 + zeta)) / n
         if (sum(counts*digamma(new.alpha + levels)) / n - digamma(new.alpha) +
            log(new.alpha) - sum(sumPhi*log(ksi)) / n > 0)
            lower.bound <- new.alpha
         else
            upper.bound <- new.alpha
      }

      # Note: the rounding prevents the oscillation of the estimates.
      alpha <- round(new.alpha, dig)

      beta.plus.1 <- sum((sumPhi.x + sumPhi.y + alpha*sumPhi) / (1 + zeta)) / (n*alpha)
      beta.gamma <- beta.plus.1 * zeta
      theta.x <- round(beta.gamma / (beta.plus.1 + beta.gamma), dig)
      theta.y <- round((beta.plus.1-1) / (beta.plus.1 + beta.gamma), dig)
      theta.alpha <- 1 - theta.x - theta.y

      cat(paste(c(alpha, theta.x), "\n"))

      if (abs(logLikelihood - oldLogLikelihood) < tol)
         break

      oldLogLikelihood <- logLikelihood

   } # for (iter in 1:maxiter)

   cat("\n")

   if (adjustInitialProb) 
      initialProb <- initial.steady.state.probabilities()
      
   vPath <- Viterbi(Q, initialProb, emissionProb, blockSizes)

   return(list(logL = logLikelihood, Q = Q, alpha = alpha, theta.x =
      theta.x, theta.y = theta.y, vPath = vPath, iterations = iter))
}


