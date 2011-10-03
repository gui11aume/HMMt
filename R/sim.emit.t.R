.sim.emit.t <- function(path, mu, S, nu) {
# Does something

    # The number of states and tracks (dimensions).
    nstates <- nrow(mu)
    ntracks <- ncol(mu)

    stopifnot (nrow(mu) == dim(S)[3])
    stopifnot (ncol(mu) == dim(S)[1])
    stopifnot (is.numeric(nu) && (length(nu) == 1))
    stopifnot (nstates <= nrow(mu))

    emissions <- matrix(, nrow = length(path), ncol = ncol(mu))
    # Compute the square root matrix of the variance matrix.
    sqrtS <- S 
    for (i in 1:nstates) {
        D <- eigen(S[,,i], symmetric = TRUE)
        sqrtS[,,i] <- D$vectors %*% diag(sqrt(D$values),
            ncol = ntracks) %*% solve(D$vectors)
        # Note: the argument ncol = ntracks is necessary in case
        # dimension is 1 (otherwise the function diag returns
        # an empty matrix.
    }

    # We simulate multidimensional t emissions by simulating
    # the 'weights' tau with a chi-square distribution and then
    # simulating a multidimensional Gaussian X with scaled
    # variance matrix.

    # Simulate tau.
    tau <- rgamma(length(path), shape = nu/2, rate = nu/2)

    # Simulate X.
    X <- matrix(rnorm(length(path)*ncol(mu)), ncol=ncol(mu))

    # Fill Y state-wise.
    for (state in 1:nstates) {
        X[path == state,] <- X[path == state,, drop = FALSE] %*%
            sqrtS[,,state] / sqrt(tau[path == state])
        X[path == state,] <- sweep(X[path == state,, drop = FALSE],
            MARGIN = 2, STATS = mu[state,], FUN = "+")
    }

    return(X)

}


