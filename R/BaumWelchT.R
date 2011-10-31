BaumWelchT <- function (x, series.length, m = 2, Q, mu, S, nu = TRUE,
     model, initial.prob, maxiter = 500, overflow = 1e-9, num.inst = 1e-9,
     tol = 1e-05, dig = 3) {
# Author: Guillaume Filion.
# Date: June 17, 2011.
# Return the estimate of the HMM with emissions distributed
# as multi dimensional correlated T variables.
# x: matrix of observations.
# series.length: length of independent blocks in matrix x.
# m: number of states.
# Q: transition matrix.
# mu: means of the emissions in different states.
# S: correlation matrix of the emissions.
# nu: degrees of freedom (in all states).
# model: a model matrix for parameter estimation.
# initial.prob: state probability at the start of series.
# maxiter: maximum number of iterations before returning.
# overflow: the probability ratio in case of numeric overflow.
# num.inst: threshold for no update (see E.step()).
# tol: threshold to difference in log-likelihood before returning.
# dig: numeric precision for parameter estimation.


###############################################
#              OPTION PROCESSING              #
###############################################

    if (is.vector(x)) {
        x <- as.matrix(x);
    }
    p <- ncol(x)
    n <- nrow(x)

    if (missing(series.length)) {
        series.length <- n;
    }
    stopifnot (sum(series.length) == n);

    if (!missing(Q)) {
        if (!all.equal(dim(Q), c(m,m)))
            stop("Q must be an m by m matrix");
    }
    else {
        Q <- matrix(0.1/(m-1), ncol = m, nrow = m);
        diag(Q) <- 0.9;
    }

    # Assign initial values to mu by k-means or quantiles.
    if (missing(mu)) {
        if ((p == 1) && (m == 2)) {
            mu <- matrix(quantile(x = x, probs = c(0.5, 0.9),
                na.rm = TRUE), nrow = m);
        }
        else {
            clusters <- kmeans(x[complete.cases(x),], centers = m);
            mu <- clusters$centers;
        }
    }
    if (!all.equal(dim(mu), c(m, p))) {
        stop ("mu must be a m x p matrix");
    }

    # Assign initial values to covariance matrices.
    if (missing(S)) {
        S <- array(NA_real_, dim = c(p, p, m));
        for (i in 1:m) {
            S[,,i] <- var(sweep(x = x, MARGIN = 2,
                STATS = mu[i,]), na.rm = TRUE);
        }
    }     

    if (!is.infinite(nu)) {
        if (!length(nu) == 1)
            stop ("nu must be a vector of length 1");
        if (is.logical(nu)) {
            if (nu) {
                nu <- 6;
            }
            else {
                nu <- Inf;
            }
        }
    }

        if (missing(model)) {
        model <- array(1:m, dim=c(m,p));
    }
    else {
        stopifnot (all(dim(model) == c(m,p)));
    }

    adjust.initial.prob <- missing(initial.prob);




###############################################
#            VARIABLE DEFINITIONS             #
###############################################


    old.loglik <- -Inf;
    iter <- 0;

    # Updated in E.step (Q as well).
    emission.prob <- matrix(NA_real_, nrow = n, ncol = m)
    phi <- matrix(NA_real_, nrow = n, ncol = m)
    phi.weights <- matrix(NA_real_, nrow = n, ncol = m)
    loglik <- NA_real_


###############################################
#           FUNCTION DEFINITIONS              #
###############################################


    initial.steady.state.probabilities <- function () {
    # Compute the steady-state initial probabilities.

        spectre <- eigen(t(Q))

        if (is.complex(spectre$values[1]))
            return(rep(1/m,m))
        if (spectre$values[1] > 0.99)
            return(spectre$vectors[,1] / sum(spectre$vectors[,1]))

        return(rep(1/m,m))

    }

    # Perform the E-step of the modified Baum-Welch algorithm.
    # Note: the double arrow assignment <<- looks for the 
    # variable on the enclosing scopes. The design prevents
    # passing parameters by value for objects that can be
    # very large, such as 'emission.prob' and 'phi.weights'.
    # The M-step can yield invalid variance estimates because
    # of numeric underflow (see update.S() for example). In that
    # case no update is performed for the given state and the
    # algorithm 'waits' for valid parameter values.
    E.step <- function () {

        weights.i <- matrix(1, nrow = n, ncol = m);

        # Emission probabilities are computed up to constant terms.
        # The following piece of code is stolen from the function
        # mahalanobis.
        for (i in 1:m) {
            # do not update upon numerical instability.
            if (any(is.na(S[,,i])) ||
                   det(as.matrix(S[,,i]) < num.inst)) {
                next;
            }
            x.minus.mu <- sweep(x = x, MARGIN = 2, STATS = mu[i,],
                check.margin = FALSE);
            psi.inv <- solve(S[,,i]);
            sqrt.det.psi <- 1/sqrt(abs(det(as.matrix(S[,,i]))));

            mahalanobis.i <- rowSums((x.minus.mu %*% psi.inv) *
                x.minus.mu);

            if (is.infinite(nu)) {
               emission.prob[,i] <<- sqrt.det.psi*exp(-mahalanobis.i/2);
            }
            else {
               emission.prob[,i] <<- sqrt.det.psi/
                  (1 + mahalanobis.i / nu)^((nu + p) / 2);
            }

            # Use Mahalanobis distances to compute the weights.
            if (!is.infinite(nu)) {
               weights.i[,i] <- (nu + p) / (nu + mahalanobis.i);
            }
        }

        # Emission probabilities of NAs are set to 1 for every states.
        emission.prob[is.na(emission.prob)] <<- 1;

        # Prevent outlier overflow.
        # If 'emission.prob' contains 'Inf', the values are scaled
        # to to the 'overflow' parameter.
        infid <- which(is.infinite(emission.prob), arr.ind = TRUE);
        emission.prob[infid[,1],] <- overflow;
        emission.prob[infid] <- 1;


        if (adjust.initial.prob) {
            initial.prob <- initial.steady.state.probabilities();
        }

        loglik <<- 0;
        cumulative.n <- 0;
        transitions <- matrix(double(m * m), nrow = m);

        counter = 0;

        for (n.i in series.length) {

            counter = counter+1;

            forwardback <- .Fortran("fwdb",
                as.integer(m),
                as.integer(n.i),
                initial.prob,
                emission.prob[(cumulative.n + 1):(cumulative.n + n.i),], Q,
                double(n.i),
                matrix(double(n.i * m), nrow = n.i),
                matrix(double(m^2), nrow = m), double(1),
                PACKAGE = "HMMt");

            loglik <<- loglik + forwardback[[9]];
            phi[(cumulative.n + 1):(cumulative.n + n.i),] <<-
                forwardback[[7]];
            transitions <- transitions + forwardback[[8]];

            cumulative.n <- cumulative.n + n.i;

        }

        Q <<- transitions / rowSums(transitions);
        phi.weights <<- phi*weights.i;

    }

    # Update the means according to the model matrix.
    update.mu <- function() {
        new.mu <- mu * NA
        for (i in 1:p) {
            sumPhi.weights.x <- colSums(phi.weights*x[,i], na.rm = TRUE);
            est <- tapply(X = sumPhi.weights.x, INDEX = model[,i],
                FUN = sum) / tapply(X = sumPhi.weights,
                INDEX = model[,i], FUN = sum);
            new.mu[,i] <- round(est[as.character(model[,i])], dig);
        }
        return (new.mu)
    }

    # Update the variances-covariances.
    update.S <- function() {
        new.S <- S * NA
        for (i in 1:p) {
            for (j in 1:i) {
                sumPhi.weights.xij <- 
                    colSums(phi.weights*x[,i]*x[,j], na.rm = TRUE);
                sumPhi.weightsmuij <- sumPhi.weights*mu[,i]*mu[,j];
                dev <- sumPhi.weights.xij - sumPhi.weightsmuij
                est <- tapply(X = dev, INDEX = model[,i], FUN = sum) / 
                    tapply(X = sumPhi, INDEX = model[,i], FUN = sum);
                new.S[i,j,] <- round(est[as.character(model[,i])], dig)
                if (i != j) {
                    new.S[j,i,] <- new.S[i,j,];
                } 
            }
        }
        # Numerical instability can lead to non-posititivity
        # of the variance-covariance matrix. In that case it is
        # set to NA (caught in E.step()).
        for (i in 1:m) {
            eig <- eigen(new.S[,,i]);
            if (any(eig$values < 0)) {
                new.S[,,i] <- NA
                #                eig$values[eig$values < 0] <- 0;
                #new.S[,,i] <- eig$vectors %*% diag(eig$values) %*%
                #    solve(eig$vectors);
            }
        }
        return (new.S)
    }

    update.nu <- function() {
        weights = rowSums(phi.weights)
        RHS <- 1 + mean((log(weights) - weights), na.rm = TRUE)
        # Find an upper bound.
        no.upper.bound <- TRUE
        new.nu <- nu
        lower.bound <- 0
        while (no.upper.bound) {
            if (digamma(new.nu/2) - digamma((new.nu + p)/2) -
                log(new.nu/2) + log((new.nu + p)/2) < RHS) {
                    lower.bound <- new.nu
                    new.nu <- 2 * new.nu
            }
            else {
                no.upper.bound <- FALSE
                upper.bound <- new.nu
            }
        }
        # The degree of freedom new.nu is estimated with a
        # precision of 0.01
        while (upper.bound - lower.bound > 0.01) {
            new.nu <- (upper.bound + lower.bound) / 2
            if (digamma(new.nu/2) - digamma((new.nu + p)/2) -
                log(new.nu/2) + log((new.nu + p)/2) < RHS)
                    lower.bound <- new.nu
            else {
                upper.bound <- new.nu
            }
        }

        # Note: rounding prevents oscillation of the estimates.
        return (round(new.nu, 1))
    }



###############################################
#                 MAIN LOOP                   #
###############################################

    for (iter in 1:maxiter) {

        cat(paste("\riteration:", iter))

        # Update Q, phi.weights and emission.prob.
        E.step()

        # CM-step 1: updating mu and S.
        sumPhi <- colSums(phi, na.rm = TRUE);
        sumPhi.weights <- colSums(phi.weights, na.rm = TRUE);

        mu <- update.mu();
        S <- update.S();

        if (!is.infinite(nu)) {

            # Update Q, phi.weights and emission.prob.
            E.step()

            # CM-step 2: updating nu.
            nu <- update.nu()
            if (nu > 250) {
                nu <- Inf;
            }

        } # if (!is.infinite(nu))


        if (abs(loglik - old.loglik) < tol)
            break

        old.loglik <- loglik

    } # for (iter in 1:maxiter)

    cat("\n")

    if (adjust.initial.prob) {
        initial.prob <- initial.steady.state.probabilities()
    }
        
    vPath <- Viterbi(Q, initial.prob, emission.prob, series.length)

    # Returns an object of thethe class BaumWelchTfit.
    return(new("fittedHiddenMarkovModel.t", Q = Q, mu = mu, S = S,
        nu = nu, model = model, ViterbiPath = vPath, phi = phi,
        logL = loglik, iterations = iter))

}

