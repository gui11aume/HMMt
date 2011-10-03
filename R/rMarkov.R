rMarkov <- function (length., Q, first = "ergodic") {
# Simulates a Markov chain of length n with transitions Q.
# The first state can be user-specified, taken from a probability
# vector, or 'ergodic', ie assuming that the chain was running
# for inifinite time before sampling the first state.
# Q must be a transition matrix, rows summing up to 1.

# Note: the dot in the paramter name 'length.' is important to
# avoid confusion with the built-in function length(). Other names
# such as 'n' turned out to create collisions within the method
# 'simulate()'.
    
    if (length(dim(Q)) != 2) {
        stop ("Q must be a matrix.")
    }

    if (nrow(Q) != ncol(Q)) {
        stop ("Q must be an mxm matrix.")
    }
    else {
        m <- nrow(Q)
    }

    if (!all.equal(rowSums(Q), rep(1, m), tolerance = 1e-6) || any(Q < 0)) {
        stop ("Q must be a transition matrix (rows must sum up to 1).")
    }

    if ((length(first) == m) && (sum(first) == 1)) {
        first = sample(1:m, size = 1, prob = first)
    }
    if (!is.integer(first) && !(first %in% 1:m) && (first != "ergodic")) {
        stopmessage <- paste(
            "The argument 'first', must be 'ergodic', a vector of
            probabilities or an initial state.")
        stop (stopmessage)
    }

    # Compute the cumulative probabilities matrix.
    sumQ <- t(apply(Q, MARGIN = 1, function(x) diffinv(x)[-1]))

    if (first == "ergodic") {
        warningmessage <- paste(
            "Matrix Q does not allow proper initiation of the sequence.",
            "Using uniform sampling instead.", sep = "\n")
        spectre <- eigen(t(Q))
        eig <- spectre$values[1]
        if (is.complex(eig) || !all.equal(eig, 1, tol=.01)) {
            warning(warningmessage)
            weights <- rep(1/m, m)
        }
        else {
            weights <- spectre$vectors[,1] / sum(spectre$vectors[,1])
        }
        first <- sample(1:m, size = 1, prob = weights)
    }
        
    
    x <- c(first, rep(0, length.-1))
    simulated <- .Fortran("simark",
        as.integer(length.),
        as.integer(m),
        array(as.double(sumQ), dim = dim(Q)),
        as.double(runif(length., min = 0, max = 1)),
        as.integer(first),
        as.integer(x),
        PACKAGE = "DmTools")

    return(simulated[[6]])
}
