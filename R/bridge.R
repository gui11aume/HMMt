bridge <- function (scores, is.sorted = FALSE) {
# Author: Guillaume Filion.
# Date: June 17, 2011.
# Return an 'bridged' series suitable for HMM.
# scores: [chromosome, start, end, x1, x2, x3, ...]

    if (missing(scores))
        stop ("Argument scores is missing")

    D <- ncol(scores)-3
    n <- nrow(scores)

    # Sort the data.
    if (!is.sorted)
        scores <- scores[order(scores[,1], scores[,2]),]

    # Compute the mean distance between consecutive start sites.
    distances <- diff(scores[-1,2])
    mean.dist <- mean(distances[distances > 0], na.rm = TRUE)

    # Segments are delimited by chromosome ends.
    breaks <- scores[-1,1] != scores[-n,1]

    # Virtual bridging: gaps are filled with NAs.
    gaps <- scores[-1,2] - scores[-n,3]
    bridges <- which(gaps > 2*mean.dist)
    virtuals <- as.integer(floor(gaps[bridges] / mean.dist) - 1)

    # Make an index vector by 'zipping' (courtesy W. Meuleman).
    nonvirtuals <- rep(rep(c(TRUE, FALSE), length(bridges)+1),
        times = rbind(diff(c(0, bridges, n)), c(virtuals, 0)))

    # Update the length of the bridged series.
    n <- length(nonvirtuals)

    # Build the virtually bridged table.
    bridged.scores <- matrix(NA_real_, nrow = n, ncol = D)
    bridged.scores[nonvirtuals,] <- as.matrix(scores[,4:ncol(scores)])

    # Update the break positions.
    updated.breaks <- rep(FALSE, n)
    updated.breaks[nonvirtuals] <- c(breaks, FALSE)
        breaks <- which(updated.breaks)

    return(list(x = bridged.scores, series.length = diff(c(0, breaks, n)),
        nonvirtuals = nonvirtuals))

}
