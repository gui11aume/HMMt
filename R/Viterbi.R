Viterbi <- function (Q, initialProb, emissionProb, series.length) {
# Author: G. Filion.
# Date: December 17, 2009.
# Return the Viterbi path of the given series.

   n <- nrow(emissionProb)
   m <- as.integer(nrow(Q))
   ViterbiPath <- integer(n)
   if(missing(series.length))
      series.length <- n

   emissionProb <- log(emissionProb)
   Q <- log(Q)

   # Prevent underflow.
   emissionProb[is.infinite(emissionProb)] <- -325
   Q[is.infinite(Q)] <- -325

   cumulative.n <- 0

   for (n.i in as.integer(series.length)) {

      viterbi <- .Fortran("vit", n.i, m, log(initialProb),
         emissionProb[(cumulative.n + 1):(cumulative.n + n.i),],
         Q, integer(n.i), matrix(double(n.i * m), nrow = n.i),
         PACKAGE = "HMMt")

      ViterbiPath[(cumulative.n + 1):(cumulative.n + n.i)] <- viterbi[[6]]
      cumulative.n <- cumulative.n + n.i

   }

   return(ViterbiPath)

}
