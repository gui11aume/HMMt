setClass ("HiddenMarkovModel",
    representation = representation(
        Q = "matrix"))

setClass ("GeneratedHiddenMarkovModel",
    representation = representation(
        path = "numeric",
        emissions = "matrix",
        first = "character"),
    contains = "HiddenMarkovModel")

setClass ("HiddenMarkovModel.t",
    representation = representation(
        mu = "matrix",
        S = "array",
        nu = "numeric",
        model = "matrix"),
    contains = "HiddenMarkovModel")

setClass ("fittedHiddenMarkovModel.t",
    representation = representation(
        ViterbiPath = "numeric",
        logL = "numeric",
        iterations = "numeric",
        method = "character"),
    contains = "HiddenMarkovModel.t",
    prototype(method = "BaumWelchT"))

setClass ("GeneratedHiddenMarkovModel.t",
    contains = c(
        "GeneratedHiddenMarkovModel",
        "HiddenMarkovModel.t"))

# Overloads $ for the class HiddenMarkovModel
getGeneric ("$")
setMethod ("$",
    signature(x = "HiddenMarkovModel"),
    definition = function(x, name) {
        return(slot(x, name))
    })

setMethod (f = "simulate",
    signature (
        object = "HiddenMarkovModel.t",
        nsim = "missing",
        seed = "missing"),
    definition = function(object, 
            length. = length(object@ViterbiPath), ...) {
        if (length. == 0) {
            stop ("Argument length is missing or 0")
        }
        path = rMarkov(length. = length., Q = object@Q, ...)
        new("GeneratedHiddenMarkovModel.t",
            path = path,
            mu = object@mu,
            S = object@S,
            nu = object@nu,
            first = ifelse(is.null(list(...)[["first"]]),
                "ergodic", list(...)[["first"]]),
            emissions = .sim.emit.t(path = path,
                mu = object@mu,
                S = object@S,
                nu = object@nu)
        )
    })
