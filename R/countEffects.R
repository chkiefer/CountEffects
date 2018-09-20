#' Average treatment effects for count outcomes
#'
#' A general function to estimate average treatment effects for a count outcome.
#'
#' @param y Dependent variable (character string). Has to be the name of a count variable (i.e., non-negative integer).
#' @param x Treatment variable (character string) currently treated as binary variable.
#' @param z Covariate variable (character string).
#'
#' @param data A data frame.
#' @param method By default a negative binomial regression model is fitted (i.e., \code{"nb"}). Alternatively, with \code{"poisson"}
#' a standard Poisson regression is fitted.
#'
#' @param distribution Distribution of the covariate. By default, normal distribution is assumed (i.e., \code{"normal"}).
#' Currently, a uniform distribution (i.e., \code{"uniform"}), a Poisson distribution (i.e., \code{"poisson"}),
#' a chisquare distribution (i.e., \code{"chisquare"}), and a negative binomial distribution (i.e., \code{"negbin"}) are possible.
#'
#' @return Object of class countEffects.
#'
#' @examples
#' ## Example with normally distribtuion covariate:
#' m1 <- countEffects(y="dv", x="treat", z="pre", data=example01,
#'                    method="poisson", distribution="normal")
#' summary(m1)
#'
#' @export

countEffects <- function(y,
                         x,
                         z,
                         data,
                         method = "nb",
                         distribution = "normal",
                         ...){
  #### Maybe some checking later on (stopifnot etc.)
  if (!is.count(data[,y])){
    stop("The dependent variable needs to consist of non-negative integers only.")
  }

  if (length(z) > 1 & distribution != "condNormal"){
    stop("Currently, just a single covariate is allowed for.")
  }


  ####
  object <- new_cATE()
  object@y <- y
  object@x <- x
  object@z <- z

  fml <- getFormula(object)

  object@model <- estimateGLM(fml, data, method)
  object@distTest <- testDistribution(data[,object@z], distribution)

  object@ate <- computeATE(x = x,
                           z = z,
                           mod = object@model,
                           data = data,
                           distribution = distribution)
  return(object)
}
