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
#' ## Example with normally distributed covariate:
#' m1 <- countEffects(y="dv", x="treat", z="pre", data=example01,
#'                    method="poisson", distribution="condNormal")
#' summary(m1)
#'
#' @export

# Result of example: -4.10494, SE=0.2972011

countEffects <- function(y,
                         x,
                         k = NULL,
                         z,
                         data,
                         method = "nb",
                         distribution = "condNormal",
                         fixed.cell = "default",
                         fixed.z = "default",
                         homoscedasticity = "default",
                         control = "default",
                         measurement = character(),
                         ...){
  #### Maybe some checking later on (stopifnot etc.)
  if (!is.count(data[,y])){
    stop("The dependent variable needs to consist of non-negative integers only.")
  }

  if (!is.null(k) & distribution != "condNormal"){
    stop("Currently, categorical covariates are not allowed for.")
  }

  if (length(z) > 1 & distribution != "condNormal"){
    stop("Currently, just a single covariate is allowed for.",
         call. = FALSE)
  }

  if (length(measurement) != 0){
    stop("Measurement models are currently under development.")
  }

  ####
  object <- new("countEffects")
  object@input <- createInput(y, x, k, z, data, measurement, method, distribution, control,
                              fixed.cell, fixed.z, homoscedasticity)
  object@parnames <- createParNames(object)

  object@syntax <- createSyntax(object)

  object@results <- computeResults(object)
  #object@results@glmresults <- estimateGLM(object)

  object@ate <- computeATE(x = x,
                           z = z,
                           mod = object@results@glmresults,
                           data = data,
                           distribution = distribution)


  return(object)
}
