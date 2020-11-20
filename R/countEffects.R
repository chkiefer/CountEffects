#' Average and conditional treatment effects for count outcomes
#'
#' A general function to estimate average and conditional treatment effects for a count outcome.
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
#'                    method="poisson", distribution="normal")
#' summary(m1)
#'
#' @export

# Result of example: -4.10494, SE=0.2972011

countEffects <- function(y,
                         x,
                         k = NULL,
                         z = NULL,
                         data,
                         method = "nb",
                         distribution = "condNormal",
                         control = "default",
                         measurement = list(),
                         na.rm = TRUE){
  #### Maybe some checking later on (stopifnot etc.)


  if (na.rm){
    z0 <- z[z %in% names(data)]
    z1 <- unlist(measurement)
    completeVec <- complete.cases(data[c(y,x,z0,z1,k)])
    data <- data[completeVec,]
  } else {
    stop("CountEffects error: missing data handling not available")
  }

  ####
  object <- new("countEffects")
  object@input <- ceff_create_input(y, x, k, z, data, measurement, method, distribution, control)

  object@fit <- ceff_estimate_creg(object)

  object@results <- ceff_compute_effects(object)


  #object@parnames <- ceff_create_parnames(object)

  #object@syntax <- ceff_create_syntax(object)



  cat("Finished.\n")
  return(object)
}
