#' Estimation of average treatment effect
#'
#' Get point and standard error estimate for the average treatment effect given a fitted regression model.
#'
#' @param x Treatment variable (character string) currently treated as binary variable.
#' @param z Covariate variable (character string).
#'
#' @param data A data frame.
#' @param mod A \code{"glm"} object with logarithmic link function.
#' a standard Poisson regression is fitted.
#'
#' @param distribution Distribution of the covariate.
#' Currently, normal distribution (i.e., \code{"normal"}),
#' uniform distribution (i.e., \code{"uniform"}), Poisson distribution (i.e., \code{"poisson"}),
#' chisquare distribution (i.e., \code{"chisquare"}), and negative binomial distribution (i.e., \code{"negbin"}) are possible.
#'
#' @return List of average treatment effect and standard error.
#'
#' @examples
#' ## Currently, no examples.
#'
#'
#' @export

computeATE <- function(x, z, mod, data, distribution){
  coefs <- coef(mod)
  vcovs <- vcov(mod)

  pnames <- c("g000", "g100", "g001", "g101")
  names(coefs) <- row.names(vcovs) <- colnames(vcovs) <- pnames

  vnames <- names(data)
  vnames[vnames == z] <- "z"
  names(data) <- vnames

  modelz <- '
  z ~ c(meanz001,meanz101)*1
  z ~~ c(varz001,varz101)*z
  group % c(groupw0,groupw1)*w
  gw0 := groupw0
  gw1 := groupw1
  mz001 := meanz001
  mz101 := meanz101
  vz001 := varz001
  vz101 := varz101
  N := exp(gw0) + exp(gw1)
  relfreq0 := exp(gw0)/N
  relfreq1 := exp(gw1)/N
  Ez1 := mz001*relfreq0 + mz101*relfreq1
  Vz1 := vz001*relfreq0 + vz101*relfreq1 + relfreq0*(mz001-Ez1)^2 + relfreq1*(mz101-Ez1)^2
  Px0 := relfreq0
  Px1 := relfreq1
  Pk0gx0 := relfreq0/Px0
  Pk0gx1 := relfreq1/Px1
  Ez1gx0 := mz001*Pk0gx0
  Ez1gx1 := mz101*Pk0gx1
  Vz1gx0 := vz001*Pk0gx0
  Vz1gx1 := vz101*Pk0gx1
  '
  require(lavaan)
  mz <- sem(modelz, data=data, group=x, group.label=c(0,1), group.w.free=TRUE)

  ## augment coefs and vcovs
  acoefs <- c(coefs, coef(mz, type="user")[-c(1:6)])
  avcovs <- lav_matrix_bdiag(vcovs, lavInspect(mz, "vcov.def", add.class = FALSE))
  row.names(avcovs) <- colnames(avcovs) <- names(acoefs)

  if (distribution == "normal"){
    Eg1 <- deltaMethod(acoefs, "exp(g000+g100)*exp((g001+g101)*Ez1+((g001+g101)^2*Vz1/2))-exp(g000)*exp(g001*Ez1+(g001^2*Vz1/2))", avcovs, func="Eg1")
  } else if (distribution == "uniform"){
    Eg1 <- deltaMethod(acoefs,
                       "exp(g000+g100)*(exp((g001+g101)*(Ez1+sqrt(3*Vz1)))-exp((g001+g101)*(Ez1-sqrt(3*Vz1))))/((g001+g101)*(2*sqrt(3*Vz1)))-exp(g000)*(exp(g001*(Ez1+sqrt(3*Vz1)))-exp(g001*(Ez1-sqrt(3*Vz1))))/(g001*(2*sqrt(3*Vz1)))",
                       avcovs, func="Eg1")
  } else if (distribution == "poisson"){
    Eg1 <- deltaMethod(acoefs, "exp(g000+g100)*exp(Ez1*(exp(g001+g101)-1))-exp(g000)*exp(Ez1*(exp(g001)-1))", avcovs, func="Eg1")
  } else if (distribution == "chisquare"){
    Eg1 <- deltaMethod(acoefs, "(exp(g000+g100)/(1-2*(g001+g101))^(Ez1/2))-(exp(g000)/(1-2*g001)^(Ez1/2))", avcovs, func="Eg1")
  } else if (distribution == "negbin"){
    Eg1 <- deltaMethod(acoefs, "exp(g000+g100)*(Ez1/(Ez1+Vz1)*exp(g001+g101)/(1-(1-Ez1/(Ez1+Vz1))*exp(g001+g101)))^(Ez1^2/(Ez1+Vz1))-exp(g000)*(Ez1/(Ez1+Vz1)*exp(g001)/(1-(1-Ez1/(Ez1+Vz1))*exp(g001)))^(Ez1^2/(Ez1+Vz1))", avcovs, func="Eg1")
  }


  res <- list(Ave=Eg1$Estimate,
              se_Ave=Eg1$SE)

  return(res)

}
