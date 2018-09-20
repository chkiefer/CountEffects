getData <- function(N){
  x <- rbinom(N, 1, .5)

  muZgx0 <- c(-1,-.5)
  muZgx1 <- c(1,.5)
  Sigmagx0 <- matrix(c(1,.5,.5,1), nrow=2)
  Sigmagx1 <- matrix(c(1,.25,.25,1), nrow=2)
  z <- MASS::mvrnorm(N, muZgx0, Sigmagx0)*(1-x) + mvrnorm(N, muZgx1, Sigmagx1)*(x)
  z1 <- z[,1]
  z2 <- z[,2]

  muY <- exp(0.1 + 0.2*x + 0.1*z1 - 0.3*x*z2)
  y <- rpois(N, muY)
  data <- data.frame(y,x,z1,z2)
}

getCoefNames <- function(nz){
  pnames <- c("g000", "g100")
  for (i in 1:nz){
    tmp <- paste0("g00",i)
    pnames <- c(pnames, tmp)
  }
  for (i in 1:nz){
    tmp <- paste0("g10",i)
    pnames <- c(pnames, tmp)
  }
  return(pnames)
}

is.count <-
  function(x, tol = .Machine$double.eps^0.5){
    x <- na.omit(x)
    tmp0 <- abs(x - round(x)) < tol
    tmp1 <- sign(x) == -1

    if (sum(tmp1) > 0){
      FALSE
      } else if (sum(tmp0)/length(x) == 1){
      TRUE
    } else  {FALSE}
}


getFormula <- function(object){
  nz <- length(object@z)
  fml <- paste0(object@y,
                "~",
                object@x,
                "*",
                object@z[1])
  if (nz > 1){
    for (i in 2:nz) {
      fml <- paste0(fml,
                    "+",
                    object@x,
                    "*",
                    object@z[i])
    }
  }
  fml <- as.formula(fml)
  return(fml)
}

estimateGLM <- function(fml, data, method){
  if (method == "nb"){
    mod <- MASS::glm.nb(fml, data)
    return(mod)
  } else if (method == "poisson"){
    mod <- glm(fml, data, family = poisson())
    return(mod)
  } else if (method == "gamma"){
    mod <- glm(fml, data, family = Gamma(link = "log"))
    return(mod)
  }
}


testDistribution <- function(z, distribution){
  if (distribution == "normal"){
    return(shapiro.test(z))
  } else if (distribution == "negbin"){
    N <- length(z)
    maxZ <- max(z)

    tmp <- tabulate(z)
    observed <- c(N - sum(tmp),
                  tmp,
                  rep(0, 1))
    names(observed) <- sapply(c(0:(maxZ + 1)), paste)

    muZ <- var(z)
    varZ <- var(z)

    p.negbin <- c(dnbinom(c(0:maxZ + 0), mu = muZ, size = muZ^2/(varZ - muZ)), 1-pnbinom((maxZ + 0), mu = muZ, size = muZ^2/(varZ - muZ)))
    chiTest <- chisq.test(observed, p=p.negbin, rescale.p = TRUE)

    return(chiTest)
  } else if (distribution == "poisson"){
    N <- length(z)
    maxZ <- max(z)

    tmp <- tabulate(z)
    observed <- c(N - sum(tmp),
                  tmp,
                  rep(0, 1))
    names(observed) <- sapply(c(0:(maxZ + 1)), paste)


    muZ <- mean(z)
    p.pois <- c(dpois(c(0:(maxZ + 0)), muZ), 1-ppois((maxZ + 0), muZ))
    chiTest <- chisq.test(observed, p=p.pois)

    return(chiTest)
  }
}
