getData <- function(N){
  z <- rnorm(N)
  x <- rbinom(N, 1, plogis(z))
  muY <- exp(0.1 + 0.2*x + 0.1*z - 0.2*x*z)
  y <- rpois(N, muY)
  data <- data.frame(y,x,z)
}


estimateGLM <- function(fml, data, method){
  if (method == "nb"){
    require(MASS)
    mod <- glm.nb(fml, data)
    return(mod)
  } else if (method == "poisson"){
    mod <- glm(fml, data, family = poisson())
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
