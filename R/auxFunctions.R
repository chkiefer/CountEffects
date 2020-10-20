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
  function(x, tol = .Machine$double.eps^0.5, na.rm = TRUE){
    if (na.rm) x <- na.omit(x)
    tmp0 <- abs(x - round(x)) < tol
    tmp1 <- sign(x) == -1

    if (sum(tmp1) > 0){
      FALSE
      } else if (sum(tmp0)/length(x) == 1){
      TRUE
    } else  {FALSE}
}


getFormula <- function(object){
  inp <- object@input

  nz <- inp@nz
  fml <- paste0(inp@vnames$y,
                "~",
                inp@vnames$x,
                "*",
                inp@vnames$z[1])
  if (nz > 1){
    for (i in 2:nz) {
      fml <- paste0(fml,
                    "+",
                    inp@vnames$x,
                    "*",
                    inp@vnames$z[i])
    }
  }
  fml <- as.formula(fml)
  return(fml)
}

estimateGLM <- function(object){
  fml <- object@syntax@formula
  data <- object@input@data
  method <- object@input@method

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


computeAlphas <- function(object, glmresults){
  nz <- object@input@nz
  coefs <- coef(glmresults)
  vcovs <- vcov(glmresults)
  pnames <- getCoefNames(nz)

  names(coefs) <- row.names(vcovs) <- colnames(vcovs) <- pnames

  tmp <- paste0("\n y ~ ", pnames[1],"*1")
  i <- 0
  for (name in pnames[-1]){
    tmp <- paste0(tmp, "+", name, "*x",i)
    i <- i +1
  }
  tmp <- paste0(tmp, "\n #Compute Alphas from Gammas\n")
  tmp <- paste0(tmp, "a000 := g000\n")
  tmp <- paste0(tmp, "a100 := g100 + g000\n")
  for (i in 1:nz){
    tmp <- paste0(tmp, "a00",i," := g00",i,"\n")
  }
  for (i in 1:nz){
    tmp <- paste0(tmp, "a10",i," := g10",i," + g00",i, "\n")
  }

  pt <- lavaanify(tmp)
  def.function <- lavaan::lav_partable_constraints_def(pt)
  JAC <- lav_func_jacobian_complex(func=def.function, x=coefs)
  info.r <- JAC %*% vcovs %*% t(JAC)
  est <- def.function(coefs)

  row.names(info.r) <- colnames(info.r) <- names(est)

  res <- list(coefs = est,
              vcovs = info.r)
  return(res)
}

