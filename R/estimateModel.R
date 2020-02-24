estimateModel <- function(object){
  objsem <- new("nbmgsem")
  objsem@ngroups <- length(object@input@vlevels$cell)
  objsem@nz <- object@input@nz
  objsem@vnames <- object@input@vnames
  objsem@method <- object@input@method


  tmp2 <- objsem@ngroups * (objsem@nz + 1)  # regression coefficients
  #tmp1 <- tmp0 + (2*objsem@nz + objsem@nz*(objsem@nz - 1)/2)*objsem@ngroups   # means, variances, covariances within one group
  #tmp2 <- tmp1 + objsem@ngroups

  objsem@npar <- as.integer(tmp2)
  #rm(tmp0); rm(tmp1); rm(tmp2)

  objsem@data <- object@input@data
  par <- numeric(objsem@npar)

  ### optimization
  objsem@nlminb <- nlminb(start = par,
                          objective = loglik,
                          objsem = objsem)
  info.m <- pracma::hessian(f = loglik, x = objsem@nlminb$par, objsem = objsem)
  objsem@varcov <- solve(info.m)/nrow(objsem@data)
  return(objsem)
}


loglik <- function(par, objsem){
    nz <- objsem@nz
    method <- objsem@method

    par_matrix <- matrix(par, nrow = objsem@ngroups, ncol = (nz + 1), byrow = T)
    df <- objsem@data
    vnames.z <- objsem@vnames$z
    vnames.y <- objsem@vnames$y

    lin_predict <- (as.matrix(cbind(1, df[,vnames.z])) * par_matrix[df$cell,]) %*% matrix(rep(1, (nz + 1)), ncol = 1)
    lambdaY <- exp(lin_predict)

    if (method == "poisson"){
      LL <- sum(dpois(x = df[,vnames.y], lambda = lambdaY, log = T))
    } else if (method == "nb"){
      LL <- sum(dnbinom(x = df[,vnames.y], mu = lambdaY, size = 1, log = T))
    }


    LL <- - LL/nrow(df)
    LL
}
