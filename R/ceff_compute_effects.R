ceff_compute_effects <- function(object){
  cat("Computing effects...")

  # Extract parameters estimates and vcov from countreg
  # without measurement model (irrelevant for effects)
  # Get CountReg partable of model
  creg_obj <- object@fit
  n_par <- max(creg_obj@fit$pt$par_free)
  creg_pt <- subset(creg_obj@fit$pt, dest != "mm" & dest != "sigmaw")
  par <- creg_pt$par
  est <- ceff_par_to_effects(par, object)

  # Get variance-covariance matrix of model
  creg_vcov_raw <- object@fit@fit$vcov_fit
  free.idx <- creg_pt$par_free
  creg_vcov <- creg_vcov_raw[free.idx, free.idx]
  JAC <- ceff_func_jacobian_simple(func=ceff_par_to_effects, x=par, object=object)
  vcov_def <- JAC %*% creg_vcov %*% t(JAC)

  se <- c(sqrt(diag(vcov_def)))
  names(se) <- names(est)
  colnames(vcov_def) <- rownames(vcov_def) <- names(est)

  tval <- est/se
  rdf <- nrow(object@input@data) - n_par


  pval <- 2*(1-pt(abs(tval),df=rdf))

  ng <- object@input@ng
  nz <- object@input@nz
  nk <- object@input@nk

  sdyx0 <- function(object){
    # TODO: change to control group value
    y <- object@input@data[, object@input@vnames$y]
    x <- object@input@data[, object@input@vnames$x]
    return(sd(y[x==0]))
  }

  sdyx0 <- sdyx0(object)


  selector <- attr(est, "type") == "Eg"
  Egx <- data.frame(est[selector],
                    se[selector],
                    tval[selector],
                    pval[selector],
                    est[selector]/sdyx0)

  names(Egx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given a treatment condition
  selector <- attr(est, "type") == "Egxgx"
  Egxgx <- data.frame(est[selector],
                      se[selector],
                      tval[selector],
                      pval[selector],
                      est[selector]/sdyx0)
  names(Egxgx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")


  cat("done\n")
  res <- new("results",
      est=est,
      se=se,
      vcov_def=vcov_def,
      Egx=Egx,
      Egxgx=Egxgx
      )
}


ceff_par_to_effects <- function(par, object){
  # join par and pt (necessary for Jacobian)
  creg_obj <- object@fit
  creg_pt <- subset(creg_obj@fit$pt, dest != "mm" & dest != "sigmaw")
  creg_pt$par <- par

  # import information from object
  input <- object@input
  ngroups <- input@ngroups
  nz <- input@nz
  nk <- input@nk
  ng <- input@ng
  vnames <- input@vnames

  celleffgrid <- expand.grid(g = 1:(ng-1), k = 0:(nk-1), x = 0:(ng-1))
  cellgrid <- expand.grid(k = 0:(nk-1), x = 0:(ng-1))

  creg_pt$group.x <- cellgrid$x[creg_pt$group]
  creg_pt$group.k <- cellgrid$k[creg_pt$group]

  groupw <- subset(creg_pt, dest == "groupw")
  regcoef <- subset(creg_pt, dest == "regcoef")
  meanz <- subset(creg_pt, !is.na(type) & type == "mean")
  varz <- subset(creg_pt, !is.na(type) & type == "var")
  covz <- subset(creg_pt, !is.na(type) & (type == "cov" | type == "cov_z_lv"))

  # Estimates and type attributes
  est <- numeric()
  type <- numeric()

  # Relative frequencies of XK-groups
  relfreq <- exp(groupw$par)/sum(exp(groupw$par))
  names(relfreq) <- paste0("relfreq_gx",cellgrid$x,"k",cellgrid$k)
  est <- c(est, relfreq)
  type <- c(type, rep("relfreq", length(relfreq)))

  # Unconditional probabilities P(K=k) and P(X=x)
  relfreq_array <- array(relfreq, dim=c(nk, ng))
  Pk <- apply(relfreq_array, 1, sum)
  names(Pk) <- paste0("Pk",0:(nk-1))
  Px <- apply(relfreq_array, 2, sum)
  names(Px) <- paste0("Px",0:(ng-1))
  est <- c(est, Pk, Px)
  type <- c(type, rep("Pk", length(Pk)), rep("Px", length(Px)))

  # Conditional probabilities P(X=x|K=k)
  Pxgk <- sapply(1:nk, function(k){relfreq_array[k, ]/Pk[k]})
  Pkgx <- sapply(1:ng, function(x){relfreq_array[ ,x]/Px[x]})
  Pxgk <- array(Pxgk, dim=c(nk, ng))
  Pkgx <- array(Pkgx, dim=c(nk, ng))

  tmp <- expand.grid(0:(nk-1), 0:(ng-1))
  Pxgk_vec <- as.numeric(Pxgk)
  Pkgx_vec <- as.numeric(Pkgx)
  names(Pxgk_vec) <- paste0("Px",tmp[,2],"gk",tmp[,1])
  names(Pkgx_vec) <- paste0("Pk",tmp[,1],"gx",tmp[,2])
  est <- c(est, Pxgk_vec, Pkgx_vec)
  type <- c(type, rep("Pxgk", length(Pxgk_vec)), rep("Pkgx", length(Pkgx_vec)))

  celleffs <- numeric()
  for (i in 1:nrow(celleffgrid)){
   # Combination of g-function and cell parameters
   cell <- celleffgrid[i,]

   # Extract relevant parameters for effect computation in this cell
   regcoef1 <- subset(regcoef, group.x == cell$g & group.k == cell$k)
   regcoef0 <- subset(regcoef, group.x == 0 & group.k == cell$k)
   meanz_xk <- subset(meanz, group.x == cell$x & group.k == cell$k)
   varz_xk <- subset(varz, group.x == cell$x & group.k == cell$k)
   covz_xk <- subset(covz, group.x == cell$x & group.k == cell$k)

   # Sort parameter tables for convenient effect computation
   regcoef0 <- regcoef0[match(c(1, vnames$z), regcoef0$rhs),]
   regcoef1 <- regcoef1[match(c(1, vnames$z), regcoef1$rhs),]
   meanz_xk <- meanz_xk[match(vnames$z, meanz_xk$lhs), ]
   varz_xk <- varz_xk[match(vnames$z, varz_xk$lhs), ]

   REGCOEF1 <- matrix(regcoef1$par, nrow = 1)
   REGCOEF0 <- matrix(regcoef0$par, nrow = 1)
   MEANZ <- matrix(c(1, meanz_xk$par), ncol = 1)

   if (nz == 1){
     SIGMAZ <- varz_xk$par
   } else if (nz >= 2){
     covz_combn <- t(combn(vnames$z, 2))
     covz_order <- integer()
     for (j in 1:nrow(covz_combn)){
       cov_lhs <- covz_combn[j,1]
       cov_rhs <- covz_combn[j,2]
       cov_pos <- which((grepl(cov_lhs, covz_xk$lhs) | grepl(cov_lhs, covz_xk$rhs)) & (grepl(cov_rhs, covz_xk$lhs) | grepl(cov_rhs, covz_xk$rhs)))
       covz_order <- c(covz_order, cov_pos)
     }
     covz_xk <- covz_xk[covz_order, ]

     SIGMAZ <- lavacreg:::creg_matrix_vech_reverse(covz_xk$par, diagonal = FALSE)
     diag(SIGMAZ) <- varz_xk$par
   } else {
     SIGMAZ <- 1
   }

   eff1 <- exp(REGCOEF1 %*% MEANZ + 0.5*(t(REGCOEF1[ ,-1]) %*% SIGMAZ %*% REGCOEF1[ ,-1]))
   eff0 <- exp(REGCOEF0 %*% MEANZ + 0.5*(t(REGCOEF0[ ,-1]) %*% SIGMAZ %*% REGCOEF0[ ,-1]))
   eff <- eff1 - eff0
   names(eff) <- paste0("eff_g",cell$g,"_gx",cell$x,"k",cell$k)
   celleffs <- c(celleffs, eff)
  }
  est <- c(est, celleffs)
  type <- c(type, rep("celleff", length(celleffs)))

  # Average Effects
  celleffs_array <- array(celleffs, dim=c(ng-1,nk,ng))
  Eg <- numeric(ng-1)
  for (i in 1:(ng-1)){
    Eg[i] <- sum(celleffs_array[i, ,]*relfreq)
  }
  names(Eg) <- paste0("Eg",1:(ng-1))
  est <- c(est, Eg)
  type <- c(type, rep("Eg", ng-1))

  # Effects given X=x
  Egxgx <- array(NA, dim = c(ng-1, ng))
  for (i in 1:(ng-1)){
    Egxgx[i,] <- apply(celleffs_array[i, ,] * Pkgx, 2, sum)
  }
  Egxgx <- as.numeric(Egxgx)
  tmp <- expand.grid(1:(ng-1),0:(ng-1))
  names(Egxgx) <- paste0("Eg",tmp[,1],"gx",tmp[,2])
  est <- c(est, Egxgx)
  type <- c(type, rep("Egxgx", length(Egxgx)))

  attr(est, which="type") <- type
  return(est)
}


# Taken from lavaan
ceff_func_jacobian_complex <- function (func, x, h = .Machine$double.eps, ..., fallback.simple = TRUE)
{
  f0 <- try(func(x * (0 + (0+1i)), ...), silent = TRUE)
  if (inherits(f0, "try-error")) {
    if (fallback.simple) {
      dx <- ceff_func_jacobian_simple(func = func, x = x,
                                     h = sqrt(h), ...)
      return(dx)
    }
    else {
      stop("function does not support non-numeric (complex) argument")
    }
  }
  nres <- length(f0)
  nvar <- length(x)
  h <- pmax(h, abs(h * x))
  tmp <- x + h
  h <- (tmp - x)
  dx <- matrix(as.numeric(NA), nres, nvar)
  for (p in seq_len(nvar)) {
    dx[, p] <- Im(func(x + h * (0+1i) * (seq.int(nvar) ==
                                           p), ...))/h[p]
  }
  dx
}

ceff_func_jacobian_simple <- function (func, x, h = sqrt(.Machine$double.eps), ...)
{
  f0 <- func(x, ...)
  nres <- length(f0)
  nvar <- length(x)
  h <- pmax(h, abs(h * x))
  tmp <- x + h
  h <- (tmp - x)
  dx <- matrix(as.numeric(NA), nres, nvar)
  for (p in seq_len(nvar)) {
    dx[, p] <- (func(x + h * (seq.int(nvar) == p), ...) -
                  func(x, ...))/h[p]
  }
  dx
}
