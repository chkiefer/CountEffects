ceff_case_independence <- function(data = dNA, intPoints = 21L, silent = FALSE){
  msg <- "This function is implemented for illustration only. It accompanies a paper by Kiefer & Mayer in the British Journal of Mathematical and Statistical Psychology. Its functionality is meant to be implemented within the main countEffects()-function in the future.\n Please stay tuned for further updates.\n"
  cat(msg)

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package \"dplyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  require(dplyr)

  ######################################
  # readData
  #####################################
  # data <- readRDS("dACTIVE.rds")
  dNA <- na.omit(data)
  N <- nrow(dNA)

  object <- new("countEffects")
  object@input <- ceff_create_input(y="dv", x="treat", k="gender", z=c("depression", "pretest"), dNA,
                                    measurement=list(depression=c("ind1", "ind2", "ind3")), method="poisson", distribution="condNormal", control="default")



  #########################################
  #### Pre-Building extended data structure for GH quadrature
  #######################################
  intPoints <- intPoints    # number of integration points
  ghPoints <- pracma::gaussHermite(intPoints)    # get GH weights and integration points

  # preparing additional variables for memory allocation
  mu.y1 <- mu.y2 <- mu.y3 <- rep(0, intPoints)
  f.y1 <- f.y2 <- f.y3 <- rep(0, intPoints)
  mu.z1 <- f.z1 <- LAMBDA <- f.counts <- rep(0, intPoints)

  # gather GH quadrature points and weights and additional variables
  extenDF <- data.frame(ghEta = ghPoints$x,
                        ghEtaFix = ghPoints$x,
                        ghWeights = ghPoints$w,
                        mu.y1, mu.y2, mu.y3,
                        f.y1, f.y2, f.y3,
                        mu.z1, f.z1,
                        LAMBDA, f.counts)
  id <- rep(c(1:N), each = intPoints)

  ## extend original data frame with GH quadrature  information
  Data <- cbind(ID = as.factor(id), dNA[id,], extenDF, out = NA)


  #######################################
  # Start values
  #####################################
  # x.start.mm <- c(1, 1, # lambda2,3
  #                 1, 1, 1, # theta1,2,3
  #                 0, 0) # nu2,3
  #
  # x.start.xk <- c(0, 1, #mu, psi in cell
  #                 3, 1, # parameters of count nb
  #                 0, 0.1, 0.1) # regression coefs y on count and latent

  x.start.kappa <- c(rep(log(N/8),8)) # 7 kappas for 8 cells

  # x.start <- c(x.start.mm,
  #              rep(x.start.xk, 8))

  # Results from previous optimization as starting values (to speed up optimization)
  x.start <- c(1.34687503867981, 1.31242266349582, 1.49916752743854, 1.35711495550324, 1.40995804145432, -0.162248625260641, -0.375687799952998, 1.37765554613443, 1.22790456993734, 6.31115184793967, 40.3536872491549, 1.31332383549179, -0.0187208071623842, 0.100391675122415, 1.63387529364496, 2.02685961860419, 5.95673076549388, 33.4550613315908, 1.21540365297411, -0.0315267899857264, 0.103776116096153, 1.36271296319571, 1.06459187616107, 7.42444508952016, 2.0876782417495, 1.32223187426111, -0.0135389408265504, 0.0952334303513592, 1.51857552169516, 1.81105825045145, 6.1624861586066, 27.1739656313697, 1.26725106693413, -0.015067899787885, 0.100767018631414, 1.62144157200804, 2.2476073479695, 6.74643341389871, 59.939553798775, 1.42012332278999, -0.0494706816481251, 0.100238064528339, 1.67337600978999, 1.84302077600215, 5.9704132613798, 34.9366516387958, 1.39764061323346, -0.0211885425635281, 0.095427117166944, 1.47304911107148, 1.87753743601293, 6.31813425036463, 49.8081258165067, 1.31454570669987, 0.00980088827365719, 0.0930780208373433, 1.58219280079135, 1.55473950422222, 5.96081130148063, 34.8675632548237, 1.26176850039523, -0.0103103003160813, 0.102758974627076)



  ########################################################################################
  # Optimization of log-likelihood
  #####################################################################################
  #####################################################
  ### Log-likelihood function for factorization case
  ####################################################

  ceff_loglik_independence <- function(x) {

    # No missing values in the parameter vector
    if(anyNA(x)){return(+Inf)}

    ## map 'x' to parameter matrices
    # invariant parameters
    # factor loadings
    lambda1 <- 1; lambda2 <- x[1]; lambda3 <- x[2]
    # residual variances
    theta1 <- x[3]; theta2 <- x[4]; theta3 <- x[5]
    # factor intercepts
    nu1 <- 0; nu2 <- x[6]; nu3 <- x[7]

    # XK-specific parameters (mu, sigma, gammas, alphas)
    xk.mat <- x[8:63] %>% matrix(ncol = 7, byrow = T)
    colnames(xk.mat) <- c("mu", "psi", "muc", "sizec", "a0", "a1", "a2")
    xk.mat <- as.data.frame(xk.mat)

    # no negative variances...
    if(theta1 <= 0 || theta2 <= 0 || theta3 <= 0 || any(xk.mat$psi <= 0) || any(xk.mat$sizec <= 0) ) return(+Inf)

    # compute xi_1* values of the Gauss-Hermite quadrature
    Data$ghEta <- sqrt(2*xk.mat$psi[Data$cell])*Data$ghEtaFix + xk.mat$mu[Data$cell]

    # expectation of z1/z2/z3 for this given value of xi_1*
    Data$mu.y1 <- nu1 + lambda1*Data$ghEta
    Data$mu.y2 <- nu2 + lambda2*Data$ghEta
    Data$mu.y3 <- nu3 + lambda3*Data$ghEta

    # f(z|\xi); likelihood of observed indicators
    Data$f.y1 <- dnorm(Data$ind1, mean = Data$mu.y1, sd = sqrt(theta1), log = FALSE)
    Data$f.y2 <- dnorm(Data$ind2, mean = Data$mu.y2, sd = sqrt(theta2), log = FALSE)
    Data$f.y3 <- dnorm(Data$ind3, mean = Data$mu.y3, sd = sqrt(theta3), log = FALSE)

    # f(xi_2) as negative binomial distributed
    # Data$f.z1 <- dnbinom(Data$pretest, mu = xk.mat$muc[Data$cell], size = xk.mat$sizec[Data$cell], log = FALSE)

    # lik counts Y for this given value of xi
    Data$LAMBDA <- exp(xk.mat$a0[Data$cell] + xk.mat$a1[Data$cell]*Data$ghEta + xk.mat$a2[Data$cell]*Data$pretest)
    Data$f.counts <- dpois(Data$dv, lambda = Data$LAMBDA, log = FALSE)

    # summands of Gauss-Hermite quadrature
    Data$out <- Data$ghWeights*Data$f.y1*Data$f.y2*Data$f.y3*Data$f.counts/sqrt(pi)

    # summing of Gauss-Hermite quadrature summands
    tmp <- setDT(Data)[, .(x = sum(out)), by=ID]
    if(anyNA(tmp)) return(+Inf)

    f.z1 <- dnbinom(dNA$pretest, mu = xk.mat$muc[dNA$cell], size = xk.mat$sizec[dNA$cell], log = FALSE)

    # computing and summing over the individual log-likelihood
    out <- tmp$x %>% log() %>% sum()
    out <- out + sum(log(f.z1))

    # rescaling for nlminb
    obj <- -1 * out

    # rescale to make nlminb happy
    obj <- obj / (nrow(Data)/intPoints)

    #cat("obj = ", obj, "\n")

    obj
  }



  ceff_loglik_independence_group <- function(x){
    # log-likelihood for group sizes

    # stochastic group size
    expkappas <- exp(x)
    n.cell <- table(Data$cell)/intPoints

    # Poisson group model
    obj_group <- sum(-expkappas + n.cell*log(expkappas) - lgamma(n.cell + 1))
    obj <- - obj_group/(nrow(Data)/intPoints)

    #cat("obj = ", obj, "\n")

    obj
  }


  ## Optimize log-likelihood for non-group part
  if (!silent){
    cat("Fitting the model...")
    time_start <- Sys.time()
  }
  time_start <- Sys.time()
  out <- nlminb(start = x.start,
                objective = ceff_loglik_independence,
                control = list(rel.tol = 1e-6)
  )

  ## Optimize log-likelihood for group part
  out_group <- nlminb(start = x.start.kappa,
                      objective = ceff_loglik_independence_group,
                      control = list(rel.tol = 1e-6)
  )

  par_final <- c(out$par, out_group$par)

  if (!silent){
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff,1), "s\n")
  }


  if (!silent){
    cat("Computing standard errors...")
    time_start <- Sys.time()
  }
  ## Standard Errors for Model (via numerical derivation)
  # for non-group part
  H <- pracma::hessian(ceff_loglik_independence, out$par)   # may take a while....
  I <- H   # Note: rescaling (-1* and /N) has already taken place in the loglik-function
  varcov <- solve(I)/N

  # for group part
  H_group <- pracma::hessian(ceff_loglik_independence_group, out_group$par)
  I_group <- H_group
  varcov_group <- solve(I_group)/N

  # join both vcovs (parameters from group and non-group part are independent)
  vcov_final <- lavaan:::lav_matrix_bdiag(varcov, varcov_group)
  # vcov_final <- readRDS("nacht2.rds")

  if (!silent){
    time_diff <-  Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff,1), "s\n")
  }

  #### Effect Estimation ####
  cat("Computing effects...")

  computeCondEffect <- function(x){
    ## This function computes conditional treatment effects within the XK-groups
    ## i.e., only the integration part given any combination of x and k
    ## aggregation to ATE and CE is done in another function

    # XK-specific parameters
    xk.mat <- x[8:63] %>% matrix(ncol = 7, byrow = T)
    xk.mat <- cbind(xk.mat, c(0,1), c(0,0,1,1,2,2,3,3))
    colnames(xk.mat) <- c("mu", "psi", "muc", "sizec", "a0", "a1", "a2", "k", "x")

    # stochastic group size
    kappas <- x[64:71]

    # re-arranging the parameters for  effect estimation
    gFuncs <- data.frame(xk.mat[,c("a0", "a1", "a2", "k", "x")])
    covFuncs <- data.frame(xk.mat[,c("mu", "psi", "muc", "sizec", "x", "k")])
    conditions <- expand.grid(k=c(0,1),x=c(0,1,2,3),gx=c(1,2,3), eff =NA)


    ## Effect for treatments t in any combination of x and k
    for (i in 1:nrow(conditions)){
      k <- conditions$k[i]
      x <- conditions$x[i]
      gx <- conditions$gx[i]

      # parameters of g-function
      A1 <- gFuncs[gFuncs$x == gx & gFuncs$k == k, c("a0", "a1", "a2")]
      A0 <- gFuncs[gFuncs$x == 0 & gFuncs$k == k, c("a0", "a1", "a2")]

      # parameters of covariates
      covPar <- covFuncs[covFuncs$x == x & covFuncs$k == k, ]

      # mgfs for of xi_1 (eta)
      mgf0_xi1 <- exp(A0$a1*covPar$mu + A0$a1^2*covPar$psi/2)
      mgf1_xi1 <- exp(A1$a1*covPar$mu + A1$a1^2*covPar$psi/2)

      # mgfs for xi_2 (pretest)
      muc <- covPar$muc
      s_sq <- covPar$muc + 1/covPar$sizec*covPar$muc^2
      r <- muc^2/(muc+s_sq)
      p <- muc/(muc+s_sq)
      mgf0_xi2 <- (p*exp(A0$a2)/(1-(1-p)*exp(A0$a2)))^r
      mgf1_xi2 <- (p*exp(A1$a2)/(1-(1-p)*exp(A1$a2)))^r

      # summands and summation of GH quadrature
      out <- exp(A1$a0)*mgf1_xi1*mgf1_xi2 - exp(A0$a0)*mgf0_xi1*mgf0_xi2
      conditions$eff[i] <- sum(out)
    }



    # computing the probabilities P(X=x,K=k)
    Pxk <- exp(kappas)/(sum(exp(kappas)))

    return(c(conditions$eff, Pxk))

  }
  # apply basic effect estimation function to estimated parameters
  ce <- computeCondEffect(par_final)

  # standard error for CondEffect function
  # computes variance of xk-conditional effects and group probabilities (transformations of model parameters)
  # delta method
  delta.method.ck <- function(func, x, varcov){
    func.jacobian <- lav_func_jacobian_complex(func, x)
    func.jacobian %*% varcov %*% t(func.jacobian)
  }
  vcov.ce <- delta.method.ck(computeCondEffect, par_final, varcov = vcov_final)
  vcov.ce <- lav_matrix_bdiag(vcov.ce[1:24,1:24], vcov.ce[25:32,25:32])
  ce.se <- sqrt(diag(vcov.ce))


  #################################################
  ####### Lavaan syntax for aggregated effects
  lsyn <- '
# Definition of free parameters re. conditional effects
## Note: These "free" parameters are plugged in from our own model estimation
g1 ~ CE100*xk00 + CE101*xk01 + CE110*xk10 + CE111*xk11 + CE120*xk20 + CE121*xk21 + CE130*xk30 + CE131*xk31
g2 ~ CE200*xk00 + CE201*xk01 + CE210*xk10 + CE211*xk11 + CE220*xk20 + CE221*xk21 + CE230*xk30 + CE231*xk31
g3 ~ CE300*xk00 + CE301*xk01 + CE310*xk10 + CE311*xk11 + CE320*xk20 + CE321*xk21 + CE330*xk30 + CE331*xk31
p1 ~ relfreq00*xk00 + relfreq01*xk01 + relfreq10*xk10 + relfreq11*xk11 + relfreq20*xk20 + relfreq21*xk21 + relfreq30*xk30 + relfreq31*xk31


## Unconditional Probabilities P(K=k)
Pk0 := relfreq00 + relfreq10 + relfreq20 + relfreq30
Pk1 := relfreq01 + relfreq11 + relfreq21 + relfreq31

## Unconditional Probabilities P(X=x)
Px0 := relfreq00 + relfreq01
Px1 := relfreq10 + relfreq11
Px2 := relfreq20 + relfreq21
Px3 := relfreq30 + relfreq31

## Conditional Probabilities P(K=k|X=x)
Pk0gx0 := relfreq00/Px0
Pk1gx0 := relfreq01/Px0
Pk0gx1 := relfreq10/Px1
Pk1gx1 := relfreq11/Px1
Pk0gx2 := relfreq20/Px2
Pk1gx2 := relfreq21/Px2
Pk0gx3 := relfreq30/Px3
Pk1gx3 := relfreq31/Px3

## Conditional Probabilities P(X=x|K=k)
Px0gk0 := relfreq00/Pk0
Px1gk0 := relfreq10/Pk0
Px2gk0 := relfreq20/Pk0
Px3gk0 := relfreq30/Pk0
Px0gk1 := relfreq01/Pk1
Px1gk1 := relfreq11/Pk1
Px2gk1 := relfreq21/Pk1
Px3gk1 := relfreq31/Pk1

# Average treatment effects
ave10 := CE100*relfreq00 + CE101*relfreq01 + CE110*relfreq10 + CE111*relfreq11 + CE120*relfreq20 + CE121*relfreq21 + CE130*relfreq30 + CE131*relfreq31
ave20 := CE200*relfreq00 + CE201*relfreq01 + CE210*relfreq10 + CE211*relfreq11 + CE220*relfreq20 + CE221*relfreq21 + CE230*relfreq30 + CE231*relfreq31
ave30 := CE300*relfreq00 + CE301*relfreq01 + CE310*relfreq10 + CE311*relfreq11 + CE320*relfreq20 + CE321*relfreq21 + CE330*relfreq30 + CE331*relfreq31

# K=k conditional effects
CE10gk0 := Px0gk0*CE100 + Px1gk0*CE110 + Px2gk0*CE120 + Px3gk0*CE130
CE20gk0 := Px0gk0*CE200 + Px1gk0*CE210 + Px2gk0*CE220 + Px3gk0*CE230
CE30gk0 := Px0gk0*CE300 + Px1gk0*CE310 + Px2gk0*CE320 + Px3gk0*CE330
CE10gk1 := Px0gk1*CE101 + Px1gk1*CE111 + Px2gk1*CE121 + Px3gk1*CE131
CE20gk1 := Px0gk1*CE201 + Px1gk1*CE211 + Px2gk1*CE221 + Px3gk1*CE231
CE30gk1 := Px0gk1*CE301 + Px1gk1*CE311 + Px2gk1*CE321 + Px3gk1*CE331

# X=x conditional effects
CE10gx0 := Pk0gx0*CE100 + Pk1gx0*CE101
CE20gx0 := Pk0gx0*CE200 + Pk1gx0*CE201
CE30gx0 := Pk0gx0*CE300 + Pk1gx0*CE301

CE10gx1 := Pk0gx1*CE110 + Pk1gx1*CE111
CE20gx1 := Pk0gx1*CE210 + Pk1gx1*CE211
CE30gx1 := Pk0gx1*CE310 + Pk1gx1*CE311

CE10gx2 := Pk0gx2*CE120 + Pk1gx1*CE121
CE20gx2 := Pk0gx2*CE220 + Pk1gx1*CE221
CE30gx2 := Pk0gx2*CE320 + Pk1gx1*CE321

CE10gx3 := Pk0gx3*CE130 + Pk1gx3*CE131
CE20gx3 := Pk0gx3*CE230 + Pk1gx3*CE231
CE30gx3 := Pk0gx3*CE330 + Pk1gx3*CE331
'
  # laavanify aggregation syntax (create partable)
  ## computes aggregated effect and corresponding standard errors
  pt <- lavaanify(lsyn)
  def.function <- lav_partable_constraints_def(pt)
  JAC <- lav_func_jacobian_complex(func=def.function, x=ce)
  info.r <- JAC %*% vcov.ce %*% t(JAC)

  # Extract the point estimates and standard errors
  est <- def.function(ce)
  se <- sqrt(diag(info.r))

  sdyx0 <- sd(dNA$dv[dNA$treat==0])

  ## show effects, standard errors and effect sizes
  tval <- est/se
  n_par <- length(par_final)
  rdf <- nrow(object@input@data) - n_par

  pval <- 2*(1-pt(abs(tval),df=rdf))

  # conditional effects given xk
  ce.tval <- ce/ce.se
  ce.pval <- 2*(1-pt(abs(ce.tval),df=rdf))

  ## Average effects
  selector <- c(23:25)
  Egx <- data.frame(est[selector],
                    se[selector],
                    tval[selector],
                    pval[selector],
                    est[selector]/sdyx0)

  names(Egx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given a treatment condition
  selector <- c(32:43)
  Egxgx <- data.frame(est[selector],
                      se[selector],
                      tval[selector],
                      pval[selector],
                      est[selector]/sdyx0)
  names(Egxgx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given K=k
  selector <- c(26:31)
  Egxgk <- data.frame(est[selector],
                      se[selector],
                      tval[selector],
                      pval[selector],
                      est[selector]/sdyx0)
  names(Egxgk) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given (X=x, K=k)
  selector <- c(1:(length(ce)-8))
  Egxgxk <- data.frame(ce[selector],
                       ce.se[selector],
                       ce.tval[selector],
                       ce.pval[selector],
                       ce[selector]/sdyx0)
  names(Egxgxk) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  object@results <- new("results",
                        est=est,
                        se=se,
                        vcov_def=info.r,
                        Egx=Egx,
                        Egxgx=Egxgx,
                        Egxgk=Egxgk,
                        Egxgxk=Egxgxk
  )
  cat("done\n")

  cat("Finished.\n")
  return(object)
}





###########################################
###################################
## CASE 2: FACTORIZATION
###################################
###########################################
ceff_case_factorization <- function(data = dNA, intPoints = 21L, silent = FALSE){
  msg <- "This function is implemented for illustration only. It accompanies a paper by Kiefer & Mayer in the British Journal of Mathematical and Statistical Psychology. Its functionality is meant to be implemented within the main countEffects()-function in the future.\n Please stay tuned for further updates.\n"
  cat(msg)

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package \"dplyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  require(dplyr)

  ######################################
  # readData
  #####################################
  # data <- readRDS("dACTIVE.rds")
  dNA <- na.omit(data)
  N <- nrow(dNA)

  object <- new("countEffects")
  object@input <- ceff_create_input(y="dv", x="treat", k="gender", z=c("depression", "pretest"), dNA,
                                    measurement=list(depression=c("ind1", "ind2", "ind3")), method="poisson", distribution="condNormal", control="default")



  #########################################
  #### Pre-Building extended data structure for GH quadrature
  #######################################
  intPoints <- intPoints    # number of integration points
  ghPoints <- pracma::gaussHermite(intPoints)    # get GH weights and integration points

  # preparing additional variables for memory allocation
  mu.y1 <- mu.y2 <- mu.y3 <- rep(0, intPoints)
  f.y1 <- f.y2 <- f.y3 <- rep(0, intPoints)
  mu.z1 <- f.z1 <- LAMBDA <- f.counts <- rep(0, intPoints)

  # gather GH quadrature points and weights and additional variables
  extenDF <- data.frame(ghEta = ghPoints$x,
                        ghEtaFix = ghPoints$x,
                        ghWeights = ghPoints$w,
                        mu.y1, mu.y2, mu.y3,
                        f.y1, f.y2, f.y3,
                        mu.z1, f.z1,
                        LAMBDA, f.counts)
  id <- rep(c(1:N), each = intPoints)

  ## extend original data frame with GH quadrature  information
  Data <- cbind(ID = as.factor(id), dNA[id,], extenDF, out = NA)


  #######################################
  # Start values
  #####################################
  # x.start.mm <- c(1, 1, # lambda2,3
  #                 1, 1, 1, # theta1,2,3
  #                 0, 0) # nu2,3
  #
  # x.start.xk <- c(0, 1, #mu, psi in cell
  #                 0, 0.1, # regression coefs count on latent
  #                 0, 0.1, 0.1) # regression coefs y on count and latent

  x.start.kappa <- c(rep(log(N/8),8)) # 7 kappas for 8 cells

  # x.start <- c(x.start.mm,
  #              rep(x.start.xk, 8))

  # Results from previous optimization as starting values (to speed up optimization)
  x.start <- c(1.32943124768397, 1.30018177817337, 1.47639768977471, 1.3856367212827, 1.41948083663654,
              -0.131964781512757, -0.354817080098272, 1.36823295333983, 1.28961536497266, 1.92093232433974,
              -0.0569845310911216, 1.30693890928587, -0.0174992075995661, 0.100876995545554, 1.64510144437603,
               2.0718692777033, 1.95603809461613, -0.110331425847837, 1.22513956132144, -0.0326942990872267,
               0.102529817650243, 1.3539235859802, 1.1331087671473, 1.97439951366681, -0.0517725892773255,
               1.30834442775628, -0.0141570410923926, 0.096989573777932, 1.51200115968336, 1.86425091408452,
               1.87771827938006, -0.0380857468066551, 1.2610770971638, -0.0139993172855461, 0.101302391253577,
               1.62335899345804, 2.01810850262393, 2.11206950009893, -0.135091767833275, 1.43082982075035,
              -0.0502360222099185, 0.0987651364804679, 1.67265655433518, 1.89435176186378, 1.78083809216978,
               0.00374384564206068, 1.39983631580143, -0.0213913816024062, 0.0951945638912023, 1.44445099970489,
               1.78927922558093, 1.95989703496735, -0.0843238847461746, 1.30490467937624, 0.0105282209468848,
               0.0942448741160402, 1.57665617453619, 1.59346857268003, 1.84202366873565, -0.0371712529208222,
               1.26606577189979, -0.0110893247908247, 0.102306186078173)

  ########################################################################################
  # Optimization of log-likelihood
  #####################################################################################
  #####################################################
  ### Log-likelihood function for factorization case
  ####################################################

  ceff_loglik_factorization <- function(x) {

    # No missing values in the parameter vector
    if(anyNA(x)){return(+Inf)}

    ## map 'x' to parameter matrices
    # invariant parameters
    # factor loadings
    lambda1 <- 1; lambda2 <- x[1]; lambda3 <- x[2]
    # residual variances
    theta1 <- x[3]; theta2 <- x[4]; theta3 <- x[5]
    # factor intercepts
    nu1 <- 0; nu2 <- x[6]; nu3 <- x[7]

    # XK-specific parameters (mu, sigma, gammas, alphas)
    xk.mat <- x[8:63] %>% matrix(ncol = 7, byrow = T)
    colnames(xk.mat) <- c("mu", "psi", "g0", "g1", "a0", "a1", "a2")
    xk.mat <- as.data.frame(xk.mat)

    # no negative variances...
    if(theta1 <= 0 || theta2 <= 0 || theta3 <= 0 || any(xk.mat$psi <= 0) ) return(+Inf)

    # compute xi_1* values of the Gauss-Hermite quadrature
    Data$ghEta <- sqrt(2*xk.mat$psi[Data$cell])*Data$ghEtaFix + xk.mat$mu[Data$cell]

    # expectation of z1/z2/z3 for this given value of xi_1*
    Data$mu.y1 <- nu1 + lambda1*Data$ghEta
    Data$mu.y2 <- nu2 + lambda2*Data$ghEta
    Data$mu.y3 <- nu3 + lambda3*Data$ghEta

    # f(z|\xi); likelihood of observed indicators
    Data$f.y1 <- dnorm(Data$ind1, mean = Data$mu.y1, sd = sqrt(theta1), log = FALSE)
    Data$f.y2 <- dnorm(Data$ind2, mean = Data$mu.y2, sd = sqrt(theta2), log = FALSE)
    Data$f.y3 <- dnorm(Data$ind3, mean = Data$mu.y3, sd = sqrt(theta3), log = FALSE)

    # f(xi_2|x_1) as Poisson regression
    Data$mu.z1 <- exp(xk.mat$g0[Data$cell] + xk.mat$g1[Data$cell]*Data$ghEta)
    Data$f.z1 <- dpois(Data$pretest, lambda = Data$mu.z1, log = FALSE)

    # lik counts Y for this given value of xi
    Data$LAMBDA <- exp(xk.mat$a0[Data$cell] + xk.mat$a1[Data$cell]*Data$ghEta + xk.mat$a2[Data$cell]*Data$pretest)
    Data$f.counts <- dpois(Data$dv, lambda = Data$LAMBDA, log = FALSE)

    # summands of Gauss-Hermite quadrature
    Data$out <- Data$ghWeights*Data$f.z1*Data$f.y1*Data$f.y2*Data$f.y3*Data$f.counts/sqrt(pi)

    # summing of Gauss-Hermite quadrature summands
    tmp <- setDT(Data)[, .(x = sum(out)), by=ID]
    if(anyNA(tmp)) return(+Inf)

    # computing and summing over the individual log-likelihood
    out <- tmp$x %>% log() %>% sum()

    # rescaling for nlminb
    obj <- -1 * out

    # rescale to make nlminb happy
    obj <- obj / (nrow(Data)/intPoints)

    #cat("obj = ", obj, "\n")

    obj
  }



  ceff_loglik_factorization_group <- function(x){
    # log-likelihood for group sizes

    # stochastic group size
    expkappas <- exp(x)
    n.cell <- table(Data$cell)/intPoints

    # Poisson group model
    obj_group <- sum(-expkappas + n.cell*log(expkappas) - lgamma(n.cell + 1))
    obj <- - obj_group/(nrow(Data)/intPoints)

    #cat("obj = ", obj, "\n")

    obj
  }


  ## Optimize log-likelihood for non-group part
  if (!silent){
    cat("Fitting the model...")
    time_start <- Sys.time()
  }
  time_start <- Sys.time()
  out <- nlminb(start = x.start,
                objective = ceff_loglik_factorization,
                control = list(rel.tol = 1e-6)
  )
  # round(out$par, 3)
  #[1]  1.329  1.300  1.476  1.386  1.419 -0.132 -0.355  1.368  1.290  1.921 -0.057  1.307 -0.017  0.101  1.645  2.072
  #[17]  1.956 -0.110  1.225 -0.033  0.103  1.354  1.133  1.974 -0.052  1.308 -0.014  0.097  1.512  1.864  1.878 -0.038
  #[33]  1.261 -0.014  0.101  1.623  2.018  2.112 -0.135  1.431 -0.050  0.099  1.673  1.894  1.781  0.004  1.400 -0.021
  #[49]  0.095  1.444  1.789  1.960 -0.084  1.305  0.011  0.094  1.577  1.593  1.842 -0.037  1.266 -0.011  0.102

  ## Optimize log-likelihood for group part
  out_group <- nlminb(start = x.start.kappa,
                      objective = ceff_loglik_factorization_group,
                      control = list(rel.tol = 1e-6)
  )
  # round(out_group$par, 3)
  #4.710 5.746 4.605 5.838 4.605 5.823 4.605 5.855

  par_final <- c(out$par, out_group$par)
  if (!silent){
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff,1), "s\n")
  }


  if (!silent){
    cat("Computing standard errors...")
    time_start <- Sys.time()
  }
  ## Standard Errors for Model (via numerical derivation)
  # for non-group part
  H <- pracma::hessian(ceff_loglik_factorization, out$par)   # may take a while....
  I <- H   # Note: rescaling (-1* and /N) has already taken place in the loglik-function
  varcov <- solve(I)/N

  # for group part
  H_group <- pracma::hessian(ceff_loglik_factorization_group, out_group$par)
  I_group <- H_group
  varcov_group <- solve(I_group)/N

  # join both vcovs (parameters from group and non-group part are independent)
  vcov_final <- lavaan:::lav_matrix_bdiag(varcov, varcov_group)
  if (!silent){
    time_diff <-  Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff,1), "s\n")
  }


  object@fit@fit <- list(par_final = par_final, vcov_final = vcov_final)

  #### Effect Estimation ####
  cat("Computing effects...")

  computeCondEffect <- function(x){
    ## This function computes conditional treatment effects within the XK-groups
    ## i.e., only the integration part given any combination of x and k
    ## aggregation to ATE and CE is done in another function

    # XK-specific parameters
    xk.mat <- x[8:63] %>% matrix(ncol = 7, byrow = T)
    xk.mat <- cbind(xk.mat, c(0,1), c(0,0,1,1,2,2,3,3))
    colnames(xk.mat) <- c("mu", "psi", "g0", "g1", "a0", "a1", "a2", "k", "x")

    # stochastic group size
    kappas <- x[64:71]

    # re-arranging the parameters for  effect estimation
    gFuncs <- data.frame(xk.mat[,c("a0", "a1", "a2", "k", "x")])
    covFuncs <- data.frame(xk.mat[,c("mu", "psi", "g0", "g1", "x", "k")])
    conditions <- expand.grid(k=c(0,1),x=c(0,1,2,3),gx=c(1,2,3), eff =NA)

    ## GH quadrature init
    intPoints <- 20
    ghPoints <- pracma::gaussHermite(intPoints)

    ## Effect for treatments t in any combination of x and k
    for (i in 1:nrow(conditions)){
      k <- conditions$k[i]
      x <- conditions$x[i]
      gx <- conditions$gx[i]

      # parameters of g-function
      A1 <- gFuncs[gFuncs$x == gx & gFuncs$k == k, c("a0", "a1", "a2")]
      A0 <- gFuncs[gFuncs$x == 0 & gFuncs$k == k, c("a0", "a1", "a2")]

      # parameters of covariates
      covPar <- covFuncs[covFuncs$x == x & covFuncs$k == k, ]

      # group-specific values of xi_1*
      eta <- sqrt(2*covPar$psi)*ghPoints$x + covPar$mu

      # Poisson regression of xi_2 on xi_1
      muPre <- exp(covPar$g0 + covPar$g1*eta)

      # moment-generating functions for xi_2
      mgf1 <- exp((muPre)*(exp(A1$a2) - 1))
      mgf0 <- exp((muPre)*(exp(A0$a2) - 1))

      # summands and summation of GH quadrature
      out <- ( exp(A1$a0 + A1$a1*eta)*mgf1 - exp(A0$a0 + A0$a1*eta)*mgf0 )*ghPoints$w/sqrt(pi)
      conditions$eff[i] <- sum(out)
    }



    # computing the probabilities P(X=x,K=k)
    Pxk <- exp(kappas)/(sum(exp(kappas)))

    return(c(conditions$eff, Pxk))

  }
  # apply basic effect estimation function to estimated parameters
  ce <- computeCondEffect(par_final)

  # standard error for CondEffect function
  # computes variance of xk-conditional effects and group probabilities (transformations of model parameters)
  # delta method
  delta.method.ck <- function(func, x, varcov){
    func.jacobian <- lav_func_jacobian_complex(func, x)
    func.jacobian %*% varcov %*% t(func.jacobian)
  }
  vcov.ce <- delta.method.ck(computeCondEffect, par_final, varcov = vcov_final)
  vcov.ce <- lav_matrix_bdiag(vcov.ce[1:24,1:24], vcov.ce[25:32,25:32])
  ce.se <- sqrt(diag(vcov.ce))


  #################################################
  ####### Lavaan syntax for aggregated effects
  lsyn <- '
# Definition of free parameters re. conditional effects
## Note: These "free" parameters are plugged in from our own model estimation
g1 ~ CE100*xk00 + CE101*xk01 + CE110*xk10 + CE111*xk11 + CE120*xk20 + CE121*xk21 + CE130*xk30 + CE131*xk31
g2 ~ CE200*xk00 + CE201*xk01 + CE210*xk10 + CE211*xk11 + CE220*xk20 + CE221*xk21 + CE230*xk30 + CE231*xk31
g3 ~ CE300*xk00 + CE301*xk01 + CE310*xk10 + CE311*xk11 + CE320*xk20 + CE321*xk21 + CE330*xk30 + CE331*xk31
p1 ~ relfreq00*xk00 + relfreq01*xk01 + relfreq10*xk10 + relfreq11*xk11 + relfreq20*xk20 + relfreq21*xk21 + relfreq30*xk30 + relfreq31*xk31


## Unconditional Probabilities P(K=k)
Pk0 := relfreq00 + relfreq10 + relfreq20 + relfreq30
Pk1 := relfreq01 + relfreq11 + relfreq21 + relfreq31

## Unconditional Probabilities P(X=x)
Px0 := relfreq00 + relfreq01
Px1 := relfreq10 + relfreq11
Px2 := relfreq20 + relfreq21
Px3 := relfreq30 + relfreq31

## Conditional Probabilities P(K=k|X=x)
Pk0gx0 := relfreq00/Px0
Pk1gx0 := relfreq01/Px0
Pk0gx1 := relfreq10/Px1
Pk1gx1 := relfreq11/Px1
Pk0gx2 := relfreq20/Px2
Pk1gx2 := relfreq21/Px2
Pk0gx3 := relfreq30/Px3
Pk1gx3 := relfreq31/Px3

## Conditional Probabilities P(X=x|K=k)
Px0gk0 := relfreq00/Pk0
Px1gk0 := relfreq10/Pk0
Px2gk0 := relfreq20/Pk0
Px3gk0 := relfreq30/Pk0
Px0gk1 := relfreq01/Pk1
Px1gk1 := relfreq11/Pk1
Px2gk1 := relfreq21/Pk1
Px3gk1 := relfreq31/Pk1

# Average treatment effects
ave10 := CE100*relfreq00 + CE101*relfreq01 + CE110*relfreq10 + CE111*relfreq11 + CE120*relfreq20 + CE121*relfreq21 + CE130*relfreq30 + CE131*relfreq31
ave20 := CE200*relfreq00 + CE201*relfreq01 + CE210*relfreq10 + CE211*relfreq11 + CE220*relfreq20 + CE221*relfreq21 + CE230*relfreq30 + CE231*relfreq31
ave30 := CE300*relfreq00 + CE301*relfreq01 + CE310*relfreq10 + CE311*relfreq11 + CE320*relfreq20 + CE321*relfreq21 + CE330*relfreq30 + CE331*relfreq31

# K=k conditional effects
CE10gk0 := Px0gk0*CE100 + Px1gk0*CE110 + Px2gk0*CE120 + Px3gk0*CE130
CE20gk0 := Px0gk0*CE200 + Px1gk0*CE210 + Px2gk0*CE220 + Px3gk0*CE230
CE30gk0 := Px0gk0*CE300 + Px1gk0*CE310 + Px2gk0*CE320 + Px3gk0*CE330
CE10gk1 := Px0gk1*CE101 + Px1gk1*CE111 + Px2gk1*CE121 + Px3gk1*CE131
CE20gk1 := Px0gk1*CE201 + Px1gk1*CE211 + Px2gk1*CE221 + Px3gk1*CE231
CE30gk1 := Px0gk1*CE301 + Px1gk1*CE311 + Px2gk1*CE321 + Px3gk1*CE331

# X=x conditional effects
CE10gx0 := Pk0gx0*CE100 + Pk1gx0*CE101
CE20gx0 := Pk0gx0*CE200 + Pk1gx0*CE201
CE30gx0 := Pk0gx0*CE300 + Pk1gx0*CE301

CE10gx1 := Pk0gx1*CE110 + Pk1gx1*CE111
CE20gx1 := Pk0gx1*CE210 + Pk1gx1*CE211
CE30gx1 := Pk0gx1*CE310 + Pk1gx1*CE311

CE10gx2 := Pk0gx2*CE120 + Pk1gx1*CE121
CE20gx2 := Pk0gx2*CE220 + Pk1gx1*CE221
CE30gx2 := Pk0gx2*CE320 + Pk1gx1*CE321

CE10gx3 := Pk0gx3*CE130 + Pk1gx3*CE131
CE20gx3 := Pk0gx3*CE230 + Pk1gx3*CE231
CE30gx3 := Pk0gx3*CE330 + Pk1gx3*CE331
'
  # laavanify aggregation syntax (create partable)
  ## computes aggregated effect and corresponding standard errors
  pt <- lavaanify(lsyn)
  def.function <- lav_partable_constraints_def(pt)
  JAC <- lav_func_jacobian_complex(func=def.function, x=ce)
  info.r <- JAC %*% vcov.ce %*% t(JAC)

  # Extract the point estimates and standard errors
  est <- def.function(ce)
  se <- sqrt(diag(info.r))

  sdyx0 <- sd(dNA$dv[dNA$treat==0])

  ## show effects, standard errors and effect sizes
  tval <- est/se
  n_par <- length(par_final)
  rdf <- nrow(object@input@data) - n_par
  pval <- 2*(1-pt(abs(tval),df=rdf))

  ce.tval <- ce/ce.se
  ce.pval <- 2*(1-pt(abs(ce.tval),df=rdf))

  ## Average effects
  selector <- c(23:25)
  Egx <- data.frame(est[selector],
                    se[selector],
                    tval[selector],
                    pval[selector],
                    est[selector]/sdyx0)

  names(Egx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given a treatment condition
  selector <- c(32:43)
  Egxgx <- data.frame(est[selector],
                      se[selector],
                      tval[selector],
                      pval[selector],
                      est[selector]/sdyx0)
  names(Egxgx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given K=k
  selector <- c(26:31)
  Egxgk <- data.frame(est[selector],
                      se[selector],
                      tval[selector],
                      pval[selector],
                      est[selector]/sdyx0)
  names(Egxgk) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given (X=x, K=k)
  selector <- c(1:(length(ce)-8))
  Egxgxk <- data.frame(ce[selector],
                      ce.se[selector],
                      ce.tval[selector],
                      ce.pval[selector],
                      ce[selector]/sdyx0)
  names(Egxgxk) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  object@results <- new("results",
             est=est,
             se=se,
             vcov_def=info.r,
             Egx=Egx,
             Egxgx=Egxgx,
             Egxgk=Egxgk,
             Egxgxk=Egxgxk
  )
  cat("done\n")

 cat("Finished.\n")
 return(object)
}







ceff_case_mixture <- function(data = dNA, intPoints = 21L, silent = FALSE){
  msg <- "This function is implemented for illustration only. It accompanies a paper by Kiefer & Mayer in the British Journal of Mathematical and Statistical Psychology. Its functionality is meant to be implemented within the main countEffects()-function in the future.\n Please stay tuned for further updates.\n"
  cat(msg)

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package \"dplyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  require(dplyr)

  ######################################
  # readData
  #####################################
  # data <- readRDS("dACTIVE.rds")
  dNA <- na.omit(data)
  N <- nrow(dNA)

  object <- new("countEffects")
  object@input <- ceff_create_input(y="dv", x="treat", k="gender", z=c("depression", "pretest"), dNA,
                                    measurement=list(depression=c("ind1", "ind2", "ind3")), method="poisson", distribution="condNormal", control="default")



  #########################################
  #### Pre-Building extended data structure for GH quadrature
  #######################################
  intPoints <- intPoints    # number of integration points
  ghPoints <- pracma::gaussHermite(intPoints)    # get GH weights and integration points

  # preparing additional variables for memory allocation
  mu.y1 <- mu.y2 <- mu.y3 <- rep(0, intPoints)
  f.y1 <- f.y2 <- f.y3 <- rep(0, intPoints)
  mu.z1 <- f.z1 <- LAMBDA <- f.counts <- rep(0, intPoints)

  # gather GH quadrature points and weights and additional variables
  extenDF <- data.frame(ghEta = ghPoints$x,
                        ghEtaFix = ghPoints$x,
                        ghWeights = ghPoints$w,
                        mu.y1, mu.y2, mu.y3,
                        f.y1, f.y2, f.y3,
                        mu.z1, f.z1,
                        LAMBDA, f.counts)
  id <- rep(c(1:N), each = intPoints)

  ## extend original data frame with GH quadrature  information
  Data <- cbind(ID = as.factor(id), dNA[id,], extenDF, out = NA)


  #######################################
  # Starting values
  # (here: some results from Mplus)
  # optimization still takes some time though
  #####################################
  # x.start.mm <- c(1.336, 1.309, # lambda2,3
  #                 1.484, 1.381, 1.400, # theta1,2,3
  #                 -0.143, -0.369) # nu2,3
  #
  # x.start.xk <- c(1.380, 1.274, #mu, psi for eta1 in cell
  #                 4.410, 7.921, 3.666, 5.686, #mus, psis for pretest
  #                 -0.482, -0.500,   # covs
  #                 1.310, -0.019, 0.101) # regression coefs y on count and latent
  #
  #
  x.start.kappa <- c(rep(log(N/8),8)) # 7 kappas for 8 cells
  #
  # x.start <- c(x.start.mm,
  #              rep(x.start.xk, 8))


  # Results from a previous optimization (to speed up the optimization here)
  x.start <- c(1.34128804620556, 1.30411102630817, 1.48400502631461, 1.36049815235144, 1.42898342676552,
              -0.15333337206284, -0.361719576044592, 1.37888042838175, 1.23033225620449, 4.42440793765966,
               7.98138674649311, 4.2465263525931, 5.3454226018001, 0.493795968805661, -0.930433988924702,
               1.30802210163957, -0.0173407285181311, 0.100734578871963, 1.64828188256694, 2.04576016795302,
               4.12899893490098, 7.85575768585768, 2.97922880864956, 6.61507568613766, -0.822776394687059,
              -1.17748072639865, 1.22433371020993, -0.0323884228106366, 0.10271652314774, 1.36781805093858,
               1.07846343500452, 5.76519655892046, 7.64146618444128, 4.13462109018567, 5.06603598479627,
              -0.0695724357900625, -0.72162201303733, 1.31212642931104, -0.013322122312099, 0.09641372475006,
               1.50404815579453, 1.83378956876989, 4.21458647919949, 8.13363725805329, 3.9319608974395,
               5.17923940449724, -0.741810853043313, -0.175088485530148, 1.2588583472501, -0.0138803691098876,
               0.101552059374981, 1.62466695456078, 2.23616533644154, 5.11265585431597, 8.31533671474113,
               3.96360577511955, 5.40516759032496, -1.00180734790174, -2.22224806666234, 1.43094502711645,
              -0.049677891523422, 0.0986701975004455, 1.67587206342194, 1.86996406954488, 4.27780388767066,
               7.63803585419686, 4.02248615155709, 5.6220398027653, -0.0883428606650102, 0.0616394226501068,
               1.39516053699646, -0.0207606804334102, 0.0956549341696078, 1.44670901628594, 1.86187148542641,
               4.51267681153719, 8.18202883558604, 3.00840235461959, 6.16899544195779, -0.281341241397972,
              -0.928230743460014, 1.31013656321926, 0.0102329749698341, 0.0935812776264315, 1.57339281502786,
               1.55950813778349, 4.45032769690392, 7.45980065581257, 2.80776822719836, 6.08223407783239,
              -0.207584491092768, -0.452144491763978, 1.27274496896602, -0.0115459896679967, 0.101578173257789)




  ########################################################################################
  # Optimization of log-likelihood
  #####################################################################################
  #####################################################
  ### Log-likelihood function for factorization case
  ####################################################

  ceff_loglik_mixture <- function(x) {

    # No missing values in the parameter vector
    if(anyNA(x)){return(+Inf)}

    ## map 'x' to parameter matrices
    # invariant parameters
    # factor loadings
    lambda1 <- 1; lambda2 <- x[1]; lambda3 <- x[2]
    # residual variances
    theta1 <- x[3]; theta2 <- x[4]; theta3 <- x[5]
    # factor intercepts
    nu1 <- 0; nu2 <- x[6]; nu3 <- x[7]

    # XK-specific parameters (mu, sigma, gammas, alphas)
    xk.mat <- x[8:95] %>% matrix(ncol = 11, byrow = T)
    colnames(xk.mat) <- c("mu1", "psi1", "mu21", "mu22", "psi21", "psi22", "cov1", "cov2", "a0", "a1", "a2")
    xk.mat <- as.data.frame(xk.mat)

    # no negative variances...
    if(theta1 <= 0 || theta2 <= 0 || theta3 <= 0 || any(xk.mat$psi1 <= 0) || any(xk.mat$psi21 <= 0) || any(xk.mat$psi22 <= 0) ) return(+Inf)

    # compute xi_1* values of the Gauss-Hermite quadrature
    Data$ghEta <- sqrt(2*xk.mat$psi1[Data$cell])*Data$ghEtaFix + xk.mat$mu1[Data$cell]

    # expectation of z1/z2/z3 for this given value of xi_1*
    Data$mu.y1 <- nu1 + lambda1*Data$ghEta
    Data$mu.y2 <- nu2 + lambda2*Data$ghEta
    Data$mu.y3 <- nu3 + lambda3*Data$ghEta

    # f(z|\xi); likelihood of observed indicators
    Data$f.y1 <- dnorm(Data$ind1, mean = Data$mu.y1, sd = sqrt(theta1), log = FALSE)
    Data$f.y2 <- dnorm(Data$ind2, mean = Data$mu.y2, sd = sqrt(theta2), log = FALSE)
    Data$f.y3 <- dnorm(Data$ind3, mean = Data$mu.y3, sd = sqrt(theta3), log = FALSE)

    # f(xi_2|x_1) as Poisson regression
    Data$mu.z11 <- xk.mat$mu21[Data$cell] + xk.mat$cov1[Data$cell]/xk.mat$psi1[Data$cell]*(Data$ghEta - xk.mat$mu1[Data$cell])
    Data$mu.z12 <- xk.mat$mu22[Data$cell] + xk.mat$cov2[Data$cell]/xk.mat$psi1[Data$cell]*(Data$ghEta - xk.mat$mu1[Data$cell])
    Data$sig.z11 <- xk.mat$psi21[Data$cell] - xk.mat$cov1[Data$cell]^2/xk.mat$psi1[Data$cell]
    Data$sig.z12 <- xk.mat$psi22[Data$cell] - xk.mat$cov2[Data$cell]^2/xk.mat$psi1[Data$cell]

    Data$f.z11 <- dnorm(Data$pretest, mean = Data$mu.z11, sd = sqrt(Data$sig.z11), log = FALSE)
    Data$f.z12 <- dnorm(Data$pretest, mean = Data$mu.z12, sd = sqrt(Data$sig.z12), log = FALSE)
    Data$f.z1 <- (Data$f.z11 + Data$f.z12)/2

    # lik counts Y for this given value of xi
    Data$LAMBDA <- exp(xk.mat$a0[Data$cell] + xk.mat$a1[Data$cell]*Data$ghEta + xk.mat$a2[Data$cell]*Data$pretest)
    Data$f.counts <- dpois(Data$dv, lambda = Data$LAMBDA, log = FALSE)


    # summands of Gauss-Hermite quadrature
    Data$out <- Data$ghWeights*Data$f.z1*Data$f.y1*Data$f.y2*Data$f.y3*Data$f.counts/sqrt(pi)

    # summing of Gauss-Hermite quadrature summands
    tmp <- setDT(Data)[, .(x = sum(out)), by=ID]
    if(anyNA(tmp)) return(+Inf)

    # computing and summing over the individual log-likelihood
    out <- tmp$x %>% log() %>% sum()

    # rescaling for nlminb
    obj <- -1 * out

    # rescale to make nlminb happy
    obj <- obj / (nrow(Data)/intPoints)

    #cat("obj = ", obj, "\n")

    obj
  }



  ceff_loglik_mixture_group <- function(x){
    # log-likelihood for group sizes

    # stochastic group size
    expkappas <- exp(x)
    n.cell <- table(Data$cell)/intPoints

    # Poisson group model
    obj_group <- sum(-expkappas + n.cell*log(expkappas) - lgamma(n.cell + 1))
    obj <- - obj_group/(nrow(Data)/intPoints)

    #cat("obj = ", obj, "\n")

    obj
  }


  ## Optimize log-likelihood for non-group part
  if (!silent){
    cat("Fitting the model...")
    time_start <- Sys.time()
  }
  time_start <- Sys.time()
  out <- nlminb(start = x.start,
                objective = ceff_loglik_mixture,
                control = list(rel.tol = 1e-6)
  )
  # round(out$par, 3)
  #[1]  1.329  1.300  1.476  1.386  1.419 -0.132 -0.355  1.368  1.290  1.921 -0.057  1.307 -0.017  0.101  1.645  2.072
  #[17]  1.956 -0.110  1.225 -0.033  0.103  1.354  1.133  1.974 -0.052  1.308 -0.014  0.097  1.512  1.864  1.878 -0.038
  #[33]  1.261 -0.014  0.101  1.623  2.018  2.112 -0.135  1.431 -0.050  0.099  1.673  1.894  1.781  0.004  1.400 -0.021
  #[49]  0.095  1.444  1.789  1.960 -0.084  1.305  0.011  0.094  1.577  1.593  1.842 -0.037  1.266 -0.011  0.102

  ## Optimize log-likelihood for group part
  out_group <- nlminb(start = x.start.kappa,
                      objective = ceff_loglik_mixture_group,
                      control = list(rel.tol = 1e-6)
  )
  # round(out_group$par, 3)
  #4.710 5.746 4.605 5.838 4.605 5.823 4.605 5.855

  par_final <- c(out$par, out_group$par)
  if (!silent){
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff,1), "s\n")
  }


  if (!silent){
    cat("Computing standard errors...")
    time_start <- Sys.time()
  }
  ## Standard Errors for Model (via numerical derivation)
  # for non-group part
  H <- pracma::hessian(ceff_loglik_mixture, out$par)   # may take a while....
  I <- H   # Note: rescaling (-1* and /N) has already taken place in the loglik-function
  varcov <- solve(I)/N

  # for group part
  H_group <- pracma::hessian(ceff_loglik_mixture_group, out_group$par)
  I_group <- H_group
  varcov_group <- solve(I_group)/N

  # join both vcovs (parameters from group and non-group part are independent)
  vcov_final <- lavaan:::lav_matrix_bdiag(varcov, varcov_group)
  if (!silent){
    time_diff <-  Sys.time() - time_start
    units(time_diff) <- "secs"
    cat("done. Took:", round(time_diff,1), "s\n")
  }


  #### Effect Estimation ####
  cat("Computing effects...")

  computeCondEffect <- function(x){
    ## This function computes conditional treatment effects within the XK-groups
    ## i.e., only the integration part given any combination of x and k
    ## aggregation to ATE and CE is done in another function

    # XK-specific parameters
    xk.mat <- x[8:95] %>% matrix(ncol = 11, byrow = T)
    xk.mat <- cbind(xk.mat, c(0,1), c(0,0,1,1,2,2,3,3))
    colnames(xk.mat) <- c("mu1", "psi1", "mu21", "mu22", "psi21", "psi22", "cov1", "cov2", "a0", "a1", "a2", "k", "x")

    # stochastic group size
    kappas <- x[96:103]

    # re-arranging the parameters for  effect estimation
    gFuncs <- data.frame(xk.mat[,c("a0", "a1", "a2", "k", "x")])
    covFuncs <- data.frame(xk.mat[,c("mu1", "psi1", "mu21", "mu22", "psi21", "psi22", "cov1", "cov2", "x", "k")])
    conditions <- expand.grid(k=c(0,1),x=c(0,1,2,3),gx=c(1,2,3), eff =NA)

    ## Effect for treatments t in any combination of x and k
    for (i in 1:nrow(conditions)){
      k <- conditions$k[i]
      x <- conditions$x[i]
      gx <- conditions$gx[i]

      # parameters of g-function
      A1 <- gFuncs[gFuncs$x == gx & gFuncs$k == k, c("a0", "a1", "a2")]
      A0 <- gFuncs[gFuncs$x == 0 & gFuncs$k == k, c("a0", "a1", "a2")]

      # parameters of covariates
      covPar <- covFuncs[covFuncs$x == x & covFuncs$k == k, ]

      effectpart1 <- exp(A1$a0 + A1$a1*covPar$mu1 + A1$a2*covPar$mu21 +
                          .5*A1$a1^2*covPar$psi1 + .5*A1$a2^2*covPar$psi21 + .5*A1$a1*A1$a2*covPar$cov1) -
                      exp(A0$a0 + A0$a1*covPar$mu1 + A0$a2*covPar$mu21 +
                          .5*A0$a1^2*covPar$psi1 + .5*A0$a2^2*covPar$psi21 + .5*A0$a1*A1$a2*covPar$cov1)
      effectpart2 <- exp(A1$a0 + A1$a1*covPar$mu1 + A1$a2*covPar$mu22 +
                          .5*A1$a1^2*covPar$psi1 + .5*A1$a2^2*covPar$psi22 + .5*A1$a1*A1$a2*covPar$cov2) -
                      exp(A0$a0 + A0$a1*covPar$mu1 + A0$a2*covPar$mu22 +
                          .5*A0$a1^2*covPar$psi1 + .5*A0$a2^2*covPar$psi22 + .5*A0$a1*A1$a2*covPar$cov2)


      conditions$eff[i] <- (effectpart1 + effectpart2)/2
    }



    # computing the probabilities P(X=x,K=k)
    Pxk <- exp(kappas)/(sum(exp(kappas)))

    return(c(conditions$eff, Pxk))

  }
  # apply basic effect estimation function to estimated parameters
  ce <- computeCondEffect(par_final)

  # standard error for CondEffect function
  # computes variance of xk-conditional effects and group probabilities (transformations of model parameters)
  # delta method
  delta.method.ck <- function(func, x, varcov){
    func.jacobian <- lav_func_jacobian_complex(func, x)
    func.jacobian %*% varcov %*% t(func.jacobian)
  }
  vcov.ce <- delta.method.ck(computeCondEffect, par_final, varcov = vcov_final)
  vcov.ce <- lav_matrix_bdiag(vcov.ce[1:24,1:24], vcov.ce[25:32,25:32])
  ce.se <- sqrt(diag(vcov.ce))


  #################################################
  ####### Lavaan syntax for aggregated effects
  lsyn <- '
# Definition of free parameters re. conditional effects
## Note: These "free" parameters are plugged in from our own model estimation
g1 ~ CE100*xk00 + CE101*xk01 + CE110*xk10 + CE111*xk11 + CE120*xk20 + CE121*xk21 + CE130*xk30 + CE131*xk31
g2 ~ CE200*xk00 + CE201*xk01 + CE210*xk10 + CE211*xk11 + CE220*xk20 + CE221*xk21 + CE230*xk30 + CE231*xk31
g3 ~ CE300*xk00 + CE301*xk01 + CE310*xk10 + CE311*xk11 + CE320*xk20 + CE321*xk21 + CE330*xk30 + CE331*xk31
p1 ~ relfreq00*xk00 + relfreq01*xk01 + relfreq10*xk10 + relfreq11*xk11 + relfreq20*xk20 + relfreq21*xk21 + relfreq30*xk30 + relfreq31*xk31


## Unconditional Probabilities P(K=k)
Pk0 := relfreq00 + relfreq10 + relfreq20 + relfreq30
Pk1 := relfreq01 + relfreq11 + relfreq21 + relfreq31

## Unconditional Probabilities P(X=x)
Px0 := relfreq00 + relfreq01
Px1 := relfreq10 + relfreq11
Px2 := relfreq20 + relfreq21
Px3 := relfreq30 + relfreq31

## Conditional Probabilities P(K=k|X=x)
Pk0gx0 := relfreq00/Px0
Pk1gx0 := relfreq01/Px0
Pk0gx1 := relfreq10/Px1
Pk1gx1 := relfreq11/Px1
Pk0gx2 := relfreq20/Px2
Pk1gx2 := relfreq21/Px2
Pk0gx3 := relfreq30/Px3
Pk1gx3 := relfreq31/Px3

## Conditional Probabilities P(X=x|K=k)
Px0gk0 := relfreq00/Pk0
Px1gk0 := relfreq10/Pk0
Px2gk0 := relfreq20/Pk0
Px3gk0 := relfreq30/Pk0
Px0gk1 := relfreq01/Pk1
Px1gk1 := relfreq11/Pk1
Px2gk1 := relfreq21/Pk1
Px3gk1 := relfreq31/Pk1

# Average treatment effects
ave10 := CE100*relfreq00 + CE101*relfreq01 + CE110*relfreq10 + CE111*relfreq11 + CE120*relfreq20 + CE121*relfreq21 + CE130*relfreq30 + CE131*relfreq31
ave20 := CE200*relfreq00 + CE201*relfreq01 + CE210*relfreq10 + CE211*relfreq11 + CE220*relfreq20 + CE221*relfreq21 + CE230*relfreq30 + CE231*relfreq31
ave30 := CE300*relfreq00 + CE301*relfreq01 + CE310*relfreq10 + CE311*relfreq11 + CE320*relfreq20 + CE321*relfreq21 + CE330*relfreq30 + CE331*relfreq31

# K=k conditional effects
CE10gk0 := Px0gk0*CE100 + Px1gk0*CE110 + Px2gk0*CE120 + Px3gk0*CE130
CE20gk0 := Px0gk0*CE200 + Px1gk0*CE210 + Px2gk0*CE220 + Px3gk0*CE230
CE30gk0 := Px0gk0*CE300 + Px1gk0*CE310 + Px2gk0*CE320 + Px3gk0*CE330
CE10gk1 := Px0gk1*CE101 + Px1gk1*CE111 + Px2gk1*CE121 + Px3gk1*CE131
CE20gk1 := Px0gk1*CE201 + Px1gk1*CE211 + Px2gk1*CE221 + Px3gk1*CE231
CE30gk1 := Px0gk1*CE301 + Px1gk1*CE311 + Px2gk1*CE321 + Px3gk1*CE331

# X=x conditional effects
CE10gx0 := Pk0gx0*CE100 + Pk1gx0*CE101
CE20gx0 := Pk0gx0*CE200 + Pk1gx0*CE201
CE30gx0 := Pk0gx0*CE300 + Pk1gx0*CE301

CE10gx1 := Pk0gx1*CE110 + Pk1gx1*CE111
CE20gx1 := Pk0gx1*CE210 + Pk1gx1*CE211
CE30gx1 := Pk0gx1*CE310 + Pk1gx1*CE311

CE10gx2 := Pk0gx2*CE120 + Pk1gx1*CE121
CE20gx2 := Pk0gx2*CE220 + Pk1gx1*CE221
CE30gx2 := Pk0gx2*CE320 + Pk1gx1*CE321

CE10gx3 := Pk0gx3*CE130 + Pk1gx3*CE131
CE20gx3 := Pk0gx3*CE230 + Pk1gx3*CE231
CE30gx3 := Pk0gx3*CE330 + Pk1gx3*CE331
'
  # laavanify aggregation syntax (create partable)
  ## computes aggregated effect and corresponding standard errors
  pt <- lavaanify(lsyn)
  def.function <- lav_partable_constraints_def(pt)
  JAC <- lav_func_jacobian_complex(func=def.function, x=ce)
  info.r <- JAC %*% vcov.ce %*% t(JAC)

  # Extract the point estimates and standard errors
  est <- def.function(ce)
  se <- sqrt(diag(info.r))

  sdyx0 <- sd(dNA$dv[dNA$treat==0])

  ## show effects, standard errors and effect sizes
  tval <- est/se
  n_par <- length(par_final)
  rdf <- nrow(object@input@data) - n_par

  pval <- 2*(1-pt(abs(tval),df=rdf))

  ce.tval <- ce/ce.se
  ce.pval <- 2*(1-pt(abs(ce.tval),df=rdf))

  ## Average effects
  selector <- c(23:25)
  Egx <- data.frame(est[selector],
                    se[selector],
                    tval[selector],
                    pval[selector],
                    est[selector]/sdyx0)

  names(Egx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given a treatment condition
  selector <- c(32:43)
  Egxgx <- data.frame(est[selector],
                      se[selector],
                      tval[selector],
                      pval[selector],
                      est[selector]/sdyx0)
  names(Egxgx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given K=k
  selector <- c(26:31)
  Egxgk <- data.frame(est[selector],
                      se[selector],
                      tval[selector],
                      pval[selector],
                      est[selector]/sdyx0)
  names(Egxgk) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given (X=x, K=k)
  selector <- c(1:(length(ce)-8))
  Egxgxk <- data.frame(ce[selector],
                       ce.se[selector],
                       ce.tval[selector],
                       ce.pval[selector],
                       ce[selector]/sdyx0)
  names(Egxgxk) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  object@results <- new("results",
                        est=est,
                        se=se,
                        vcov_def=info.r,
                        Egx=Egx,
                        Egxgx=Egxgx,
                        Egxgk=Egxgk,
                        Egxgxk=Egxgxk
  )
  cat("done\n")

  cat("Finished.\n")
  return(object)
}







