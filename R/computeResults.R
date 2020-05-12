computeResults <- function(object){
  library(lavaan)

  if (object@input@distribution == "condNormal"){
    glmresults <- estimateGLM(object)
    semresults <- estimateModel(object)
  } else {
    glmresults <- estimateGLM(object)
    semresults <- NULL
  }



  allcoefs <- computeAdditionalCoefficients(object, glmresults, semresults)
  lavresults <- allcoefs$m1_sem
  est <- allcoefs$est
  se <- allcoefs$se
  vcov.def <- allcoefs$vcov.def ## vcov of defined parameters
  tval <- est/se

  if (!is.null(semresults)){
    rdf <- nrow(object@input@data) - semresults@npar
  } else {
    rdf <- glmresults$df.residual
  }

  pval <- 2*(1-pt(abs(tval),df=rdf))

  ng <- object@input@ng
  nz <- object@input@nz
  nk <- object@input@nk

  sdyx0 <- function(object){
    y <- object@input@data[, object@input@vnames$y]
    x <- object@input@data[, object@input@vnames$x]
    return(sd(y[x==0]))
  }

  sdyx0 <- sdyx0(object)



  Egx <- data.frame(est[object@parnames@Egx],
                    se[object@parnames@Egx],
                    tval[object@parnames@Egx],
                    pval[object@parnames@Egx],
                    est[object@parnames@Egx]/sdyx0)

  names(Egx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")

  ## Effects given a treatment condition
  Egxgx <- data.frame(est[object@parnames@Egxgx],
                      se[object@parnames@Egxgx],
                      tval[object@parnames@Egxgx],
                      pval[object@parnames@Egxgx],
                      est[object@parnames@Egxgx]/sdyx0)
  names(Egxgx) <- c("Estimate", "SE", "Est./SE", "p-value", "Effect Size")


  ## g functions
  gammas <- matrix(object@parnames@alphas, ncol=ng)
  gammalabels <- matrix(object@parnames@gammalabels, ncol=ng)
  gx <- vector("list",ng)

  for(i in 1:ng){
    tmp <- data.frame(gammas[,i],
                      est[gammas[,i]],
                      se[gammas[,i]],
                      tval[gammas[,i]],
                      pval[gammas[,i]])
    names(tmp) <- c("Coefficient", "Estimate", "SE", "Est./SE", "p-value")
    rownames(tmp) <- gammalabels[,i]
    gx[[i]] <- tmp
  }

  ## adjusted means
  adjmeans <- data.frame(est[object@parnames@adjmeans],
                         se[object@parnames@adjmeans],
                         tval[object@parnames@adjmeans])
  names(adjmeans) <- c("Estimate", "SE", "Est./SE")


  res <- new("results",
      lavresults=lavresults,
      glmresults=glmresults,
      semresults=new("nbmgsem"),
      est=est,
      se=se,
      vcov.def=vcov.def,
      Egx=Egx,
      Egxgx=Egxgx,
      gx=gx,
      adjmeans=adjmeans#,
      # AveEffZ="data.frame",
      # condeffects="data.frame"
      )



}



## additional parameters based on glm() results

computeAdditionalCoefficients <- function(object, glmresults, semresults){
  sem.args <- list(model=object@syntax@model,
                   group="cell",
                   fixed.x=object@input@fixed.z,
                   group.label=object@input@vlevels$cell,
                   data=object@input@data,
                   group.w.free = !object@input@fixed.cell)
  m1_sem <- do.call("sem", sem.args)

  nz <- object@input@nz

  if (!is.null(semresults)){
    coefs <- semresults@nlminb$par
    vcovs <- semresults@varcov
    pnames <- matrix(object@parnames@alphas, ncol = 1)
    names(coefs) <- rownames(vcovs) <- colnames(vcovs) <- pnames

  } else {
    m1_glm <- computeAlphas(object, glmresults)
    coefs <- m1_glm$coefs
    vcovs <- m1_glm$vcovs
    pnames <- names(coefs)
  }




  pt <- parTable(m1_sem)
  vcovs_tmp <- m1_sem@vcov$vcov
  row.names(vcovs_tmp) <- colnames(vcovs_tmp) <- m1_sem@ParTable$label[1:nrow(vcovs_tmp)]
  vcovs_tmp[pnames, pnames] <- vcovs

  coefs_tmp <- pt$est[1:nrow(vcovs_tmp)]
  names(coefs_tmp) <- m1_sem@ParTable$label[1:nrow(vcovs_tmp)]
  coefs_tmp[pnames] <- coefs

  ## compute effects
  def.function <- lav_partable_constraints_def(pt)
  JAC <- lav_func_jacobian_complex(func=def.function, x=coefs_tmp)
  info.r <- JAC %*% vcovs_tmp %*% t(JAC)

  se <- sqrt(diag(info.r))
  est <- def.function(coefs_tmp)
  vcov.def <- info.r

  row.names(vcov.def) <- colnames(vcov.def) <- names(se) <- names(est)

  est <- c(coefs, est)
  se <- c(sqrt(diag(vcovs)), se)


  return(list(est=est, se=se, vcov.def=vcov.def, m1_sem=m1_sem))

}

