ceff_create_syntax <- function(object){

  ## lavaan syntax
  model <- character()
  model <- createLavaanSyntax(object)

  ## formula for PMET case (only computed if necessary)
  formula <- formula()

  res <- new("syntax",
             model=model,
             formula=formula
  )
}



createLavaanSyntax <- function(obj) {

  inp <- obj@input
  parnames <- obj@parnames

  ## input information
  y <- inp@vnames$y
  z <- inp@vnames$z
  ng <- inp@ng
  nz <- inp@nz
  nk <- inp@nk
  fixed.cell <- inp@fixed.cell
  fixed.z <- inp@fixed.z
  sampmeanz <- inp@sampmeanz
  homoscedasticity <- inp@homoscedasticity
  observed.freq <- inp@observed.freq

  ## parnames
  alphas <- parnames@alphas
  cellmeanz <- parnames@cellmeanz
  cellvarz <- parnames@cellvarz
  cellcovz <- parnames@cellcovz
  relfreq <- parnames@relfreq
  groupw <- parnames@groupw
  cellexpc <- parnames@cellexpc
  Eygx <- parnames@Eygx
  Eygxgx <- parnames@Eygxgx
  meanz <- parnames@meanz
  pk <- parnames@pk
  px <- parnames@px
  Ezk <- parnames@Ezk
  Egx <- parnames@Egx
  adjmeans <- parnames@adjmeans
  Pkgx <- parnames@Pkgx
  Pxgk <- parnames@Pxgk
  Ezgx <- parnames@Ezgx
  Ezgk <- parnames@Ezgk
  Ezkgx <- parnames@Ezkgx
  Egxgx <- parnames@Egxgx
  Egxgk <- parnames@Egxgk
  Egxgxk <- parnames@Egxgxk
  AveEffZ <- parnames@AveEffZ

  model <- "#### lavaan Syntax for CountEffects Model ####"


  model <- paste0(model, "\n\n## Structural Model \n")

  ## syntax intercepts
  model <- paste0(model, create_syntax_intercepts(y,alphas))

  ## syntax regression coefficients in each cell
  model <- paste0(model, create_syntax_regcoef(y,z,nz,alphas))

  ## mean z in each cell
  model <- paste0(model, create_syntax_cellmeanz(z, nz, fixed.z, cellmeanz,
                                                 sampmeanz))

  ## variance z in each cell
  model <- paste0(model, create_syntax_cellvarz(z, nz, fixed.z, cellvarz,
                                                 sampvarz))

  ## covariances between stochastic z
  model <- paste0(model, create_syntax_covz(z, nz, nk, ng, fixed.z, cellcovz))

  ## homoscedastic residual variances
  model <- paste0(model, create_syntax_homoscedasticity(y,ng,nk,homoscedasticity))

  ## compute relative group frequencies
  model <- paste0(model, create_syntax_group_freq(fixed.cell, relfreq,
                                                  observed.freq, groupw))

  ## compute cell-conditional expectation given cell E[E(Y|X=x,K=k,Z)|X=x,K=k]
  model <- paste0(model, create_syntax_cellexpc(nz, ng, nk, cellexpc, cellmeanz, cellvarz, cellcovz, alphas))

  ## compute cell-conditional expectation E[E(Y|X=x,K=k,Z)]
  model <- paste0(model, create_syntax_eygx(nz, ng, nk, Eygx, cellexpc, relfreq))

  ## compute unconditional means of z
  model <- paste0(model, create_syntax_ez(nz, meanz, cellmeanz, relfreq))

  ## compute unconditional probabilities of K*=k
  model <- paste0(model, create_syntax_pk(nk, pk, relfreq))

  ## compute unconditional probabilities of X=x
  model <- paste0(model, create_syntax_ex(px, ng, nk, relfreq))

  ## compute average effects
  model <- paste0(model, create_syntax_Egx(ng, Egx, Eygx))

  ## compute adjusted means
  model <- paste0(model, create_syntax_adjmeans(ng, adjmeans, Eygx))

  ## conditional probabilities of K=k given X=x (Pkgx)
  model <- paste0(model, create_syntax_Pkgx(ng, nk, relfreq, Pkgx, px))

  ## conditional probabilities of X=x given K=k (Pxgk)
  model <- paste0(model, create_syntax_Pxgk(ng, nk, relfreq, Pxgk, pk))

  ## compute cell-conditional expectation E[E(Y|X=x,K=k,Z)|X=x0]
  model <- paste0(model, create_syntax_eygxgx(nz, ng, nk, Eygxgx, cellexpc, Pkgx))

  ## conditional expectations of Z given X=x (Ezgx)
  model <- paste0(model, create_syntax_Ezgx(ng, nk, nz, Ezgx, Pkgx, cellmeanz))

  ## effects given a treatment condition
  model <- paste0(model, create_syntax_Egxgx(ng,Egxgx,Eygxgx))




  return(model)
}




## functions that create parts of the lavaan/methods syntax


create_syntax_intercepts <- function(y, alphas){

  res <- paste0(y," ~ c(", paste(alphas[1,,],collapse=","), ")*1")
  return(res)

}

create_syntax_regcoef <- function(y, z, nz, alphas){

  res <- NULL
  if (nz>0) {
    for (i in 1:nz) {
      tmp <- paste0(y," ~ c(", paste(alphas[i+1,,],collapse=","), ")*",z[i])
      res <- paste0(res, "\n", tmp)
    }
  }
  return(res)

}

create_syntax_cellmeanz <- function(z, nz, fixed.z, cellmeanz, sampmeanz){

  res <- NULL
  if(!fixed.z){
    ## stochastic z
    if (nz>0) {
      cellmeanz <- matrix(cellmeanz, nrow=nz)
      for (i in 1:nz) {
        tmp <- paste0(z[i]," ~ c(", paste(cellmeanz[i,],collapse=","), ")*1")
        res <- paste0(res, "\n", tmp)
      }
    }

  }else if(fixed.z){
    ## fixed z
    if (nz>0) {
      res <- paste0(res, "\n\n## Fixed Means of Z")
      cellmeanz <- matrix(cellmeanz, nrow=nz)
      tmp <- paste0(cellmeanz, " := ",  sampmeanz, collapse="\n")
      res <- paste0(res, "\n", tmp)
    }
  }

  return(res)
}

create_syntax_cellvarz <- function(z, nz, fixed.z, cellvarz, sampvarz){

  res <- NULL
  if(!fixed.z){
    ## stochastic z
    if (nz>0) {
      cellvarz <- matrix(cellvarz, nrow=nz)
      for (i in 1:nz) {
        tmp <- paste0(z[i]," ~~ c(", paste(cellvarz[i,],collapse=","), ")*",z[i])
        res <- paste0(res, "\n", tmp)
      }
    }

  }else if(fixed.z){
    ## fixed z
    if (nz>0) {
      res <- paste0(res, "\n\n## Fixed Means of Z")
      cellvarz <- matrix(cellvarz, nrow=nz)
      tmp <- paste0(cellmeanz, " := ",  sampvarz, collapse="\n")
      res <- paste0(res, "\n", tmp)
    }
  }

  return(res)
}


create_syntax_covz <- function(z, nz, nk, ng, fixed.z, cellcovz){

  if (nz <= 1){return(character())}

  nzz <- nz*(nz-1)/2
  cellcovz <- array(cellcovz, dim=c(nk, ng, nzz))
  zz <- t(combn(z, 2))

  ## covariances between stochastic z
  res <- NULL
  if(!fixed.z){
    ## stochastic z
    if (nz>0) {
      for (i in 1:nzz) {
        tmp <- paste0(zz[i,1]," ~~ c(", paste(cellcovz[, ,i],collapse=","), ")*",zz[i,2])
        res <- paste0(res, "\n", tmp)
      }
    }

  return(res)
  }
}

create_syntax_homoscedasticity <- function(y, ng, nk, homoscedasticity){

  res <- NULL
  if(homoscedasticity){
    tmp <- paste0(y, " ~~ c(",
                  paste(rep("veps", times=ng*nk),collapse=","),
                  ")*", y)
    res <- paste0(res, "\n", tmp)
  }

  return(res)
}

create_syntax_group_freq <- function(fixed.cell, relfreq, observed.freq, groupw){


  res <- "\n\n## Relative Group Frequencies \n"

  if(fixed.cell){
    tmp <- paste(paste0(relfreq, " := ", observed.freq), collapse="\n")
    res <- paste0(res, tmp)

  }else if(!fixed.cell){

    ## syntax group weights
    tmp <- paste0("group % c(", paste(groupw, collapse=","), ")*w")
    res <- paste0(res, tmp)

    tmp <- paste(paste0("exp(", groupw, ")"), collapse=" + ")
    tmp <- paste0("N := ",tmp)
    res <- paste0(res, "\n", tmp)

    tmp <- paste(paste0(relfreq, " := exp(", groupw, ")/N"), collapse="\n")
    res <- paste0(res, "\n", tmp)
  }

  return(res)
}


create_syntax_cellexpc <- function(nz, ng, nk, cellexpc, cellmeanz, cellvarz, cellcovz, alphas){

  res <- NULL
  res <- paste0(res, "\n\n## Cell-conditional Expectations given cell E[E(Y|X=x, K=k, Z)|X=x0,K=k]")
  cellexpc <- array(cellexpc, dim = c(ng, nk, ng))
  alphas <- array(alphas, dim = c((nz+1), nk, ng))

  if (nz>0){
    nzz <- nz*(nz-1)/2

    cellmeanz <- array(cellmeanz, dim=c(ng, nk, nz))
    cellvarz <- array(cellvarz, dim=c(ng, nk, nz))
    cellcovz <- array(cellcovz, dim=c(nk, ng, nzz))

    tmp0 <- expand.grid(g=1:(ng-1), x=0:(ng-1), k=0:(nk-1), gx=0:(ng-1))

    for (j in 1:nrow(tmp0)){
      x <- tmp0$x[j] + 1
      k <- tmp0$k[j] + 1
      gx <- tmp0$gx[j] + 1

      tmp <- paste0(cellexpc[x, k, gx]," := exp(", alphas[1,k,x])
        for (i in 1:nz) {
          tmp <- paste0(tmp, paste(" + ", cellmeanz[gx,k,i], "*", alphas[(i+1),k,x]))
          tmp <- paste0(tmp, paste(" + .5*", cellvarz[gx,k,i], "*", alphas[(i+1),k,x],"^2"))
        }
        if(nzz){
          counter <- 1L
          for (i in 2:nz){
            for (j in 1:(i-1)){
              tmp <- paste0(tmp, paste(" + .5*", cellcovz[k, x, counter], "*", alphas[(i+1),k,x], "*", alphas[(j+1),k,x]))
              counter <- counter + 1L
            }
          }
        }
        tmp <- paste0(tmp, ")")
        res <- paste0(res, "\n", tmp)
      }
  } else {

  }


  #cat(res)
  return(res)
}


create_syntax_eygx <- function(nz, ng, nk, Eygx, cellexpc, relfreq){
 ### TODO anpassen
  res <- NULL
  res <- paste0(res, "\n\n## Cell-conditional Expectations E[E(Y|X=x, K, Z)]")
  cellexpc <- array(cellexpc, dim = c(ng, nk, ng))
  Eygx <- array(Eygx, dim = c(ng))
  relfreq <- array(relfreq, dim = c(ng, nk))

  for (x in 1:ng){
      tmp <- paste0(Eygx[x], " := ", paste(cellexpc[x, , ], relfreq, sep="*", collapse = "+"))
      res <- paste0(res, "\n", tmp)
  }

  #cat(res)
  return(res)
}

create_syntax_eygxgx <- function(nz, ng, nk, Eygxgx, cellexpc, Pkgx){
  ### TODO anpassen
  res <- NULL
  res <- paste0(res, "\n\n## Cell-conditional Expectations E[E(Y|X=x, K, Z)|X=x0]")
  cellexpc <- array(cellexpc, dim = c(ng, nk, ng))
  Eygxgx <- array(Eygxgx, dim = c(ng, ng))
  Pkgx <- array(Pkgx, dim = c(nk))

  for (x in 1:ng){
    for (x0 in 1:ng){
      tmp <- paste0(Eygxgx[x,x0], " := ", paste(cellexpc[x, ,x0], Pkgx, sep="*", collapse = "+"))
      res <- paste0(res, "\n", tmp)
    }
  }

  #cat(res)
  return(res)
}

create_syntax_ez <- function(nz, meanz, cellmeanz, relfreq){

  res <- NULL
  if (nz>0) {
    res <- paste0(res, "\n\n## Unconditional Expectations E(Z)")
    cellmeanz <- matrix(cellmeanz, nrow=nz)
    for (i in 1:nz) {
      tmp <- paste0(meanz[i]," := ", paste(cellmeanz[i,], relfreq, sep="*", collapse=" + "))
      res <- paste0(res, "\n", tmp)
    }
  }

  return(res)
}


create_syntax_pk <- function(nk, pk, relfreq){

  res <- NULL
  if (nk>1) {
    res <- paste0(res, "\n\n## Unconditional Probabilities P(K=k)")
    relfreq <- matrix(relfreq, nrow=nk)
    for (i in 1:nk) {
      tmp <- paste0(pk[i], " := ", paste(relfreq[i,], collapse=" + "))
      res <- paste0(res, "\n", tmp)
    }
  }

  return(res)
}


create_syntax_ex <- function(px, ng, nk, relfreq){

  res <- "\n\n## Unconditional Probabilities P(X=x)"

  relfreq <- matrix(relfreq, nrow=nk)
  for (i in 1:ng) {
    tmp <- paste0(px[i], " := ", paste(relfreq[,i], collapse=" + "))
    res <- paste0(res, "\n", tmp)
  }

  return(res)
}

create_syntax_Egx <- function(ng, Egx, Eygx){

  Eygx <- array(Eygx, dim = c(ng))
  Egx <- array(Egx, dim = c(ng))

  res <- "\n\n## Average Effects"


  ## average total effects
  for(i in 2:ng){
    tmp <- paste0(Egx[i-1]," := ",
                  paste(Eygx[i], "-", Eygx[1]))
    res <- paste0(res, "\n", tmp)
  }
  #cat(res)
  return(res)
}

create_syntax_adjmeans <- function(ng, adjmeans, Eygx){

  ##TODO compute adjusted means based on gammas so that it can be used in
  ## single group models as well

  res <- "\n\n## Adjusted Means"


  ## adjusted means
  tmp <- paste0(adjmeans[1]," := ",
                paste(Eygx[1]))
  res <- paste0(res, "\n", tmp)

  for(i in 2:ng){
    tmp <- paste(Eygx[i])
    tmp <- paste0(adjmeans[i]," := ", tmp)
    res <- paste0(res, "\n", tmp)
  }
  #cat(res)
  return(res)
}


create_syntax_Pkgx <- function(ng, nk, relfreq, Pkgx, px){

  res <- "\n\n## Conditional Probabilities P(K=k|X=x)"

  relfreq <- matrix(relfreq, nrow=nk)
  Pkgx <- matrix(Pkgx, nrow=nk)

  for(i in 1:ng){
    for(k in 1:nk){
      Pkgx[k,i] <- paste0(Pkgx[k,i], " := ", relfreq[k,i], "/", px[i])
    }
  }
  res <- paste0(res, "\n", paste(Pkgx, collapse="\n"))

  return(res)
}


create_syntax_Pxgk <- function(ng, nk, relfreq, Pxgk, pk){

  res <- NULL
  if(nk>1){
    res <- paste0(res, "\n\n## Conditional Probabilities P(X=x|K=k)")
    relfreq <- matrix(relfreq, nrow=nk)
    Pxgk <- matrix(Pxgk, nrow=ng)

    for(i in 1:ng){
      for(k in 1:nk){
        Pxgk[i,k] <- paste0(Pxgk[i,k], " := ", relfreq[k,i], "/", pk[k])
      }
    }
    res <- paste0(res, "\n", paste(Pxgk, collapse="\n"))
  }

  return(res)
}



create_syntax_Ezgx <- function(ng,nk,nz,Ezgx,Pkgx,cellmeanz){

  res <- NULL

  if(nz!=0){
    res <- paste0(res, "\n\n## Conditional Expectations E(Z|X=x)")
    cellmeanz <- array(cellmeanz, dim=c(nz,nk,ng))
    Ezgx <- matrix(Ezgx, nrow=nz)
    Pkgx <- matrix(Pkgx, nrow=nk)
    for(i in 1:ng){
      for(q in 1:nz){
        Ezgx[q,i] <- paste0(Ezgx[q,i], " := ",
                            paste(cellmeanz[q,1:nk,i], "1", sep="*", collapse=" + "))
      }
    }
    res <- paste0(res, "\n", paste(Ezgx, collapse="\n"))
  }

  return(res)
}

create_syntax_Egxgx <- function(ng,Egxgx,Eygxgx){

  #TODO: ist noch nicht richtig
  ## Effects given a treatment condition
  res <-  "\n\n## Effects given X=x"

  Eygxgx <- array(Eygxgx, dim = c(ng, ng))
  Egxgx <- matrix(Egxgx, ncol=ng)

  #
  # ## average total effects
  # for(i in 2:ng){
  #   tmp <- paste0(Egx[i-1]," := ",
  #                 paste(Eygx[i], "-", Eygx[1]))
  #   res <- paste0(res, "\n", tmp)
  # }
  # #cat(res)
  # return(res)

  for(gx in 1:ng){
    for(i in 2:ng){
      tmp <- paste0(Egxgx[i-1,gx], " := ",
                    paste(Eygxgx[i,gx],"-",Eygxgx[1,gx]))
      res <- paste0(res, "\n", tmp)

    }
  }
 #cat(res)
  return(res)





}






