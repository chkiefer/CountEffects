computeATE <- function(x, z, mod, data, distribution){
  nz <- length(z)

  coefs <- coef(mod)
  vcovs <- vcov(mod)

  pnames <- getCoefNames(nz)
  names(coefs) <- row.names(vcovs) <- colnames(vcovs) <- pnames



  if (distribution == "sn"){
    fml <- as.formula(paste0(z,"~1"))

    mz <- sn::selm(fml, data = data, method = "MPLE")
    acoefs <- c(coefs, mz@param$dp)

    if (is.null(mz@param.var$dp)){mz@param.var$dp <- diag(0,3)}
    avcovs <- lavaan::lav_matrix_bdiag(vcovs, mz@param.var$dp)
  } else if (distribution == "condSN"){
    d0 <- data[data[,names(data)==x] == 0,]
    d1 <- data[data[,names(data)==x] == 1,]

    mz0 <- sn::selm(z ~ 1, data = d0, method = "MPLE")
    if (is.null(mz0@param.var$dp)){mz0@param.var$dp <- diag(0,3)}
    mz1 <- sn::selm(z ~ 1, data = d1, method = "MPLE")
    if (is.null(mz1@param.var$dp)){mz1@param.var$dp <- diag(0,3)}

    gw <- '
    y ~ c(a, b)*1
    group % c(gw0, gw1)*w
    N := exp(gw0) + exp(gw1)
    relfreq0 := exp(gw0)/N
    relfreq1 := exp(gw1)/N
    Px0 := relfreq0
    Px1 := relfreq1
    '
    groupmodel <- lavaan::sem(gw, data=data, group = x,group.label=c(0,1), group.w.free=TRUE)

    acoefs <- c(coefs, mz0@param$dp, mz1@param$dp, lavaan::coef(groupmodel, type="user")[c(10:11)])
    avcovs <- lavaan::lav_matrix_bdiag(vcovs, mz0@param.var$dp, mz1@param.var$dp, lavaan::lavInspect(groupmodel, "vcov.def", add.class = FALSE)[4:5,4:5])
    row.names(avcovs) <- colnames(avcovs) <- names(acoefs) <- c(pnames, "xi0", "omega0", "alpha0", "xi1", "omega1", "alpha1", "Px0", "Px1")
  } else {
    modelz <- getModelSyntax(z, nz)
    mz <- lavaan::sem(modelz, data=data, group=x, group.label=c(0,1), group.w.free=TRUE)

    tmp <- 4*nz + nz*(nz - 1) + 2
    acoefs <- c(coefs, lavaan::coef(mz, type="user")[-c(1:tmp)])
    avcovs <- lavaan::lav_matrix_bdiag(vcovs, lavaan::lavInspect(mz, "vcov.def", add.class = FALSE))
  }

  ## augment coefs and vcovs
  ## TODO: DIMENSIONEN ANPASSEN

  row.names(avcovs) <- colnames(avcovs) <- names(acoefs)

  if (distribution == "normal"){
    Eg1 <- car::deltaMethod(acoefs,
                       "exp(g000+g100)*exp((g001+g101)*Ez1+((g001+g101)^2*Vz1/2))-exp(g000)*exp(g001*Ez1+(g001^2*Vz1/2))",
                       avcovs,
                       func="Eg1")
  } else if (distribution == "uniform"){
    Eg1 <- car::deltaMethod(acoefs,
                       "exp(g000+g100)*(exp((g001+g101)*(Ez1+sqrt(3*Vz1)))-exp((g001+g101)*(Ez1-sqrt(3*Vz1))))/((g001+g101)*(2*sqrt(3*Vz1)))-exp(g000)*(exp(g001*(Ez1+sqrt(3*Vz1)))-exp(g001*(Ez1-sqrt(3*Vz1))))/(g001*(2*sqrt(3*Vz1)))",
                       avcovs,
                       func="Eg1")
  } else if (distribution == "poisson"){
    Eg1 <- car::deltaMethod(acoefs,
                       "exp(g000+g100)*exp(Ez1*(exp(g001+g101)-1))-exp(g000)*exp(Ez1*(exp(g001)-1))",
                       avcovs,
                       func="Eg1")
  } else if (distribution == "chisquare"){
    Eg1 <- car::deltaMethod(acoefs,
                       "(exp(g000+g100)/(1-2*(g001+g101))^(Ez1/2))-(exp(g000)/(1-2*g001)^(Ez1/2))",
                       avcovs,
                       func="Eg1")
  } else if (distribution == "negbin"){
    Eg1 <- car::deltaMethod(acoefs,
                       "exp(g000+g100)*(Ez1/(Ez1+Vz1)*exp(g001+g101)/(1-(1-Ez1/(Ez1+Vz1))*exp(g001+g101)))^(Ez1^2/(Ez1+Vz1))-exp(g000)*(Ez1/(Ez1+Vz1)*exp(g001)/(1-(1-Ez1/(Ez1+Vz1))*exp(g001)))^(Ez1^2/(Ez1+Vz1))",
                       avcovs,
                       func="Eg1")
  } else if (distribution == "sn"){
    Eg1 <- car::deltaMethod(acoefs,
                            "2*pnorm(alpha/sqrt(1+alpha^2)*sqrt(omega)*(g001+g101))*
                            exp(g000+g100)*exp((g001+g101)*xi+((g001+g101)^2*omega/2))-
                            2*pnorm(alpha/sqrt(1+alpha^2)*sqrt(omega)*g001)*
                            exp(g000)*exp(g001*xi+(g001^2*omega/2))",
                            avcovs,
                            func="Eg1")
  } else if (distribution == "condNormal"){
    # parFunc <- "Px0*exp(g000+g100)*exp((g001+g101)*Ez1gx0+((g001+g101)^2*Vz1gx0/2))-Px0*exp(g000)*exp(g001*Ez1gx0+(g001^2*Vz1gx0/2))+
      #          Px1*exp(g000+g100)*exp((g001+g101)*Ez1gx1+((g001+g101)^2*Vz1gx1/2))-Px1*exp(g000)*exp(g001*Ez1gx1+(g001^2*Vz1gx1/2))"
    parFunc <- getParFunction(nz)
    Eg1 <- car::deltaMethod(acoefs,
                       parFunc$ParFunction,
                       avcovs,
                       func="Eg1")
    Eg1x0 <- car::deltaMethod(acoefs,
                         #"exp(g000+g100)*exp((g001+g101)*Ez1gx0+((g001+g101)^2*Vz1gx0/2))-exp(g000)*exp(g001*Ez1gx0+(g001^2*Vz1gx0/2))",
                         parFunc$ParFunctiongx0,
                         avcovs,
                         func="Eg1x0")
    Eg1x1 <- car::deltaMethod(acoefs,
                         #"exp(g000+g100)*exp((g001+g101)*Ez1gx1+((g001+g101)^2*Vz1gx1/2))-exp(g000)*exp(g001*Ez1gx1+(g001^2*Vz1gx1/2))",
                         parFunc$ParFunctiongx1,
                         avcovs,
                         func="Eg1x1")
  } else if (distribution == "condSN"){
    Eg1 <- car::deltaMethod(acoefs,
                            "Px0*2*pnorm(alpha0/sqrt(1+alpha0^2)*sqrt(omega0)*(g001+g101))*
                            exp(g000+g100)*exp((g001+g101)*xi0+((g001+g101)^2*omega0/2))-
                            Px0*2*pnorm(alpha0/sqrt(1+alpha0^2)*sqrt(omega0)*g001)*
                            exp(g000)*exp(g001*xi0+(g001^2*omega0/2))+
                            Px1*2*pnorm(alpha1/sqrt(1+alpha1^2)*sqrt(omega1)*(g001+g101))*
                            exp(g000+g100)*exp((g001+g101)*xi1+((g001+g101)^2*omega1/2))-
                            Px1*2*pnorm(alpha1/sqrt(1+alpha1^2)*sqrt(omega1)*g001)*
                            exp(g000)*exp(g001*xi1+(g001^2*omega1/2))",
                            avcovs,
                            func="Eg1")
    Eg1x0 <- car::deltaMethod(acoefs,
                              "2*pnorm(alpha0/sqrt(1+alpha0^2)*omega0*(g001+g101))*
                              exp(g000+g100)*exp((g001+g101)*xi0+((g001+g101)^2*omega0^2/2))-
                              2*pnorm(alpha0/sqrt(1+alpha0^2)*omega0*g001)*
                              exp(g000)*exp(g001*xi0+(g001^2*omega0^2/2))",
                              avcovs,
                              func="Eg1x0")
    Eg1x1 <- car::deltaMethod(acoefs,
                              "2*pnorm(alpha1/sqrt(1+alpha1^2)*omega1*(g001+g101))*
                              exp(g000+g100)*exp((g001+g101)*xi1+((g001+g101)^2*omega1^2/2))-
                              2*pnorm(alpha1/sqrt(1+alpha1^2)*omega1*g001)*
                              exp(g000)*exp(g001*xi1+(g001^2*omega1^2/2))",
                              avcovs,
                              func="Eg1x1")
  } else {stop("Distribution not supported!")}

  if (distribution != "condNormal"){
    res <- list(Ave=Eg1$Estimate,
                se_Ave=Eg1$SE)
  } else {
    res <- list(Ave=Eg1$Estimate,
                se_Ave=Eg1$SE,
                AveX0=Eg1x0$Estimate,
                se_AveX0=Eg1x0$SE,
                AveX1=Eg1x1$Estimate,
                se_AveX1=Eg1x1$SE)
  }


  return(res)

}
