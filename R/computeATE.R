computeATE <- function(x, z, mod, data, distribution){
  nz <- length(z)

  coefs <- coef(mod)
  vcovs <- vcov(mod)

  pnames <- getCoefNames(nz)
  names(coefs) <- row.names(vcovs) <- colnames(vcovs) <- pnames

  modelz <- getModelSyntax(z, nz)
  mz <- lavaan::sem(modelz, data=data, group=x, group.label=c(0,1), group.w.free=TRUE)

  ## augment coefs and vcovs
  ## TODO: DIMENSIONEN ANPASSEN
  tmp <- 4*nz + nz*(nz - 1) + 2
  acoefs <- c(coefs, coef(mz, type="user")[-c(1:tmp)])
  avcovs <- lav_matrix_bdiag(vcovs, lavInspect(mz, "vcov.def", add.class = FALSE))
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
