getParFunction <- function(nz){
  tmp <- "exp(g000+g100"

  for (i in 1:nz){
    tmp <- paste0(tmp,"+(g00",i,"+g10",i,")*Ez",i,"gx0")
    tmp <- paste0(tmp, "+(g00",i,"+g10",i,")^2*Vz",i,"gx0/2")
  }

  tmp <- paste0(tmp, ")")
  tmp <- paste0(tmp, "-")

  tmp <- paste0(tmp,"exp(g000")
  for (i in 1:nz) {
    tmp <- paste0(tmp,"+g00",i,"*Ez",i,"gx0")
    tmp <- paste0(tmp, "+g00",i,"^2*Vz",i,"gx0/2")
  }

  tmp <- paste0(tmp, ")")
  ParFunctiongx0 <- tmp

  tmp <- "exp(g000+g100"
  for (i in 1:nz){
    tmp <- paste0(tmp,"+(g00",i,"+g10",i,")*Ez",i,"gx1")
    tmp <- paste0(tmp, "+(g00",i,"+g10",i,")^2*Vz",i,"gx1/2")
  }

  tmp <- paste0(tmp, ")")
  tmp <- paste0(tmp, "-")

  tmp <- paste0(tmp,"exp(g000")
  for (i in 1:nz){
    tmp <- paste0(tmp,"+g00",i,"*Ez",i,"gx1")
    tmp <- paste0(tmp, "+g00",i,"^2*Vz",i,"gx1/2")
  }

  tmp <- paste0(tmp, ")")
  ParFunctiongx1 <- tmp

  tmp <- paste0("Px0*(", ParFunctiongx0, ")+","Px1*(", ParFunctiongx1, ")")
  ParFunction <- tmp

  res <- list(ParFunction = ParFunction,
              ParFunctiongx0 = ParFunctiongx0,
              ParFunctiongx1 = ParFunctiongx1)
  return(res)
}

