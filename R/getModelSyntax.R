getModelSyntax <- function(z, nz){
  tmp <- '## Covariate Moments\n'
  for (i in 1:nz){
    tmp <- paste0(tmp, z[i]," ~ ","c(mz00",i,", mz10",i,")*1\n")
    tmp <- paste0(tmp, z[i]," ~~ ","c(vz00",i,", vz10",i,")*",z[i],"\n")
  }

  if (nz > 1){
    for (i in 1:(nz-1)){
      for (j in (i+1):nz){
        tmp <- paste0(tmp, z[i]," ~~ ", "c(cvz00",i,j,", cvz10",i,j,")*",z[j],"\n")
      }
    }
  }


  tmp <- paste0(tmp,"\n","## Relative Group Frequencies\n")
  tmp <- paste0(tmp, "group"," % ","c(gw0, gw1)*w\n")
  tmp <- paste0(tmp, "N := exp(gw0) + exp(gw1)\n")
  tmp <- paste0(tmp, "relfreq0 := exp(gw0)/N\n")
  tmp <- paste0(tmp, "relfreq1 := exp(gw1)/N\n")

  tmp <- paste0(tmp, "\n", "## Defnition of Required Parameters\n")

  for (i in 1:nz){
    tmp <- paste0(tmp, "Ez",i," := mz00",i,"*relfreq0 + mz10",i,"*relfreq1\n")
    tmp <- paste0(tmp, "Vz",i," := vz00",i,"*relfreq0 + vz10",i,"*relfreq1 + relfreq0*(mz00",i,"-Ez",i,")^2 + relfreq1*(mz10",i,"-Ez",i,")^2\n")
  }

  tmp <- paste0(tmp, "Px0 := relfreq0\n")
  tmp <- paste0(tmp, "Px1 := relfreq1\n")
  tmp <- paste0(tmp, "Pk0gx0 := relfreq0/Px0\n")
  tmp <- paste0(tmp, "Pk0gx1 := relfreq1/Px1\n")

  for (i in 1:nz){
    tmp <- paste0(tmp, "Ez",i,"gx0 := mz00",i,"*Pk0gx0\n")
    tmp <- paste0(tmp, "Ez",i,"gx1 := mz10",i,"*Pk0gx1\n")
    tmp <- paste0(tmp, "Vz",i,"gx0 := vz00",i,"*Pk0gx0\n")
    tmp <- paste0(tmp, "Vz",i,"gx1 := vz10",i,"*Pk0gx1\n")
  }
  if (nz > 1){
    for (i in 1:(nz-1)){
      for (j in (i+1):nz){
        tmp <- paste0(tmp, "CVz",i,"z",j, "gx0 := ", "cvz00",i,j,"*Pk0gx0\n")
        tmp <- paste0(tmp, "CVz",i,"z",j, "gx1 := ", "cvz10",i,j,"*Pk0gx1\n")
      }
    }
  }


  return(tmp)
}
