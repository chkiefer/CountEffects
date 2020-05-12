
#' @export
setMethod("summary", signature(object = "countEffects"),
          function(object, ...) {
            if (object@input@distribution == "condNormal"){
              show(object)
            } else {
              cat("Average Treatment Effect\n")
              cat("--------------------------------------------------\n")
              ate <- unlist(object@ate)
              cat("Estimate:", ate[1],"\n")
              cat("SE:", ate[2], "\n\n")
              cat("Goodness-of-Fit\n")
              cat("--------------------------------------------------\n")
              pseudoR <- 1-object@results@glmresults$deviance/object@results@glmresults$null.deviance
              cat("Pseudo-R-squared (GLM):", pseudoR, "\n")
            }
          }
)


setMethod("show", "countEffects", function(object){

  ng <- object@input@ng
  nk <- object@input@nk
  nz <- object@input@nz
  vnames <- object@input@vnames
  vlevels <- object@input@vlevels
  gammas <- object@parnames@alphas    # ! Caution
  gammalabels <- object@parnames@gammalabels
  label.g.function <- object@parnames@label.g.function
  label.covs <- object@parnames@label.covs
  label.Egx <- object@parnames@label.Egx


  cat("\n\n--------------------- Variables  --------------------- \n\n")
  cat("Outcome variable Y: ", paste0(vnames$y), "\n")
  cat("Treatment variable X: ", paste0(vnames$x), "  (Reference group: ",
      paste0(object@input@control, ")\n"))
  if(!is.null(vnames$k)){
    cat("Categorical covariates K: ", paste0(vnames$k), "\n")
  }
  if(!is.null(vnames$z)){
    tmp <- "Continuous covariates in Z=("
    tmp <- paste0(tmp, paste0("Z",1:nz, collapse=","), "): ")
    tmp <- paste0(tmp, paste0(paste0("Z", 1:nz, "="), vnames$z, collapse=" "))
    cat(tmp, "\n")
  }
  cat("\n\n --------------------- Regression Model --------------------- \n")

  tmp <- paste0("E(Y|X",label.covs,") = ")
  tmp <- paste0(tmp, "g0",label.g.function," + ")
  tmp <- paste0(tmp, paste0("g",1:(ng-1),label.g.function,"*I_X=",1:(ng-1),
                            collapse=" + "))
  cat("\n",tmp, "\n")

  gammalabels2 <- gammalabels[,,1]
  gammalabels2[1] <- ""

  for(i in 1:ng){
    tmp <- paste0("  g",i-1,label.g.function," = exp(")
    tmp <- paste0(tmp, paste(gammas[,,i], gammalabels2, sep=" * ", collapse=" + "))
    tmp <- paste0(tmp, ")")
    tmp <- gsub("*  ", "", tmp, fixed=TRUE)

    ## Different Effect Function for Poisson Regression
    if (i > 1){
      tmp <- paste0(tmp," - exp(")
      tmp <- paste0(tmp, paste(gammas[,,1], gammalabels2, sep=" * ", collapse=" + "))
      tmp <- paste0(tmp, ")")
      tmp <- gsub("*  ", "", tmp, fixed=TRUE)
    }


    if(length(gammalabels2)==1){tmp <- gsub("*", "", tmp, fixed=TRUE)}

    if(nchar(tmp) > 80){
      ## split g function over several lines
      tmp <- unlist(strsplit(tmp, " + ", fixed=TRUE))
      tmp <- capture.output(cat(tmp, sep=" + ", fill=80))
      tmp[2:length(tmp)] <- paste0("            + ",tmp[2:length(tmp)])
      cat(tmp, sep="\n")
    } else{
      cat(tmp, "\n")
    }
  }
  # print coefficients of g-Functions
  for(i in 1:ng){
    if(i==1){
      tmp <- paste0("Intercept Function g",i-1,label.g.function)
      cat("\n",tmp, "\n\n")
    }else{
      tmp <- paste0("Effect Function g",i-1,label.g.function)
      tmp <- paste0(tmp, "   [", object@input@vnames$x,
                    ": ", object@input@vlevels$levels.x.original[i],
                    " vs. ",object@input@vlevels$levels.x.original[1], "]")
      cat("\n",tmp, "\n\n")
    }
    tmp <- object@results@gx[[i]]
    tmp[,2:5] <- round(tmp[,2:5], digits=3)
    print(tmp, print.gap=3, row.names=FALSE)
  }
  cat("\n\n--------------------- Cell Counts  --------------------- \n\n")


  if(nk>1){
    cat("\nCells \n")
    tmp <- expand.grid(K=vlevels$kstar, X=vlevels$levels.x.original)[,2:1]
    tmp$Cell <- vlevels$cell
    print(tmp, print.gap=3)

  }

  if(nk==1){
    cat("Cells \n")
    tmp <- data.frame(X=vlevels$levels.x.original)
    print(tmp, row.names=F, print.gap=3)

  }

  cat("\n")
  cat("Cell Counts \n\n")
  cat("This table shows cell counts including missings. \n")
  cat("See also output under lavaan results for number of observations \n")
  cat("actually used in the analysis. \n\n")

  if(nk==1){
    print(ftable(object@input@data[vnames$x]), print.gap=3)
  }else{
    cellcounts <- as.formula(paste0(paste(vnames$k, collapse="+"),
                                    "~", vnames$x))
    print(ftable(cellcounts, data=object@input@data), print.gap=3)
  }


  cat("\n\n --------------------- Adjusted Means --------------------- \n\n")
  namesadjmeans <- paste0("Adj.Mean",0:(ng-1))
  adjmeans <- object@results@adjmeans
  row.names(adjmeans) <- namesadjmeans
  print(adjmeans, digits=3, print.gap=3)


  cat("\n\n --------------------- Average Effects --------------------- \n\n")
  Egx <- object@results@Egx
  row.names(Egx) <- label.Egx
  print(Egx, digits=3, print.gap=3)

  if(!(nz==0 & nk==1)){
    cat("\n\n --------------------- Effects given a Treatment Condition --------------------- \n\n")
    tmp <- expand.grid(g=1:(ng-1), x=0:(ng-1))
    namesEgxgx <- paste0("E[g",tmp$g,label.g.function,"|X=",tmp$x, "]")
    Egxgx <- object@results@Egxgx
    row.names(Egxgx) <- namesEgxgx
    print(Egxgx, digits=3, print.gap=3)

  }

  }

)

