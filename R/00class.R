#' @export
setClass("input",
         representation(
           method="character",
           vnames="list", ## variable names
           vlevels="list", ## variable levels (for x, k, kstar and cell)
           control="character",
           ng="integer", ## number of treatment groups
           nz="integer", ## number of z
           nk="integer", ## number of unfolded categories of K
           data="data.frame",
           fixed.cell ="logical",
           fixed.z ="logical",
           observed.freq="numeric", ## observed group frequencies (fixed.cell only)
           sampmomentz="array" ## manifest sample means for continuous covariates
         )
)

#' @export
setClass("parnames", representation(
  alphas="array",
  betas="array",
  gammas="array",
  gammalabels="array",
  label.g.function="character",
  label.covs="character",
  label.Egx="character",
  cellmeanz="character",
  cellvarz="character",
  cellskewz="character",
  meanz="character",
  pk="character",
  px="character",
  Ezk="character",
  Pkgx="character", ## P(K=k|X=x)
  Pxgk="character", ## P(X=x|K=k)
  Ezgx="character", ## E(Z|X=x)
  Ezgk="character", ## E(Z|K=k)
  Ezkgx="character", ## E(Z*K|X=x)
  groupw="character",
  relfreq="character",
  Egx="character",
  Egxgx="character", ## E(gx|X=x)
  Egxgk="character", ## E(gx|K=k)
  Egxgxk="character", ## E(gx|X=x,K=k)
  adjmeans="character",
  AveEffZ="character" ## average effect continuous covariate Z
)
)

#' @export
setClass("syntax",
         representation(
           model="character",
           formula="formula"
         ))

#' @export
setClass("countEffects",
         representation(
           input              = "input",
           parnames           = "parnames",
           syntax             = "syntax",
           model              = "ANY",  # a glm model
           ate                = "list"
         )
)





