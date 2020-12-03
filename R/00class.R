#' @export
setClass("input",
         representation(
           method="character",
           vnames="list", ## variable names
           vlevels="list", ## variable levels (for x, k, kstar and cell)
           ngroups="integer",
           control="character",
           ng="integer", ## number of treatment groups
           nz="integer", ## number of z
           nk="integer", ## number of unfolded categories of K
           data="data.frame",
           measurement="list",
           distribution="character",
           forml="formula"
         )
)

#' @export
setClass("parnames", representation(
  alphas="array",
  label.g.function="character",
  label.covs="character",
  label.Egx="character",
  cellmeanz="character",
  cellvarz="character",
  cellcovz="character",
  cellskewz="character",
  cellexpc="character",
  Eygx="character", ## E[E(Y|X=x,K,Z)]
  Eygxgx="character", ## E[E(Y|X=x,K,Z)|X=x0]
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
setClass("results", representation(
  est="numeric",
  se="numeric",
  vcov_def="matrix",
  Egx="data.frame",
  Egxgx="data.frame",
  Egxgk="data.frame",
  Egxgxk="data.frame"
)
)


#' @export
setClass("countEffects",
         representation(
           input              = "input",
           fit                = "lavacreg",
           parnames           = "parnames",
           syntax             = "syntax",
           results            = "results",
           ate                = "list"
         )
)





