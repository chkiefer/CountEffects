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
           measurement="character",
           distribution="character",
           fixed.cell ="logical",
           fixed.z ="logical",
           observed.freq="numeric", ## observed group frequencies (fixed.cell only)
           sampmeanz="array", ## manifest sample means for continuous covariates
           sampvarz="array",
           sampcovz="array",
           homoscedasticity="logical"   ### only relevant for NB regression - equal overdispersion in groups
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
setClass("nbmgsem", representation = representation(
  ngroups = "integer",
  nz = "integer",
  npar = "integer",
  vnames = "list",
  method = "character",
  distribution = "character",
  data = "data.frame",
  par = "numeric",
  nlminb = "list",
  varcov = "matrix"
)
)


#' @export
setClass("results", representation(
  semresults="nbmgsem",
  lavresults="lavaan",
  glmresults="ANY", ## glm or negbin
  est="numeric",
  se="numeric",
  vcov.def="matrix",
  hypotheses="data.frame",
  hypothesesk="data.frame",
  Egx="data.frame",
  AdditionalEffects="data.frame",
  Egxgx="data.frame",
  Egxgk="data.frame",
  Egxgxk="data.frame",
  gx="list",
  adjmeans="data.frame",
  AveEffZ="data.frame",
  condeffects="data.frame"
)
)


#' @export
setClass("countEffects",
         representation(
           input              = "input",
           parnames           = "parnames",
           syntax             = "syntax",
           results            = "results",
           ate                = "list"
         )
)





