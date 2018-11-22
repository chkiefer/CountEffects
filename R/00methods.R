
#' @export
setMethod("summary", signature(object = "countEffects"),
          function(object, ...) {
            cat("Average Treatment Effect\n")
            cat("--------------------------------------------------\n")
            ate <- unlist(object@ate)
            cat("Estimate:", ate[1],"\n")
            cat("SE:", ate[2], "\n\n")
            cat("Goodness-of-Fit\n")
            cat("--------------------------------------------------\n")
            pseudoR <- 1-object@model$deviance/object@model$null.deviance
            cat("Pseudo-R-squared (GLM):", pseudoR, "\n")
          }
)


setMethod("show", "countEffects",
          function(object){
            summary(object)
          }
)

