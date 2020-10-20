ceff_estimate_creg <- function(object){
  input <- object@input

  forml <- input@forml
  measurement <- input@measurement
  if (!length(measurement)){
    measurement <- NULL
  }
  data <- input@data
  method <- input@method

  fit <- countreg(forml = forml,
           lv = measurement,
           group = "cell",
           data = data,
           family = method)
  return(fit)

}


