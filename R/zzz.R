.onLoad <- function(libname, pkgname) {

  if (!requireNamespace("CountReg", quietly = TRUE)) {
    stop("Package \"CountReg\" needed for this package to work. Please install it.\n \"CountReg\" is available via github.com/chkiefer/CountReg",
      call. = FALSE)
  } else {
    require(CountReg)
    require(lavaan)
  }
}
