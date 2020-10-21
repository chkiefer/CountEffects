.onLoad <- function(libname, pkgname) {

  if (!requireNamespace("lavacreg", quietly = TRUE)) {
    stop("Package \"lavacreg\" needed for this package to work. Please install it.\n \"lavacreg\" is available via github.com/chkiefer/lavacreg",
      call. = FALSE)
  } else {
    require(lavacreg)
    require(lavaan)
  }
}
