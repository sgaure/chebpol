.onAttach <- function(libname,pkgname) {
  if(!havefftw()) {
    packageStartupMessage("*** ",pkgname,": FFTW not used.\n*** You should install it from http://fftw.org\n*** or check if your OS-distribution provides it, and recompile.")
  }
}

.onLoad <- function(libname,pkgname) {
  if(is.na(thr <- as.integer(Sys.getenv('CHEBPOL_THREADS')))) thr <- 1L
  options(chebpol.threads=thr)
}
