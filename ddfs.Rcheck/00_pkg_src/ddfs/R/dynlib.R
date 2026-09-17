#' @useDynLib ddfs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils packageDescription
NULL

.onAttach <- function(libname, pkgname) {
  #nice basic info
  vers <- packageDescription("ddfs")[["Version"]]
  packageStartupMessage("##################################################################\n",
                        "\n",
                        "This is ddfs version ", vers, ". \n",
                        "See ",
                        sQuote("package?ddfs"), " for a brief introduction.\n",
                        "Type ", sQuote("citation('ddfs')"), " to learn how to cite this package.\n",
                        "\n",
                        "Please report any issue or bug to the authors (See the description\n",
                        "file)\n",
                        "\n",
                        "##################################################################\n",
                        appendLF = TRUE)
  return(TRUE)
}

