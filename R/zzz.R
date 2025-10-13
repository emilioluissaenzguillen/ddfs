# #' @importFrom reticulate py_module_available py_install import
# .onLoad <- function(libname, pkgname) {
#   if (!reticulate::py_module_available("scipy")) {
#     reticulate::py_install("scipy", pip = TRUE)
#   }
#
#   # Make scipy available globally
#   scipy <<- reticulate::import("scipy.interpolate",
#                                delay_load = TRUE, convert = FALSE)
# }

# #' @importFrom reticulate py_module_available py_install import
# .onLoad <- function(libname, pkgname) {
#   if (!reticulate::py_module_available("scipy")) {
#     reticulate::py_install("scipy", pip = TRUE)
#   }
#
#
#   # Import scipy.interpolate and assign to package-level variable
#   scipy_fitpack <<- reticulate::import("scipy.interpolate.dfitpack",
#                                        delay_load = TRUE,
#                                        convert = FALSE)
# }
