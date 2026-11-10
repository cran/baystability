#' Deprecated names of baystability 0.2.0
#'
#' @description
#' These functions of baystability 0.2.0 still work in 1.0.0, with a warning, and will be removed in a later version. Use
#' the new names, which follow those of the package \pkg{stability} and its Bayesian version \pkg{bayesstability}:
#' \tabular{ll}{
#'   old name \tab new name \cr
#'   \code{ge_mean()} \tab \code{\link{bys_ls_means}()} \cr
#'   \code{ge_var()} \tab \code{\link{bys_ls_var}()} \cr
#'   \code{g_eff()} \tab \code{\link{bys_ls_gen_effects}()} \cr
#'   \code{e_eff()} \tab \code{\link{bys_ls_env_effects}()} \cr
#'   \code{ge_eff()} \tab \code{\link{bys_ls_ge_effects}()} \cr
#'   \code{ge_model()} \tab \code{\link{bys_ls_model}()} \cr
#'   \code{ge_ammi()} \tab \code{\link{bys_ls_ammi}()} \cr
#'   \code{matrix_k()} \tab \code{\link{bys_contrasts}()} \cr
#'   \code{orthnorm()} \tab \code{\link{bys_orthonormal}()} \cr
#' }
#' The arguments and results are unchanged.
#'
#' @param .data,.y,.gen,.env,.rep,n,u,basis,norm
#'   the arguments of the earlier version, passed on to the new function.
#' @return The result of the new function.
#' @name bys-deprecated
#' @keywords internal
NULL

# call the new function in place of a deprecated one, from the caller's frame
.bys_deprecated <- function(old, new, cl, env, rename = NULL, fixed = list()) {
  .Deprecated(new = new, package = "baystability", old = old)
  cl[[1L]] <- call("::", as.name("baystability"), as.name(new))
  for (a in names(rename)) names(cl)[names(cl) == a] <- rename[[a]]
  for (a in names(fixed)) cl[[a]] <- fixed[[a]]
  eval(cl, env)
}

#' @rdname bys-deprecated
#' @export
ge_mean <- function(.data, .y, .gen, .env) UseMethod("ge_mean")

#' @rdname bys-deprecated
#' @export
ge_mean.default <- function(.data, .y, .gen, .env) {
  .bys_deprecated("ge_mean", "bys_ls_means", match.call(), parent.frame())
}

#' @rdname bys-deprecated
#' @export
ge_var <- function(.data, .y, .gen, .env) UseMethod("ge_var")

#' @rdname bys-deprecated
#' @export
ge_var.default <- function(.data, .y, .gen, .env) {
  .bys_deprecated("ge_var", "bys_ls_var", match.call(), parent.frame())
}

#' @rdname bys-deprecated
#' @export
g_eff <- function(.data, .y, .gen, .env) UseMethod("g_eff")

#' @rdname bys-deprecated
#' @export
g_eff.default <- function(.data, .y, .gen, .env) {
  .bys_deprecated("g_eff", "bys_ls_gen_effects", match.call(), parent.frame())
}

#' @rdname bys-deprecated
#' @export
e_eff <- function(.data, .y, .gen, .env) UseMethod("e_eff")

#' @rdname bys-deprecated
#' @export
e_eff.default <- function(.data, .y, .gen, .env) {
  .bys_deprecated("e_eff", "bys_ls_env_effects", match.call(), parent.frame())
}

#' @rdname bys-deprecated
#' @export
ge_eff <- function(.data, .y, .gen, .env) UseMethod("ge_eff")

#' @rdname bys-deprecated
#' @export
ge_eff.default <- function(.data, .y, .gen, .env) {
  .bys_deprecated("ge_eff", "bys_ls_ge_effects", match.call(), parent.frame())
}

#' @rdname bys-deprecated
#' @export
ge_model <- function(.data, .y, .gen, .env, .rep) UseMethod("ge_model")

#' @rdname bys-deprecated
#' @export
ge_model.default <- function(.data, .y, .gen, .env, .rep) {
  .bys_deprecated("ge_model", "bys_ls_model", match.call(), parent.frame())
}

#' @rdname bys-deprecated
#' @export
ge_ammi <- function(.data, .y, .gen, .env, .rep) UseMethod("ge_ammi")

#' @rdname bys-deprecated
#' @export
ge_ammi.default <- function(.data, .y, .gen, .env, .rep) {
  .bys_deprecated("ge_ammi", "bys_ls_ammi", match.call(), parent.frame())
}

#' @rdname bys-deprecated
#' @export
matrix_k <- function(n) UseMethod("matrix_k")

#' @rdname bys-deprecated
#' @export
matrix_k.default <- function(n) {
  .bys_deprecated("matrix_k", "bys_contrasts", match.call(), parent.frame())
}

#' @rdname bys-deprecated
#' @export
orthnorm <- function(u = NULL, basis = TRUE, norm = TRUE) UseMethod("orthnorm")

#' @rdname bys-deprecated
#' @export
orthnorm.default <- function(u = NULL, basis = TRUE, norm = TRUE) {
  .bys_deprecated("orthnorm", "bys_orthonormal", match.call(), parent.frame())
}
