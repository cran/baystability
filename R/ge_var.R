#' @name    ge_var
#' @aliases ge_var
#' @title   Genotype by Environment Interaction Variances
#' @description Calcuates Genotype by Environment Interaction Variances
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Genotype by Environment Interaction Variances
#'
#' @author
#' \enumerate{
#'     \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'    }
#'
#' @references
#'  Perez-Elizalde, S., Jarquin, D., and Crossa, J. (2011)
#'  A General Bayesian Estimation Method of Linear–Bilinear Models
#'  Applied to Plant Breeding Trials With Genotype × Environment Interaction.
#'  \emph{Journal of Agricultural, Biological, and Environmental Statistics},
#'   17, 15–37.  (\href{https://link.springer.com/article/10.1007/s13253-011-0063-9}{doi:10.1007/s13253-011-0063-9})
#'
#' @import tidyverse
#' @import tidyr
#' @import rlang
#' @importFrom magrittr %>%
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#'
#' data(cultivo2008)
#' ge_var(
#'     .data  = cultivo2008
#'    , .y    = y
#'    , .gen  = entry
#'    , .env  = site
#'    )
#'
#'


if(getRversion() >= "2.15.1"){
  utils::globalVariables(
    c(
      "Gen"
      , "Env"
      , "Y"
      , "GEVar"
      , "."
    )
  )
}

ge_var <- function(.data, .y, .gen, .env) {
  UseMethod("ge_var")
}


#' @export
#' @rdname ge_var

ge_var.default <-
  function(.data, .y, .gen, .env){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
      Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))


    ge_variances <-
      df1 %>%
      dplyr::group_by(Gen, Env) %>%
      dplyr::summarise(GEVar = var(Y)) %>%
      tidyr::spread(key = Env, value = GEVar) %>%
     # magrittr::set_rownames(.$Gen) %>%
      dplyr::ungroup() %>%
      dplyr::select(- Gen) %>%
      as.matrix()

    return(list(
      ge_variances   = ge_variances
      ))
}
