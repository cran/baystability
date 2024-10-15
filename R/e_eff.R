#' @name    e_eff
#' @aliases e_eff
#' @title   Environment Effects
#' @description Calcuates Environment Effects
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Environment Effects
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
#' @import dplyr ggplot2 lme4 rlang rstiefel tibble tidyr scales
#' @importFrom MASS Null mvrnorm
#' @importFrom stats rexp rgamma rnorm runif sigma var
#'
#' @export
#'
#' @examples
#'
#' data(cultivo2008)
#' e_eff(
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
      "Env"
      , "Y"
      , "EMean"
      , "."
      , "EEffs"
    )
  )
}


e_eff <- function(.data, .y, .gen, .env) {
  UseMethod("e_eff")
}


#' @export
#' @rdname e_eff

e_eff.default <-
  function(.data, .y, .gen, .env){


    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
      Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))


    e_effects <-
      df1 %>%
      dplyr::group_by(Env) %>%
      dplyr::summarise(EMean = mean(Y)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(EEffs = EMean - mean(EMean)) %>%
    #  magrittr::set_rownames(.$Env) %>%
      dplyr::ungroup() %>%
      dplyr::select(EEffs) %>%
      as.matrix()

    return(e_effects)
}
