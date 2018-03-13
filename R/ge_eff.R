#' @name    ge_eff
#' @aliases ge_eff
#' @title   Genotype by Environment Interaction Effects
#' @description Calcuates Genotype by Environment Interaction Effects
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Genotype by Environment Interaction Effects
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
#'
#' @export
#'
#' @examples
#'
#' data(cultivo2008)
#' ge_eff(
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
      , "GEMean"
      , "GMean"
      , "EMean"
      , "Mean"
      , "GEEffs"
      , "."
    )
  )
}


ge_eff <- function(.data, .y, .gen, .env) {
  UseMethod("ge_eff")
}


#' @export
#' @rdname ge_eff

ge_eff.default <-
  function(.data, .y, .gen, .env){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
      Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))


    ge_effects <-
      df1 %>%
      dplyr::group_by(Gen, Env) %>%
      dplyr::mutate(GEMean = mean(Y)) %>%
      dplyr::group_by(Gen) %>%
      dplyr::mutate(GMean = mean(Y)) %>%
      dplyr::group_by(Env) %>%
      dplyr::mutate(EMean = mean(Y)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Mean = mean(Y)) %>%
      dplyr::group_by(Gen, Env) %>%
      dplyr::summarize(GEEffs = mean(GEMean - GMean - EMean + Mean)) %>%
      tidyr::spread(key = Env, value = GEEffs) %>%
  #    magrittr::set_rownames(.$Gen) %>%
      dplyr::ungroup() %>%
      dplyr::select(- Gen) %>%
      as.matrix()

    ge_svd <- svd(ge_effects)

    return(list(
        ge_effects = ge_effects
      , ge_svd     = ge_svd
      ))
}
