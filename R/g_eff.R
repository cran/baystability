#' @name    g_eff
#' @aliases g_eff
#' @title   Genotype Effects
#' @description Calcuates Genotype Effects
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Genotype Effects
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
#' g_eff(
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
      , "Y"
      , "GMean"
      , "."
      , "GEffs"
    )
  )
}
g_eff <- function(.data, .y, .gen, .env) {
  UseMethod("g_eff")
}


#' @export
#' @rdname g_eff

g_eff.default <-
  function(.data, .y, .gen, .env){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
      Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))

    g_effects <-
      df1 %>%
        dplyr::group_by(Gen) %>%
        dplyr::summarise(GMean = mean(Y)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(GEffs = GMean - mean(GMean)) %>%
      #  magrittr::set_rownames(.$Gen) %>%
        dplyr::ungroup() %>%
        dplyr::select(GEffs) %>%
        as.matrix()

    return(g_effects)
}

