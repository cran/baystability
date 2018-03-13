#' @name    ge_mean
#' @aliases ge_mean
#' @title   Genotype by Environment Interaction Means
#' @description Calcuates Genotype by Environment Interaction Means
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Genotype by Environment Interaction Means
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
#' ge_mean(
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
      , "."
    )
  )
}

ge_mean <- function(.data, .y, .gen, .env) {
  UseMethod("ge_mean")
}


#' @export
#' @rdname ge_mean

ge_mean.default <-
  function(.data, .y, .gen, .env){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))

    df1 <- tibble::as_tibble(data.frame(
        Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Y   = .data[[.y]]
    ))

    ge_means <-
      df1 %>%
      dplyr::group_by(Gen, Env) %>%
      dplyr::summarise(GEMean = mean(Y)) %>%
      tidyr::spread(key = Env, value = GEMean) %>%
    #  magrittr::set_rownames(.$Gen) %>%
      dplyr::ungroup() %>%
      dplyr::select(- Gen) %>%
      as.matrix()

    grand_mean <- mean(ge_means)

    return(list(
        ge_means   = ge_means
      , grand_mean = grand_mean
      ))
}
