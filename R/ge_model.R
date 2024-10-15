#' @name    ge_model
#' @aliases ge_model
#' @title   Genotype by Environment Interaction Model
#' @description Calcuates Genotype by Environment Interaction Model
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#' @param .rep   Replication Factor
#'
#' @return Genotype by Environment Interaction Model
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
#' fm1 <-
#'    ge_model(
#'       .data  = cultivo2008
#'      , .y    = y
#'      , .gen  = entry
#'      , .env  = site
#'      , .rep  = rep
#'      )
#'
#'

ge_model <- function(.data, .y, .gen, .env, .rep) {
  UseMethod("ge_model")
}


#' @export
#' @rdname ge_model

ge_model.default <-
  function(.data, .y, .gen, .env, .rep){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))
    .rep  <- deparse(substitute(.rep))

    df1 <- data.frame(
            Y = .data[[.y]]
          , Env = factor(.data[[.env]])
          , Gen = factor(.data[[.gen]])
          , Rep = factor(.data[[.rep]])
          )

     ge_fm <-
       lme4::lmer(Y ~ Env + Gen + Env:Gen + (1|Env:Rep), data = df1)
      return(ge_fm)
}
