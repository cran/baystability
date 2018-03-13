#' @name    ge_ammi
#' @aliases ge_ammi
#' @title   AMMI of Genotype by Environment Interaction Model
#' @description Performs Additive Main Effects and Multiplication Interaction Analysis of Genotype by Environment Interaction Model
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
#' @import tidyverse
#' @import tidyr
#' @import lme4
#' @import rlang
#' @importFrom magrittr %>%
#' @importFrom stats sigma
#'
#' @export
#'
#' @examples
#'
#' data(cultivo2008)
#' fm1 <-
#'    ge_ammi(
#'       .data  = cultivo2008
#'      , .y    = y
#'      , .gen  = entry
#'      , .env  = site
#'      , .rep  = rep
#'      )
#'
#'
#' data(cultivo2009)
#' fm2 <-
#'    ge_ammi(
#'       .data  = cultivo2009
#'      , .y    = y
#'      , .gen  = entry
#'      , .env  = site
#'      , .rep  = rep
#'      )
#'
#'

ge_ammi <- function(.data, .y, .gen, .env, .rep) {
  UseMethod("ge_ammi")
}


#' @export
#' @rdname ge_ammi

ge_ammi.default <-
  function(.data, .y, .gen, .env, .rep){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))
    .rep  <- deparse(substitute(.rep))

    df1 <- tibble::as_tibble(data.frame(
        Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Rep = factor(.data[[.rep]])
      , Y   = .data[[.y]]
    ))

    g   <- length(levels(df1$Gen))
    e   <- length(levels(df1$Env))
    Rep <- length(levels(df1$Rep))
    k   <- min(g, e) - 1

    ge_means <-
      ge_mean(
         .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
        )

    ge_effects <-
      ge_eff(
          .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
      )


    g_effects <-
      g_eff(
        .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
      )

    e_effects <-
      e_eff(
        .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
      )

    ge_model <-
      ge_model(
         .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
        , .rep  = Rep
      )

    return(list(
      g        = g
    , e        = e
    , Rep      = Rep
    , k        = k
    , ge_means = ge_means
    , mu       = ge_means$grand_mean
    , sigma2   = sigma(ge_model)^2
    , tau      = 1/sigma(ge_model)^2
    , tao      = g_effects
    , delta    = e_effects
    , lambdas  = ge_effects$ge_svd$d[1:k]
    , alphas   = ge_effects$ge_svd$u[,1:k]
    , gammas   = ge_effects$ge_svd$v[,1:k]
      ))
}
