#' @name    bayes_ammi
#' @aliases bayes_ammi
#' @title   Bayesian Estimation of Genotype by Environment Interaction Model
#' @description Bayesian estimation method of linear–bilinear models for Genotype by Environment Interaction Model
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#' @param .rep   Replication Factor
#' @param .nIter Number of Iterations
#'
#' @return Genotype by Environment Interaction Model
#'
#' @author
#' \enumerate{
#'     \item Muhammad Yaseen (\email{myaseen208@gmail.com})
#'     \item Diego Jarquin (\email{diego.jarquin@gmail.com})
#'     \item Sergio Perez-Elizalde (\email{sergiop@colpos.mx})
#'     \item Juan Burgueño (\email{j.burgueno@cgiar.org})
#'     \item Jose Crossa (\email{j.crossa@cgiar.org})
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
#' @importFrom stats rexp rgamma rnorm runif setNames sigma var
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(baystability)
#' library(dplyr)
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
#' r0 <- fm1$g
#' c0 <- fm1$e
#' n0 <- fm1$Rep
#' k0 <- fm1$k
#'
#' mu0      <- fm1$mu
#' sigma20  <- fm1$sigma2
#' tau0     <- fm1$tau
#' tao0     <- fm1$tao
#' delta0   <- fm1$delta
#' lambdas0 <- fm1$lambdas
#' alphas0  <- fm1$alphas
#' gammas0  <- fm1$gammas
#'
#' ge_means0 <- fm1$ge_means$ge_means
#'
#' data(cultivo2008)
#'
#' fm2 <-
#'  ge_ammi(
#'    .data = cultivo2009
#'    , .y    = y
#'    , .gen  = entry
#'    , .env  = site
#'    , .rep  = rep
#'   )
#'
#' k        <- fm2$k
#' alphasa  <- fm2$alphas
#' gammasa  <- fm2$gammas
#'
#' alphas1  <-
#'            tibble::as_tibble(fm2$alphas, .name_repair = "unique") %>%
#'            setNames(paste0("V", 1:ncol(.)))
#' gammas1  <-
#'            tibble::as_tibble(fm2$gammas, .name_repair = "unique") %>%
#'            setNames(paste0("V", 1:ncol(.)))
#'
#'
#' # Biplots OLS
#' library(ggplot2)
#'     BiplotOLS1 <-
#'       ggplot(data = alphas1, mapping = aes(x = V1, y = V2)) +
#'       geom_point() +
#'       geom_hline(yintercept = 0) +
#'       geom_vline(xintercept = 0) +
#'       geom_text(aes(label = 1:nrow(alphas1)), vjust = "inward", hjust = "inward") +
#'       scale_x_continuous(
#'                limits = c(-max(abs(c(range(alphas1[, 1:2]))))
#'                          , max(abs(c(range(alphas1[, 1:2])))))) +
#'       scale_y_continuous(
#'                limits = c(-max(abs(c(range(alphas1[, 1:2]))))
#'                          , max(abs(c(range(alphas1[, 1:2])))))) +
#'       labs(title = "OLS", x = expression(u[1]), y = expression(u[2])) +
#'       theme_bw() +
#'       theme(plot.title = element_text(hjust = 0.5))
#'       print(BiplotOLS1)
#'
#'
#'     BiplotOLS2 <-
#'       ggplot(data = gammas1, mapping = aes(x = V1, y = V2)) +
#'       geom_point() +
#'       geom_hline(yintercept = 0) +
#'       geom_vline(xintercept = 0) +
#'       geom_text(aes(label = 1:nrow(gammas1)), vjust = "inward", hjust = "inward") +
#'       scale_x_continuous(
#'                 limits = c(-max(abs(c(range(gammas1[, 1:2]))))
#'                          , max(abs(c(range(gammas1[, 1:2])))))) +
#'       scale_y_continuous(
#'                 limits = c(-max(abs(c(range(gammas1[, 1:2]))))
#'                           , max(abs(c(range(gammas1[, 1:2])))))) +
#'       labs(title = "OLS", x = expression(v[1]), y = expression(v[2])) +
#'       theme_bw() +
#'       theme(plot.title = element_text(hjust = 0.5))
#'       print(BiplotOLS2)
#'
#'
#'     BiplotOLS3 <-
#'       ggplot(data = alphas1, mapping = aes(x = V1, y = V2)) +
#'       geom_point() +
#'       geom_hline(yintercept = 0) +
#'       geom_vline(xintercept = 0) +
#'       geom_text(aes(label = 1:nrow(alphas1)), vjust = "inward", hjust = "inward") +
#'       geom_point(data = gammas1, mapping = aes(x = V1, y = V2)) +
#'       geom_segment(data = gammas1, aes(x = 0, y = 0, xend = V1, yend = V2),
#'                     arrow = arrow(length = unit(0.2, "cm")), alpha = 0.75, color = "red") +
#'       geom_text(data = gammas1,
#'               aes(x = V1, y = V2, label = paste0("E", 1:nrow(gammasa)))
#'              , vjust = "inward", hjust = "inward") +
#'       scale_x_continuous(
#'               limits = c(-max(abs(c(range(alphas1[, 1:2], gammas1[, 1:2]))))
#'                         , max(abs(c(range(alphas1[, 1:2], gammas1[, 1:2])))))) +
#'       scale_y_continuous(
#'               limits = c(-max(abs(c(range(alphas1[, 1:2], gammas1[, 1:2]))))
#'                        , max(abs(c(range(alphas1[, 1:2], gammas1[, 1:2])))))) +
#'       labs(title = "OLS", x = expression(PC[1]), y = expression(PC[2])) +
#'       theme_bw() +
#'       theme(plot.title = element_text(hjust = 0.5))
#'       print(BiplotOLS3)
#'
#' data(cultivo2009)
#' fm3 <-
#'   bayes_ammi(
#'     .data = cultivo2009
#'     , .y     = y
#'     , .gen   = entry
#'     , .env   = site
#'     , .rep   = rep
#'     , .nIter = 200
#'   )
#'
#'  Mean_Alphas <- fm3$Mean_Alphas %>% setNames(paste0("V", 1:ncol(.)))
#'  Mean_Gammas <- fm3$Mean_Gammas %>% setNames(paste0("V", 1:ncol(.)))
#'
#'
#' # Biplots Bayesian
#' BiplotBayes1 <-
#'   ggplot(data = Mean_Alphas, mapping = aes(x = V1, y = V2)) +
#'   geom_point() +
#'   geom_hline(yintercept = 0) +
#'   geom_vline(xintercept = 0) +
#'   geom_text(aes(label = 1:nrow(Mean_Alphas)),
#'              vjust = "inward"
#'            , hjust = "inward") +
#'   scale_x_continuous(
#'      limits = c(-max(abs(c(range(Mean_Alphas[, 1:2]))))
#'                , max(abs(c(range(Mean_Alphas[, 1:2])))))) +
#'   scale_y_continuous(
#'       limits = c(-max(abs(c(range(Mean_Alphas[, 1:2]))))
#'                 , max(abs(c(range(Mean_Alphas[, 1:2])))))) +
#'   labs(title = "Bayes", x = expression(u[1]), y = expression(u[2])) +
#'   theme_bw() +
#'   theme(plot.title = element_text(hjust = 0.5))
#'
#' print(BiplotBayes1)
#'
#'
#' BiplotBayes2 <-
#'   ggplot(data = Mean_Gammas, mapping = aes(x = V1, y = V2)) +
#'   geom_point() +
#'   geom_hline(yintercept = 0) +
#'   geom_vline(xintercept = 0) +
#'   geom_text(aes(label = 1:nrow(Mean_Gammas)), vjust = "inward", hjust = "inward") +
#'   scale_x_continuous(
#'             limits = c(-max(abs(c(range(Mean_Gammas[, 1:2]))))
#'                       , max(abs(c(range(Mean_Gammas[, 1:2])))))) +
#'   scale_y_continuous(
#'             limits = c(-max(abs(c(range(Mean_Gammas[, 1:2]))))
#'                       , max(abs(c(range(Mean_Gammas[, 1:2])))))) +
#'   labs(title = "Bayes", x = expression(v[1]), y = expression(v[2])) +
#'   theme_bw() +
#'   theme(plot.title = element_text(hjust = 0.5))
#'
#' print(BiplotBayes2)
#'
#'
#' BiplotBayes3 <-
#'   ggplot(data = Mean_Alphas, mapping = aes(x = V1, y = V2)) +
#'   geom_point() +
#'   geom_hline(yintercept = 0) +
#'   geom_vline(xintercept = 0) +
#'   geom_text(aes(label = 1:nrow(Mean_Alphas)),
#'              vjust = "inward", hjust = "inward") +
#'   geom_point(data = Mean_Gammas, mapping = aes(x = V1, y = V2)) +
#'   geom_segment(data = Mean_Gammas,
#'                 aes(x = 0, y = 0, xend = V1, yend = V2),
#'                arrow = arrow(length = unit(0.2, "cm"))
#'                , alpha = 0.75, color = "red") +
#'   geom_text(data = Mean_Gammas,
#'             aes(x = V1, y = V2,
#'             label = paste0("E", 1:nrow(Mean_Gammas))),
#'             vjust = "inward", hjust = "inward") +
#'   scale_x_continuous(
#'             limits = c(-max(abs(c(range(Mean_Alphas[, 1:2], Mean_Gammas[, 1:2]))))
#'                       , max(abs(c(range(Mean_Alphas[, 1:2], Mean_Gammas[, 1:2])))))) +
#'   scale_y_continuous(
#'            limits = c(-max(abs(c(range(Mean_Alphas[, 1:2], Mean_Gammas[, 1:2]))))
#'                    , max(abs(c(range(Mean_Alphas[, 1:2], Mean_Gammas[, 1:2])))))) +
#'   labs(title = "Bayes", x = expression(PC[1]), y = expression(PC[2])) +
#'   theme_bw() +
#'   theme(plot.title = element_text(hjust = 0.5))
#'
#' print(BiplotBayes3)
#'
#' Plot1Mu <-
#'   ggplot(data = fm3$mu1, mapping = aes(x = 1:nrow(fm3$mu1), y = mu)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(mu), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Mu)
#'
#' Plot2Mu <-
#'   ggplot(data = fm3$mu1, mapping = aes(mu)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(mu)) +
#'   theme_bw()
#' print(Plot2Mu)
#'
#'
#' Plot1Sigma2 <-
#'   ggplot(data = fm3$tau1, mapping = aes(x = 1:nrow(fm3$tau1), y = tau)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(sigma^2), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Sigma2)
#'
#'
#' Plot2Sigma2 <-
#'   ggplot(data = fm3$tau1, mapping = aes(tau)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(sigma^2)) +
#'   theme_bw()
#' print(Plot2Sigma2)
#'
#'
#' # Plot of Alphas
#' Plot1Alpha1 <-
#'   ggplot(data = fm3$tao1, mapping = aes(x = 1:nrow(fm3$tao1), y = tao1)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(alpha[1]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Alpha1)
#'
#' Plot2Alpha1 <-
#'   ggplot(data = fm3$tao1, mapping = aes(tao1)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(alpha[1])) +
#'   theme_bw()
#' print(Plot2Alpha1)
#'
#' Plot1Alpha2 <-
#'   ggplot(data = fm3$tao1, mapping = aes(x = 1:nrow(fm3$tao1), y = tao2)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(alpha[2]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Alpha2)
#'
#' Plot2Alpha2 <-
#'   ggplot(data = fm3$tao1, mapping = aes(tao2)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(alpha[2])) +
#'   theme_bw()
#' print(Plot2Alpha2)
#'
#' # Plot of Betas
#' Plot1Beta1 <-
#'   ggplot(data = fm3$delta1, mapping = aes(x = 1:nrow(fm3$delta1), y = delta1)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(beta[1]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Beta1)
#'
#' Plot2Beta1 <-
#'   ggplot(data = fm3$delta1, mapping = aes(delta1)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(beta[1])) +
#'   theme_bw()
#' print(Plot2Beta1)
#'
#'
#' Plot1Beta2 <-
#'   ggplot(data = fm3$delta1, mapping = aes(x = 1:nrow(fm3$delta1), y = delta2)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(beta[2]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Beta2)
#'
#' Plot2Beta2 <-
#'   ggplot(data = fm3$delta1, mapping = aes(delta2)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(beta[2])) +
#'   theme_bw()
#' print(Plot2Beta2)
#'
#'
#' Plot1Beta3 <-
#'   ggplot(data = fm3$delta1, mapping = aes(x = 1:nrow(fm3$delta1), y = delta3)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(beta[3]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Beta3)
#'
#' Plot2Beta3 <-
#'   ggplot(data = fm3$delta1, mapping = aes(delta3)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(beta[3])) +
#'   theme_bw()
#' print(Plot2Beta3)
#' }
#'

if(getRversion() >= "2.15.1"){
  utils::globalVariables(
    c(
      "alphas0"
      , "c0"
      , "delta0"
      , "gammas0"
      , "lambdas0"
      , "ge_means0"
      , "mu0"
      , "n0"
      , "tao0"
      , "tau0"
    )
  )
}



bayes_ammi <- function(.data, .y, .gen, .env, .rep, .nIter){
  UseMethod("bayes_ammi")
}


#' @export
#' @rdname bayes_ammi


bayes_ammi.default <- function(.data, .y, .gen, .env, .rep, .nIter){

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

fm1   <-
    ge_ammi(
          .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
        , .rep  = Rep
      )

  Rep          <- fm1$Rep
  alphas       <- fm1$alphas
  gammas       <- fm1$gammas
  ge_means     <- fm1$ge_means$ge_means
  tau          <- fm1$tau
  tao          <- fm1$tao
  lambdas      <- fm1$lambdas
  delta        <- fm1$delta
  mu           <- fm1$ge_means$grand_mean
  g            <- fm1$g
  e            <- fm1$e
  k            <- min(c(g, e)) - 1

  ge_variances <- ge_var(.data = df1, .y = Y, .gen = Gen, .env = Env)$ge_variances

  phi          <- 0
  uni          <- 1
  med_lamb     <- (Rep*tau*(sum((alphas %*% t(gammas))*ge_means)) + n0*tau*lambdas0)/(Rep*tau + n0*tau)
  u_neg        <- - med_lamb/sqrt((Rep*tau + n0*tau)^(-1))
  u_pos        <- (lambdas[-length(lambdas)] - med_lamb[-1])/sqrt((Rep*tau + n0*tau)^(-1))
  alpha_star   <- (u_neg + sqrt(u_neg^2 + 4))/2
  z            <- c(runif(n = 1, min = u_neg, max = u_pos), rexp(n = (k-1), rate = alpha_star) + u_neg[-1])
  lambdas      <- sqrt((Rep*tau + n0*tau)^(-1))*z + med_lamb


  #===================Gibbs Sampling Part(No Change)=========
  rmf.matrix.gibbs3 <- function(M, X){
    sM <- svd(M)
    H  <- sM$u %*% diag(sM$d)
    Y  <- X %*% sM$v
    m  <- dim(H)[1]
    R  <- dim(H)[2]

    for(r in sample (1:R)){
      N      <- MASS::Null(Y[ ,-r])
      y      <- rstiefel::rmf.vector(t(N) %*% H[ ,r])
      Y[ ,r] <- N %*% y
    }

    entrada  <- Y %*% t(sM$v)
    entrada1 <- cbind(rep(1, m), entrada)

    salida2  <- orthnorm(u = entrada1, basis = TRUE, norm = TRUE)
    salida3  <- salida2[ ,2:(k+1)]
    return(salida3)
  }

#===================Sampler Function Part==================
  kg <- t(matrix_k(g))
  ke <- t(matrix_k(e))

  mu1      <- NULL
  tau1     <- NULL
  tao1     <- NULL
  delta1   <- NULL
  lambdas1 <- NULL
  alphas1  <- NULL
  gammas1  <- NULL

  r0   <- g
  Jg   <- matrix(rep(1, g), ncol = 1)
  Je   <- matrix(rep(1, e), ncol = 1)

  E    <-
      ge_means                        -
      mu*Jg %*% t(Je)                 -
      matrix(tao, ncol = 1) %*% t(Je) -
      Jg %*% matrix(delta, ncol = e)  -
      alphas %*% diag(lambdas) %*% t(gammas)


  times <- 1
  for(time in 1:.nIter) {
  tau     <-
      rgamma(
          n     =  1
        , shape = Rep*g*e/2+(n0-1)/2
        , rate  = (Rep/2) * sum(diag(E %*% t(E))) +  (1/2) * (Rep-1) * t(Jg) %*% ge_variances %*% Je + (n0-1)/2*tau0^(-1)
      )

  mu       <-
      rnorm(
          n    = 1
        , mean = (Rep* t(Jg) %*% ge_means %*% Je + n0*r0*c0*mu0)/(Rep*g*e+n0*r0*c0)
        , sd   = sqrt(tau^(-1)/(Rep*g*e + n0*r0*c0))
      )

  tao      <- kg %*%
      MASS::mvrnorm(
          n     = 1
        , mu    = (Rep*e*tau*t(kg) %*% ge_means %*% Je/e+n0*c0*tau*t(kg)%*%tao0)/(tau*(Rep*e+n0*c0))
        , Sigma = diag(1, g-1)*(tau*(Rep*e+n0*c0))^(-1)
      )

    delta    <- ke %*%
      MASS::mvrnorm(
          n     = 1
        , mu    = (Rep*g*tau*t(ke) %*% t(ge_means) %*% Jg/g + n0*r0*tau*t(ke) %*% delta0)/(Rep*g*tau+n0*r0*tau)
        , Sigma = diag(1, e-1)*(Rep*g*tau+n0*r0*tau)^(-1)
      )

    alphas   <- rmf.matrix.gibbs3(M = tau*(Rep*(ge_means) %*% gammas %*% diag(lambdas)+n0*ge_means0 %*% gammas0 %*% diag(lambdas0)), X = alphas)
    gammas   <- rmf.matrix.gibbs3(M = tau*(Rep*t(ge_means) %*% alphas %*% diag(lambdas)+n0*t(ge_means0) %*% alphas0 %*% diag(lambdas0)), X = gammas)

    alphas2  <- matrix(alphas, nrow = 1)
    gammas2  <- matrix(gammas, nrow = 1)


    mu1      <-  rbind(mu1, mu)
    tau1     <-  rbind(tau1, tau)
    tao1     <-  cbind(tao1, tao)
    delta1   <-  cbind(delta1, delta)
    lambdas1 <-  cbind(lambdas1, lambdas)
    alphas1  <-  rbind(alphas1, alphas2)
    gammas1  <-  rbind(gammas1, gammas2)

    times   <- times + 1
    }

    mu1      <-  tibble::as_tibble(mu1, .name_repair = "unique")
    tau1     <-  tibble::as_tibble(tau1, .name_repair = "unique")
    tao1     <-  tibble::as_tibble(t(tao1), .name_repair = "unique")
    delta1   <-  tibble::as_tibble(t(delta1), .name_repair = "unique")
    lambdas1 <-  tibble::as_tibble(t(lambdas1), .name_repair = "unique")
    alphas1  <-  tibble::as_tibble(alphas1, .name_repair = "unique")
    gammas1  <-  tibble::as_tibble(gammas1, .name_repair = "unique")

     colnames(mu1)      <- c("mu")
     colnames(tau1)     <- c("tau")
     colnames(tao1)     <- paste0("tao", 1:g)
     colnames(delta1)   <- paste0("delta", 1:g)
     colnames(lambdas1) <- paste0("lambdas", 1:k)
     colnames(alphas1)  <- paste0("alphas", 1:(g*k))
     colnames(gammas1)  <- paste0("gammas", 1:(e*k))

     Mean_Alphas <- tibble::as_tibble(matrix(colMeans(alphas1), ncol = k), .name_repair = "unique")
     Mean_Gammas <- tibble::as_tibble(matrix(colMeans(gammas1), ncol = k), .name_repair = "unique")



    return(list(
      mu1         = mu1
    , tau1        = tau1
    , tao1        = tao1
    , delta1      = delta1
    , lambdas1    = lambdas1
    , alphas1     = alphas1
    , gammas1     = gammas1
    , Mean_Alphas = Mean_Alphas
    , Mean_Gammas = Mean_Gammas
    ))
  }
