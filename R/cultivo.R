#' Maize trials of 2008 and 2009
#'
#' @description
#' Grain yield of 12 maize hybrids in 25 environments in two consecutive years, the multi-environment trial of
#' Perez-Elizalde et al. (2012), Section 3.2: \code{cultivo2008} (the first year, used as the prior) and
#' \code{cultivo2009} (the second year). Randomized complete blocks with 3 replicates in every environment.
#'
#' @format Data frames with 900 rows (12 genotypes x 25 environments x 3 replicates) and 4 columns:
#' \describe{
#'   \item{\code{site}}{environment, 1 to 25 (the same environments in both years).}
#'   \item{\code{rep}}{replicate, 1 to 3.}
#'   \item{\code{entry}}{genotype, 1 to 12.}
#'   \item{\code{y}}{grain yield (kg/ha).}
#' }
#' @details The paper describes two replicates; the data have three. The posterior means of its Table B.1 are, to
#'   within a posterior standard deviation of a hundredth, the averages of the least-squares estimates of the two
#'   years, as the method gives with equal weights for both years (see the vignette on the published analysis).
#'   The help page of version 0.2.0 described another data set.
#' @source Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of
#'   linear-bilinear models applied to plant breeding trials with genotype x environment interaction. \emph{Journal of
#'   Agricultural, Biological, and Environmental Statistics}, \strong{17}(1), 15--37. \doi{10.1007/s13253-011-0063-9}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @examples
#' data(cultivo2008)
#' data(cultivo2009)
#' str(cultivo2009)
#' bys_ls_means(
#'     .data = cultivo2009
#'   , .y    = y
#'   , .gen  = entry
#'   , .env  = site
#'   )$grand_mean
#' @name cultivo2008
#' @aliases cultivo2009
#' @docType data
#' @keywords datasets
NULL
