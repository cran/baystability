#-------------------------------------------------------------------------------------------------------------------
# The sources of every method: one registry gives the short citation printed under each result, the references of
# the interpretations and the application, bys_cite(), inst/CITATION and inst/REFERENCES.bib (written from it by
# tools/write-citation.R).
#-------------------------------------------------------------------------------------------------------------------

.bys_ref <- function(key, short, text, ...) {
  list(key = key, short = short, text = text, bib = utils::bibentry(key = key, ...))
}

.bys_p <- function(given, family) utils::person(given = given, family = family)

.bys_refs_list <- function() {
  p <- .bys_p
  refs <- list(
    .bys_ref("Crossa2011", "Crossa et al. (2011)",
             paste("Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L.",
                   "(2011). Bayesian estimation of the additive main effects and multiplicative interaction model.",
                   "Crop Science, 51(4), 1458-1469. doi:10.2135/cropsci2010.06.0343"),
             bibtype = "Article",
             title = "Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model",
             author = c(p("Jose", "Crossa"), p("Sergio", "Perez-Elizalde"), p("Diego", "Jarquin"),
                        p(c("Jose", "Miguel"), "Cotes"), p("Kert", "Viele"), p("Genzhou", "Liu"),
                        p(c("Paul", "L."), "Cornelius")),
             journal = "Crop Science", year = "2011", volume = "51", number = "4", pages = "1458--1469",
             doi = "10.2135/cropsci2010.06.0343"),
    .bys_ref("PerezElizalde2012", "Perez-Elizalde et al. (2012)",
             paste("Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of",
                   "linear-bilinear models applied to plant breeding trials with genotype x environment interaction.",
                   "Journal of Agricultural, Biological, and Environmental Statistics, 17(1), 15-37.",
                   "doi:10.1007/s13253-011-0063-9"),
             bibtype = "Article",
             title = paste("A General {Bayesian} Estimation Method of Linear-Bilinear Models Applied to Plant Breeding",
                           "Trials with Genotype x Environment Interaction"),
             author = c(p("Sergio", "Perez-Elizalde"), p("Diego", "Jarquin"), p("Jose", "Crossa")),
             journal = "Journal of Agricultural, Biological, and Environmental Statistics", year = "2012",
             volume = "17", number = "1", pages = "15--37", doi = "10.1007/s13253-011-0063-9"),
    .bys_ref("Viele2000", "Viele and Srinivasan (2000)",
             paste("Viele, K. and Srinivasan, C. (2000). Parsimonious estimation of multiplicative interaction in",
                   "analysis of variance using Kullback-Leibler information. Journal of Statistical Planning and",
                   "Inference, 84(1-2), 201-219."),
             bibtype = "Article",
             title = paste("Parsimonious Estimation of Multiplicative Interaction in Analysis of Variance Using",
                           "{Kullback-Leibler} Information"),
             author = c(p("Kert", "Viele"), p("C.", "Srinivasan")),
             journal = "Journal of Statistical Planning and Inference", year = "2000", volume = "84",
             number = "1-2", pages = "201--219"),
    .bys_ref("Hoff2009", "Hoff (2009)",
             paste("Hoff, P. D. (2009). Simulation of the matrix Bingham-von Mises-Fisher distribution, with",
                   "applications to multivariate and relational data. Journal of Computational and Graphical",
                   "Statistics, 18(2), 438-456. doi:10.1198/jcgs.2009.07177"),
             bibtype = "Article",
             title = paste("Simulation of the Matrix {Bingham-von Mises-Fisher} Distribution, with Applications to",
                           "Multivariate and Relational Data"),
             author = p(c("Peter", "D."), "Hoff"),
             journal = "Journal of Computational and Graphical Statistics", year = "2009", volume = "18",
             number = "2", pages = "438--456", doi = "10.1198/jcgs.2009.07177"),
    .bys_ref("Yan2000", "Yan et al. (2000)",
             paste("Yan, W., Hunt, L. A., Sheng, Q. and Szlavnics, Z. (2000). Cultivar evaluation and",
                   "mega-environment investigation based on the GGE biplot. Crop Science, 40(3), 597-605.",
                   "doi:10.2135/cropsci2000.403597x"),
             bibtype = "Article",
             title = "Cultivar Evaluation and Mega-Environment Investigation Based on the {GGE} Biplot",
             author = c(p("W.", "Yan"), p(c("L.", "A."), "Hunt"), p("Q.", "Sheng"), p("Z.", "Szlavnics")),
             journal = "Crop Science", year = "2000", volume = "40", number = "3", pages = "597--605",
             doi = "10.2135/cropsci2000.403597x"),
    .bys_ref("Spiegelhalter2002", "Spiegelhalter et al. (2002)",
             paste("Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A. (2002). Bayesian measures",
                   "of model complexity and fit. Journal of the Royal Statistical Society B, 64(4), 583-639.",
                   "doi:10.1111/1467-9868.00353"),
             bibtype = "Article", title = "{Bayesian} Measures of Model Complexity and Fit",
             author = c(p(c("David", "J."), "Spiegelhalter"), p(c("Nicola", "G."), "Best"),
                        p(c("Bradley", "P."), "Carlin"), p("Angelika", "van der Linde")),
             journal = "Journal of the Royal Statistical Society B", year = "2002", volume = "64", number = "4",
             pages = "583--639", doi = "10.1111/1467-9868.00353"),
    .bys_ref("Watanabe2010", "Watanabe (2010)",
             paste("Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable",
                   "information criterion in singular learning theory. Journal of Machine Learning Research, 11,",
                   "3571-3594."),
             bibtype = "Article",
             title = paste("Asymptotic Equivalence of {Bayes} Cross Validation and Widely Applicable Information",
                           "Criterion in Singular Learning Theory"),
             author = p("Sumio", "Watanabe"), journal = "Journal of Machine Learning Research", year = "2010",
             volume = "11", pages = "3571--3594"),
    .bys_ref("Vehtari2017", "Vehtari et al. (2017)",
             paste("Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using",
                   "leave-one-out cross-validation and WAIC. Statistics and Computing, 27(5), 1413-1432.",
                   "doi:10.1007/s11222-016-9696-4"),
             bibtype = "Article",
             title = "Practical {Bayesian} Model Evaluation Using Leave-One-Out Cross-Validation and {WAIC}",
             author = c(p("Aki", "Vehtari"), p("Andrew", "Gelman"), p("Jonah", "Gabry")),
             journal = "Statistics and Computing", year = "2017", volume = "27", number = "5",
             pages = "1413--1432", doi = "10.1007/s11222-016-9696-4"),
    .bys_ref("Gabry2019", "Gabry et al. (2019)",
             paste("Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman, A. (2019). Visualization in",
                   "Bayesian workflow. Journal of the Royal Statistical Society A, 182(2), 389-402.",
                   "doi:10.1111/rssa.12378"),
             bibtype = "Article", title = "Visualization in {Bayesian} Workflow",
             author = c(p("Jonah", "Gabry"), p("Daniel", "Simpson"), p("Aki", "Vehtari"), p("Michael", "Betancourt"),
                        p("Andrew", "Gelman")),
             journal = "Journal of the Royal Statistical Society A", year = "2019", volume = "182", number = "2",
             pages = "389--402", doi = "10.1111/rssa.12378")
    )
  stats::setNames(refs, vapply(refs, `[[`, "", "key"))
}

# the sources of a result, by its class (and model)
.bys_sources_of <- function(x) {
  cls <- .bys_cls(x)
  s <- x$settings
  switch(cls,
    bayes_ammi = c("PerezElizalde2012", "Crossa2011", if (identical(s$model, "SREG")) "Yan2000", "Hoff2009",
                   "Spiegelhalter2002", "Watanabe2010"),
    biplots = c("PerezElizalde2012", "Crossa2011"),
    ge_ammi = "PerezElizalde2012",
    compare = c("Spiegelhalter2002", "Watanabe2010", "Vehtari2017"),
    character(0))
}

# "Sources: Crossa et al. (2011); ..." for the end of a print
.bys_sources_line <- function(x) {
  k <- unique(.bys_sources_of(x))
  if (!length(k)) return(invisible(NULL))
  r <- .bys_refs_list()[k]
  txt <- strwrap(paste0("Sources: ", paste(vapply(r, `[[`, "", "short"), collapse = "; "),
                        ". Full references: bys_cite(x)."), width = getOption("width", 80L))
  cat("\n", paste(txt, collapse = "\n"), "\n", sep = "")
  invisible(NULL)
}

#' References of a result or of the package
#'
#' @description
#' \code{bys_cite()} gives the references of the methods behind a result of the package, or the references of
#' every method in the package, as text or BibTeX. \code{citation("baystability")} cites the package and its main
#' methods.
#'
#' @param x a result of a function of the package, or \code{NULL} (default) for every reference of the package.
#' @param style \code{"text"} (default) or \code{"bibtex"}.
#' @return An object of class \code{bibentry} (printed in the chosen style).
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @examples
#' bys_cite(
#'     x     = NULL
#'   , style = "text"
#'   )
#' data(cultivo2009)
#' g <-
#'   bys_ls_ammi(
#'       .data = cultivo2009
#'     , .y    = y
#'     , .gen  = entry
#'     , .env  = site
#'     , .rep  = rep
#'     )
#' bys_cite(
#'     x     = g
#'   , style = "bibtex"
#'   )
#' @export
bys_cite <- function(x = NULL, style = c("text", "bibtex")) UseMethod("bys_cite")

#' @rdname bys_cite
#' @export
bys_cite.default <- function(x = NULL, style = c("text", "bibtex")) {
  style <- match.arg(style)
  refs <- .bys_refs_list()
  keys <- if (is.null(x)) names(refs) else {
    if (!inherits(x, "bys_result")) .bys_stop("bys_cite", "`x` must be a result of a function of the package or NULL.")
    unique(.bys_sources_of(x))
  }
  out <- do.call(c, lapply(refs[keys], `[[`, "bib"))
  if (is.null(out)) out <- utils::bibentry("Misc", title = "none", author = "none", year = "none")[0L]
  structure(out, bys_style = style, class = c("bys_cite", class(out)))
}

#' @export
print.bys_cite <- function(x, ...) {
  style <- attr(x, "bys_style")
  y <- x
  class(y) <- setdiff(class(y), "bys_cite")
  if (!length(y)) {
    cat("No references.\n")
    return(invisible(x))
  }
  if (identical(style, "bibtex")) print(y, style = "Bibtex") else {
    refs <- .bys_refs_list()
    keys <- vapply(unclass(y), function(e) attr(e, "key"), "")
    cat(paste(vapply(refs[keys], `[[`, "", "text"), collapse = "\n\n"), "\n")
  }
  invisible(x)
}
