#' Bayesian AMMI and SREG analysis in a web browser
#'
#' @description
#' \code{bys_app()} starts the application of \pkg{baystability}: load a trial (the example \code{\link{cultivo2009}}
#' or a CSV file with one row per plot) and, optionally, an earlier trial of the same genotypes and environments as the
#' prior (the example \code{\link{cultivo2008}} or a CSV file), fit the Bayesian AMMI or SREG (GGE) model of
#' Perez-Elizalde et al. (2012), check convergence, draw the biplot with HPD regions and compare models of different
#' numbers of terms, with and without the earlier trial. Every page shows an interpretation computed from the results
#' (\code{\link{bys_interpret}()}), the output as R prints it, formatted tables, plots and the R code that repeats the
#' analysis, each with a Copy button.
#'
#' @details
#' The application needs the suggested packages \pkg{shiny}, \pkg{bslib} and \pkg{DT}. Plots are sent to the browser
#' as WebP images when \pkg{magick} is installed (smaller than PNG; PNG otherwise) and can be downloaded as PNG, PDF
#' or WebP. Text is set in IBM Plex Sans and output in IBM Plex Mono when the browser can reach Google Fonts, in the
#' system fonts otherwise.
#'
#' @param launch.browser passed to \code{shiny::runApp()}: \code{TRUE} (default in an interactive session) opens the
#'   application in the web browser.
#' @param ... further arguments to \code{shiny::runApp()}, such as \code{port}.
#' @return \code{bys_app()} is called for its side effect and returns the value of \code{shiny::runApp()} when the
#'   application is closed; \code{bys_app_object()} returns the application object (for \code{shiny::runApp()},
#'   \code{shiny::testServer()} or deployment on a server).
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @references
#' Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A General Bayesian Estimation Method of Linear-Bilinear
#' Models Applied to Plant Breeding Trials With Genotype x Environment Interaction. \emph{Journal of Agricultural,
#' Biological, and Environmental Statistics}, \bold{17}(1), 15--37. (\doi{10.1007/s13253-011-0063-9})
#' @seealso \code{\link{bayes_ammi}}, \code{\link{bayes_ammi_biplot}}, \code{\link{bys_compare}}, \code{\link{bys_interpret}}
#' @examples
#' if (interactive()) {
#'   bys_app(
#'       launch.browser = TRUE
#'     )
#' }
#' @export
bys_app <- function(launch.browser = interactive(), ...) {
  .bys_app_need()
  shiny::runApp(bys_app_object(), launch.browser = launch.browser, ...)
}

#' @rdname bys_app
#' @export
bys_app_object <- function() {
  .bys_app_need()
  shiny::shinyApp(ui = .bys_app_ui(), server = .bys_app_server)
}

.bys_app_need <- function() {
  miss <- c("shiny", "bslib", "DT")[!vapply(c("shiny", "bslib", "DT"), requireNamespace, NA, quietly = TRUE)]
  if (length(miss)) .bys_stop("bys_app", "the application needs the package(s) ", paste(miss, collapse = ", "),
                              ": install.packages(c(", paste0("\"", miss, "\"", collapse = ", "), ")).")
}

#-------------------------------------------------------------------------------------------------------------------
# look: theme, fonts, styles, Copy buttons
#-------------------------------------------------------------------------------------------------------------------
.bys_app_font <- function(family, kind = c("sans", "mono")) {
  fb <- switch(match.arg(kind),
               sans = c("system-ui", "-apple-system", "Segoe UI", "Helvetica Neue", "Arial", "sans-serif"),
               mono = c("SFMono-Regular", "Menlo", "Consolas", "Liberation Mono", "DejaVu Sans Mono", "monospace"))
  do.call(bslib::font_collection, c(list(bslib::font_google(family, local = FALSE)), as.list(fb)))
}

.bys_app_theme <- function() {
  bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#1f4e9e",
                  base_font = .bys_app_font("IBM Plex Sans", "sans"),
                  heading_font = .bys_app_font("IBM Plex Sans", "sans"),
                  code_font = .bys_app_font("IBM Plex Mono", "mono"))
}

.bys_app_head <- function() {
  shiny::tags$head(
      shiny::tags$style(shiny::HTML(paste(
          "pre.shiny-text-output{font-family:var(--bs-font-monospace);font-size:.85em;line-height:1.5;white-space:pre;"
        , "  overflow-x:auto;background:#f6f8fa;border:1px solid #dde3ea;border-radius:6px;padding:10px 12px;}"
        , "table.dataTable thead th{background:#1f4e9e;color:#fff;font-weight:600;}"
        , ".bys-bar{display:flex;gap:.4rem;justify-content:flex-end;margin-bottom:.4rem}"
        , ".dataTables_wrapper .dt-buttons{float:left;margin:0 1rem .5rem 0}"
        , ".dataTables_wrapper .dataTables_length{float:left;margin-bottom:.5rem}"
        , ".dataTables_wrapper .dataTables_filter{float:right;margin-bottom:.5rem}"
        , ".dataTables_wrapper .dataTables_length label,.dataTables_wrapper .dataTables_filter label{display:inline-flex;"
        , "  align-items:center;gap:.4rem;font-weight:400;margin:0}"
        , ".dataTables_wrapper select,.dataTables_wrapper input[type=search]{display:inline-block;width:auto}"
        , ".dataTables_wrapper .dataTables_info{float:left;padding-top:.6rem;font-size:.9em}"
        , ".dataTables_wrapper .dataTables_paginate{float:right;padding-top:.4rem}"
        , ".dataTables_wrapper table.dataTable{clear:both}"
        , ".dataTables_wrapper::after{content:'';display:table;clear:both}"
        , ".bys-plot img{max-width:100%;height:auto}")))
    , shiny::tags$script(shiny::HTML(paste(
          "$(document).on('click', '.bys-copy', function() {"
        , "  var el = document.getElementById($(this).data('target')); if (!el) return;"
        , "  var txt = el.innerText || el.textContent, b = $(this);"
        , "  var done = function() { b.find('.bys-copy-label').text('Copied'); setTimeout(function() { b.find('.bys-copy-label').text('Copy'); }, 1500); };"
        , "  if (navigator.clipboard && window.isSecureContext) { navigator.clipboard.writeText(txt).then(done); return; }"
        , "  var ta = document.createElement('textarea'); ta.value = txt; document.body.appendChild(ta); ta.select();"
        , "  document.execCommand('copy'); document.body.removeChild(ta); done();"
        , "});", sep = "\n")))
    )
}

.bys_copy_bar <- function(target) {
  shiny::div(class = "bys-bar",
             shiny::tags$button(type = "button", class = "btn btn-outline-secondary btn-sm bys-copy",
                                `data-target` = target, shiny::icon("copy"),
                                shiny::span(class = "bys-copy-label", "Copy")))
}

.bys_btn <- "btn btn-sm btn-outline-secondary me-1"

# a table for the application: Copy, CSV and Excel buttons, all rows on request, numbers rounded
.bys_dt <- function(d, title = "table") {
  d <- as.data.frame(d, check.names = FALSE)
  num <- which(vapply(d, is.numeric, NA))
  t <- DT::datatable(d, rownames = FALSE, extensions = "Buttons", style = "auto",
                     class = "compact stripe hover cell-border",
                     options = list(pageLength = 15L, scrollX = TRUE, dom = "Blfrtip",
                                    lengthMenu = list(c(15L, 50L, -1L), c("15", "50", "All")),
                                    buttons = list(list(extend = "copy", title = NULL, className = .bys_btn),
                                                   list(extend = "csv", title = title, className = .bys_btn),
                                                   list(extend = "excel", title = title, className = .bys_btn))))
  pcol  <- num[grepl("^p_|^prob|p_value|^Probability", names(d)[num])]
  whole <- setdiff(num[vapply(d[num], function(v) all(is.na(v) | v == round(v)), NA)], pcol)
  dec   <- setdiff(num, c(whole, pcol))
  if (length(pcol)) t <- DT::formatRound(t, pcol, digits = 3L)
  if (length(dec)) t <- DT::formatSignif(t, dec, digits = 5L)
  if (length(whole)) t <- DT::formatRound(t, whole, digits = 0L)
  t
}

# a ggplot as an image for renderImage(): WebP through magick when installed, PNG otherwise
.bys_plot_file <- function(p, width = 900, height = 600, res = 110) {
  png <- tempfile(fileext = ".png")
  if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(png, width = width, height = height, res = res) else
    grDevices::png(png, width = width, height = height, res = res)
  print(p)
  grDevices::dev.off()
  if (.bys_webp_ok()) {
    webp <- sub("\\.png$", ".webp", png)
    ok <- tryCatch({
      magick::image_write(magick::image_read(png), webp, format = "webp", quality = 100,
                          defines = c(`webp:lossless` = "true"))
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(list(src = webp, contentType = "image/webp", width = "100%", alt = "plot"))
  }
  list(src = png, contentType = "image/png", width = "100%", alt = "plot")
}

# a card with tabs: Interpretation, Output, Table, Plot (each with a Copy button where it is text)
.bys_result_tabs <- function(id, plot = TRUE, table = TRUE, code = FALSE, extra = NULL) {
  ns <- shiny::NS(id)
  tabs <- list(
      bslib::nav_panel("Interpretation", .bys_copy_bar(ns("int_text")), shiny::uiOutput(ns("int")),
                       shiny::tags$pre(id = ns("int_text"), class = "shiny-text-output", style = "display:none"))
    , if (plot) bslib::nav_panel("Plot", extra, shiny::div(class = "bys-plot", shiny::imageOutput(ns("plot"), height = "auto")),
                                 shiny::div(shiny::downloadButton(ns("png"), "PNG", class = "btn-sm"),
                                            shiny::downloadButton(ns("pdf"), "PDF", class = "btn-sm"),
                                            shiny::downloadButton(ns("webp"), "WebP", class = "btn-sm")))
    , if (table) bslib::nav_panel("Table", DT::DTOutput(ns("table")))
    , bslib::nav_panel("Output", .bys_copy_bar(ns("out")), shiny::verbatimTextOutput(ns("out")))
    , if (code) bslib::nav_panel("R code", .bys_copy_bar(ns("code")), shiny::verbatimTextOutput(ns("code")))
    )
  do.call(bslib::navset_card_tab, c(list(id = ns("tabs")), Filter(Negate(is.null), tabs)))
}

.bys_int_html <- function(it) {
  li <- function(s) if (length(s)) shiny::tags$ul(lapply(s, shiny::tags$li)) else NULL
  shiny::tagList(
      shiny::div(class = "alert alert-primary", shiny::strong(it$headline))
    , if (length(it$findings)) shiny::tagList(shiny::h6("Findings"), li(it$findings))
    , if (length(it$recommendations)) shiny::tagList(shiny::h6("Recommendations"), li(it$recommendations))
    , if (length(it$cautions)) shiny::tagList(shiny::h6("Cautions"), li(it$cautions))
    , shiny::tags$details(shiny::tags$summary("How this is read"), shiny::p(it$guide))
    )
}

# the server side of a result card: res() a bys_result, table_fun(res) a data frame, plot_fun(res) a ggplot
.bys_result_server <- function(id, res, table_fun = NULL, plot_fun = NULL, code_fun = NULL, hb = shiny::reactive(TRUE)) {
  force(res)
  force(table_fun)
  force(plot_fun)
  force(code_fun)
  shiny::moduleServer(id, function(input, output, session) {
    it <- shiny::reactive({
      shiny::req(res())
      bys_interpret(res(), higher_better = hb())
    })
    output$int <- shiny::renderUI(.bys_int_html(it()))
    output$int_text <- shiny::renderText(paste(format(it(), style = "text", width = 100L), collapse = "\n"))
    output$out <- shiny::renderPrint({
      shiny::req(res())
      op <- options(width = 120L)
      on.exit(options(op))
      print(res())
    })
    if (!is.null(table_fun)) output$table <- DT::renderDT({
      shiny::req(res())
      .bys_dt(table_fun(res()), id)
    })
    if (!is.null(plot_fun)) {
      pl <- shiny::reactive({
        shiny::req(res())
        plot_fun(res(), input)
      })
      output$plot <- shiny::renderImage(.bys_plot_file(pl()), deleteFile = TRUE)
      output$png <- shiny::downloadHandler(paste0(id, ".png"), function(file)
        ggplot2::ggsave(file, pl(), width = 10, height = 7, dpi = 200))
      output$pdf <- shiny::downloadHandler(paste0(id, ".pdf"), function(file)
        ggplot2::ggsave(file, pl(), width = 10, height = 7))
      output$webp <- shiny::downloadHandler(paste0(id, ".webp"), function(file) bys_save(pl(), file, width = 10, height = 7))
    }
    if (!is.null(code_fun)) output$code <- shiny::renderText(code_fun())
  })
}

#-------------------------------------------------------------------------------------------------------------------
# pages
#-------------------------------------------------------------------------------------------------------------------
.bys_app_ui <- function() {
  side <- function(...) bslib::sidebar(width = 300, ...)
  bslib::page_navbar(
      title = "baystability", id = "nav", theme = .bys_app_theme(), fillable = FALSE, header = .bys_app_head()
    , bslib::nav_panel("Home", .bys_app_home())
    , bslib::nav_panel("Data", bslib::layout_sidebar(
          sidebar = side(
              shiny::radioButtons("source", "Trial", c(`Example: cultivo2009 (12 genotypes, 25 sites)` = "cultivo",
                                                       `CSV file` = "csv"))
            , shiny::conditionalPanel("input.source == 'csv'", shiny::fileInput("file", "CSV file (one row per plot)", accept = ".csv"))
            , shiny::selectInput("col_y", "Response", NULL), shiny::selectInput("col_gen", "Genotype", NULL)
            , shiny::selectInput("col_env", "Environment", NULL), shiny::selectInput("col_rep", "Replicate", NULL)
            , shiny::checkboxInput("hb", "Higher values are better", TRUE)
            , shiny::radioButtons("prior_source", "Earlier trial as prior",
                                  c(none = "none", `Example: cultivo2008` = "cultivo", `CSV file (same columns)` = "csv"))
            , shiny::conditionalPanel("input.prior_source == 'csv'",
                                      shiny::fileInput("prior_file", "CSV file of the earlier trial", accept = ".csv")))
        , bslib::card(bslib::card_header("The trial"), .bys_copy_bar("data_summary"),
                      shiny::verbatimTextOutput("data_summary"), DT::DTOutput("data_table"))))
    , bslib::nav_panel("Model", bslib::layout_sidebar(
          sidebar = side(
              shiny::selectInput("model", "Model", c(AMMI = "AMMI", `SREG (GGE)` = "SREG"))
            , shiny::numericInput("k", "Multiplicative terms (0: all)", 2, min = 0, step = 1)
            , shiny::numericInput("n0", "Weight of the earlier trial (empty: its replicates)", NA, min = 0.01)
            , shiny::numericInput("iter", "Iterations kept per chain", 3000, min = 100, step = 500)
            , shiny::numericInput("burnin", "Burn-in", 1000, min = 0, step = 500)
            , shiny::numericInput("thin", "Thinning", 1, min = 1, step = 1)
            , shiny::numericInput("chains", "Chains", 2, min = 1, max = 6, step = 1)
            , shiny::numericInput("seed", "Seed", 1, step = 1)
            , shiny::numericInput("level", "HPD probability", 0.95, min = 0.5, max = 0.999, step = 0.01)
            , shiny::tags$details(shiny::tags$summary("Vague priors (without an earlier trial)"),
                                  shiny::numericInput("lambda_var", "Prior variance of the singular values", 1e8),
                                  shiny::numericInput("nu", "Degrees of freedom of the error variance prior", 1, min = 0.01))
            , shiny::actionButton("run", "Fit the model", class = "btn-primary w-100", icon = shiny::icon("play")))
        , .bys_result_tabs("fit", code = TRUE, extra = shiny::selectInput("fit_plot", NULL,
                           c(`traces` = "trace", `posterior histograms` = "density", `biplot of posterior means` = "biplot",
                             `main effects` = "effects")))))
    , bslib::nav_panel("Biplot", bslib::layout_sidebar(
          sidebar = side(
              shiny::checkboxGroupInput("pb", "HPD regions", c(`95%` = 0.05, `99%` = 0.01), selected = c(0.05, 0.01))
            , shiny::radioButtons("which", "Show", c(both = "both", genotypes = "genotype", environments = "environment"))
            , shiny::checkboxInput("stable", "Regions containing the origin", TRUE)
            , shiny::sliderInput("bp_thin", "Share of the draws used", 0.05, 1, 0.25, step = 0.05))
        , .bys_result_tabs("biplot")))
    , bslib::nav_panel("Model choice", bslib::layout_sidebar(
          sidebar = side(shiny::checkboxGroupInput("cmp_k", "Terms", 1:4, selected = 1:3, inline = TRUE),
                         shiny::checkboxInput("cmp_prior", "Each also without the earlier trial", FALSE),
                         shiny::actionButton("compare", "Compare", class = "btn-primary w-100"),
                         shiny::helpText("Uses the settings of the Data and Model pages; each model is fitted once."))
        , .bys_result_tabs("compare")))
    , bslib::nav_panel("References", bslib::card(bslib::card_header("Cite the methods you use"),
                                                 .bys_copy_bar("refs"), shiny::verbatimTextOutput("refs")))
    )
}

.bys_app_home <- function() {
  card <- function(title, text) bslib::card(bslib::card_header(title), shiny::p(text))
  shiny::tagList(
      shiny::h3("Bayesian AMMI and SREG analysis with an earlier trial as prior")
    , shiny::p("The Gibbs sampler of Perez-Elizalde, Jarquin and Crossa (2012): the estimates of an earlier trial of",
               "the same genotypes and environments become conjugate priors for the new trial, weighted by the",
               "earlier trial's replicates. Without an earlier trial the priors are vague. Every parameter has a",
               "posterior distribution, so every statement about genotypes and environments comes with a probability.")
    , bslib::layout_column_wrap(width = 1 / 3,
          card("1. Data", "The example cultivo2009 or your CSV file, one row per plot, and optionally the earlier trial.")
        , card("2. Model", "AMMI or SREG (GGE), number of terms, weight of the earlier trial, chains; traces and diagnostics.")
        , card("3. Biplot", "Highest posterior density regions: stable genotypes and groups of environments.")
        , card("4. Model choice", "WAIC and DIC of several numbers of terms, with and without the earlier trial.")
        , card("5. References", "The articles behind every method, to cite with the results.")
        , card("Every page", "Interpretation, output as R prints it, tables and plots, all with Copy or download."))
    )
}

#-------------------------------------------------------------------------------------------------------------------
# server
#-------------------------------------------------------------------------------------------------------------------
.bys_app_example <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "baystability", envir = e)
  e[[name]]
}

.bys_app_server <- function(input, output, session) {
  raw <- shiny::reactive({
    if (identical(input$source, "csv")) {
      shiny::req(input$file)
      data.table::fread(input$file$datapath, data.table = FALSE)
    } else .bys_app_example("cultivo2009")
  })
  prior <- shiny::reactive({
    switch(if (is.null(input$prior_source)) "none" else input$prior_source,
           none = NULL,
           cultivo = .bys_app_example("cultivo2008"),
           csv = {
             shiny::req(input$prior_file)
             data.table::fread(input$prior_file$datapath, data.table = FALSE)
           })
  })
  shiny::observeEvent(raw(), {
    nm <- names(raw())
    pick <- function(pat, default) { h <- grep(pat, nm, ignore.case = TRUE, value = TRUE); if (length(h)) h[1L] else default }
    shiny::updateSelectInput(session, "col_y", choices = nm, selected = pick("^y$|yield", nm[length(nm)]))
    shiny::updateSelectInput(session, "col_gen", choices = nm, selected = pick("gen|entry|variety", nm[1L]))
    shiny::updateSelectInput(session, "col_env", choices = nm, selected = pick("env|site|loc", nm[2L]))
    shiny::updateSelectInput(session, "col_rep", choices = nm, selected = pick("rep|block", nm[3L]))
  })
  output$data_summary <- shiny::renderPrint({
    d <- raw()
    shiny::req(input$col_y %in% names(d), input$col_gen %in% names(d), input$col_env %in% names(d))
    N <- table(d[[input$col_gen]], d[[input$col_env]])
    cat(sprintf("%d plots; %d genotypes; %d environments; replicates per cell from %d to %d\n",
                nrow(d), nrow(N), ncol(N), min(N), max(N)))
    p <- prior()
    if (is.null(p)) cat("No earlier trial: vague priors.\n\n") else
      cat(sprintf("Earlier trial as prior: %d plots\n\n", nrow(p)))
    print(summary(d[[input$col_y]]))
  })
  output$data_table <- DT::renderDT(.bys_dt(raw(), "data"))

  fit <- shiny::reactiveVal(NULL)
  args <- shiny::reactive(list(
      .data = raw(), .y = input$col_y, .gen = input$col_gen, .env = input$col_env, .rep = input$col_rep
    , .nIter = as.integer(input$iter), .prior_data = prior()
    , .n0 = if (is.null(input$n0) || is.na(input$n0)) NULL else input$n0
    , .model = input$model, .k = if (is.null(input$k) || is.na(input$k) || input$k < 1) NULL else as.integer(input$k)
    , .burnin = as.integer(input$burnin), .thin = as.integer(input$thin), .chains = as.integer(input$chains)
    , .priors = bys_priors(lambda_var = input$lambda_var, nu = input$nu), .level = input$level
    , .seed = as.integer(input$seed), .align = TRUE, .waic = TRUE))
  shiny::observeEvent(input$run, {
    r <- tryCatch(shiny::withProgress(message = "Sampling", value = 0.5, do.call(bayes_ammi, args())),
                  error = function(e) e)
    if (inherits(r, "error")) {
      shiny::showNotification(conditionMessage(r), type = "error", duration = 10)
      return()
    }
    fit(r)
  })
  hb <- shiny::reactive(isTRUE(input$hb))
  code <- shiny::reactive(.bys_app_code(args(), input$source, input$prior_source))
  .bys_result_server("fit", fit, function(x) x$summary, function(x, i) plot(x, which = input$fit_plot), code, hb)
  bp <- shiny::reactive({
    shiny::req(fit(), length(input$pb) > 0L)
    if (fit()$settings$k < 2L) return(NULL)
    .bys_biplots(fit(), fn = "bayes_ammi_biplot", burnin = 0, thin = input$bp_thin, pb = as.numeric(input$pb), plot_stable = isTRUE(input$stable),
            plot_unstable = TRUE, .which = input$which)
  })
  .bys_result_server("biplot", bp, function(x) x$biplot_data, function(x, i) x$biplot, hb = hb)
  cmp <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$compare, {
    a <- args()
    ks <- as.integer(input$cmp_k)
    with_prior <- if (!is.null(a$.prior_data) && isTRUE(input$cmp_prior)) c(TRUE, FALSE) else !is.null(a$.prior_data)
    grid <- expand.grid(k = ks, prior = with_prior)
    p0 <- a$.prior_data
    fits <- list()
    shiny::withProgress(message = "Fitting models", value = 0, {
      for (r in seq_len(nrow(grid))) {
        a$.k <- grid$k[r]
        a[".prior_data"] <- list(if (grid$prior[r]) p0 else NULL)
        nm <- sprintf("%s%d%s", a$.model, a$.k, if (grid$prior[r]) "-prior" else "")
        fits[[nm]] <- tryCatch(do.call(bayes_ammi, a), error = function(e) NULL)
        shiny::incProgress(1 / nrow(grid))
      }
    })
    fits <- Filter(Negate(is.null), fits)
    if (length(fits)) cmp(do.call(bys_compare, fits))
  })
  .bys_result_server("compare", cmp, function(x) x$table, function(x, i) plot(x), hb = hb)
  output$refs <- shiny::renderPrint({
    op <- options(width = 110L)
    on.exit(options(op))
    print(bys_cite())
  })
}

# the R code of the current fit, every argument on its own line
.bys_app_code <- function(a, source, prior_source) {
  csv <- identical(source, "csv")
  pdat <- switch(if (is.null(prior_source)) "none" else prior_source, none = "NULL", cultivo = "cultivo2008",
                 csv = "earlier")
  vals <- c(.data = if (csv) "trial" else "cultivo2009", .y = a$.y, .gen = a$.gen, .env = a$.env, .rep = a$.rep,
            .nIter = a$.nIter, .prior_data = pdat, .n0 = if (is.null(a$.n0)) "NULL" else format(a$.n0),
            .model = sprintf("\"%s\"", a$.model), .k = if (is.null(a$.k)) "NULL" else a$.k, .burnin = a$.burnin,
            .thin = a$.thin, .chains = a$.chains,
            .priors = sprintf("bys_priors(lambda_var = %s, nu = %s)", format(a$.priors$lambda_var), format(a$.priors$nu)),
            .level = a$.level, .seed = a$.seed, .align = "TRUE", .waic = "TRUE")
  w <- max(nchar(names(vals)))
  body <- paste0("    ", c(" ", rep(",", length(vals) - 1L)), " ", formatC(names(vals), width = -w), " = ", vals,
                 collapse = "\n")
  load <- c(if (csv) "trial <- data.table::fread(\"trial.csv\")" else "data(cultivo2009)",
            switch(pdat, cultivo2008 = "data(cultivo2008)", earlier = "earlier <- data.table::fread(\"earlier.csv\")",
                   NULL))
  paste0("library(baystability)\n", paste0(load, "\n", collapse = ""), "fit <-\n  bayes_ammi(\n", body,
         "\n    )\nfit\nbys_interpret(\n    x             = fit\n  , higher_better = TRUE\n  , top           = 5\n  )\n")
}
