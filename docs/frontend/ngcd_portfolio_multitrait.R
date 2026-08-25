# ngcd_portfolio_multitrait.R ---------------------------------------------------------------
#
# Drop-in support for MULTI-TRAIT "Portfolio & risk" in the NextGenCrossDesign Shiny frontend.
# Backend contract: nextgenCrossDesign >= 0.19.0 (main @ f3444a2).
#
# WHAT CHANGED IN THE BACKEND
#   Multi-trait runs used to return NONE of the portfolio columns -- the runner annotated them
#   only for single-trait runs. They are now emitted for every run, with the two axes resolved
#   on the SELECTION INDEX instead of on one trait:
#       cross_level  = w'm          index of mid-parent GEBVs
#       cross_upside = sqrt(w'Sw)   index SD within the family, exact cross-trait covariance
#   Same column names, same types, same quadrant labels, same tertile risk bins as single-trait.
#
# IF YOU CHANGE NOTHING, the existing panel will already render for multi-trait runs. The three
# things below are what stop it being MISLEADING. (1) is not optional.
#
#   (1) BADGE `portfolio_basis == "linearized_rank_index"`. For the rank-based index methods
#       (auto / weighted / threshold -- and `auto` promotes to `weighted` whenever trait weights
#       are present, so this is the COMMON case) there is no linear index in genetic units; the
#       axes come from reinterpreting the trait weights as standardized-unit coefficients. The
#       quadrant is internally consistent but is NOT a decomposition of `multi_trait_score`.
#       Never present it as "why this cross was chosen".
#   (2) AXIS LABELS lose their units. Single-trait level is a mid-parent GEBV in trait units;
#       index level is unitless. Label by role and put the weights in a drawer.
#   (3) CONFIDENCE TOOLTIP. `confidence_method` is now `midparent_pev_index[_partial]`. The index
#       PEV is a block-diagonal approximation (per-trait ridge fits are independent, so there is
#       no cross-trait estimation-error covariance). Ranking input, not a calibrated interval.
#
# Everything in section 1 is plain R with no Shiny dependency, so it can be unit-tested directly.
# -------------------------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a


## 1. Pure adapters ==========================================================================

#' Is there a portfolio to draw at all?
#' Columns are absent when the backend could not resolve an index (no coefficients, missing
#' per-trait columns). Treat that as the empty state, not as an error.
ngfe_has_portfolio <- function(crosses) {
  need <- c("cross_level", "cross_upside", "portfolio_profile")
  is.data.frame(crosses) && nrow(crosses) > 0L &&
    all(need %in% names(crosses)) && any(is.finite(crosses$cross_level))
}

ngfe_is_multitrait <- function(diag) identical(diag$basis %||% "", "multi_trait_index")

#' Axis labels. Single-trait axes carry the trait's units; index axes do not.
ngfe_portfolio_axes <- function(diag, trait_label = NULL) {
  if (ngfe_is_multitrait(diag)) {
    list(level  = "Index level (weighted mid-parent GEBV)",
         upside = "Index spread within family (SD)")
  } else {
    lab <- trait_label %||% "trait"
    list(level  = sprintf("Mid-parent GEBV (%s)", lab),
         upside = sprintf("Within-family SD (%s)", lab))
  }
}

#' The required badge. Returns NULL when no caveat applies, so callers can just test for it.
ngfe_portfolio_badge <- function(diag) {
  if (!identical(diag$portfolio_basis %||% "", "linearized_rank_index")) return(NULL)
  list(tone   = "warning",
       label  = "Indicative axes",
       detail = diag$portfolio_basis_note %||%
         paste("This index combines rank-normalized traits, so the level and spread axes are",
               "indicative rather than a breakdown of the score the crosses were ranked on."))
}

#' All run-level warnings the panel must surface, in priority order.
#' Each element is list(tone, label, detail).
ngfe_portfolio_warnings <- function(diag) {
  out <- list()
  b <- ngfe_portfolio_badge(diag)
  if (!is.null(b)) out[[length(out) + 1L]] <- b

  # Index confidence dominated by one trait -> risk_bin is a single-trait statement.
  note <- diag$pev_concentration_note %||% NA_character_
  if (!is.na(note)) {
    out[[length(out) + 1L]] <- list(
      tone = "warning", label = "Confidence is driven by one trait", detail = note)
  }

  dis <- unlist(diag$risk_disproportionate_traits %||% list())
  if (length(dis)) {
    out[[length(out) + 1L]] <- list(
      tone  = "info",
      label = sprintf("More risk than gain: %s", paste(dis, collapse = ", ")),
      detail = paste("These traits contribute more of the index prediction error than of the",
                     "index spread. Consider dropping them from the index, or improving how",
                     "they are phenotyped."))
  }
  out
}

#' Per-trait index weights for the "how was this scored" drawer.
#'
#' TWO THINGS TO GET RIGHT WHEN DISPLAYING THESE.
#'   - A minimized trait carries a NEGATIVE weight. Render it as direction, not a minus sign.
#'   - MAGNITUDE IS NOT COMPARABLE ACROSS TRAITS. These are raw-unit weights,
#'     w_k = coef_k * sign_k / scale_k, so a trait whose GEBVs barely vary gets a huge weight
#'     purely because its scale is small. A real example from the shipped contract file:
#'     yield scale 3.62 -> weight 0.138, disease scale 0.000141 -> weight -3540.7. The ratio of
#'     weights (25667) is just the inverse ratio of scales (25670) -- it does NOT mean disease is
#'     25000x more important. Never render these as a bar chart or "importance" ranking.
#'     The comparable quantities are the shares in ngfe_index_traits(); lead with those.
ngfe_index_weights <- function(diag) {
  w <- diag$index_weights
  if (is.null(w) || !length(w)) return(NULL)
  v <- unlist(w)
  data.frame(Trait  = names(v),
             Effect = ifelse(unname(v) >= 0, "raises index", "lowers index"),
             `Raw weight` = signif(unname(v), 4),   # diagnostic only -- see note above
             check.names = FALSE, stringsAsFactors = FALSE)
}

#' Per-trait contribution table: what each trait actually buys, and what it costs in certainty.
#' Spread share may be negative -- an antagonistic trait genuinely REMOVES spread from the index.
ngfe_index_traits <- function(diag) {
  it <- diag$index_traits
  if (is.null(it) || !length(it)) return(NULL)
  pct <- function(x) if (is.null(x) || !is.finite(x)) NA_real_ else round(100 * x, 1)
  do.call(rbind, lapply(it, function(r) data.frame(
    Trait            = r$trait,
    Direction        = r$direction,
    Weight           = round(r$weight, 4),
    Reliability      = round(r$marker_effect_reliability %||% NA_real_, 3),
    `Spread %`       = pct(r$mean_variance_share),
    `Error %`        = pct(r$mean_pev_share),
    check.names = FALSE, stringsAsFactors = FALSE)))
}

#' Confidence tooltip copy, keyed on how confidence was actually resolved.
ngfe_confidence_note <- function(diag) {
  switch(diag$confidence_method %||% "",
    midparent_pev_index_partial = paste(
      "Ranked on the index prediction-error variance, combined across traits as if the",
      "per-trait marker-effect fits were independent. The merit also carries a variance term",
      "this does not cover, so read it as a ranking, not an interval."),
    midparent_pev_index = paste(
      "Ranked on the index prediction-error variance, combined across traits as if the",
      "per-trait marker-effect fits were independent."),
    midparent_pev_partial = paste(
      "Ranked on the mid-parent prediction-error variance. The merit also carries a variance",
      "term this does not cover."),
    midparent_pev = "Ranked on the mid-parent prediction-error variance.",
    reliability   = "Per-cross confidence is unavailable for this run, so risk bins are not shown.",
    NULL)
}

#' Vectorized per-cross risk driver label: "driven by protein (62% of index error)".
#' Returns NA for single-trait runs, which have no components to attribute to.
ngfe_risk_driver_label <- function(crosses) {
  if (!all(c("risk_driver_trait", "risk_driver_share") %in% names(crosses)))
    return(rep(NA_character_, nrow(crosses)))
  ifelse(is.na(crosses$risk_driver_trait), NA_character_,
         sprintf("driven by %s (%.0f%% of index error)",
                 crosses$risk_driver_trait, 100 * crosses$risk_driver_share))
}

#' Human-readable quadrant labels. Keep the backend factor levels as the source of truth.
ngfe_profile_labels <- function() c(
  breakthrough = "Breakthrough", workhorse = "Workhorse",
  long_shot    = "Long shot",    deprioritize = "Deprioritize")


## 2. Shiny module ===========================================================================
## Wire with:  ngcdPortfolioUI("portfolio")
##             ngcdPortfolioServer("portfolio", run = reactive(rv$result))
## where `run` returns the list from ng_run_cross_prediction().

ngcdPortfolioUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("warnings")),
    shiny::fluidRow(
      shiny::column(8, plotly::plotlyOutput(ns("scatter"), height = "460px")),
      shiny::column(4,
        shiny::h5("Index weights"),
        shiny::tableOutput(ns("weights")),
        shiny::h5("Trait contributions"),
        shiny::tableOutput(ns("traits")),
        shiny::helpText(shiny::textOutput(ns("confnote")))
      )
    )
  )
}

#' @param axis_x which quantity goes on x. The existing single-trait tab plots upside on x;
#'   keep whichever your app already uses so breeders are not re-oriented mid-release.
ngcdPortfolioServer <- function(id, run, axis_x = c("upside", "level")) {
  axis_x <- match.arg(axis_x)
  shiny::moduleServer(id, function(input, output, session) {

    crosses <- shiny::reactive({
      r <- run(); shiny::req(r)
      r$selected_crosses %||% NULL
    })
    diagnostics <- shiny::reactive({
      r <- run(); shiny::req(r)
      r$priority_risk_diagnostics %||% list()
    })

    output$warnings <- shiny::renderUI({
      w <- ngfe_portfolio_warnings(diagnostics())
      if (!length(w)) return(NULL)
      shiny::tagList(lapply(w, function(x) shiny::div(
        class = if (identical(x$tone, "warning")) "alert alert-warning" else "alert alert-info",
        shiny::tags$strong(x$label), shiny::tags$br(), shiny::tags$small(x$detail))))
    })

    output$scatter <- plotly::renderPlotly({
      cr <- crosses(); d <- diagnostics()
      shiny::validate(shiny::need(ngfe_has_portfolio(cr),
        "No portfolio available for this run."))

      ax  <- ngfe_portfolio_axes(d, trait_label = d$trait_label)
      lab <- ngfe_profile_labels()
      drv <- ngfe_risk_driver_label(cr)

      hover <- paste0(
        cr$parent1, " x ", cr$parent2,
        "<br>", lab[as.character(cr$portfolio_profile)],
        "<br>level: ",  signif(cr$cross_level, 4),
        "<br>spread: ", signif(cr$cross_upside, 4),
        ifelse(is.na(cr$cross_confidence), "",
               paste0("<br>confidence: ", sprintf("%.2f", cr$cross_confidence))),
        ifelse(is.na(drv), "", paste0("<br>", drv)))

      xs <- if (axis_x == "upside") cr$cross_upside else cr$cross_level
      ys <- if (axis_x == "upside") cr$cross_level  else cr$cross_upside
      xlab <- if (axis_x == "upside") ax$upside else ax$level
      ylab <- if (axis_x == "upside") ax$level  else ax$upside

      p <- plotly::plot_ly(
        x = xs, y = ys, type = "scatter", mode = "markers",
        color = cr$risk_bin, colors = c(low = "#2E7D5B", med = "#C9922F", high = "#B4522F"),
        text = hover, hoverinfo = "text",
        marker = list(size = 10, line = list(width = 1, color = "rgba(0,0,0,.25)")))

      # Median cuts: these ARE the quadrant boundaries the backend used.
      plotly::layout(p,
        xaxis = list(title = xlab), yaxis = list(title = ylab),
        legend = list(title = list(text = "Risk")),
        shapes = list(
          list(type = "line", x0 = stats::median(xs, na.rm = TRUE),
               x1 = stats::median(xs, na.rm = TRUE), y0 = min(ys, na.rm = TRUE),
               y1 = max(ys, na.rm = TRUE), line = list(dash = "dot", width = 1)),
          list(type = "line", y0 = stats::median(ys, na.rm = TRUE),
               y1 = stats::median(ys, na.rm = TRUE), x0 = min(xs, na.rm = TRUE),
               x1 = max(xs, na.rm = TRUE), line = list(dash = "dot", width = 1))))
    })

    output$weights  <- shiny::renderTable(ngfe_index_weights(diagnostics()), na = "-")
    output$traits   <- shiny::renderTable(ngfe_index_traits(diagnostics()),  na = "-")
    output$confnote <- shiny::renderText(ngfe_confidence_note(diagnostics()) %||% "")
  })
}


## 3. Rules to enforce in the UI =============================================================
#
# WITHIN-RUN ONLY. `cross_confidence` is a min-max normalization and `risk_bin` is tertiles of
# the full post-filter candidate pool, copied unchanged to selected rows. Never carry these
# across runs, and never phrase a count of high-risk crosses as an absolute statement about plan
# quality. Ties are kept together, so bins may deliberately be unequal.
#
# NEGATIVE IS MEANINGFUL, twice over:
#   - a minimized trait has a negative index weight (higher raw value lowers the index);
#   - a trait antagonistic to the rest of the index has a negative spread share, because it
#     genuinely removes spread from the index. "Protein narrows the index spread by 15%."
#
# SINGLE-TRAIT RUNS still work through all of the above: `portfolio_basis` is "single_trait",
# no badge is produced, the risk-driver columns are absent and ngfe_risk_driver_label() returns
# NA, and ngfe_index_weights()/ngfe_index_traits() return NULL. Any exhaustive switch on
# `portfolio_basis` must handle "single_trait".
#
# VERIFY AGAINST:
#   1 trait                     -> basis "single_trait",  no badge, trait-unit axes
#   2 traits + weight           -> index_method "weighted",       basis "linearized_rank_index",
#                                  BADGE SHOWN
#   2 traits + economic_weight  -> index_method "economic_index", basis "linear_index",
#                                  no badge, portfolio_basis_note is NA
#   any trait set to decrease   -> negative index weight, drawer reads "lowers index"
#   a poorly phenotyped trait   -> named by risk_driver_trait; listed in
#                                  risk_disproportionate_traits; if it takes >90% of the index
#                                  error, pev_concentration_note is populated
