# Per-trait check-threshold veto: flag (and optionally exclude) crosses whose mid-parent for a
# trait is on the wrong side of a breeder-chosen check line, on a GEBV or phenotype basis. A
# candidate-eligibility filter upstream of the optimizer (index objective untouched). See
# docs/design/2026-07-25-trait-check-threshold-design.md.

# Build a validated per-trait check spec (one check per trait). `direction` NA is filled from
# `trait_direction` (increase -> reject if the mid-parent is BELOW the check; decrease -> ABOVE).
ng_trait_check_spec <- function(trait, check, direction = NA, basis = "gebv",
                                trait_direction = NULL) {
  trait <- as.character(trait); check <- as.character(check)
  n <- length(trait)
  if (!n) ng_stop("ng_trait_check_spec needs at least one trait")
  if (length(check) != n) ng_stop("check must be one per trait")
  direction <- tolower(as.character(rep_len(direction, n)))
  basis <- tolower(as.character(rep_len(basis, n)))
  if (any(!basis %in% c("gebv", "phenotype")))
    ng_stop("basis must be 'gebv' or 'phenotype'")
  # resolve NA directions from trait_direction
  need <- is.na(direction) | !nzchar(direction) | direction == "na"
  if (any(need)) {
    if (is.null(trait_direction))
      ng_stop("direction is NA and no trait_direction supplied to resolve it for: ",
              paste(trait[need], collapse = ", "))
    td <- tolower(as.character(trait_direction[trait[need]]))
    if (anyNA(td)) ng_stop("trait_direction has no entry for: ",
                           paste(trait[need][is.na(td)], collapse = ", "))
    direction[need] <- ifelse(td == "increase", "below",
                       ifelse(td == "decrease", "above", NA_character_))
  }
  if (any(!direction %in% c("above", "below")))
    ng_stop("direction must resolve to 'above' or 'below'")
  dup <- unique(trait[duplicated(trait)])
  if (length(dup))
    ng_stop("ng_trait_check_spec: one check per trait (v1); duplicate trait(s): ",
            paste(dup, collapse = ", "))
  data.frame(trait = trait, check = check, reject_if = direction, basis = basis,
             stringsAsFactors = FALSE)
}

# Apply trait checks to a scoring data.frame: compute mid-parent trait values, test against
# check thresholds, flag violations, optionally exclude them, emit diagnostics.
# 'trait_values' is a named list keyed by trait: list(yield = list(gebv = <named numeric>, phenotype = NULL), ...)
ng_apply_trait_checks <- function(scores, spec, trait_values, exclude = FALSE) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(scores)))
    ng_stop("scores must have parent1 and parent2 columns")
  n <- nrow(scores)
  p1 <- as.character(scores$parent1); p2 <- as.character(scores$parent2)
  violation <- matrix(NA, nrow = n, ncol = nrow(spec),
                      dimnames = list(NULL, spec$trait))
  n_not_evaluable <- 0L
  for (k in seq_len(nrow(spec))) {
    tr <- spec$trait[[k]]; ck <- spec$check[[k]]; b <- spec$basis[[k]]; dir <- spec$reject_if[[k]]
    vals <- trait_values[[tr]][[b]]
    if (is.null(vals)) {                      # basis unavailable for this trait
      scores[[paste0(tr, "_check_value")]] <- NA_real_
      scores[[paste0(tr, "_check_violation")]] <- NA
      n_not_evaluable <- n_not_evaluable + n
      next
    }
    thr <- suppressWarnings(as.numeric(vals[ck]))
    mp <- 0.5 * (suppressWarnings(as.numeric(vals[p1])) + suppressWarnings(as.numeric(vals[p2])))
    v <- if (identical(dir, "above")) mp > thr else mp < thr   # NA-safe: NA stays NA
    if (!is.finite(thr)) v <- rep(NA, n)
    violation[, k] <- v
    n_not_evaluable <- n_not_evaluable + sum(is.na(v))
    scores[[paste0(tr, "_check_value")]] <- thr
    scores[[paste0(tr, "_check_violation")]] <- v
  }
  fail <- apply(violation, 1L, function(r) any(r %in% TRUE))       # TRUE only where a real violation (NA never counts)
  scores$threshold_violation <- vapply(seq_len(n), function(i) {
    hit <- spec$trait[which(violation[i, ] %in% TRUE)]
    paste(hit, collapse = ",")
  }, character(1))
  scores$threshold_ok <- !fail
  diag <- list(
    active = spec,
    n_by_trait = as.list(colSums(violation == TRUE, na.rm = TRUE)),
    n_flagged = sum(fail),
    n_not_evaluable = as.integer(n_not_evaluable),
    n_excluded = if (isTRUE(exclude)) sum(fail) else 0L)
  if (isTRUE(exclude)) scores <- scores[!fail, , drop = FALSE]
  rownames(scores) <- NULL
  attr(scores, "trait_check_diagnostics") <- diag
  scores
}
