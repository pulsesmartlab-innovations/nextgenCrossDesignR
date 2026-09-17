ng_stop <- function(..., call. = FALSE) stop(paste0(...), call. = call.)

# Null-coalescing operator. `%||%` only entered base R in 4.4.0, but DESCRIPTION declares
# R (>= 4.1.0), and the operator is used across the package -- including R/03_metrics.R, which
# sits in ng_score_crosses()'s core scoring path. Without this definition the package fails
# immediately on R 4.1-4.3 and is silently fine on 4.4+, so the breakage is invisible to anyone
# developing on a current R. Defining it here shadows the base version identically on 4.4+.
`%||%` <- function(a, b) if (is.null(a)) b else a


ng_as_numeric_matrix <- function(x, name = deparse(substitute(x))) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(rownames(x))) ng_stop(name, " must have row names")
  if (is.null(colnames(x))) colnames(x) <- paste0("m", seq_len(ncol(x)))
  x
}

ng_check_same_ids <- function(x, ids, object_name = "matrix") {
  if (is.null(rownames(x))) ng_stop(object_name, " must have row names")
  miss <- setdiff(ids, rownames(x))
  if (length(miss)) ng_stop(object_name, " is missing ", length(miss), " requested ids")
  x[ids, , drop = FALSE]
}

ng_make_pairs <- function(ids, include_self = FALSE) {
  ids <- as.character(ids)
  out_n <- if (include_self) length(ids) * (length(ids) + 1) / 2 else length(ids) * (length(ids) - 1) / 2
  p1 <- character(out_n)
  p2 <- character(out_n)
  k <- 1L
  for (i in seq_along(ids)) {
    j0 <- if (include_self) i else i + 1L
    if (j0 <= length(ids)) {
      for (j in j0:length(ids)) {
        p1[k] <- ids[i]
        p2[k] <- ids[j]
        k <- k + 1L
      }
    }
  }
  data.frame(parent1 = p1, parent2 = p2, stringsAsFactors = FALSE)
}

ng_selection_intensity <- function(prop_selected, n_progeny = NULL) {
  p <- as.numeric(prop_selected)
  if (!is.finite(p) || p <= 0 || p >= 1) ng_stop("prop_selected must be in (0, 1)")
  if (is.null(n_progeny) || !is.finite(n_progeny[[1L]]) || n_progeny[[1L]] <= 0) {
    z <- stats::qnorm(1 - p)
    return(stats::dnorm(z) / p)
  }
  # Finite-population correction. Use the expected mean of the top k order
  # statistics of an N(0,1) sample of size n, via Blom (1958) plotting
  # positions: E[Z_{(j:n)}] ~= Phi^{-1}((j - 3/8) / (n + 1/4)). This is
  # deterministic, fast, and converges to phi(z)/p as n -> Inf.
  n <- max(1L, as.integer(round(n_progeny[[1L]])))
  k <- max(1L, as.integer(round(n * p)))
  k <- min(k, n)
  j <- seq.int(n - k + 1L, n)
  mean(stats::qnorm((j - 0.375) / (n + 0.25)))
}

ng_scale01 <- function(x, bigger_is_better = TRUE) {
  x <- as.numeric(x)
  r <- range(x[is.finite(x)], na.rm = TRUE)
  if (!all(is.finite(r)) || abs(diff(r)) < .Machine$double.eps) return(rep(0.5, length(x)))
  y <- (x - r[1]) / diff(r)
  if (!bigger_is_better) y <- 1 - y
  y
}

ng_standardize <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

ng_match_vector <- function(x, ids, name = deparse(substitute(x))) {
  if (is.null(names(x))) {
    if (length(x) != length(ids)) ng_stop(name, " must be named or have length ids")
    names(x) <- ids
  }
  miss <- setdiff(ids, names(x))
  if (length(miss)) ng_stop(name, " is missing ", length(miss), " ids")
  as.numeric(x[ids])
}

ng_parent_kinship <- function(geno, method = c("vanraden", "yang")) {
  method <- match.arg(method)
  geno <- ng_as_numeric_matrix(geno, "geno")
  p <- colMeans(geno, na.rm = TRUE) / 2
  X <- sweep(geno, 2, 2 * p, "-")
  X[!is.finite(X)] <- 0
  if (identical(method, "yang")) {
    # Yang / GCTA: standardize each marker to unit variance (equal weight), then average.
    s <- sqrt(2 * p * (1 - p)); s[!is.finite(s) | s <= 0] <- Inf   # monomorphic -> X/Inf = 0
    Z <- sweep(X, 2, s, "/"); Z[!is.finite(Z)] <- 0
    denom <- sum(is.finite(s)); if (denom <= 0) denom <- ncol(geno)
    K <- tcrossprod(Z) / denom
  } else {
    # VanRaden (default): single overall allele-frequency scaling.
    denom <- sum(2 * p * (1 - p), na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- ncol(geno)
    K <- tcrossprod(X) / denom
  }
  rownames(K) <- rownames(geno)
  colnames(K) <- rownames(geno)
  attr(K, "method") <- method
  attr(K, "relationship_scale") <- "additive_relationship"
  attr(K, "coancestry_divisor") <- 2
  K
}

ng_pair_key <- function(parent1, parent2) paste(parent1, parent2, sep = "||")

# Evaluate `expr` under a fixed RNG seed, restoring the caller's global RNG state
# afterwards so seeded internal randomness (optimizers, benchmark fixtures) is
# reproducible without perturbing the surrounding stream. `expr` is a promise and is
# only forced AFTER the seed is set. A non-finite/NULL seed evaluates `expr` as-is.
ng_with_rng_seed <- function(seed, expr) {
  s <- suppressWarnings(as.numeric(seed[[1L]]))
  if (!length(s) || !is.finite(s)) return(force(expr))
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  set.seed(as.integer(s))
  force(expr)
}

# Per-parent diagnostics for the inbred-parent assumption used by the
# DH/RIL recombination-variance kernel. Dosage is considered inbred at a
# marker when min(d, ploidy - d) <= tolerance (e.g. 0 or 2 +/- noise for
# diploids). Returns the row indices and IDs of parents whose fraction of
# heterozygous-looking markers exceeds `fraction_tolerance`, plus the worst
# observed fraction.
ng_audit_inbred_dosage <- function(geno, ploidy = 2,
                                   tolerance = 0.05,
                                   fraction_tolerance = 0.02) {
  if (!is.matrix(geno)) geno <- as.matrix(geno)
  d <- as.numeric(geno)
  dim(d) <- dim(geno)
  het_dist <- pmin(abs(d), abs(ploidy - d))
  het_dist[!is.finite(het_dist)] <- 0
  het_marker <- het_dist > tolerance
  frac <- rowMeans(het_marker, na.rm = TRUE)
  frac[!is.finite(frac)] <- 0
  bad <- which(frac > fraction_tolerance)
  ids <- if (!is.null(rownames(geno))) rownames(geno) else as.character(seq_len(nrow(geno)))
  list(
    violators = ids[bad],
    violators_index = bad,
    fraction = frac,
    max_fraction = if (length(frac)) max(frac, na.rm = TRUE) else 0,
    tolerance = tolerance,
    fraction_tolerance = fraction_tolerance
  )
}

# Optional external packages (SimpleMating, genomicMateSelectR) are used only for
# opt-in benchmark comparisons and are intentionally NOT declared dependencies:
# they are not on CRAN, and the package is fully functional without them (callers
# fall back to native methods). These helpers reach them through variable
# arguments so R CMD check does not treat them as undeclared imports, and so pak
# never tries to resolve them during dependency setup.
ng_has_optional_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

ng_optional_pkg_fun <- function(pkg, fun) {
  get(fun, envir = asNamespace(pkg), mode = "function")
}

# Set ctx fields WITHOUT losing NULLs. `ctx$key <- NULL` deletes the key, so a stage that later
# does list2env(ctx, environment()) and reads `key` by bare name gets an unbound-variable error
# instead of NULL. Bracket assignment with a length-1 list stores the NULL and keeps the key.
#
# The first formal is named `.ctx` (not `ctx`) even though every call site passes it
# positionally as `ng_ctx_put(ctx, ...)`. A bare `ctx` formal placed before `...` is a real
# partial-matching trap here: R matches named `...` args against un-consumed formals BEFORE
# positional matching, so a value literally named `c` (a one-letter prefix of "ctx") would bind
# to the `ctx` formal instead of falling into `...`, silently displacing the real ctx list. The
# leading dot makes that collision impossible for any plausible field name.
ng_ctx_put <- function(.ctx, ...) {
  vals <- list(...)
  nms <- names(vals)
  if (is.null(nms) || any(!nzchar(nms))) ng_stop("ng_ctx_put requires named values")
  ctx <- .ctx
  for (nm in nms) ctx[nm] <- list(vals[[nm]])
  ctx
}


# ---- Order-independent per-trait RNG seeds ---------------------------------------------------
#
# A trait's scientific result must never depend on WHERE the breeder happened to list it in the
# direction file. Before 0.28.0 the per-trait loops derived their seeds from the trait's ROW
# POSITION (`seed + i - 1L`, `seed + j`), so permuting the direction file gave a trait a
# different CV fold split -> a different ridge lambda -> different marker effects -> a different
# progeny variance. The observed spread on a 24-parent / 40-marker panel was seven orders of
# magnitude in `<trait>_vpm`, and it changed the selected crossing plan.
#
# `ng_name_hash32()` maps a trait NAME to a stable non-negative integer. It is written out in
# plain R arithmetic on the string's UTF-8 bytes precisely so that it does NOT depend on any
# hashing internal that R, a platform or a locale is free to change: `enc2utf8()` fixes the byte
# sequence, and every intermediate (max 131 * (2^31 - 2) + 255 ~= 2.8e11) is exactly
# representable in a double, so the result is bit-identical on every platform and R version.
ng_name_hash32 <- function(x) {
  modulus <- 2147483647            # 2^31 - 1, keeps the result inside integer range
  vapply(as.character(x), function(s) {
    if (is.na(s)) s <- "NA"
    bytes <- as.integer(charToRaw(enc2utf8(s)))
    h <- 0
    for (b in bytes) h <- (h * 131 + b) %% modulus
    h
  }, numeric(1L), USE.NAMES = FALSE)
}

# Seed for a per-trait RNG stream: identity-derived, never position-derived. `salt` separates
# independent uses of the same trait within one run (e.g. the posterior draw stream vs anything
# else keyed off the same base seed). Returns a value in [0, 2^31 - 2], safe for set.seed().
ng_trait_rng_seed <- function(base_seed, trait, salt = 0L) {
  modulus <- 2147483647
  base <- suppressWarnings(as.numeric(base_seed))
  if (length(base) != 1L || !is.finite(base)) base <- 0
  salt <- suppressWarnings(as.numeric(salt))
  if (length(salt) != 1L || !is.finite(salt)) salt <- 0
  h <- ng_name_hash32(trait)
  as.integer((base %% modulus + salt %% modulus + h) %% modulus)
}

# ---------------------------------------------------------------------------------------
# Sizing a batch by MEMORY rather than by cores.
#
# A batch of single-trait analyses is bounded by RAM, not CPU. Each concurrent job holds
# its own genotype copy and, when the marker count is small enough for the dense beta
# covariance, its own Sigma_beta -- ~288 MB at 6000 markers, per worker. Sizing a fan-out
# from detectCores() ignores that, and in a container it is wrong twice over:
# detectCores() reports the HOST's cpus, and the cgroup memory cap is invisible to base R.
# The failure mode is an OOM kill deep into a multi-hour run, losing everything unwritten.
# ---------------------------------------------------------------------------------------

# Bytes of memory this process can actually expect to use, with an attribute naming how it
# was determined so a surprising worker count can be explained afterwards. Ordered most- to
# least-authoritative: a cgroup cap binds regardless of what the host reports.
ng_available_memory_bytes <- function() {
  out <- function(bytes, source) {
    bytes <- suppressWarnings(as.numeric(bytes))
    if (!length(bytes) || !is.finite(bytes) || bytes <= 0) return(NULL)
    structure(bytes, source = source)
  }
  read1 <- function(path) {
    if (!file.exists(path)) return(NULL)
    tryCatch(trimws(readLines(path, n = 1L, warn = FALSE)), error = function(e) NULL)
  }
  shell1 <- function(cmd, args) {
    tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
             error = function(e) NULL)
  }

  # cgroup v2 -- the container case. "max" means uncapped, so fall through rather than
  # treating the literal string as a number.
  v2 <- read1("/sys/fs/cgroup/memory.max")
  if (!is.null(v2) && !identical(v2, "max")) {
    used <- suppressWarnings(as.numeric(read1("/sys/fs/cgroup/memory.current") %||% NA))
    cap <- suppressWarnings(as.numeric(v2))
    free <- if (is.finite(used)) cap - used else cap
    r <- out(free, "cgroup v2 (memory.max)")
    if (!is.null(r)) return(r)
  }
  # cgroup v1 -- "unlimited" is encoded as a near-INT64_MAX sentinel, not as a word.
  v1 <- read1("/sys/fs/cgroup/memory/memory.limit_in_bytes")
  if (!is.null(v1)) {
    cap <- suppressWarnings(as.numeric(v1))
    if (is.finite(cap) && cap < 2^62) {
      used <- suppressWarnings(as.numeric(
        read1("/sys/fs/cgroup/memory/memory.usage_in_bytes") %||% NA))
      free <- if (is.finite(used)) cap - used else cap
      r <- out(free, "cgroup v1 (memory.limit_in_bytes)")
      if (!is.null(r)) return(r)
    }
  }
  # Linux host: MemAvailable is the kernel's own estimate of what a new workload can get.
  if (file.exists("/proc/meminfo")) {
    mi <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character(0))
    line <- grep("^MemAvailable:", mi, value = TRUE)
    if (length(line)) {
      kb <- suppressWarnings(as.numeric(sub("^MemAvailable:\\s*([0-9]+).*$", "\\1", line[[1L]])))
      r <- out(kb * 1024, "/proc/meminfo MemAvailable")
      if (!is.null(r)) return(r)
    }
  }
  # macOS: free + inactive + speculative pages. Inactive pages are reclaimable, so counting
  # only "free" understates what is available by a large margin on a warm machine.
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    vm <- shell1("vm_stat", character(0))
    if (length(vm)) {
      psize <- suppressWarnings(as.numeric(sub(".*page size of ([0-9]+) bytes.*", "\\1", vm[[1L]])))
      if (!is.finite(psize)) psize <- 4096
      pages <- function(label) {
        ln <- grep(label, vm, value = TRUE, fixed = TRUE)
        if (!length(ln)) return(0)
        suppressWarnings(as.numeric(gsub("[^0-9]", "", ln[[1L]])))
      }
      free <- sum(vapply(c("Pages free:", "Pages inactive:", "Pages speculative:"),
                         function(l) { p <- pages(l); if (is.finite(p)) p else 0 }, numeric(1)))
      r <- out(free * psize, "vm_stat (free + inactive + speculative)")
      if (!is.null(r)) return(r)
    }
    r <- out(suppressWarnings(as.numeric(shell1("sysctl", c("-n", "hw.memsize"))[[1L]])) / 2,
             "sysctl hw.memsize / 2")
    if (!is.null(r)) return(r)
  }
  if (identical(.Platform$OS.type, "windows")) {
    ps <- shell1("powershell", c("-NoProfile", "-Command",
      "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"))
    if (length(ps)) {
      kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", ps[[length(ps)]])))
      r <- out(kb * 1024, "Win32_OperatingSystem FreePhysicalMemory")
      if (!is.null(r)) return(r)
    }
  }
  # Nothing could be measured. Claim little rather than much: under-reporting costs speed,
  # over-reporting costs the whole run.
  structure(2 * 1024^3, source = "fallback (could not measure; assuming 2 GB)")
}

# Peak resident bytes one single-trait job should be expected to need. Deliberately a
# coarse upper-ish estimate built from the terms that actually dominate -- an estimate that
# is precise about small terms and wrong about Sigma_beta would be worse than useless.
ng_batch_job_bytes <- function(n_parents, n_markers, dense_beta_cov = TRUE) {
  n <- as.numeric(n_parents); m <- as.numeric(n_markers)
  dbl <- 8
  geno <- n * m * dbl
  # The ridge fit rbinds the training augmentation onto the parents, so a second copy of the
  # genotype matrix is live at the same time as the first.
  fit_copy <- geno
  # Sigma_beta is m x m and is the single largest allocation whenever it is built at all.
  sigma <- if (isTRUE(dense_beta_cov)) m * m * dbl else 0
  # The scored cross table: every unordered pair, with of the order of a dozen numeric
  # columns per trait.
  pairs <- (n * (n - 1) / 2) * 12 * dbl
  # A daemon with R and this package loaded, before any data.
  base_process <- 150 * 1024^2
  geno + fit_copy + sigma + pairs + base_process
}

# How many jobs to run at once. Memory is a second ceiling on top of cores, never a
# replacement for it, and the answer is never zero: running slowly beats refusing to start.
ng_batch_worker_count <- function(n_jobs, per_job_bytes, memory_budget_bytes = NULL,
                                  cores = NULL) {
  n_jobs <- max(1L, as.integer(n_jobs))
  if (is.null(cores)) {
    cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (!is.finite(cores) || cores < 1) cores <- 1L
    cores <- max(1L, as.integer(cores) - 1L)
  }
  cores <- max(1L, as.integer(cores))
  if (is.null(memory_budget_bytes)) memory_budget_bytes <- ng_available_memory_bytes()
  by_mem <- suppressWarnings(
    as.integer(floor(as.numeric(memory_budget_bytes) / max(1, as.numeric(per_job_bytes)))))
  if (!length(by_mem) || is.na(by_mem)) by_mem <- cores
  n <- max(1L, min(cores, by_mem, n_jobs))
  basis <- if (n == 1L && by_mem < 1L) {
    "memory: the budget does not fit even one job, so jobs run one at a time"
  } else if (identical(n, by_mem) && by_mem <= min(cores, n_jobs)) {
    sprintf("memory: %s budget / %s per job", ng_bytes_label(memory_budget_bytes),
            ng_bytes_label(per_job_bytes))
  } else if (identical(n, n_jobs)) {
    "jobs: fewer jobs than the machine could run at once"
  } else {
    sprintf("cores: %d usable", cores)
  }
  structure(n, basis = basis)
}

ng_bytes_label <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return("unknown")
  u <- c("B", "KB", "MB", "GB", "TB"); i <- 1L
  while (x >= 1024 && i < length(u)) { x <- x / 1024; i <- i + 1L }
  sprintf("%.1f %s", x, u[[i]])
}
