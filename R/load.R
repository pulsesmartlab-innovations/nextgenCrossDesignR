ng_load <- function(path = NULL, use_cpp = TRUE, verbose = FALSE) {
  if (is.null(path)) {
    candidates <- c(
      Sys.getenv("NG_REPO_ROOT", unset = NA_character_),
      getwd(),
      file.path(getwd(), "nextgen_cross_design"),
      file.path(getwd(), ".."),
      file.path(getwd(), "..", "nextgen_cross_design")
    )
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
    candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
    is_root <- vapply(candidates, function(candidate) {
      desc <- file.path(candidate, "DESCRIPTION")
      if (!file.exists(file.path(candidate, "R", "load.R")) || !file.exists(desc)) {
        return(FALSE)
      }
      any(grepl("^Package:\\s*nextgenCrossDesign\\s*$", readLines(desc, n = 20L, warn = FALSE)))
    }, logical(1L))
    roots <- candidates[is_root]
    if (!length(roots)) stop("Could not locate nextgenCrossDesign package root", call. = FALSE)
    path <- roots[[1L]]
  }
  r_dir <- file.path(path, "R")
  files <- list.files(r_dir, pattern = "^[0-9].*[.]R$", full.names = TRUE)
  files <- files[basename(files) != "load.R"]
  files <- sort(files)
  # 0.30.0: source exactly what `R CMD build` would put in the tarball.
  #
  # This glob is a directory listing, so it picked up ANY file in R/ whose name
  # starts with a digit -- including files that are not part of the package. The
  # working tree carries eleven untracked legacy prototype scripts in R/
  # (01_relationships.R, 03_ld.R, 05_pmv.R, 07_optimizers.R, 10_marker_effects.R
  # and the rest), and their names collide with the package's own numbering, so
  # every ng_load() session sourced a DIFFERENT file set from the installed
  # package -- silently, and after the real sources, so a colliding definition
  # would have won. That makes every locally run test a test of something other
  # than the package.
  #
  # THE RULE: exclude what .Rbuildignore excludes. Dev-load should match the
  # BUILT package, and .Rbuildignore is the single declaration of what the built
  # package contains, so deriving the load set from it makes the two definitions
  # the same definition rather than two lists that can drift apart. It is also
  # the only one of the candidate rules that needs no external tooling: git may
  # be absent (an unpacked tarball, a CI image, a copied directory), and an
  # explicit manifest in this file would be a third list to keep in sync. A file
  # added to R/ and not build-ignored is loaded, exactly as it would be shipped.
  #
  # Two files the tarball DOES contain are still not sourced here, for reasons
  # that have nothing to do with .Rbuildignore: `load.R` itself (this file, the
  # bootstrap), and `RcppExports.R`, whose .Call() stubs dispatch into a compiled
  # DLL that a source-loaded session does not have -- ng_load_cpp() supplies
  # those symbols via Rcpp::sourceCpp() instead. Both are outside the `^[0-9]`
  # glob, so the rule below never sees them.
  files <- ng_load_apply_rbuildignore(files, path)
  # A file we DO intend to load must be readable. On a OneDrive/iCloud/DFS
  # working copy a cloud-only placeholder opens successfully and reads as ZERO
  # bytes while file.size() still reports the real size, so sys.source() on one
  # is a silent no-op: the package would appear to load and then be missing
  # whatever that file defined. Refuse loudly instead, naming the file.
  for (f in files) ng_load_require_materialised(f)
  for (f in files) sys.source(f, envir = .GlobalEnv)
  if (isTRUE(use_cpp)) ng_load_cpp(file.path(path, "src", "ng_kernels.cpp"), verbose = verbose)
  invisible(TRUE)
}

# Apply the package's own .Rbuildignore to a vector of absolute paths, using the
# same rule `R CMD build` does: each non-blank, non-comment line is a Perl
# regular expression, matched UNANCHORED against the path relative to the
# package root (see tools:::.Rbuildignore / R-exts "Building package tarballs").
# A missing .Rbuildignore excludes nothing.
#
# One deliberate deviation: `R CMD build` treats EVERY non-empty line as a
# pattern, `#` comments included, so a comment that happened to match a path
# would silently drop that file from the tarball. Comment lines are skipped here
# instead, because a prose comment accidentally matching a source path is a bug
# in either direction and this is the safe one. Verified on this package's own
# .Rbuildignore that the two readings give an IDENTICAL verdict on every file in
# R/ (11 dropped either way).
ng_load_apply_rbuildignore <- function(files, path) {
  if (!length(files)) return(files)
  ignore_file <- file.path(path, ".Rbuildignore")
  if (!file.exists(ignore_file)) return(files)
  patterns <- tryCatch(readLines(ignore_file, warn = FALSE), error = function(e) character(0))
  patterns <- trimws(patterns)
  patterns <- patterns[nzchar(patterns) & !startsWith(patterns, "#")]
  if (!length(patterns)) return(files)
  root <- normalizePath(path, winslash = "/", mustWork = FALSE)
  abs <- normalizePath(files, winslash = "/", mustWork = FALSE)
  prefix <- paste0(sub("/$", "", root), "/")
  # Paths are always inside <root>/R here; the basename fallback keeps the
  # function total if a caller ever hands it something else.
  rel <- ifelse(startsWith(abs, prefix),
                substring(abs, nchar(prefix) + 1L),
                paste0("R/", basename(abs)))
  drop <- rep(FALSE, length(rel))
  for (p in patterns) {
    hit <- tryCatch(grepl(p, rel, perl = TRUE), error = function(e) rep(FALSE, length(rel)))
    drop <- drop | hit
  }
  files[!drop]
}

# Refuse a source file that cannot actually be read. `file.size()` is not enough:
# a cloud-storage placeholder reports its real size and then yields no bytes.
ng_load_require_materialised <- function(f) {
  size <- file.size(f)
  if (is.na(size)) {
    stop("ng_load(): cannot stat the package source file ", f,
         " -- it is listed in R/ and is not build-ignored, so it must be loadable.",
         call. = FALSE)
  }
  if (size == 0) return(invisible(TRUE))   # genuinely empty file: nothing to read
  con <- tryCatch(file(f, "rb"), error = function(e) e)
  if (inherits(con, "error")) {
    stop("ng_load(): cannot open the package source file ", f, " -- ",
         conditionMessage(con), call. = FALSE)
  }
  on.exit(close(con), add = TRUE)
  got <- tryCatch(length(readBin(con, "raw", n = 1L)), error = function(e) 0L)
  if (got == 0L) {
    stop("ng_load(): the package source file ", f, " reports ", size,
         " bytes but reads as EMPTY. This is the signature of a cloud-storage ",
         "placeholder (OneDrive/iCloud 'Files On-Demand') whose contents are not ",
         "materialised on this machine. Sourcing it would silently define ",
         "nothing and leave the loaded package incomplete. Download/pin the ",
         "file locally, or remove it from R/, and re-run.",
         call. = FALSE)
  }
  invisible(TRUE)
}

ng_load_cpp <- function(cpp_file, verbose = FALSE) {
  skip <- tolower(trimws(Sys.getenv("NGCD_SKIP_CPP", unset = "0")))
  if (skip %in% c("1", "true", "yes", "y")) return(FALSE)
  if (!requireNamespace("Rcpp", quietly = TRUE)) return(FALSE)
  if (!file.exists(cpp_file)) return(FALSE)
  ok <- tryCatch({
    Rcpp::sourceCpp(cpp_file, verbose = verbose, rebuild = FALSE)
    TRUE
  }, error = function(e) {
    if (isTRUE(verbose)) message("C++ acceleration unavailable: ", conditionMessage(e))
    FALSE
  })
  ok
}
