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
  for (f in sort(files)) sys.source(f, envir = .GlobalEnv)
  if (isTRUE(use_cpp)) ng_load_cpp(file.path(path, "src", "ng_kernels.cpp"), verbose = verbose)
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
