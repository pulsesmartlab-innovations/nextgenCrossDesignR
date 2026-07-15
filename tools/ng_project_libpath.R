# Standalone helper (no package dependencies): prepend a project-local library
# directory to .libPaths() ONLY when it was built for the current platform.
#
# Why this exists: the project-local `.Rlib` holds binary packages. If a `.Rlib`
# built on one OS is moved to another (e.g. a Windows library with libs/x64/*.dll
# copied onto macOS/Linux), blindly prepending it makes R resolve packages like
# lpSolve/Rcpp from that library and then fail to load their native code
# ("shared object 'lpSolve.so' not found"), SHADOWING the correct system/user
# packages. `requireNamespace()` then returns FALSE and callers behave as if the
# dependency were missing (e.g. the OCS MIP silently degrades / errors).
#
# The platform check is intentionally CHEAP: it reads the `Built:` platform of a
# single installed package's DESCRIPTION rather than scanning the whole tree
# (a deep scan of a cloud-synced library, e.g. OneDrive, can stall for minutes
# while files are hydrated).

# Return TRUE if `lib` looks built for the current OS (or if that can't be
# determined, in which case we do not block).
ng_project_lib_platform_ok <- function(lib) {
  pkgs <- list.dirs(lib, recursive = FALSE, full.names = TRUE)
  if (!length(pkgs)) return(TRUE)  # empty library: harmless to prepend
  cur_os <- tolower(R.version$os)         # e.g. "darwin23.0", "mingw32", "linux-gnu"
  cur_win <- .Platform$OS.type == "windows"
  for (p in pkgs) {
    desc <- file.path(p, "DESCRIPTION")
    if (!file.exists(desc)) next
    built <- grep("^Built:", readLines(desc, warn = FALSE), value = TRUE)
    if (!length(built)) next
    b <- tolower(built[[1L]])
    lib_win <- grepl("mingw|w64|windows", b)
    lib_mac <- grepl("darwin|apple", b)
    lib_lin <- grepl("linux", b)
    if (!lib_win && !lib_mac && !lib_lin) return(TRUE)  # unknown token: don't block
    if (cur_win) return(lib_win)
    if (grepl("darwin", cur_os)) return(lib_mac)
    if (grepl("linux", cur_os)) return(lib_lin)
    return(TRUE)  # unrecognized current OS: don't block
  }
  TRUE  # no readable DESCRIPTION with a Built field: treat as safe
}

ng_prepend_project_lib <- function(lib) {
  lib <- normalizePath(lib, mustWork = FALSE)
  if (!dir.exists(lib)) return(invisible(FALSE))
  compatible <- isTRUE(tryCatch(ng_project_lib_platform_ok(lib), error = function(e) TRUE))
  if (compatible) {
    .libPaths(c(lib, .libPaths()))
  } else {
    message(sprintf(
      "ng_prepend_project_lib: skipping project library '%s' (built for a different platform); using system/user libraries instead.",
      lib))
  }
  invisible(compatible)
}
