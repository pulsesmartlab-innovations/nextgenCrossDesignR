# Never hand the shell a binary the platform cannot run.
#
# ng_alphamate_default_executable() probed only two names, "AlphaMate.exe" and
# "AlphaMate", and on non-Windows it listed the bare name first but kept the .exe as a
# fallback. The bundled binaries are actually named AlphaMate.exe (Windows PE32+) and
# AlphaMate_Unix (ELF), so on macOS and Linux the bare name never matched, the probe
# fell through to the Windows executable, and the run died with
#
#     AlphaMate failed with exit code 126: ... cannot execute binary file
#
# 126 means "found, but not executable" -- the least informative way to report "there is
# no binary for your platform". Worse, on Linux, where AlphaMate_Unix DOES run, the
# picker still chose the Windows .exe because AlphaMate_Unix was never a candidate. The
# package is required to behave identically on Windows, macOS and Linux; silently
# selecting a foreign-architecture binary is the opposite of that.
#
# A Windows PE cannot execute on Unix under any circumstance, so offering it there is
# never useful: it converts a clear "not found" into a confusing exec failure.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
old_env <- Sys.getenv("NG_ALPHAMATE_EXE", unset = NA_character_)
Sys.unsetenv("NG_ALPHAMATE_EXE")
on.exit(if (!is.na(old_env)) Sys.setenv(NG_ALPHAMATE_EXE = old_env), add = TRUE)

mk <- function(files) {
  d <- tempfile("am_"); b <- file.path(d, "external", "AlphaMate", "binaries")
  dir.create(b, recursive = TRUE)
  for (f in files) writeLines("stub", file.path(b, f))
  d
}

is_windows <- identical(.Platform$OS.type, "windows")

# ---- 1. the Unix-named binary is discovered at all -------------------------
setwd(mk("AlphaMate_Unix"))
got <- ng_alphamate_default_executable()
if (!is_windows) {
  stopifnot(grepl("AlphaMate_Unix", got, fixed = TRUE))
  stopifnot(file.exists(got))
}

# ---- 2. a Windows .exe is never chosen on a Unix platform ------------------
# This is the case that produced exit 126: the ONLY file present is the Windows
# binary. The honest outcome is "no usable binary", not a path that cannot exec.
setwd(mk("AlphaMate.exe"))
got_exe_only <- ng_alphamate_default_executable()
if (!is_windows) {
  if (file.exists(got_exe_only) && grepl("[.]exe$", got_exe_only)) {
    stop("ng_alphamate_default_executable() returned a Windows .exe on a Unix platform: ",
         got_exe_only, call. = FALSE)
  }
}

# ---- 3. with both present, each platform gets its own -----------------------
setwd(mk(c("AlphaMate.exe", "AlphaMate_Unix")))
both <- ng_alphamate_default_executable()
if (is_windows) {
  stopifnot(grepl("AlphaMate[.]exe$", both))
} else {
  stopifnot(grepl("AlphaMate_Unix", both, fixed = TRUE))
}

# ---- 4. an explicit override still wins -------------------------------------
# A user who built their own binary must not be second-guessed by name matching.
setwd(mk("AlphaMate.exe"))
Sys.setenv(NG_ALPHAMATE_EXE = "/somewhere/custom/AlphaMate")
stopifnot(identical(ng_alphamate_default_executable(), "/somewhere/custom/AlphaMate"))
Sys.unsetenv("NG_ALPHAMATE_EXE")

setwd(old_wd)
cat("alphamate_binary_matches_the_platform: PASS\n")
