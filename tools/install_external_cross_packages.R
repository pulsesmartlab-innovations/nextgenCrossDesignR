.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
dir.create(".Rlib", showWarnings = FALSE)

repos <- getOption("repos")
repos["CRAN"] <- "https://cloud.r-project.org"
options(repos = repos)

cran_pkgs <- c("PopVar", "GenomicMating", "remotes")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = ".Rlib", dependencies = TRUE)
  }
}

if (!requireNamespace("SimpleMating", quietly = TRUE) &&
    requireNamespace("remotes", quietly = TRUE)) {
  remotes::install_github(
    "Resende-Lab/SimpleMating",
    lib = ".Rlib",
    dependencies = TRUE,
    upgrade = "never"
  )
}

pkgs <- c("PopVar", "genomicMateSelectR", "GenomicMating", "SimpleMating")
print(data.frame(pkg = pkgs, installed = sapply(pkgs, requireNamespace, quietly = TRUE)))
