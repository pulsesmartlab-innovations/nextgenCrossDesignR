\name{nextgenCrossDesign-package}
\alias{nextgenCrossDesign-package}
\alias{nextgenCrossDesign}
\title{NextGen Cross Design Backend}
\description{
Backend package for genomic cross prediction, mate allocation, multi-trait
selection, calibration, external baseline comparison, crop validation, and
polyploid model-family experiments.
}
\details{
The primary supported production path is diploid DH/RIL cross design. Polyploid
functions are exposed with explicit model-family scope and validation caveats.
The native \code{var_complex}, SimpleMating-style, and AlphaMate-style methods
are implemented as package components, with optional exact external calls
available only when those tools are installed or configured. Data preflight is
part of the package
workflow: \code{ng_preflight_input_tables()} can flag ordinary import issues,
screen putative duplicate genotype profiles, remove redundant duplicate parents
from the analysis tables, and return \code{qc$cleaned_tables} for downstream
scoring and allocation.
}
\keyword{package}
