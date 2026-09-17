# Worker count must be derived from MEMORY, not from core count.
#
# Each concurrent job holds its own copy of the genotype matrix and, when the marker
# count is tractable enough for the dense beta covariance (R/39: want_beta_cov, markers
# <= 6000), its own Sigma_beta. At 6000 markers that single matrix is 6000^2 * 8 bytes
# = ~288 MB per worker, before anything else. Ten workers is ~3 GB of Sigma_beta alone.
#
# ng_run_cp_parallel_cores() (R/39:434) takes detectCores(logical = FALSE) - 1 and knows
# nothing about memory. In a container that is doubly wrong: detectCores() reports the
# HOST's cpus, and the cgroup memory cap is invisible to base R entirely. A batch sized
# that way is how a 30-hour run gets OOM-killed at hour 29 with nothing written.
#
# So the sizing must read the real budget and divide. These are the properties that make
# it safe; the exact worker count on any given machine is not one of them.
# helper_load.R resolves the package root and prepends the project library ONLY when it was
# built for this platform -- a blunt .libPaths(".Rlib") prepend silently shadows working
# packages with ones built for another OS, and a mirai daemon inherits the parent's library
# paths, so it fails inside a worker with an opaque "shared object not found".
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- the budget itself -----------------------------------------------------------------
avail <- ng_available_memory_bytes()
stopifnot(is.numeric(avail), length(avail) == 1L, is.finite(avail), avail > 0)
# Sanity: any machine that can run this package has more than 64 MB and less than 100 TB.
# Wide on purpose -- this asserts the helper returns BYTES, not KB (a units slip is the
# realistic failure, and it would silently oversubscribe by 1024x).
stopifnot(avail > 64 * 1024^2, avail < 100 * 1024^4)
stopifnot(is.character(attr(avail, "source")), nzchar(attr(avail, "source")))

# --- per-job footprint ------------------------------------------------------------------
small <- ng_batch_job_bytes(n_parents = 200L, n_markers = 1000L, dense_beta_cov = FALSE)
big   <- ng_batch_job_bytes(n_parents = 200L, n_markers = 6000L, dense_beta_cov = TRUE)
stopifnot(is.numeric(small), is.numeric(big), small > 0, big > small)
# Sigma_beta must dominate when it is dense: 6000^2 * 8 is ~288 MB, and no other term in a
# 200 x 6000 problem comes close. If this fails, the estimate has stopped modelling the
# thing that actually drives the memory.
stopifnot(big > 6000^2 * 8)
# ... and must NOT be charged when the run does not build it.
stopifnot(ng_batch_job_bytes(200L, 6000L, dense_beta_cov = FALSE) < big / 2)

# --- the sizing decision ----------------------------------------------------------------
# A budget that fits exactly three jobs must yield three workers, never more, however many
# cores the box has.
w <- ng_batch_worker_count(n_jobs = 17L, per_job_bytes = 100 * 1024^2,
                           memory_budget_bytes = 350 * 1024^2, cores = 64L)
# as.integer() strips the `basis` attribute the count deliberately carries.
stopifnot(identical(as.integer(w), 3L))

# Never more workers than there is work.
stopifnot(identical(as.integer(ng_batch_worker_count(2L, 1024, 1024^4, cores = 64L)), 2L))

# Cores still cap it -- memory is a second ceiling, not a replacement for the first.
stopifnot(identical(as.integer(ng_batch_worker_count(17L, 1024, 1024^4, cores = 4L)), 4L))

# A budget too small for even one job must still run, serially. Refusing to start is worse
# than running slowly: the breeder would be left with no result at all.
stopifnot(identical(as.integer(ng_batch_worker_count(17L, 10 * 1024^3, 1024^2, cores = 64L)), 1L))

# The decision has to be explicable afterwards -- a batch that quietly ran 2-wide on a
# 32-core box looks like a bug unless it says why.
basis <- attr(ng_batch_worker_count(17L, 100 * 1024^2, 350 * 1024^2, cores = 64L), "basis")
stopifnot(is.character(basis), grepl("memory", basis))

cat("PASS: batch_workers_fit_the_memory_budget\n")
