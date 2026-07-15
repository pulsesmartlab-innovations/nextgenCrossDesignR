#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

class NgDisjointSet {
 public:
  explicit NgDisjointSet(int n) : parent_(n), rank_(n, 0) {
    for (int i = 0; i < n; ++i) parent_[static_cast<std::size_t>(i)] = i;
  }

  int find(int x) {
    int px = parent_[static_cast<std::size_t>(x)];
    if (px != x) {
      parent_[static_cast<std::size_t>(x)] = find(px);
    }
    return parent_[static_cast<std::size_t>(x)];
  }

  void unite(int a, int b) {
    int ra = find(a);
    int rb = find(b);
    if (ra == rb) return;
    unsigned char rank_a = rank_[static_cast<std::size_t>(ra)];
    unsigned char rank_b = rank_[static_cast<std::size_t>(rb)];
    if (rank_a < rank_b) std::swap(ra, rb);
    parent_[static_cast<std::size_t>(rb)] = ra;
    if (rank_a == rank_b) ++rank_[static_cast<std::size_t>(ra)];
  }

 private:
  std::vector<int> parent_;
  std::vector<unsigned char> rank_;
};

static double ng_ld_pair_r2_cpp(NumericMatrix geno, int col_a, int col_b) {
  const int n = geno.nrow();
  double sum_a = 0.0;
  double sum_b = 0.0;
  double sum_aa = 0.0;
  double sum_bb = 0.0;
  double sum_ab = 0.0;
  int count = 0;
  for (int row = 0; row < n; ++row) {
    const double a = geno(row, col_a);
    const double b = geno(row, col_b);
    if (!R_finite(a) || !R_finite(b)) continue;
    sum_a += a;
    sum_b += b;
    sum_aa += a * a;
    sum_bb += b * b;
    sum_ab += a * b;
    ++count;
  }
  if (count < 2) return NA_REAL;
  const double dn = static_cast<double>(count);
  const double denom = (dn * sum_aa - sum_a * sum_a) * (dn * sum_bb - sum_b * sum_b);
  if (!R_finite(denom) || denom <= 0.0) return NA_REAL;
  const double cov_num = dn * sum_ab - sum_a * sum_b;
  double r2 = (cov_num * cov_num) / denom;
  if (!R_finite(r2)) return NA_REAL;
  if (r2 < 0.0) r2 = 0.0;
  if (r2 > 1.0) r2 = 1.0;
  return r2;
}

// [[Rcpp::export]]
IntegerVector ng_ld_prune_graph_cpp(NumericMatrix geno,
                                    int window = 100,
                                    double r2_threshold = 0.9,
                                    double maf_threshold = 0.01,
                                    double ploidy = 2.0) {
  const int n = geno.nrow();
  const int m = geno.ncol();
  if (n < 1 || m < 1) return IntegerVector(0);
  if (window < 1) window = 100;
  if (!R_finite(r2_threshold) || r2_threshold < 0.0) r2_threshold = 0.9;
  if (r2_threshold > 1.0) r2_threshold = 1.0;
  if (!R_finite(maf_threshold) || maf_threshold < 0.0) maf_threshold = 0.01;
  if (maf_threshold > 0.5) maf_threshold = 0.5;
  if (!R_finite(ploidy) || ploidy <= 0.0) ploidy = 2.0;

  std::vector<double> maf(static_cast<std::size_t>(m), NA_REAL);
  std::vector<unsigned char> pass_maf(static_cast<std::size_t>(m), 0);
  int pass_count = 0;
  for (int col = 0; col < m; ++col) {
    double sum = 0.0;
    int count = 0;
    for (int row = 0; row < n; ++row) {
      const double value = geno(row, col);
      if (!R_finite(value)) continue;
      sum += value;
      ++count;
    }
    if (count > 0) {
      const double allele_frequency = (sum / static_cast<double>(count)) / ploidy;
      double marker_maf = std::min(allele_frequency, 1.0 - allele_frequency);
      if (!R_finite(marker_maf)) marker_maf = NA_REAL;
      maf[static_cast<std::size_t>(col)] = marker_maf;
      if (R_finite(marker_maf) && marker_maf >= maf_threshold) {
        pass_maf[static_cast<std::size_t>(col)] = 1;
        ++pass_count;
      }
    }
  }
  if (pass_count == 0) return IntegerVector(0);
  if (m == 1) return pass_maf[0] ? IntegerVector::create(1) : IntegerVector(0);

  NgDisjointSet dsu(m);
  for (int i = 0; i < m - 1; ++i) {
    if (!pass_maf[static_cast<std::size_t>(i)]) continue;
    const int j_end = std::min(m - 1, i + window);
    for (int j = i + 1; j <= j_end; ++j) {
      if (!pass_maf[static_cast<std::size_t>(j)]) continue;
      const double r2 = ng_ld_pair_r2_cpp(geno, i, j);
      if (R_finite(r2) && r2 > r2_threshold) dsu.unite(i, j);
    }
    if ((i & 127) == 0) Rcpp::checkUserInterrupt();
  }

  std::vector<int> best(static_cast<std::size_t>(m), -1);
  for (int col = 0; col < m; ++col) {
    if (!pass_maf[static_cast<std::size_t>(col)]) continue;
    const int root = dsu.find(col);
    const int current = best[static_cast<std::size_t>(root)];
    if (current < 0 ||
        maf[static_cast<std::size_t>(col)] > maf[static_cast<std::size_t>(current)] ||
        (maf[static_cast<std::size_t>(col)] == maf[static_cast<std::size_t>(current)] && col < current)) {
      best[static_cast<std::size_t>(root)] = col;
    }
  }

  std::vector<int> keep;
  keep.reserve(static_cast<std::size_t>(pass_count));
  for (int col = 0; col < m; ++col) {
    if (!pass_maf[static_cast<std::size_t>(col)]) continue;
    const int root = dsu.find(col);
    if (best[static_cast<std::size_t>(root)] == col) keep.push_back(col + 1);
  }
  return wrap(keep);
}

// [[Rcpp::export]]
DataFrame ng_dh_recomb_pairs_cpp(NumericMatrix geno,
                                 NumericVector beta,
                                 NumericVector beta_var,
                                 IntegerVector chr,
                                 NumericVector pos_cm,
                                 CharacterVector ids,
                                 CharacterVector pair_parent1,
                                 CharacterVector pair_parent2,
                                 double window_cm = -1.0) {
  const int n_pairs = pair_parent1.size();
  const int m = geno.ncol();
  std::unordered_map<std::string, int> id_index;
  for (int i = 0; i < ids.size(); ++i) id_index[as<std::string>(ids[i])] = i;

  NumericVector dh_recomb_var(n_pairs);
  NumericVector dh_pmv_var(n_pairs);

  for (int r = 0; r < n_pairs; ++r) {
    const int i = id_index[as<std::string>(pair_parent1[r])];
    const int j = id_index[as<std::string>(pair_parent2[r])];
    double v = 0.0;
    double pmv = 0.0;

    double carry = 0.0;
    int prev_chr = -2147483647;
    double prev_pos = 0.0;
    for (int k = 0; k < m; ++k) {
      if (k == 0 || chr[k] != prev_chr) {
        carry = 0.0;
      } else {
        const double dist_prev = std::fabs(pos_cm[k] - prev_pos);
        carry *= std::exp(-2.0 * dist_prev / 100.0);
      }
      const double xi = geno(i, k);
      const double xj = geno(j, k);
      if (R_finite(xi) && R_finite(xj)) {
        const double dk = 0.5 * (xi - xj);
        const double ak = dk * beta[k];
        const double e_beta2 = beta[k] * beta[k] + std::max(0.0, beta_var[k]);
        v += ak * ak + 2.0 * ak * carry;
        pmv += dk * dk * e_beta2 + 2.0 * ak * carry;
        carry += ak;
      }
      prev_chr = chr[k];
      prev_pos = pos_cm[k];
    }
    dh_recomb_var[r] = std::max(0.0, v);
    dh_pmv_var[r] = std::max(0.0, pmv);
  }

  return DataFrame::create(
    Named("dh_recomb_var") = dh_recomb_var,
    Named("dh_pmv_var") = dh_pmv_var
  );
}

// ---- Full off-diagonal posterior PMV (genomicMateSelectR formulation) -----
//
// Per pair: a = d .* beta, where d_k = 0.5 (geno[p1, k] - geno[p2, k]).
//   VPM      = a' R a
//   PMV_diag = VPM + sum_k d_k^2 * Var(beta_k)         (diagonal posterior)
//   PMV_full = VPM + d' (R o Sigma_beta) d             (full off-diag posterior)
//
// R_decay and R_had_Sigma (= R .* Sigma_beta) are precomputed once on the R
// side and passed in. The inner matrix-vector multiplies are done with
// cache-friendly outer-over-columns, inner-over-rows loops (NumericMatrix is
// column-major in R, matching this access pattern). On a 5K-marker fixture
// this runs in seconds vs minutes for the R reference.
//
// [[Rcpp::export]]
DataFrame ng_dh_recomb_pairs_full_posterior_cpp(NumericMatrix geno,
                                                NumericVector beta,
                                                NumericVector beta_var_diag,
                                                NumericMatrix R_decay,
                                                NumericMatrix R_had_Sigma,
                                                IntegerVector pair_p1_zero,
                                                IntegerVector pair_p2_zero) {
  const int m = beta.size();
  const int n_pairs = pair_p1_zero.size();
  if (R_decay.nrow() != m || R_decay.ncol() != m) {
    Rcpp::stop("R_decay must be m x m matching length(beta)");
  }
  if (R_had_Sigma.nrow() != m || R_had_Sigma.ncol() != m) {
    Rcpp::stop("R_had_Sigma must be m x m matching length(beta)");
  }
  if (beta_var_diag.size() != m) {
    Rcpp::stop("beta_var_diag must have length m");
  }

  NumericVector out_v(n_pairs);
  NumericVector out_pmv_diag(n_pairs);
  NumericVector out_pmv_full(n_pairs);

  std::vector<double> d(static_cast<std::size_t>(m));
  std::vector<double> a(static_cast<std::size_t>(m));
  std::vector<double> Ra(static_cast<std::size_t>(m));
  std::vector<double> Sd(static_cast<std::size_t>(m));

  const double* R_ptr = R_decay.begin();
  const double* S_ptr = R_had_Sigma.begin();

  for (int r = 0; r < n_pairs; ++r) {
    const int p1 = pair_p1_zero[r];
    const int p2 = pair_p2_zero[r];

    // Build d, a and the diagonal PMV correction in one pass.
    double pmv_diag_extra = 0.0;
    for (int k = 0; k < m; ++k) {
      double dk = 0.5 * (geno(p1, k) - geno(p2, k));
      if (!R_finite(dk)) dk = 0.0;
      d[static_cast<std::size_t>(k)] = dk;
      double ak = dk * beta[k];
      if (!R_finite(ak)) ak = 0.0;
      a[static_cast<std::size_t>(k)] = ak;
      double bv = beta_var_diag[k];
      if (!R_finite(bv) || bv < 0.0) bv = 0.0;
      pmv_diag_extra += dk * dk * bv;
    }

    // Ra[i] = sum_k R_decay(i, k) * a[k]; Sd[i] = sum_k R_had_Sigma(i, k) * d[k].
    // Outer loop over column k (stride-1 in column-major), inner over row i.
    std::fill(Ra.begin(), Ra.end(), 0.0);
    std::fill(Sd.begin(), Sd.end(), 0.0);
    for (int k = 0; k < m; ++k) {
      const double ak = a[static_cast<std::size_t>(k)];
      const double dk = d[static_cast<std::size_t>(k)];
      const std::size_t col_off = static_cast<std::size_t>(k) * static_cast<std::size_t>(m);
      const double* R_col = R_ptr + col_off;
      const double* S_col = S_ptr + col_off;
      // Skip zero columns when ak == 0 AND dk == 0 (banded R has many).
      if (ak == 0.0 && dk == 0.0) continue;
      for (int i = 0; i < m; ++i) {
        Ra[static_cast<std::size_t>(i)] += R_col[i] * ak;
        Sd[static_cast<std::size_t>(i)] += S_col[i] * dk;
      }
    }

    // VPM and full-posterior PMV via final dot products.
    double v = 0.0;
    double pmv_off_extra = 0.0;
    for (int i = 0; i < m; ++i) {
      v += a[static_cast<std::size_t>(i)] * Ra[static_cast<std::size_t>(i)];
      pmv_off_extra += d[static_cast<std::size_t>(i)] * Sd[static_cast<std::size_t>(i)];
    }

    out_v[r]        = std::max(0.0, v);
    out_pmv_diag[r] = std::max(0.0, v + pmv_diag_extra);
    out_pmv_full[r] = std::max(0.0, v + pmv_off_extra);
  }

  return DataFrame::create(
    Named("dh_recomb_var") = out_v,
    Named("dh_pmv_var") = out_pmv_diag,
    Named("dh_pmv_var_full_posterior") = out_pmv_full
  );
}

// ---- Banded sparse-COO recombination kernel (Kosambi / RIL / windowed DH) -
//
// Consumes the upper-triangle COO band built by R's ng_recomb_decay_banded():
// triplets (i_idx, j_idx, x_val) with i_idx <= j_idx, both 0-based. Diagonal
// entries (i == j) contribute once; off-diagonals contribute twice
// (representing R symmetry). The R reference uses ifelse(is_diag, 1, 2);
// we pre-compute the weight via the i == j comparison.
//
// VPM = sum over triplets of a[i] * a[j] * x * weight
// PMV = VPM + sum_k d_k^2 * Var(beta_k)
//
// [[Rcpp::export]]
DataFrame ng_dh_recomb_pairs_banded_cpp(NumericMatrix geno,
                                        NumericVector beta,
                                        NumericVector beta_var,
                                        IntegerVector band_i_zero,
                                        IntegerVector band_j_zero,
                                        NumericVector band_x,
                                        IntegerVector pair_p1_zero,
                                        IntegerVector pair_p2_zero) {
  const int m = beta.size();
  const int n_pairs = pair_p1_zero.size();
  const int n_band = band_i_zero.size();
  if (band_j_zero.size() != n_band || band_x.size() != n_band) {
    Rcpp::stop("band_i, band_j, band_x must have the same length");
  }
  if (beta_var.size() != m) {
    Rcpp::stop("beta_var must have length m");
  }

  NumericVector out_v(n_pairs);
  NumericVector out_pmv(n_pairs);

  std::vector<double> d(static_cast<std::size_t>(m));
  std::vector<double> a(static_cast<std::size_t>(m));

  for (int r = 0; r < n_pairs; ++r) {
    const int p1 = pair_p1_zero[r];
    const int p2 = pair_p2_zero[r];

    double pmv_extra = 0.0;
    for (int k = 0; k < m; ++k) {
      double dk = 0.5 * (geno(p1, k) - geno(p2, k));
      if (!R_finite(dk)) dk = 0.0;
      d[static_cast<std::size_t>(k)] = dk;
      double ak = dk * beta[k];
      if (!R_finite(ak)) ak = 0.0;
      a[static_cast<std::size_t>(k)] = ak;
      double bv = beta_var[k];
      if (!R_finite(bv) || bv < 0.0) bv = 0.0;
      pmv_extra += dk * dk * bv;
    }

    double v = 0.0;
    for (int t = 0; t < n_band; ++t) {
      const int i = band_i_zero[t];
      const int j = band_j_zero[t];
      const double x = band_x[t];
      const double contribution = a[static_cast<std::size_t>(i)] *
                                  a[static_cast<std::size_t>(j)] * x;
      // Off-diagonal entries contribute twice for symmetry; diagonal once.
      v += (i == j) ? contribution : 2.0 * contribution;
    }

    out_v[r]   = std::max(0.0, v);
    out_pmv[r] = std::max(0.0, v + pmv_extra);
  }

  return DataFrame::create(
    Named("dh_recomb_var") = out_v,
    Named("dh_pmv_var") = out_pmv
  );
}

// ---- Bhattacharya-Chakraborty-Mallick (2016) posterior beta sampler ------
//
// Draw S samples from N(beta_hat, sigma_e2 * (X'X + lambda I)^{-1}) given:
//   1. theta_s ~ N(0, (sigma_e2 / lambda) I_m)
//   2. eta_s   ~ N(0, sigma_e2 I_n)
//   3. nu_s    = X theta_s + eta_s
//   4. w_s     = A^{-1} (yc - nu_s),  A = XX' + lambda I (Cholesky-factored once)
//   5. beta_s  = theta_s + X' w_s
//
// X is centered n x m. A_chol is the lower-triangular Cholesky factor of A
// (precomputed on the R side via base::chol(); we take its transpose because
// R chol returns upper-triangular, but we accept the upper factor here and
// do the standard back-then-forward solve). Inputs are stride-1 and indexed
// directly via NumericMatrix's column-major layout.
//
// Returns an m x n_draws matrix of beta draws.
//
// [[Rcpp::export]]
NumericMatrix ng_bcm_posterior_sampler_cpp(NumericMatrix X_centered,
                                           NumericVector yc,
                                           NumericMatrix A_chol_upper,
                                           double sigma_e2,
                                           double lambda,
                                           int n_draws,
                                           int seed) {
  const int n = X_centered.nrow();
  const int m = X_centered.ncol();
  if (yc.size() != n) Rcpp::stop("yc must have length nrow(X)");
  if (A_chol_upper.nrow() != n || A_chol_upper.ncol() != n) {
    Rcpp::stop("A_chol_upper must be n x n");
  }
  if (sigma_e2 <= 0.0 || !R_finite(sigma_e2)) Rcpp::stop("sigma_e2 must be > 0");
  if (lambda <= 0.0 || !R_finite(lambda)) Rcpp::stop("lambda must be > 0");
  if (n_draws < 1) Rcpp::stop("n_draws must be >= 1");

  NumericMatrix beta_draws(m, n_draws);

  // Seed R's RNG so the draw sequence matches stats::rnorm under
  // set.seed(seed). R::rnorm pulls from the same Mersenne Twister state as
  // R's own rnorm, so this gives exact reproducibility with the R reference.
  Rcpp::Function set_seed("set.seed");
  set_seed(seed);

  const double sd_prior = std::sqrt(sigma_e2 / lambda);
  const double sd_noise = std::sqrt(sigma_e2);

  std::vector<double> theta(static_cast<std::size_t>(m));
  std::vector<double> eta(static_cast<std::size_t>(n));
  std::vector<double> rhs(static_cast<std::size_t>(n));   // yc - nu
  std::vector<double> w(static_cast<std::size_t>(n));     // A^{-1} rhs

  const double* X_ptr = X_centered.begin();
  const double* U_ptr = A_chol_upper.begin();

  for (int s = 0; s < n_draws; ++s) {
    // Step 1-2: draw theta (m) and eta (n) from N(0, sd^2) using R's RNG
    // directly via R::rnorm (no per-call function dispatch overhead).
    for (int k = 0; k < m; ++k) theta[static_cast<std::size_t>(k)] = R::rnorm(0.0, sd_prior);
    for (int i = 0; i < n; ++i) eta[static_cast<std::size_t>(i)]   = R::rnorm(0.0, sd_noise);

    // Step 3: nu = X theta + eta; rhs = yc - nu
    //         X is n x m, column-major. For each column k, accumulate
    //         rhs[i] -= X(i,k) * theta[k] starting from rhs[i] = yc[i] - eta[i].
    for (int i = 0; i < n; ++i) {
      rhs[static_cast<std::size_t>(i)] = yc[i] - eta[static_cast<std::size_t>(i)];
    }
    for (int k = 0; k < m; ++k) {
      const double tk = theta[static_cast<std::size_t>(k)];
      if (tk == 0.0) continue;
      const std::size_t col_off = static_cast<std::size_t>(k) * static_cast<std::size_t>(n);
      const double* X_col = X_ptr + col_off;
      for (int i = 0; i < n; ++i) {
        rhs[static_cast<std::size_t>(i)] -= X_col[i] * tk;
      }
    }

    // Step 4: w = A^{-1} rhs via two triangular solves on the Cholesky factor.
    //         A = U' U where U is upper-triangular. Solve U' z = rhs, then U w = z.
    //         Both solves are O(n^2).
    //
    // Forward solve: U' z = rhs  (U' is lower-triangular: U(j, i) = U'(i, j))
    for (int i = 0; i < n; ++i) {
      double sum = rhs[static_cast<std::size_t>(i)];
      for (int j = 0; j < i; ++j) {
        // U' element (i, j) is U(j, i), at U_ptr[i*n + j].
        sum -= U_ptr[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                     static_cast<std::size_t>(j)] *
               w[static_cast<std::size_t>(j)];
      }
      const double pivot = U_ptr[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                                 static_cast<std::size_t>(i)];
      w[static_cast<std::size_t>(i)] = sum / pivot;
    }
    // Back solve: U w_final = z (z is currently in w; overwrite top-down)
    for (int i = n - 1; i >= 0; --i) {
      double sum = w[static_cast<std::size_t>(i)];
      for (int j = i + 1; j < n; ++j) {
        sum -= U_ptr[static_cast<std::size_t>(j) * static_cast<std::size_t>(n) +
                     static_cast<std::size_t>(i)] *
               w[static_cast<std::size_t>(j)];
      }
      const double pivot = U_ptr[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                                 static_cast<std::size_t>(i)];
      w[static_cast<std::size_t>(i)] = sum / pivot;
    }

    // Step 5: beta = theta + X' w. X' is m x n; (X' w)[k] = sum_i X(i, k) * w[i].
    // With column-major X, X(i, k) for fixed k is stride-1 in i, perfect.
    for (int k = 0; k < m; ++k) {
      double xtw = 0.0;
      const std::size_t col_off = static_cast<std::size_t>(k) * static_cast<std::size_t>(n);
      const double* X_col = X_ptr + col_off;
      for (int i = 0; i < n; ++i) {
        xtw += X_col[i] * w[static_cast<std::size_t>(i)];
      }
      beta_draws(k, s) = theta[static_cast<std::size_t>(k)] + xtw;
    }
  }

  // Preserve marker names on rows of beta_draws if X_centered carries colnames.
  Rcpp::List dimnames(2);
  if (X_centered.hasAttribute("dimnames")) {
    Rcpp::List xdn = X_centered.attr("dimnames");
    if (xdn.size() == 2L) dimnames[0] = xdn[1];  // markers
  }
  beta_draws.attr("dimnames") = dimnames;

  return beta_draws;
}

// ---- Greedy local-swap OCS optimizer (with fused objective evaluation) ---
//
// Replaces ng_local_swap()'s outer-and-inner R loops. The inner cost is
// dominated by ng_plan_objective() recomputing parent counts and group
// coancestry per swap candidate (R-level loop + O(p^2) K-matrix multiply).
// In C++ we fuse the two operations: incrementally maintain parent counts
// across swap candidates, and compute the K-quadratic form directly from
// the current count vector without R-level conversion.
//
// Inputs:
//   linear_gain      length-n_pairs vector of scores$.linear_gain
//   pair_p1_zero     length-n_pairs 0-based row index of parent1 into parents[]
//   pair_p2_zero     length-n_pairs 0-based row index of parent2 into parents[]
//   parent_K         p x p kinship matrix (rows/cols indexed 0..p-1, same
//                    order as the parent_p1/p2 indices)
//   selected_zero    initial selection (0-based pair indices)
//   max_per_parent   parent-use cap
//   lambda_group     group-coancestry penalty weight
//   local_iter       outer iteration cap
//   off_pool_size    inner-loop cap on "off" candidates considered (R uses 1000)
//
// Returns the 0-based pair indices of the swap-optimised selection.
//
// [[Rcpp::export]]
IntegerVector ng_local_swap_cpp(NumericVector linear_gain,
                                IntegerVector pair_p1_zero,
                                IntegerVector pair_p2_zero,
                                NumericMatrix parent_K,
                                IntegerVector selected_zero,
                                int max_per_parent,
                                double lambda_group,
                                double lambda_parent_use,
                                int local_iter,
                                int off_pool_size) {
  const int n_pairs = linear_gain.size();
  const int p = parent_K.nrow();
  if (parent_K.ncol() != p) Rcpp::stop("parent_K must be square");
  if (pair_p1_zero.size() != n_pairs || pair_p2_zero.size() != n_pairs) {
    Rcpp::stop("pair index vectors must match length(linear_gain)");
  }

  // Working state.
  std::vector<char> selected_flag(static_cast<std::size_t>(n_pairs), 0);
  std::vector<int>  selected(selected_zero.begin(), selected_zero.end());
  for (int idx : selected) {
    if (idx < 0 || idx >= n_pairs) Rcpp::stop("selected index out of range");
    selected_flag[static_cast<std::size_t>(idx)] = 1;
  }
  std::vector<int> counts(static_cast<std::size_t>(p), 0);
  for (int idx : selected) {
    counts[static_cast<std::size_t>(pair_p1_zero[idx])] += 1;
    counts[static_cast<std::size_t>(pair_p2_zero[idx])] += 1;
  }

  // Helper: objective = sum(linear_gain[selected]) - lambda * cvec' K cvec
  // where cvec = counts / sum(counts). The base sum is incrementally
  // maintained as gain_sum below; only the K-quadratic part recomputes.
  auto group_coancestry = [&](const std::vector<int>& cnt) -> double {
    int total = 0;
    for (int v : cnt) total += v;
    if (total <= 0) return 0.0;
    const double inv_total = 1.0 / static_cast<double>(total);
    // cvec[i] = cnt[i] * inv_total. q = sum_{i,j} cvec[i] * K(i,j) * cvec[j].
    // Compute K %*% cvec then dot with cvec; outer over column j is stride-1
    // in column-major.
    const double* K_ptr = parent_K.begin();
    std::vector<double> kc(static_cast<std::size_t>(p), 0.0);
    for (int j = 0; j < p; ++j) {
      if (cnt[static_cast<std::size_t>(j)] == 0) continue;
      const double cj = cnt[static_cast<std::size_t>(j)] * inv_total;
      const double* K_col = K_ptr + static_cast<std::size_t>(j) * static_cast<std::size_t>(p);
      for (int i = 0; i < p; ++i) kc[static_cast<std::size_t>(i)] += K_col[i] * cj;
    }
    double q = 0.0;
    for (int i = 0; i < p; ++i) {
      if (cnt[static_cast<std::size_t>(i)] == 0) continue;
      q += cnt[static_cast<std::size_t>(i)] * inv_total * kc[static_cast<std::size_t>(i)];
    }
    return q;
  };

  // Parent-use penalty = sum_p (cnt_p / total)^2, matching
  // ng_plan_objective_contribution() so the C++ kernel optimizes the SAME
  // objective as the MIP path. This removes the slow pure-R fallback that was
  // previously required whenever lambda_parent_use > 0.
  auto parent_use_sq = [&](const std::vector<int>& cnt) -> double {
    long long total = 0;
    for (int v : cnt) total += v;
    if (total <= 0) return 0.0;
    double ss = 0.0;
    for (int v : cnt) ss += static_cast<double>(v) * static_cast<double>(v);
    return ss / (static_cast<double>(total) * static_cast<double>(total));
  };

  double gain_sum = 0.0;
  for (int idx : selected) gain_sum += linear_gain[idx];
  double best_obj = gain_sum;
  if (lambda_group > 0.0) best_obj -= lambda_group * group_coancestry(counts);
  if (lambda_parent_use > 0.0) best_obj -= lambda_parent_use * parent_use_sq(counts);

  // Stable working buffers reused across iterations.
  std::vector<int> off; off.reserve(static_cast<std::size_t>(n_pairs));
  std::vector<int> on;  on.reserve(static_cast<std::size_t>(selected.size()));
  std::vector<int> trial_counts(static_cast<std::size_t>(p), 0);

  for (int it = 0; it < local_iter; ++it) {
    bool improved = false;

    // Build off/on index lists.
    off.clear(); on.clear();
    for (int idx = 0; idx < n_pairs; ++idx) {
      if (selected_flag[static_cast<std::size_t>(idx)]) on.push_back(idx);
      else off.push_back(idx);
    }
    // Sort: off descending by linear_gain (highest-gain swap-ins first),
    //       on ascending  by linear_gain (lowest-gain swap-outs first).
    // stable_sort keeps equal-gain entries in ascending index order, so ties are
    // resolved deterministically (reproducible plans on degenerate/tied landscapes).
    std::stable_sort(off.begin(), off.end(), [&](int a, int b) {
      return linear_gain[a] > linear_gain[b];
    });
    std::stable_sort(on.begin(), on.end(), [&](int a, int b) {
      return linear_gain[a] < linear_gain[b];
    });
    const int off_limit = std::min(static_cast<int>(off.size()), off_pool_size);

    for (std::size_t oi = 0; oi < on.size() && !improved; ++oi) {
      const int drop_idx = on[oi];
      const int drop_p1 = pair_p1_zero[drop_idx];
      const int drop_p2 = pair_p2_zero[drop_idx];

      // counts after the drop.
      std::copy(counts.begin(), counts.end(), trial_counts.begin());
      trial_counts[static_cast<std::size_t>(drop_p1)] -= 1;
      trial_counts[static_cast<std::size_t>(drop_p2)] -= 1;
      const double gain_drop = linear_gain[drop_idx];

      for (int ai = 0; ai < off_limit; ++ai) {
        const int add_idx = off[static_cast<std::size_t>(ai)];
        const int add_p1 = pair_p1_zero[add_idx];
        const int add_p2 = pair_p2_zero[add_idx];

        // Parent-use feasibility check.
        // Drop+add of the same parent edge keeps that parent's count flat;
        // need to handle the case where p1 == p2 (selfing pair).
        int post_p1_use = trial_counts[static_cast<std::size_t>(add_p1)] + 1;
        int post_p2_use = trial_counts[static_cast<std::size_t>(add_p2)] +
          ((add_p1 == add_p2) ? 1 : 1);
        if (add_p1 == add_p2) {
          // Selfing: add_p1 used twice by this single cross.
          if (trial_counts[static_cast<std::size_t>(add_p1)] + 2 > max_per_parent) continue;
        } else {
          if (post_p1_use > max_per_parent) continue;
          if (post_p2_use > max_per_parent) continue;
        }

        // candidate counts = trial_counts + delta from add.
        trial_counts[static_cast<std::size_t>(add_p1)] += 1;
        trial_counts[static_cast<std::size_t>(add_p2)] += 1;

        // Candidate objective.
        const double candidate_gain = gain_sum - gain_drop + linear_gain[add_idx];
        double obj = candidate_gain;
        if (lambda_group > 0.0) {
          obj -= lambda_group * group_coancestry(trial_counts);
        }
        if (lambda_parent_use > 0.0) {
          obj -= lambda_parent_use * parent_use_sq(trial_counts);
        }

        if (obj > best_obj + 1e-10) {
          // Commit swap.
          selected_flag[static_cast<std::size_t>(drop_idx)] = 0;
          selected_flag[static_cast<std::size_t>(add_idx)]  = 1;
          // Update counts to the trial state.
          std::copy(trial_counts.begin(), trial_counts.end(), counts.begin());
          // Update gain_sum and best_obj.
          gain_sum = candidate_gain;
          best_obj = obj;
          // Replace drop_idx with add_idx in selected vector.
          for (std::size_t k = 0; k < selected.size(); ++k) {
            if (selected[k] == drop_idx) { selected[k] = add_idx; break; }
          }
          improved = true;
          break;
        }
        // Undo trial_counts add.
        trial_counts[static_cast<std::size_t>(add_p1)] -= 1;
        trial_counts[static_cast<std::size_t>(add_p2)] -= 1;
      }
    }
    if (!improved) break;
  }

  // Return as 0-based IntegerVector matching the R selected[] convention
  // used elsewhere (R caller adds 1L for R 1-based indexing).
  IntegerVector out(selected.size());
  for (std::size_t k = 0; k < selected.size(); ++k) out[k] = selected[k];
  return out;
}

// Dominance-aware cross scoring: accumulate the O(n_crosses x markers) loop over the precomputed
// progeny-moment tables. For each cross c and marker k, look up the moments at the parental dosage
// pair (di, dj) and accumulate mid-parent breeding value, heterosis (expected progeny dominance),
// additive and dominance within-family variance, and the additive x dominance covariance.
// Returns an (n_crosses x 5) matrix: [mid_bv, heterosis, add_var, dom_var, cov_ad].
// [[Rcpp::export]]
NumericMatrix ng_poly_dominance_scores_cpp(IntegerMatrix M,
                                           IntegerVector i1, IntegerVector i2,
                                           NumericMatrix mu, NumericMatrix varX,
                                           NumericMatrix EH, NumericMatrix varH, NumericMatrix covXH,
                                           NumericVector ba, NumericVector bd,
                                           NumericVector cen_a, NumericVector hbar,
                                           double intercept, bool has_dom) {
  int n = i1.size();
  int m = M.ncol();
  NumericMatrix out(n, 5);
  for (int c = 0; c < n; ++c) {
    int r1 = i1[c], r2 = i2[c];
    double bv = intercept, het = 0.0, av = 0.0, dv = 0.0, cv = 0.0;
    for (int k = 0; k < m; ++k) {
      int di = M(r1, k), dj = M(r2, k);
      double a = ba[k];
      bv += a * (mu(di, dj) - cen_a[k]);
      av += a * a * varX(di, dj);
      if (has_dom) {
        double d = bd[k];
        het += d * (EH(di, dj) - hbar[k]);
        dv += d * d * varH(di, dj);
        cv += 2.0 * a * d * covXH(di, dj);
      }
    }
    out(c, 0) = bv; out(c, 1) = het; out(c, 2) = av; out(c, 3) = dv; out(c, 4) = cv;
  }
  return out;
}
