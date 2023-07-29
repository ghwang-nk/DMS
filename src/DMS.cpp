#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double func_weight(double t, double gamma) {
  return(std::pow(t * (1 - t), -gamma));
}

//' CUSUM statistics
//'
//' @param X \eqn{n\times p}, component-wise unit variances
//' @param gamma 0 = CUSUM, 0.5 = two-sample
//'
//' @return CUSUM statistics, \eqn{(n - 1)\times p\times \code{length(gamma)}}
//'
//' @export
// [[Rcpp::export]]
arma::cube cusum(arma::mat X, arma::vec gamma) {
  int n = X.n_rows;
  arma::cube out = arma::zeros(n - 1, X.n_cols, gamma.n_rows);
  arma::rowvec Sn = arma::sum(X, 0);
  arma::rowvec Si = arma::zeros(1, X.n_cols);
  for (int i = 1; i < n; i++) {
    Si += X.row(i - 1);
    arma::rowvec temp = (Si - i * Sn / n) / std::sqrt(n);
    for (int it = 0; it < gamma.n_rows; it++) {
      if (gamma(it) == 0) {
        out(arma::span(i - 1, i - 1), arma::span::all, arma::span(it, it)) = temp;
      } else {
        double weight = func_weight(double(i) / n, gamma(it));
        out(arma::span(i - 1, i - 1), arma::span::all, arma::span(it, it)) = weight * temp;
      }
    }
  }
  return(out);
}

//' Bootstrapped CUSUM statistics (two-sample version)
//'
//' @param X \eqn{n}-dimensional, unit variance
//' @param e multipliers, \eqn{n\times B}
//'
//' @return Bootstrapped CUSUM statistics, \eqn{(n - 1)\times B}
//'
// [[Rcpp::export]]
arma::mat cusum_2sample_boots(arma::vec X, arma::mat e) {

  int n = e.n_rows;
  int B = e.n_cols;
  arma::mat out = arma::zeros(n - 1, B);

  arma::rowvec sum_e = arma::sum(e, 0);
  double sum_X = arma::sum(X);
  arma::mat eX = e.each_col() % X;
  arma::rowvec sum_eX = arma::sum(eX, 0);

  int n_tl = 0;
  int n_tr = n;
  arma::rowvec sum_e_tl(B, arma::fill::zeros);
  arma::rowvec sum_e_tr = sum_e;
  double sum_X_tl = 0;
  double sum_X_tr = sum_X;
  arma::rowvec sum_eX_tl(B, arma::fill::zeros);
  arma::rowvec sum_eX_tr = sum_eX;

  for (int t = 0; t < n - 1; t++) {  // t = 1:(n-1)
    n_tl += 1;
    n_tr -= 1;
    sum_e_tl += e.row(t);
    sum_e_tr -= e.row(t);
    sum_X_tl += X(t);
    sum_X_tr -= X(t);
    sum_eX_tl += eX.row(t);
    sum_eX_tr -= eX.row(t);
    arma::rowvec dl = sum_eX_tl / n_tl - sum_e_tl / (std::pow(n_tl, 2) / sum_X_tl);
    arma::rowvec dr = sum_eX_tr / n_tr - sum_e_tr / (std::pow(n_tr, 2) / sum_X_tr);
    out.row(t) = std::sqrt(n_tl * n_tr / double(n)) * (dl - dr);
  }
  return(out);
}

//' Bootstrapped Max(0.5) statistics
//'
//' @param X \eqn{n\times p}
//' @param e multipliers, \eqn{n\times B}
//' @param n_boundary \code{[n_boundary:(n - n_boundary)]}
//'
//' @return Bootstrapped Max(0.5) statistics
//'
//' @references
//' Yu, M. and Chen, X. (2021) Finite Sample Change Point Inference and Identification for High-Dimensional Mean Vectors. Journal of the Royal Statistical Society Series B: Statistical Methodology, 83, 247–270
//'
//' @seealso \code{\link{YC21}}
//'
// [[Rcpp::export]]
arma::rowvec cal_YC21(arma::mat X, arma::mat e, int n_boundary) {
  int t_start = n_boundary;
  int t_end = X.n_rows - n_boundary;
  arma::mat res = arma::zeros(X.n_cols, e.n_cols);
  for (int j = 0; j < X.n_cols; j++) {
    arma::mat cusum_2sample_boots_j = cusum_2sample_boots(X.col(j), e);
    res.row(j) = arma::max(arma::abs(cusum_2sample_boots_j.rows(arma::span(t_start - 1, t_end - 1))), 0);
  }
  arma::rowvec out = arma::max(res, 0);
  return(out);
}

//' Difference-based variance estimates
//'
//' @param X \eqn{n\times p}
//'
//' @return Component-wise variance estimates (\eqn{1\times p})
//'
//' @references
//' Rice, J. (1984). Bandwidth choice for nonparametric regression. The Annals of Statistics, 12, 1215-1230.
//'
//' @export
// [[Rcpp::export]]
List est_var_diff(arma::mat X) {
  List out;
  arma::rowvec s2(X.n_cols, arma::fill::zeros);
  arma::mat d2 = arma::zeros(X.n_rows - 1, X.n_cols);
  for (int i = 0; i < X.n_rows - 1; i++) {
    d2.row(i) = arma::pow(X.row(i + 1) - X.row(i), 2);
    s2 += d2.row(i);
  }
  s2 *= 1.0 / (2 * (X.n_rows - 1));
  out["s2"] = s2;
  out["d2"] = d2;
  return(out);
}

//' Normalization of Sum(0.5) statistic
//'
//' @param X \eqn{n\times p}
//'
//' @return Estimated expectation and variance of the Sum(0.5) statistic
//'
//' @references
//' Wang, Y., Zou, C., Wang, Z. and Yin, G. (2019) Multiple change-points detection in high dimension. Random Matrices: Theory and Applications, 08, 1950014.
//'
//' @seealso \code{\link{WZWY19}}
//'
// [[Rcpp::export]]
List cal_WZWY19(arma::mat X) {

  List out;

  List temp = est_var_diff(X);
  arma::rowvec s2 = temp["s2"];
  arma::rowvec sum_d2 = 2 * (X.n_rows - 1) * s2;
  arma::mat d2 = temp["d2"];

  arma::rowvec Zero(X.n_cols, arma::fill::zeros);
  arma::mat d2_ext = arma::join_cols(Zero, Zero, d2, Zero);

  double est_tr_R2 = 0;
  int K = 4;  // leave-K-out
  arma::rowvec up = arma::sum(d2_ext.rows(arma::span(0, K)), 0);
  for (int i = 1; i <= X.n_rows - K + 1; i++) {
    up += d2_ext.row(i + K);
    up -= d2_ext.row(i - 1);
    arma::rowvec s2_i = sum_d2 - up;
    if (i > 1 & i < X.n_rows - K + 1) {
      int i_start = i + K;
      int i_end = i - 1;
      s2_i = sum_d2 - up + arma::pow(X.row(i_start - 1) - X.row(i_end - 1), 2);
    }
    s2_i *= 1.0 / (2 * (X.n_rows - K - 1));
    double dT_DInv_d = arma::sum((X.row(i - 1) - X.row(i)) % (X.row(i + 1) - X.row(i + 2)) / s2_i);
    est_tr_R2 += dT_DInv_d * dT_DInv_d;
  }
  est_tr_R2 *= 1.0 / (4 * (X.n_rows - 3));

  double est_epsT_tR_eps2 = 0;
  K = 3;  // leave-K-out
  up = arma::sum(d2_ext.rows(arma::span(0, K)), 0);
  for (int i = 1; i <= X.n_rows - K + 1; i++) {
    up += d2_ext.row(i + K);
    up -= d2_ext.row(i - 1);
    arma::rowvec s2_i = sum_d2 - up;
    if (i > 1 & i < X.n_rows - K + 1) {
      int i_start = i + K;
      int i_end = i - 1;
      s2_i = sum_d2 - up + arma::pow(X.row(i_start - 1) - X.row(i_end - 1), 2);
    }
    s2_i *= 1.0 / (2 * (X.n_rows - K - 1));
    double dT_DInv_d = arma::sum((X.row(i - 1) - X.row(i)) % (X.row(i) - X.row(i + 1)) / s2_i);
    est_epsT_tR_eps2 += dT_DInv_d * dT_DInv_d;
  }
  est_epsT_tR_eps2 *= 1.0 / (X.n_rows - 2);
  est_epsT_tR_eps2 -= 3 * est_tr_R2;

  double pi = arma::datum::pi;
  out["mean"] = (X.n_rows + 2) * X.n_cols;
  out["var"] = (2 * pi * pi - 18) / 3 * X.n_rows * X.n_rows * est_tr_R2 + (15 - pi * pi) / 3 * X.n_rows * (est_epsT_tR_eps2 - X.n_cols * X.n_cols);
  out["s2"] = s2;
  out["tr_R2"] = est_tr_R2;
  out["epsT_tR_eps2"] = est_epsT_tR_eps2;
  return(out);
}

//' \eqn{(s_0, p)}-norm
//'
//' @param C \eqn{(n - 1)\times p}
//' @param n_boundary \code{[n_boundary:(n - n_boundary)]}
//' @param s0 \eqn{s_0}
//' @param p \eqn{p}
//'
//' @return Maximum and maximizer over all feasible changepoints for each \eqn{p}
//'
//' @references
//' Liu, B., Zhou, C., Zhang, X. and Liu, Y. (2020). A unified data-adaptive framework for high dimensional change point detection. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 82 933–63.
//'
//' @seealso \code{\link{YC21}} \code{\link{cal_LZZL20}}
//'
// [[Rcpp::export]]
List norm_s0_p(arma::mat C, int n_boundary, int s0, arma::vec p) {
  List out;
  arma::vec stat = arma::zeros(p.n_rows);
  arma::vec th = arma::zeros(p.n_rows);
  int n = C.n_rows + 1;
  int t_start = n_boundary;
  int t_end = n - n_boundary;
  arma::mat C_abs = arma::abs(C);
  arma::mat C_abs_sort = arma::sort(C_abs, "descend", 1);
  arma::mat C_abs_sort_s0 = C_abs_sort(arma::span(t_start - 1, t_end - 1), arma::span(0, s0 - 1));
  for (int ip = 0; ip < p.n_rows; ip++) {
    arma::vec temp = arma::zeros(C.n_rows);
    if (std::isfinite(p(ip))) {
      temp = arma::pow(arma::sum(arma::pow(C_abs_sort_s0, p(ip)), 1), 1 / p(ip));
    } else {
      temp = C_abs_sort_s0.col(0);
    }
    stat(ip) = arma::max(temp);
    th(ip) = arma::index_max(temp) + t_start;
  }
  out["stat"] = stat;
  out["th"] = th;
  return(out);
}

//' Bootstrapped LZZL statistics
//'
//' @param X \eqn{n\times p}
//' @param e multipliers, \eqn{n\times B}
//' @param n_boundary \code{[n_boundary:(n - n_boundary)]}
//' @param s0 \eqn{s_0 = p / 2}
//' @param p \eqn{p \in \{1,2,3,4,5,\infty\}}
//' @param sig estimates of component-wise standard deviations
//'
//' @return Bootstrapped LZZL statistics
//'
//' @references
//' Liu, B., Zhou, C., Zhang, X. and Liu, Y. (2020). A unified data-adaptive framework for high dimensional change point detection. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 82 933–63.
//'
//' @seealso \code{\link{YC21}} \code{\link{norm_s0_p}}
//'
// [[Rcpp::export]]
arma::mat cal_LZZL20(arma::mat X, arma::mat e, int n_boundary, int s0, arma::vec p, arma::mat sig) {
  int n = X.n_rows;
  arma::cube C_boots = arma::zeros(n - 1, X.n_cols, e.n_cols);
  arma::vec weight = arma::sqrt(arma::linspace(1, n - 1, n - 1) % arma::linspace(n - 1, 1, n - 1)) / double(n);
  for (int j = 0; j < X.n_cols; j++) {
    arma::mat cusum_2sample_boots_j = cusum_2sample_boots(X.col(j), e);
    C_boots(arma::span::all, arma::span(j, j), arma::span::all) = cusum_2sample_boots_j.each_col() % (weight / sig.col(j));
  }
  arma::mat out = arma::zeros(p.n_rows, C_boots.n_slices);
  for (int b = 0; b < C_boots.n_slices; b++) {
    List out_b = norm_s0_p(C_boots.slice(b), n_boundary, s0, p);
    arma::vec stat_b = out_b["stat"];
    out.col(b) = stat_b;
  }
  return(out);
}

//' ZWS statistics
//'
//' The following implementation is an adaptation of the Matlab code created by Zhang, Wang, and Shao (2020).
//' The original Matlab code can be accessed at the following link:
//' https://doi.org/10.1080/01621459.2021.1884562.
//'
//' @param X \eqn{n\times p}
//' @param q \eqn{q \in \{2, 6\}}, suggested
//' @param truc \eqn{\code{truc} = q}
//'
//' @return ZWS statistics
//'
//' @references
//' Zhang, Y., Wang, R. and Shao, X. (2022) Adaptive Inference for Change Points in High-Dimensional Data. Journal of the American Statistical Association, 117, 1751–1762.
//'
// [[Rcpp::export]]
List cal_ZWS22(arma::mat X, int q, int truc) {

  List out;

  int n = X.n_rows;
  int p = X.n_cols;

  arma::vec fact = arma::zeros(q);
  arma::mat choose = arma::zeros(q, n);
  for (int i = 1; i <= q; i++) {
    fact(i - 1) = R::gammafn(i + 1);
  }
  for (int j = 1; j <= q; j++) {
    for (int i = j; i <= n; i++) {
      choose(j - 1, i - 1) = R::choose(i, j);
    }
  }

  arma::uvec index = arma::conv_to<arma::uvec>::from(arma::regspace(2 * truc, n - 2 * truc));
  arma::mat num = arma::zeros(p, index.n_rows);
  arma::mat denom = arma::zeros(index.n_rows, n - 4 * truc + 2);
  arma::mat temp = arma::zeros(index.n_rows, n - 4 * truc + 2);
  arma::vec tempsum = arma::zeros(q);
  arma::rowvec curr = arma::zeros<arma::rowvec>(n - 4 * truc + 2);
  arma::mat presum = arma::zeros(q, n);
  arma::mat postsum = arma::zeros(q, n);

  for (int l = 1; l <= p; l++) {

    arma::rowvec x = X.col(l - 1).t();
    presum.row(0) = arma::cumsum(x);
    postsum.row(0) = arma::cumsum(arma::reverse(x));
    for (int c = 2; c <= q; c++) {
      for (int i = c; i <= n + c - q; i++) {
        presum(c - 1, i - 1) = presum(c - 1, i - 2) + presum(c - 2, i - 2) * x(i - 1);
        postsum(c - 1, i - 1) = postsum(c - 1, i - 2) + postsum(c - 2, i - 2) * x(n - i);
      }
    }

    arma::rowvec presum_q = presum.row(q - 1);
    arma::rowvec postsum_q = postsum.row(q - 1);
    arma::rowvec choose_q = choose.row(q - 1);
    arma::vec tmp = fact(q - 1) * presum_q(index - 1) % choose_q(n - index - 1) + std::pow(-1, q) * fact(q - 1) * postsum_q(n - index - 1) % choose_q(index - 1);
    num.row(l - 1) += tmp.t();
    for (int c = 1; c <= q - 1; c++) {
      arma::rowvec presum_c = presum.row(c - 1);
      arma::rowvec choose_c = choose.row(c - 1);
      arma::rowvec choose_qc = choose.row(q - c - 1);
      arma::rowvec postsum_qc = postsum.row(q - c - 1);
      arma::vec tmp = std::pow(-1, q - c) * fact(c - 1) * fact(q - c - 1) * presum_c(index - 1) % choose_c(n - index - q + c - 1) % choose_qc(index - c - 1) % postsum_qc(n - index - 1);
      num.row(l - 1) += tmp.t();
    }

    for (int k : index) {
      for (int t = truc; t <= k - truc; t++) {
        double dq = presum(q - 1, t - 1) * choose(q - 1, k - t - 1) * fact(q - 1);
        tempsum(0) = arma::sum(x(arma::span(t, k - 1)));
        for (int c = 2; c <= q; c++) {
          dq += tempsum(c - 2) * presum(q - c, t - 1) * std::pow(-1, c - 1) * choose(q - c, k - t - c) * choose(c - 2, t - q + c - 2) * fact(c - 2) * fact(q - c);
          tempsum(c - 1) = presum(c - 1, k - 1) - presum(c - 1, t - 1);
          for (int d = 1; d <= c - 1; d++) {
            tempsum(c - 1) = tempsum(c - 1) - presum(d - 1, t - 1) * tempsum(c - d - 1);
          }
        }
        curr(t - truc) = dq + tempsum(q - 1) * std::pow(-1, q) * choose(q - 1, t - 1) * fact(q - 1);
      }
      for (int t = truc; t <= n - truc - k; t++) {
        double dq = postsum(q - 1, t - 1) * std::pow(-1, q) * choose(q - 1, n - t - k - 1) * fact(q - 1);
        tempsum(0) = arma::sum(x(arma::span(k, n - t - 1)));
        for (int c = 2; c <= q; c++) {
          dq += tempsum(c - 2) * postsum(q - c, t - 1) * std::pow(-1, q - c + 1) * choose(q - c, n - t - k - c) * choose(c - 2, t - q + c - 2) * fact(c - 2) * fact(q - c);
          tempsum(c - 1) = postsum(c - 1, n - k - 1) - postsum(c - 1, t - 1);
          for (int d = 1; d <= c - 1; d++) {
            tempsum(c - 1) = tempsum(c - 1) - postsum(d - 1, t - 1) * tempsum(c - d - 1);
          }
        }
        curr(t + k - 3 * truc + 1) = dq + tempsum(q - 1) * choose(q - 1, t - 1) * fact(q - 1);
      }
      temp.row(k - 2 * truc) = curr;
    }
    denom += temp;
  }

  arma::rowvec ts_num = arma::pow(arma::sum(num, 0), 2);
  arma::rowvec ts_denom = (arma::sum(arma::pow(denom, 2), 1)).t();
  double ts = arma::max(ts_num / ts_denom) * n;

  out["ts"] = ts;
  out["num"] = ts_num;
  out["denom"] = ts_denom;
  return(out);
}
