#' Distribution function of the Gumbel distribution
#'
#' @param x vector of quantiles
#' @param mu location
#' @param beta scale
#'
pgumbel <- function(x, mu = 0, beta = 1) {
  exp(-exp(-(x - mu) / beta))
}

#' Normalizing factor A
#'
#' @param x x
#'
f_A <- function(x) {
  sqrt(2 * log(x))
}

#' Normalizing factor D
#'
#' @param x x
#' @param d dimension
#'
f_D <- function(x, d = 1) {
  2 * log(x) + d / 2 * log(log(x)) - log(gamma(d / 2))
}

#' Max(0) (Jirak, 2015)
#'
#' @param X \eqn{n\times p}, component-wise unit variances
#' @param CUSUM \eqn{(n - 1)\times p}, or \code{NULL}
#' @param mod \code{TRUE} = use modified normalizing factors in Wang and Feng (2023)
#'
#' @return The p-value and the test statistic
#'
#' @references
#' Jirak, M. (2015) Uniform change point tests in high dimension. The Annals of Statistics, 43, 2451–2483.
#'
#' Wang, G. and Feng, L. (2023) Computationally efficient and data-adaptive changepoint inference in high dimension. Journal of the Royal Statistical Society Series B: Statistical Methodology, 85, 936–958.
#'
#' @export
J15 <- function(X, CUSUM = NULL, mod = TRUE) {
  if (is.null(CUSUM)) {
    CUSUM <- cusum(X, 0)[, , 1]
  }
  M <- max(abs(CUSUM))
  if (mod) {
    pv <- 1 - pgumbel(2 * M ^ 2 - log(2 * ncol(X)))
  } else {
    ep <- sqrt(2 * log(2 * ncol(X)))
    fp <- ep / 2 - log(3 * log(2 * ncol(X))) / ep
    pv <- 1 - pgumbel(ep * (M - fp))
  }
  return(list(pv = pv, M = M))
}

#' Max(0.5) (Wang and Feng, 2023)
#'
#' @param X \eqn{n\times p}, component-wise unit variances
#' @param n.boundary \code{[n.boundary:(n - n.boundary)]}
#' @param CUSUM \eqn{(n - 1)\times p}, or \code{NULL}
#'
#' @return The p-value and the test statistic
#'
#' @references
#' Wang, G. and Feng, L. (2023) Computationally efficient and data-adaptive changepoint inference in high dimension. Journal of the Royal Statistical Society Series B: Statistical Methodology, 85, 936–958.
#'
#' @export
max_2sample <- function(X, n.boundary, CUSUM = NULL) {
  if (is.null(CUSUM)) {
    CUSUM <- cusum(X, 0.5)[, , 1]
  }
  M <- max(abs(CUSUM[n.boundary:(nrow(X) - n.boundary), ]))
  temp <- ncol(X) * log((nrow(X) / n.boundary - 1) ^ 2)
  pv <- 1 - pgumbel(f_A(temp) * M - f_D(temp))
  return(list(pv = pv, M = M))
}

#' Max(0.5), bootstrap calibrated (Yu and Chen, 2021)
#'
#' @param X \eqn{n\times p}
#' @param n.boundary \code{[n.boundary:(n - n.boundary)]}
#' @param B number of bootstrap replications
#' @param CUSUM \eqn{(n - 1)\times p}, or \code{NULL}
#' @param M Max(0.5) statistic, or \code{NULL}
#'
#' @return The p-value and the test statistic
#'
#' @references
#' Yu, M. and Chen, X. (2021) Finite Sample Change Point Inference and Identification for High-Dimensional Mean Vectors. Journal of the Royal Statistical Society Series B: Statistical Methodology, 83, 247–270.
#'
#' @seealso \code{\link{cal_YC21}}
#'
#' @importFrom stats rnorm
#' @export
YC21 <- function(X, n.boundary = 40, B = 200, CUSUM = NULL, M = NULL) {
  if (is.null(CUSUM) & is.null(M)) {
    CUSUM <- cusum(X, 0.5)[, , 1]
    M <- max(abs(CUSUM[n.boundary:(nrow(X) - n.boundary), ]))
  }
  e <- matrix(rnorm(nrow(X) * B), nrow = nrow(X))
  M_boots <- cal_YC21(X, e, n.boundary)[1, ]  # row-vector returned by cal_YC21
  pv <- mean(M_boots >= M)
  return(list(pv = pv, M = M))
}

#' Sum(0.5), with power-enhanced version (Wang, Zou, Wang and Yin, 2019)
#'
#' @param X \eqn{n\times p}
#' @param CUSUM \eqn{(n - 1)\times p}, or \code{NULL}
#' @param PE \code{TRUE} = include the power-enhanced version
#' @param PE.c \eqn{c}
#' @param PE.h \eqn{h}
#' @param n.boundary \code{[n.boundary:(n - n.boundary)]} (for Max(0.5))
#'
#' @return The p-value and the test statistic
#'
#' @references
#' Wang, Y., Zou, C., Wang, Z. and Yin, G. (2019) Multiple change-points detection in high dimension. Random Matrices: Theory and Applications, 08, 1950014.
#'
#' @seealso \code{\link{cal_WZWY19}}
#'
#' @importFrom stats pnorm
#' @export
WZWY19 <- function(X, CUSUM = NULL, PE = TRUE, PE.c = 100, PE.h = (2 * log(nrow(X) * ncol(X))) ^ 1.1, n.boundary = 0.1 * nrow(X)) {
  if (is.null(CUSUM)) {
    CUSUM <- cusum(X, 0.5)[, , 1]
  }
  res <- cal_WZWY19(X)
  s2 <- res$s2[1, ]
  CUSUM2_s <- sweep(CUSUM ^ 2, 2, s2, "/")
  S <- sum(CUSUM2_s)
  S <- (S - res$mean) / sqrt(res$var)
  pv_sum <- 1 - pnorm(S)
  ts <- NULL
  pv <- NULL
  if (PE) {
    M2 <- max(abs(CUSUM2_s[n.boundary:(nrow(X) - n.boundary), ]))
    ts <- S + PE.c * sqrt(res$var) * (M2 > PE.h)
    pv <- 1 - pnorm(ts)
  }
  return(list(pv = pv, ts = ts, pv_sum = pv_sum, S = S, CUSUM2_s = CUSUM2_s))
}

#' LZZL20
#'
#' @param X \eqn{n\times p}
#' @param CUSUM CUSUM(0)
#' @param CUSUM_2sample CUSUM(0.5), for \eqn{\hat{\sigma}_{t,j}}
#' @param LZZL.s0 \eqn{s_0 = p / 2}
#' @param LZZL.p \eqn{p \in \{1,2,3,4,5,\infty\}}
#' @param LZZL.B number of bootstrap replications
#' @param LZZL.boundary \code{[n.boundary:(n - n.boundary)]}
#' @param est \code{TRUE} = enable changepoint estimation
#'
#' @return The p-value
#'
#' @references
#' Liu, B., Zhou, C., Zhang, X. and Liu, Y. (2020) A Unified Data-Adaptive Framework for High Dimensional Change Point Detection. Journal of the Royal Statistical Society Series B: Statistical Methodology, 82, 933–963.
#'
#' @seealso \code{\link{norm_s0_p}} \code{\link{cal_LZZL20}}
#'
#' @importFrom stats var rnorm
#' @export
LZZL20 <- function(X, CUSUM = NULL, CUSUM_2sample = NULL, LZZL.s0 = ncol(X) / 2, LZZL.p = c(1:5, Inf), LZZL.B = 500, LZZL.boundary = 0.2 * nrow(X), est = FALSE) {
  if (is.null(CUSUM) & is.null(CUSUM_2sample)) {
    temp <- cusum(X, c(0, 0.5))
    CUSUM <- temp[, , 1]
    CUSUM_2sample <- temp[, , 2]
  }
  s2_0 <- (nrow(X) - 1) * apply(X, 2, var)
  s_tj <- sqrt(-sweep(CUSUM_2sample ^ 2, 2, s2_0, "-") / nrow(X))
  C <- CUSUM / s_tj
  temp <- norm_s0_p(C, LZZL.boundary, LZZL.s0, LZZL.p)
  C_norm <- temp$stat
  e <- matrix(rnorm(nrow(X) * LZZL.B), nrow = nrow(X))
  C_norm_boots <- cal_LZZL20(X, e, LZZL.boundary, LZZL.s0, LZZL.p, s_tj)
  pv_C_norm <- sapply(1:length(LZZL.p), function (ip) {
    sum(C_norm_boots[ip, ] > C_norm[ip]) / (LZZL.B + 1)
  })
  T_ad <- min(pv_C_norm)
  pv_LCB <- t(sapply(1:length(LZZL.p), function (ip) {
    if (length(unique(C_norm_boots[ip, ])) < ncol(C_norm_boots)) {
      stop("Ties when calculating pv_LCB!")
    }
    (ncol(C_norm_boots) - rank(C_norm_boots[ip, ])) / ncol(C_norm_boots)
  }))
  T_ad_boots <- apply(pv_LCB, 2, min)
  pv_LZZL <- sum(T_ad_boots <= T_ad) / (LZZL.B + 1)
  th <- NULL
  if (est) {
    th <- temp$th[which.min(pv_C_norm)]
  }
  return(list(pv = pv_LZZL, pvs = pv_C_norm, th = th, ths = temp$th))
}

#' P-value combination
#'
#' @param pvs p-values
#' @param type combination rule
#'
#' @return The combined p-value
#'
#' @importFrom stats pchisq pcauchy
#' @export
pvc <- function(pvs, type = "min") {
  pv <- dplyr::case_when(
    type == "min" ~ 1 - (1 - min(pvs)) ^ length(pvs),
    type == "Fisher" ~ 1 - pchisq(-2 * sum(log(pvs)), 2 * length(pvs)),
    type == "Cauchy" ~ pcauchy(mean(tan((0.5 - pvs) * pi)), lower.tail = FALSE)
  )
  return(pv)
}

#' ZWS22
#'
#' @param X \eqn{n\times p}
#'
#' @return The p-value
#'
#' @references
#' Zhang, Y., Wang, R. and Shao, X. (2022) Adaptive Inference for Change Points in High-Dimensional Data. Journal of the American Statistical Association, 117, 1751–1762.
#'
#' @seealso \code{\link{cal_ZWS22}}
#'
#' @export
ZWS22 <- function(X) {
  q <- 2
  res <- cal_ZWS22(X, q, q)
  pv_2 <- mean(cal_q2 >= res$ts)
  q <- 6
  res <- cal_ZWS22(X, q, q)
  pv_6 <- mean(cal_q6 >= res$ts)
  pv_ZWS <- pvc(c(pv_2, pv_6), "min")
  return(pv = pv_ZWS)
}
