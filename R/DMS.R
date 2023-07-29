#' DMS method for adaptive changepoint testing (for independent observations)
#'
#' @param X A data matrix of dimension \eqn{n\times p}, where each row represents an observation, and each column stands for a variable.
#' @param gam A numeric value (\eqn{0} or \eqn{0.5}) indicating which version of the DMS statistic to use.
#' @param lam An integer between \eqn{1} and \eqn{n/2}, representing a pre-specified boundary removal parameter if \code{gam = 0.5}.
#'
#' @return \code{DMS} returns a list containing p-values of the DMS, max-, and sum-based testing procedures, respectively.
#'
#' @export
#'
#' @examples
#' library(DMS)
#' set.seed(6)
#' n <- 200
#' p <- 200
#' tau <- floor(n / 2)
#' mu <- matrix(0, n, p)
#' mu[(tau + 1):n, 1:10] <- 0.5
#' eps <- matrix(rnorm(n * p), n, p)
#' X <- mu + eps
#' DMS(X, gam = 0)
#' # $pv.DMS
#' # [1] 1.528718e-07
#' # $pv.Max
#' # [1] 0.0004963364
#' # $pv.Sum
#' # [1] 1.565607e-05
#' DMS(X, gam = 0.5, lam = 20)
#' # $pv.DMS
#' # [1] 1.404542e-06
#' # $pv.Max
#' # [1] 0.005177261
#' # $pv.Sum
#' # [1] 1.565607e-05
DMS <- function(X, gam = 0, lam = NULL) {

  n <- nrow(X)
  p <- ncol(X)

  # cusum ----
  temp <- cusum(X, c(0, 0.5))
  CUSUM_2sample <- temp[, , 2]

  # var estimates ----
  se <- sqrt(drop(est_var_diff(X)$s2))
  CUSUM_2sample_s <- sweep(CUSUM_2sample, 2, se, "/")

  # Max ----
  pv_max <- NULL
  if (gam == 0) {
    CUSUM <- temp[, , 1]
    CUSUM_s <- sweep(CUSUM, 2, se, "/")
    pv_max <- J15(X, CUSUM_s)$pv
  } else if (gam == 0.5) {
    if (is.null(lam)) {
      stop("Argument lam must be specified!")
    }
    pv_max <- max_2sample(X, lam, CUSUM_2sample_s)$pv
  } else {
    stop("Check Argument gam!")
  }

  # Sum ----
  pv_sum <- WZWY19(X, CUSUM_2sample, PE = FALSE)$pv_sum

  # DMS ----
  pv_DMS <- pvc(c(pv_max, pv_sum), type = "Fisher")

  return(list(pv.DMS = pv_DMS, pv.Max = pv_max, pv.Sum = pv_sum))
}
