# R version of the C routine compute_CUSUM_and_Threshold
compute_CUSUM_and_Threshold_R <- function(cumsums, cumsum_vec, t, grid,
                                          as, nu_as, estimate_mean) {
  len <- length(grid)


  if (estimate_mean) {
    CUSUM <- sweep(cumsums, 2, sqrt(grid / (t * (t - grid))), FUN = "*") -
      sweep(cumsum_vec - cumsums, 2, sqrt((t - grid) / (t * grid)), FUN = "*")
  } else {
    CUSUM <- sweep(cumsum_vec - cumsums, 2, sqrt(grid), FUN = "/")
  }

  CUSUM_sq <- CUSUM^2
  CUSUM_abs <- abs(CUSUM)

  teststat <- matrix(0, nrow = length(as), ncol = len)
  for (h in seq_along(as)) {
    teststat[h, ] <- colSums((CUSUM_sq - nu_as[h]) * (CUSUM_abs > as[h]))
  }

  return(teststat)
}
