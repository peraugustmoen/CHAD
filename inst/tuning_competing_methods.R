# Modification of code from the R package 'ocd: High-Dimensional Multiscale
# Online Changepoint Detection' (https://CRAN.R-project.org/package=ocd) and
# Github Repository "Focused" (https://github.com/guillemr/Focused).
# Original code by: Yudong Chen, Tengyao Wang, Richard J. Samworth (ocd) and
# Gaetano Romano, Liudmila Pishchagina and Guillem Rigaill (Focused).

# Description of modifications: The functions for choosing threshold values
# via Monte Carlo simulation are adjusted to attain a given false alarm
# probability.


MC_ocd_FA <- function(
  dim, false_alarm_prob, N, beta = 1,
  sparsity = "auto", MC_reps,
  est_length = 0, seed = 123
) {
  set.seed(seed)
  N <- N - est_length
  peak_stat <- matrix(0, MC_reps, 3)
  colnames(peak_stat) <- c("diag", "off_d", "off_s")
  if (sparsity == "sparse") peak_stat <- peak_stat[, -2]
  if (sparsity == "dense") peak_stat <- peak_stat[, -3]

  # run MC_reps simulations for peak statistics of S_diag,
  # S_{off,d} and S_{off,s}
  cat("Running MC simulation for OCD\n")
  mean_est <- rep(0, dim)


  for (rep in 1:MC_reps) {
    if (rep %% 100 == 0) {
      cat("Iteration: ", rep, "\n")
    }
    if (est_length > 0) {
      y_est <- matrix(rnorm(dim * est_length), nrow = dim, ncol = est_length)
      mean_est <- rowSums(y_est) / est_length
    }
    A <- matrix(0, dim, 1)
    tail <- matrix(0, dim, floor(log2(dim)) * 2 + 4)

    for (i in 1:N) {
      x_new <- rnorm(dim)
      ret <- ocd::ocd_update(x_new - mean_est, A, tail, beta, sparsity)
      A <- ret$A
      tail <- ret$tail
      peak_stat[rep, ] <- pmax(peak_stat[rep, ], ret$stat)
    }
  }
  qq <- 1 - false_alarm_prob / dim(peak_stat)[2] # Bonferroni correction
  # compute the MC thresholds from the peak statistics
  thresh_est <- function(v) quantile(sort(v), qq)
  th_individual <- apply(peak_stat, 2, thresh_est)
  th <- th_individual
  names(th) <- colnames(peak_stat)
  return(th)
}

MC_Mei_FA <- function(
  dim, false_alarm_prob, N, beta = 1,
  b = beta / sqrt(dim), MC_reps,
  est_length = 0, seed = 123
) {
  set.seed(seed)
  N <- N - est_length
  peak_stat <- matrix(0, MC_reps, 2)
  colnames(peak_stat) <- c("max", "sum")

  # run MC_reps simulations for peak statistics
  cat("Running MC simulation for Mei\n")
  mean_est <- rep(0, dim)

  for (rep in 1:MC_reps) {
    if (rep %% 100 == 0) {
      cat("Iteration: ", rep, "\n")
    }

    if (est_length > 0) {
      y_est <- matrix(rnorm(dim * est_length), nrow = dim, ncol = est_length)
      mean_est <- rowSums(y_est) / est_length
    }
    R <- matrix(0, dim, 2)

    for (i in 1:N) {
      x_new <- rnorm(dim)
      ret <- ocd::Mei_update(x_new - mean_est, R, b)
      R <- ret$R
      peak_stat[rep, ] <- pmax(peak_stat[rep, ], ret$stat)
    }
  }

  qq <- 1 - false_alarm_prob / dim(peak_stat)[2] # bonferroni correction

  # compute the MC thresholds from the peak statistics
  thresh_est <- function(v) quantile(sort(v), qq)
  th_individual <- apply(peak_stat, 2, thresh_est)

  th <- th_individual
  names(th) <- colnames(peak_stat)
  return(th)
}

MC_XS_FA <- function(
  dim, false_alarm_prob, N, p0 = 1 / sqrt(dim),
  w = 200, MC_reps,
  est_length = 0, seed = 123
) {
  set.seed(seed)
  N <- N - est_length
  peak_stat <- rep(-Inf, MC_reps)

  cat("Running MC simulation for XS\n")
  mean_est <- rep(0, dim)

  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps) {
    if (rep %% 100 == 0) {
      cat("Iteration: ", rep, "\n")
    }
    if (est_length > 0) {
      y_est <- matrix(rnorm(dim * est_length), nrow = dim, ncol = est_length)
      mean_est <- rowSums(y_est) / est_length
    }
    X_recent <- CUSUM <- matrix(0, dim, w)

    for (i in 1:N) {
      x_new <- rnorm(dim)
      ret <- ocd::XS_update(x_new - mean_est, X_recent, CUSUM, p0, w)
      X_recent <- ret$X_recent
      CUSUM <- ret$CUSUM
      peak_stat[rep] <- pmax(peak_stat[rep], ret$stat)
    }
  }

  qq <- 1 - false_alarm_prob

  # compute the MC thresholds from the peak statistics
  th <- quantile(sort(peak_stat), qq)
  return(th)
}

MC_Chan_FA <- function(
  dim, false_alarm_prob, N,
  p0 = 1 / sqrt(dim), w = 200,
  lambda = sqrt(8) - 2, MC_reps,
  est_length = 0, seed = 123
) {
  set.seed(seed)
  N <- N - est_length
  peak_stat <- rep(-Inf, MC_reps)

  cat("Running MC simulation for CHAN\n")
  mean_est <- rep(0, dim)

  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps) {
    if (rep %% 100 == 0) {
      cat("Iteration: ", rep, "\n")
    }
    X_recent <- CUSUM <- matrix(0, dim, w)
    if (est_length > 0) {
      y_est <- matrix(rnorm(dim * est_length), nrow = dim, ncol = est_length)
      mean_est <- rowSums(y_est) / est_length
    }
    for (i in 1:N) {
      x_new <- rnorm(dim)
      ret <- ocd::Chan_update(
        x_new - mean_est, X_recent,
        CUSUM, p0, w, lambda
      )
      X_recent <- ret$X_recent
      CUSUM <- ret$CUSUM
      peak_stat[rep] <- pmax(peak_stat[rep], ret$stat)
    }
  }

  qq <- 1 - false_alarm_prob

  # compute the MC thresholds from the peak statistics
  th <- quantile(sort(peak_stat), qq)
  return(th)
}

## Source the R implementation of (sparsity-adaptive) MdFOCuS from the
## installed CHAD package:
source(system.file("MdFocus_MeanGaussian_md.R", package = "CHAD"))
MC_mdfocus_FA <- function(dim, false_alarm_prob, N, MC_reps, seed = 123) {
  cat("Running MC simulation for mdfocus\n")
  set.seed(seed)
  p <- dim

  # sparsity_levels = 2^1, 2^2, ..., 2^{floor(log_2(p))}:
  sparsity_levels <- 2^seq_len(floor(log2(p)))

  ll <- list() # storing maximum statistics for each MC sample

  for (k in 1:MC_reps) {
    if (k %% 100 == 0) {
      cat("Iteration: ", k, "\n")
    }
    data <- data.frame(matrix(rnorm(N * p), nrow = N, ncol = p))
    res <- FocusCH_HighDim(data,
      get_opt_cost = \(...) get_partial_opt(...,
        cost = cost_lr_partial0, which_par = sparsity_levels
      ),
      # dim_indexes = as.list(1:ncol(data)),
      # common_ratio_step = 1.3,
      threshold = rep(Inf, 2 + length(sparsity_levels))
    )
    res <- -(res$opt.cost |> reduce(rbind)) |> apply(2, max)
    ll[[k]] <- res
  }
  rr2 <- reduce(ll, rbind)

  quant <- false_alarm_prob / ncol(rr2) # quantile with Bonferroni correction

  thresholds_mdfocus <- apply(rr2, 2, quantile, probs = 1 - quant)

  return(thresholds_mdfocus)
}

MC_mdfocus_nonzeromean_FA <- function(dim, false_alarm_prob, N, MC_reps, seed = 123) {
  cat("Running MC simulation for mdfocus\n")
  set.seed(seed)
  p <- dim
  # sparsity_levels = 2^1, 2^2, ..., 2^{floor(log_2(p))}:
  sparsity_levels <- 2^seq_len(floor(log2(p)))
  ll <- list()
  for (k in 1:MC_reps) {
    if (k %% 100 == 0) {
      cat("Iteration: ", k, "\n")
    }
    data <- data.frame(matrix(rnorm(N * p), nrow = N, ncol = p))
    res <- FocusCH_HighDim(data,
      get_opt_cost = \(...) get_partial_opt(...,
        which_par = sparsity_levels
      ),
      # dim_indexes = as.list(1:ncol(data)),
      # common_ratio_step = 1.3,
      threshold = rep(Inf, 2 + length(sparsity_levels))
    )
    res <- -(res$opt.cost |> reduce(rbind)) |> apply(2, max)
    ll[[k]] <- res
  }
  rr2 <- reduce(ll, rbind)
  quant <- false_alarm_prob / ncol(rr2)
  thresholds_mdfocus <- apply(rr2, 2, quantile, probs = 1 - quant)
  return(thresholds_mdfocus)
}
