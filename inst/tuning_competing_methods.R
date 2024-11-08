# Modification of code from the R package 'ocd: High-Dimensional Multiscale
#   Online Changepoint Detection'
# Original code by: Yudong Chen, Tengyao Wang, Richard J. Samworth
# URL: https://CRAN.R-project.org/package=ocd
# Description of modifications: The functions for choosing threshold values
# via Monte Carlo simulation are adjusted to attain a given false alarm
# probability, as opposed to the patience
# Date: 7th of November 2024

MC_ocd_FA <- function(dim, false_alarm_prob, N, beta=1, sparsity='auto', MC_reps,
                      est_length = 0, seed = 123){
  set.seed(seed)
  N = N - est_length
  peak_stat <- matrix(0, MC_reps, 3)
  colnames(peak_stat) <- c('diag','off_d','off_s')
  if (sparsity == 'sparse') peak_stat <- peak_stat[,-2]
  if (sparsity == 'dense') peak_stat <- peak_stat[,-3]

  # run MC_reps simulations for peak statistics of S_diag, S_{off,d} and S_{off,s}
  cat("Running MC simulation for OCD\n")
  mean_est = rep(0,dim)
  if(est_length>0){
    y_est = matrix(rnorm(dim*est_length), nrow = dim, ncol = est_length)
    mean_est = rowSums(y_est)/est_length
  }

  for (rep in 1:MC_reps){
    if(rep %%100==0){
      cat("Iteration: ", rep, "\n")
    }
    A <- matrix(0, dim, 1)
    tail <- matrix(0, dim, floor(log2(dim))*2+4)

    for (i in 1:N){
      x_new <- rnorm(dim)
      ret <- ocd::ocd_update(x_new - mean_est, A, tail, beta, sparsity)
      A <- ret$A; tail <- ret$tail
      peak_stat[rep,] <- pmax(peak_stat[rep,], ret$stat)
    }
  }
  qq = 1-false_alarm_prob/ dim(peak_stat)[2] # Bonferroni correction
  # compute the MC thresholds from the peak statistics
  thresh_est <- function(v) quantile(sort(v), qq)
  th_individual <- apply(peak_stat, 2, thresh_est)
  th <- th_individual
  names(th) <- colnames(peak_stat)
  return(th)
}

MC_Mei_FA <- function(dim, false_alarm_prob, N, beta = 1, b=beta/sqrt(dim), MC_reps,
                      est_length = 0, seed = 123){
  set.seed(seed)
  N = N - est_length
  peak_stat <- matrix(0, MC_reps, 2)
  colnames(peak_stat) <- c('max','sum')

  # run MC_reps simulations for peak statistics
  cat("Running MC simulation for Mei\n")
  mean_est = rep(0,dim)
  if(est_length>0){
    y_est = matrix(rnorm(dim*est_length), nrow = dim, ncol = est_length)
    mean_est = rowSums(y_est)/est_length
  }
  for (rep in 1:MC_reps){
    if(rep %%100==0){
      cat("Iteration: ", rep, "\n")
    }
    R <- matrix(0, dim, 2)

    for (i in 1:N){
      x_new <- rnorm(dim)
      ret <- ocd::Mei_update(x_new - mean_est, R, b)
      R <- ret$R
      peak_stat[rep,] <- pmax(peak_stat[rep,], ret$stat)
    }
  }

  qq = 1-false_alarm_prob/ dim(peak_stat)[2] # bonferroni correction

  # compute the MC thresholds from the peak statistics
  thresh_est <- function(v) quantile(sort(v), qq)
  th_individual <- apply(peak_stat, 2, thresh_est)

  th <- th_individual
  names(th) <- colnames(peak_stat)
  return(th)

}

MC_XS_FA <- function(dim, false_alarm_prob, N, p0=1/sqrt(dim), w=200, MC_reps,
                     est_length = 0, seed=123){
  set.seed(seed)
  N = N - est_length
  peak_stat <- rep(-Inf, MC_reps)

  cat("Running MC simulation for XS\n")
  mean_est = rep(0,dim)
  if(est_length>0){
    y_est = matrix(rnorm(dim*est_length), nrow = dim, ncol = est_length)
    mean_est = rowSums(y_est)/est_length
  }
  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps){
    if(rep %%100==0){
      cat("Iteration: ", rep, "\n")
    }
    X_recent <- CUSUM <- matrix(0, dim, w)

    for (i in 1:N){
      x_new <- rnorm(dim)
      ret <- ocd::XS_update(x_new - mean_est, X_recent, CUSUM, p0, w)
      X_recent <- ret$X_recent; CUSUM <- ret$CUSUM
      peak_stat[rep] <- pmax(peak_stat[rep], ret$stat)
    }
  }

  qq = 1-false_alarm_prob

  # compute the MC thresholds from the peak statistics
  th <- quantile(sort(peak_stat), qq)
  return(th)
}

MC_Chan_FA <- function(dim, false_alarm_prob, N,  p0=1/sqrt(dim), w=200,
                    lambda=sqrt(8)-2, MC_reps, est_length = 0, seed = 123){
  set.seed(seed)
  N = N - est_length
  peak_stat <- rep(-Inf, MC_reps)

  cat("Running MC simulation for CHAN\n")
  mean_est = rep(0,dim)
  if(est_length>0){
    y_est = matrix(rnorm(dim*est_length), nrow = dim, ncol = est_length)
    mean_est = rowSums(y_est)/est_length
  }
  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps){
    if(rep %%100==0){
      cat("Iteration: ", rep, "\n")
    }
    X_recent <- CUSUM <- matrix(0, dim, w)

    for (i in 1:N){
      x_new <- rnorm(dim)
      ret <- ocd::Chan_update(x_new - mean_est, X_recent, CUSUM, p0, w, lambda)
      X_recent <- ret$X_recent; CUSUM <- ret$CUSUM
      peak_stat[rep] <- pmax(peak_stat[rep], ret$stat)
    }
  }

  qq = 1-false_alarm_prob

  # compute the MC thresholds from the peak statistics
  th <- quantile(sort(peak_stat), qq)
  return(th)
}
