#' @useDynLib CHAD compute_CUSUM_and_Threshold
#' @importFrom stats dnorm pnorm quantile rnorm setNames
#' @importFrom utils read.csv



#' @title Constructor method for the CHAngepoint Detector (CHAD) class
#' @param p Dimension of the data
#' @param method Test statistic to be used, either \code{"mean"} or
#' \code{"covariance"}. See details below
#' @param leading_constant A numeric vector specifying the leading constants for
#' the penalty function, or the string "MC" for Monte Carlo
#' simulation of these
#' If "MC" is given as input, the arguments \code{false_alarm_prob},
#' \code{N} and \code{MC_reps} are used in the MC simulation
#' @param baseline_sd If \code{method = "mean"}, a numeric vector specifying
#' the baseline standard deviation for each coordinate of the data
#' Only applicable when \code{method = "mean"}
#' @param baseline_operatornorm If \code{method = "covariance"}, a numeric value
#' specifying the baseline operator norm of the pre-change covariance matrix.
#' If \code{baseline_operatornorm = NA}, the operator norm is estimated
#' in real-time
#' @param estimate_mean Boolean indicating whether the data should be estimated
#' in real time, or taken to be zero
#' @param constant_penalty Boolean indicating whether critical values (and
#' tresholds) should be constant with respect to the sample size
#' @param min_prechange_obs Minimum number of of pre-change observations
#' required for a candidate changepoint location to be considered.
#' Note that if \code{method = "covariance"}, and
#' \code{baseline_operatornorm = NA}, it is reccomended to set
#' \code{min_prechange_obs = 10} to avoid highly variable estimates of the
#' baseline noise level.
#' @param false_alarm_prob False alarm probability for MC simulations
#'  (between 0 and 1). Only required when \code{leading_constant='MC'}
#' @param MC_reps Number of Monte Carlo repetitions to estimate the
#' leading constant(s) of the penalty function(s)
#' Only required when \code{leading_constant='MC'}
#' @param N Maximum sample size of which Monte Carlo simulations are performed.
#' Only required when \code{leading_constant='MC'}
#' @return An object of class 'CHAD' with the following attributes:
#' \itemize{
#' \item \code{class} - subclass and class
#' \item \code{p} - Data dimension
#' \item \code{method} - Method used for online changepoint detection
#' \item \code{leading_constant} - A vector of leading constant(s) for the
#' penalty function(s)
#' \item \code{estimate_mean} - Boolean indicating whether mean estimation
#' should be performed
#' \item \code{n_obs} - Number of currently monitored observations
#' \item \code{n_obs_tot} - Number of observations processed in total
#' \item \code{baseline_sd} - Vector of standard deviations (for mean changes)
#' \item \code{baseline_operatornorm} - Baseline operator norm of the covariance
#' matrix (for covariance changes)
#' \item \code{statistics} - A vector of test statistics
#' \item \code{overall_max_statistics} - Vector of largest statistics
#' recorded so far
#' \item \code{allstats} - A list of various auxiliary statistics
#' \item \code{constant_penalty} - Boolean indicating whether penalties and
#' thresholds are constant with respect to the sample size
#' \item \code{min_prechange_obs} - Minimum number of pre-change observations
#' \item \code{status} - Either monitoring or changepoint detected
#' }
#' @details When \code{method = 'mean'}, the test statistic of
#' Liu et al. (2021) is used to test for a changepoint. If
#' \code{method = 'covariance'}, the test statistic of Moen (2024) is used.
#' @examples
#' detector_mean <- CHAD::CHAD(
#'   p = 100, method = "mean",
#'   leading_constant = c(2, 2)
#' )
#' detector_covariance <- CHAD(
#'   p = 100, method = "covariance",
#'   leading_constant = c(2)
#' )
#' @references
#' \itemize{
#' \item Liu et al. (2021) Minimax rates in sparse, high-dimensional change
#'  point detection. Ann. Statist. 49(2): 1081-1112.
#'  \item Moen (2024) Minimax rates in variance and covariance changepoint
#'  testing. arXiv:2405.07757.
#' }
#' @export
CHAD <- function(p,
                 method = c("mean", "covariance"),
                 leading_constant,
                 baseline_sd = 1,
                 baseline_operatornorm = NA,
                 estimate_mean = TRUE,
                 constant_penalty = FALSE,
                 min_prechange_obs = 1,
                 false_alarm_prob = 0.05,
                 MC_reps = 100,
                 N = 1000) {
  if (identical(leading_constant, "MC")) {
    leading_constant <- switch(method,
      mean = MC_mean(
        p = p,
        false_alarm_prob = false_alarm_prob,
        constant_penalty = constant_penalty,
        MC_reps = MC_reps,
        N = N,
        estimate_mean = estimate_mean,
        min_prechange_obs = min_prechange_obs,
        baseline_sd = baseline_sd
      ),
      covariance = MC_covariance(
        p = p,
        false_alarm_prob = false_alarm_prob,
        constant_penalty = constant_penalty,
        MC_reps = MC_reps,
        N = N,
        estimate_mean = estimate_mean,
        baseline_operatornorm = baseline_operatornorm,
        min_prechange_obs = min_prechange_obs
      )
    )
  }

  detector <- switch(method,
    mean = new_mean(
      p = p,
      leading_constant = leading_constant,
      baseline_sd = baseline_sd,
      constant_penalty = constant_penalty,
      estimate_mean = estimate_mean,
      min_prechange_obs = min_prechange_obs
    ),
    covariance = new_covariance(
      p = p,
      leading_constant = leading_constant,
      estimate_mean = estimate_mean,
      constant_penalty = constant_penalty,
      baseline_operatornorm = baseline_operatornorm,
      min_prechange_obs = min_prechange_obs
    )
  )
  return(detector)
}

#' @title Constructor of subclass 'meanDetector' of class 'CHAD'
#' @param p Dimension of the data
#' @param leading_constant A numeric vector of length 2. The first entry is
#' the leading constant of the penalty function \eqn{r(k,p,t)} for
#' \eqn{k\geq\sqrt{p\log n}}. The second entry is the leading constant for the
#' penalty function when \eqn{k\leq\sqrt{p\log n}}
#' @param baseline_sd Baseline standard deviation for each component of the
#' data
#' @param constant_penalty Boolean indicating whether leading constants and
#' critical values should be constant with respect to the sample size or not
#' @param min_prechange_obs Minimum number of observations before a changepoint
#' @param estimate_mean Boolean indicating whether mean-centering should be
#' performed or if pre-change mean should be assumed to be zero
#' @return An object of subclass 'meanDetector' of class 'CHAD'
#' @details Uses the test statistic of Liu et al. (2021) to test for a change
#' in mean online. The penalty function is given by
#' \eqn{r(k,p,t) = log(t) + k\log(1 + \sqrt{p\log t}/k)} for
#' \eqn{k<\sqrt{p\log t}} and \eqn{r(k,p,t) = \sqrt{p\log t} + \log t} for
#' \eqn{k\geq\sqrt{p\log t}}. The critical values are given by
#' \eqn{r(k,p,t)} multiplied by the leading constant(s) specified in the
#' \code{leading_constant} argument.
#' @examples
#' detector <- CHAD::new_mean(
#'   p = 100, leading_constant = c(7, 6), baseline_sd = 1,
#'   constant_penalty = FALSE, estimate_mean = TRUE,
#'   min_prechange_obs = 1
#' )
#' @references
#' \itemize{
#' \item Liu et al. (2021) Minimax rates in sparse, high-dimensional change
#'  point detection. Ann. Statist. 49(2): 1081-1112.}
#' @export
new_mean <- function(p, leading_constant, baseline_sd, constant_penalty,
                     min_prechange_obs, estimate_mean) {
  grid <- c()

  stats <- c(0, 0) # one for sparse, one for dense regime
  overall_max_stats <- c(0, 0)
  allstats <- list() # various auxiliary statistics

  detector <- structure(list(),
    class = c("meanDetector", "CHAD"),
    p = p,
    method = "mean",
    estimate_mean = estimate_mean,
    leading_constant = leading_constant,
    n_obs = 0,
    n_obs_tot = 0,
    grid = grid,
    baseline_sd = baseline_sd,
    statistics = stats,
    overall_max_statistics = overall_max_stats,
    allstats = allstats,
    constant_penalty = constant_penalty,
    min_prechange_obs = min_prechange_obs,
    status = "monitoring"
  )
  return(detector)
}

#' @title Constructor of subclass 'covarianceDetector' of class 'CHAD'
#' @param p Data dimension
#' @param leading_constant A single numeric value acting
#'  as the leading constant for the penalty function for the covariance
#' change detector.
#' @param estimate_mean Boolean indicating whether mean-centering should be
#' performed.
#' @param constant_penalty Boolean indicating whether leading constants and
#' critical values should be constant with respect to the sample size or not
#' @param min_prechange_obs Minimum number of pre-change observations before
#' a candidate changepoint location is considered
#' @param baseline_operatornorm A numeric value
#' specifying the baseline operator norm of the covariance matrix
#' for each coordinate of the data. If \code{baseline_operatornorm = NA}, the
#' operator norm is estimated in real-time
#' @return An object of subclass 'covarianceDetector' of class 'CHAD'
#' @details Implements the test statistic of Moen (2024) for dense
#' changes in covariance. The penalty function is given by
#' \eqn{h(p,t,g) = \max\left(\frac{p + \log t}{g},
#' \sqrt{\frac{p + \log t}{g}}\right)}, and the critical values are given by
#' \eqn{h(p,t,g)} multiplied by the leading constant specified in the
#' \code{leading_constant} argument.
#' @examples
#' detector <- CHAD::new_covariance(
#'   p = 10, leading_constant = c(5),
#'   estimate_mean = TRUE, constant_penalty = FALSE,
#'   baseline_operatornorm = NA, min_prechange_obs = 1
#' )
#' @references
#' \itemize{
#'  \item Moen (2024) Minimax rates in variance and covariance changepoint
#'  testing. arXiv:2405.07757.
#' }
#' @export
new_covariance <- function(p,
                           leading_constant,
                           estimate_mean,
                           constant_penalty,
                           min_prechange_obs,
                           baseline_operatornorm) {
  grid <- c()
  stats <- c(1)
  overall_max_stats <- c(1)

  detector <- structure(list(),
    class = c("covarianceDetector", "CHAD"),
    p = p,
    method = "covariance",
    leading_constant = leading_constant,
    estimate_mean = estimate_mean,
    n_obs = 0,
    n_obs_tot = 0,
    grid = grid,
    statistics = stats,
    overall_max_statistics = overall_max_stats,
    constant_penalty = constant_penalty,
    baseline_operatornorm = baseline_operatornorm,
    min_prechange_obs = min_prechange_obs,
    status = "monitoring"
  )
  return(detector)
}




##### Methods for the 'CHAD' class #####

#' @title Accessor functions for the 'CHAD' class
#' @param detector Object of class 'CHAD'
#' @name accessor
NULL

#' @rdname accessor
#' @export
data_dim <- function(detector) attr(detector, "p")

#' @rdname accessor
#' @export
method <- function(detector) attr(detector, "method")

#' @rdname accessor
#' @export
n_obs <- function(detector) attr(detector, "n_obs")

#' @rdname accessor
#' @export
n_obs_tot <- function(detector) attr(detector, "n_obs_tot")

#' @rdname accessor
#' @export
min_prechange_obs <- function(detector) attr(detector, "min_prechange_obs")

#' @rdname accessor
#' @export
false_alarm_probability <- function(detector) {
  attr(
    detector,
    "false_alarm_probability"
  )
}

#' @rdname accessor
#' @export
leading_constant <- function(detector) attr(detector, "leading_constant")

#' @rdname accessor
#' @export
baseline_sd <- function(detector) attr(detector, "baseline_sd")

#' @rdname accessor
#' @export
statistics <- function(detector) attr(detector, "statistics")

#' @rdname accessor
#' @export
all_statistics <- function(detector) attr(detector, "allstats")

#' @rdname accessor
#' @export
status <- function(detector) attr(detector, "status")

#' @rdname accessor
#' @export
bool_estimate_mean <- function(detector) attr(detector, "estimate_mean")



#' Reset an object of the class 'CHAD'
#' @param detector Object of class 'CHAD'
#' @return Updated object `detector`
#' @export
reset <- function(detector) UseMethod("reset")


#' @title Compute maximum ratio between the test statistic(s) divided by the
#' penalty function(s) / critical values(s)
#' @param detector Object of class 'CHAD'
#' @return maximum of the ratio between the current test statistics and their
#' respective thresholds.
#' @export
normalizedStatistics <- function(detector) {
  max(statistics(detector) / leading_constant(detector))
}

#' @title Check if a changepoint has occurred
#' @param detector Object of class 'CHAD'
#' @return Updated object \code{detector}
#' @details The function \code{\link{normalizedStatistics}} is applied to check
#' if a changepoint should be declared.
#' @export
checkChange <- function(detector) {
  if (normalizedStatistics(detector) >= 1) {
    n <- n_obs(detector)
    attr(detector, "status") <- setNames(n, "declared at")
    cat("Changepoint declared at time =", n, "\n")
  }
  return(detector)
}



#' @title Function for processing a new data point
#' @param detector Object of class 'CHAD'
#' @param y_new A new data point. It must be of the same dimension as
#' specified in the \code{p} attribute of \code{detector}.
#' @return Updated object \code{detector}
#' @export
getData <- function(detector, y_new) UseMethod("getData")


#' @describeIn getData Process a new data point for subclass 'meanDetector'
#' @export
getData.meanDetector <- function(detector, y_new) {
  n_obs <- attr(detector, "n_obs") + 1
  n_obs_tot <- attr(detector, "n_obs_tot") + 1
  attr(detector, "n_obs") <- n_obs
  attr(detector, "n_obs_tot") <- n_obs_tot
  p <- data_dim(detector)

  y_new <- y_new / baseline_sd(detector)


  new_grid <- get_grid(n_obs)

  if (n_obs == 1) {
    old_cumulative_sum <- rep(0, p)
  } else {
    old_cumulative_sum <- attr(detector, "allstats")$cumsum
  }

  # Store Y_1 + ... + Y_t
  attr(detector, "allstats")$cumsum <- old_cumulative_sum + y_new


  # update cumulative sums over the grid
  if (n_obs == 1) {
    attr(detector, "allstats")$cumsums <- matrix(y_new, nrow = p, ncol = 1)
  } else {
    newinds <- (attr(detector, "grid") + 1) %in% new_grid
    len <- sum(newinds) + 1
    new_cumsums <- matrix(NA, nrow = p, ncol = len)
    new_cumsums[, 1] <- old_cumulative_sum
    if (len > 1) {
      new_cumsums[, 2:len] <- attr(detector, "allstats")$cumsums[, newinds]
    }

    attr(detector, "allstats")$cumsums <- new_cumsums
  }

  attr(detector, "grid") <- new_grid[]

  if (n_obs < 2) {
    return(detector)
  }

  if (attr(detector, "constant_penalty")) {
    mean_penalty_and_threshold <- mean_penalty_and_threshold(2, p)
  } else {
    mean_penalty_and_threshold <- mean_penalty_and_threshold(n_obs_tot, p)
  }

  as <- mean_penalty_and_threshold[[1]]
  nu_as <- mean_penalty_and_threshold[[2]]
  r_values <- mean_penalty_and_threshold[[3]]
  s_values <- mean_penalty_and_threshold[[4]]

  cumsum <- attr(detector, "allstats")$cumsum
  estimate_mean <- attr(detector, "estimate_mean")

  ## Now that the cumulative sums, thresholding values and penalty functions
  ## are computed, a C function is called, which computes CUSUMs, thresholds
  ## and sums these, and divide by the penalty function. The matrix 'Amatrix'
  ## is returned. Each column represents a candidate changepoint location as
  ## as implied by the grid g. Each row represents a candidate sparsity level.
  Amatrix <- .Call(
    compute_CUSUM_and_Threshold, new_cumsums, cumsum, as.integer(n_obs),
    as.integer(new_grid), as.integer(length(new_grid)), as.integer(p),
    as.numeric(as), as.numeric(nu_as), as.integer(length(as)),
    as.integer(estimate_mean)
  )

  Amatrix <- matrix(Amatrix,
    nrow = length(as), ncol = length(new_grid),
    byrow = TRUE
  )
  attr(detector, "allstats")$Amatrix <- Amatrix

  # divide each A by the corresponding penalty value in the r_value vector
  A_scaled <- sweep(Amatrix, 1, r_values, FUN = "/")
  attr(detector, "allstats")$A_scaled <- A_scaled[, ]

  admissable_g_inds <- rep(TRUE, length(new_grid))
  if (min_prechange_obs(detector) > 1) {
    admissable_g_inds <- (n_obs - new_grid) >= min_prechange_obs(detector)

    if (sum(admissable_g_inds) == 0) {
      return(detector)
    }
    new_grid <- new_grid[admissable_g_inds]
    A_scaled <- A_scaled[, admissable_g_inds]
  }




  attr(detector, "statistics")[1] <- max(A_scaled[1, ])
  if (length(as) == 2) {
    attr(detector, "statistics")[2] <- max(A_scaled[2, ])
  } else {
    attr(detector, "statistics")[2] <- max(A_scaled[2:length(as), ])
  }

  largest_index <- which(A_scaled == max(A_scaled), arr.ind = TRUE)

  attr(detector, "allstats")$argmax_s <- s_values[largest_index[1]]
  attr(detector, "allstats")$argmax_g <- new_grid[largest_index[2]]

  if (attr(detector, "statistics")[1] > attr(
    detector,
    "overall_max_statistics"
  )[1]) {
    largest_index <- which(A_scaled[1, ] == max(A_scaled[1, ]), arr.ind = TRUE)
    attr(detector, "overall_max_statistics")[1] <-
      attr(detector, "statistics")[1]
    attr(detector, "allstats")$overall_argmax_g_dense <-
      new_grid[largest_index[1]]
    attr(detector, "allstats")$overall_argmax_n_obs_dense <- n_obs
  }
  if (attr(detector, "statistics")[2] > attr(
    detector,
    "overall_max_statistics"
  )[2]) {
    largest_index <-
      which(A_scaled[-1, ] == max(A_scaled[-1, ]), arr.ind = TRUE)
    attr(detector, "overall_max_statistics")[2] <-
      attr(detector, "statistics")[2]
    attr(detector, "allstats")$overall_argmax_s_sparse <-
      s_values[largest_index[1] + 1]
    attr(detector, "allstats")$overall_argmax_g_sparse <-
      new_grid[largest_index[2]]
    attr(detector, "allstats")$overall_argmax_n_obs_sparse <-
      n_obs
  }







  if (status(detector) == "monitoring") detector <- checkChange(detector)
  return(detector)
}

#' @describeIn getData Process a new data point for subclass
#' "covarianceDetector"
#' @export
getData.covarianceDetector <- function(detector, y_new) {
  n_obs <- attr(detector, "n_obs") + 1
  n_obs_tot <- attr(detector, "n_obs_tot") + 1
  attr(detector, "n_obs") <- n_obs
  attr(detector, "n_obs_tot") <- n_obs_tot
  p <- data_dim(detector)

  estimate_mean <- attr(detector, "estimate_mean")

  new_grid <- get_grid(n_obs)

  if (n_obs == 1) {
    old_cumulative_sum <- matrix(0, nrow = p, ncol = p)
    old_cumulative_sum_mean <- rep(0, p)
  } else {
    old_cumulative_sum <- attr(detector, "allstats")$cumsum
    old_cumulative_sum_mean <- attr(detector, "allstats")$cumsum_mean
  }

  # Store Y_1Y_1^\top + ... + Y_tY_t^\top
  attr(detector, "allstats")$cumsum <- old_cumulative_sum + y_new %*% t(y_new)
  attr(detector, "allstats")$cumsum_mean <- old_cumulative_sum_mean +
    y_new



  # update cumulative sums over the grid
  if (n_obs == 1) {
    attr(detector, "allstats")$cumsums <- y_new %*% t(y_new)
    attr(detector, "allstats")$cumsums_mean <- y_new
    attr(detector, "allstats")$log2cumsums <- y_new %*% t(y_new)
    attr(detector, "allstats")$log2cumsums_mean <- y_new
  } else {
    newinds <- (attr(detector, "grid") + 1) %in% new_grid
    len <- sum(newinds) + 1
    new_cumsums <- array(NA, dim = c(len, p, p))
    new_cumsums[1, , ] <- old_cumulative_sum
    if (len > 1) {
      new_cumsums[2:len, , ] <- attr(detector, "allstats")$cumsums[newinds, , ]
    }
    attr(detector, "allstats")$cumsums <- new_cumsums

    # means
    new_cumsums_mean <- matrix(NA, nrow = p, ncol = len)
    new_cumsums_mean[, 1] <- old_cumulative_sum_mean
    if (len > 1) {
      new_cumsums_mean[, 2:len] <-
        attr(detector, "allstats")$cumsums_mean[, newinds]
    }
    attr(detector, "allstats")$cumsums_mean <- new_cumsums_mean

    # store first log_2(n_obs) cumulative sum as well
    if (floor(log(n_obs, base = 2)) > floor(log(n_obs - 1, base = 2))) {
      len <- floor(log(n_obs, base = 2)) + 1
      newlog2cumsum <- array(NA, dim = c(len, p, p))
      newlog2cumsum[1:(len - 1), , ] <- attr(detector, "allstats")$log2cumsums
      newlog2cumsum[len, , ] <- attr(detector, "allstats")$cumsum
      attr(detector, "allstats")$log2cumsums <- newlog2cumsum

      newlog2cumsum_mean <- matrix(NA, nrow = p, ncol = len)
      newlog2cumsum_mean[, 1:(len - 1)] <-
        attr(detector, "allstats")$log2cumsums_mean
      newlog2cumsum_mean[, len] <- attr(detector, "allstats")$cumsum_mean
      attr(detector, "allstats")$log2cumsums_mean <- newlog2cumsum_mean
    }
  }

  attr(detector, "grid") <- new_grid

  if (n_obs == 1) {
    return(detector)
  }

  ## compute test statistic
  teststats <- rep(NA, length(new_grid))
  S_t <- attr(detector, "allstats")$cumsum
  S_t_min_g <- attr(detector, "allstats")$cumsums
  for (i in 1:length(new_grid)) {
    g <- new_grid[i]
    if (g < min_prechange_obs(detector)) {
      teststats[i] <- 0
      next
    }
    if (estimate_mean) {
      if (g < 2) {
        teststats[i] <- 0
        next
      }
    }
    Sigmahat_right <- (S_t - S_t_min_g[i, , ]) / g
    Sigmahat_left <-
      attr(detector, "allstats")$log2cumsums[floor(log(g, base = 2)) + 1, , ] *
        2^(-floor(log(g, base = 2)))

    if (estimate_mean) {
      mean_right <- (attr(detector, "allstats")$cumsum_mean -
        attr(detector, "allstats")$cumsums_mean[, i]) / g
      Sigmahat_right <- Sigmahat_right - mean_right %*% t(mean_right)

      mean_left <-
        attr(detector, "allstats")$log2cumsums_mean[
          ,
          floor(log(g, base = 2)) + 1
        ] * 2^(-floor(log(g, base = 2)))

      Sigmahat_left <- Sigmahat_left - mean_left %*% t(mean_left)
    }



    pen <- 0
    if (attr(detector, "constant_penalty")) {
      pen <- max((p + log(2)) / g, sqrt((p + log(2)) / g))
    } else {
      pen <- max((p + log(n_obs_tot)) / g, sqrt((p + log(n_obs_tot)) / g))
    }
    if (is.na(attr(detector, "baseline_operatornorm"))) {
      sigmahat2 <- norm(Sigmahat_left, type = "2")
      teststats[i] <- max(norm(Sigmahat_left - Sigmahat_right, type = "2")) /
        sigmahat2 / pen
    } else {
      teststats[i] <- max(norm(Sigmahat_left - Sigmahat_right, type = "2")) /
        attr(detector, "baseline_operatornorm") / pen
    }
  }


  attr(detector, "statistics") <- max(teststats)
  attr(detector, "allstats")$teststats <- teststats
  attr(detector, "allstats")$argmax <- new_grid[which.max(teststats)]
  ## update overall max statistics
  if (attr(detector, "statistics") > attr(
    detector,
    "overall_max_statistics"
  )) {
    attr(detector, "overall_max_statistics") <- attr(detector, "statistics")
    attr(detector, "allstats")$overall_argmax <-
      attr(detector, "allstats")$argmax
    attr(detector, "allstats")$overall_argmax_n_obs <- n_obs
  }

  if (status(detector) == "monitoring") detector <- checkChange(detector)

  return(detector)
}


#' @title Reset changepoint detector
#' @param detector Object of class 'CHAD'
#' @return Updated object \code{detector}
#' @export
reset <- function(detector) UseMethod("reset")


#' @describeIn reset Reset object of subclass 'meanDetector'
#' @export
reset.meanDetector <- function(detector) {
  p <- data_dim(detector)

  attr(detector, "n_obs") <- 0
  attr(detector, "statistics") <- 0
  attr(detector, "allstats") <- list()
  attr(detector, "status") <- "monitoring"
  return(detector)
}

#' @describeIn reset Reset object of subclass 'covarianceDetector'
#' @export
reset.covarianceDetector <- function(detector) {
  p <- data_dim(detector)
  attr(detector, "n_obs") <- 0
  attr(detector, "statistics") <- c()
  attr(detector, "allstats") <- list()
  attr(detector, "status") <- "monitoring"
  return(detector)
}

#' @title MC simulation for meanDetector
#' @description Chooses the leading constants for the penalty function for
#' the 'meanDetector' subclass of 'CHAD'.
#' @param p Data dimension
#' @param false_alarm_prob Desired false alarm probability (between 0 and 1)
#' @param constant_penalty If \code{true}, the thresholding values and penalty
#' function do not grow with the sample size.
#' @param min_prechange_obs Minimum number of pre-change observations required
#' for a candidate changepoint location to be considered
#' @param MC_reps number of MC iterations
#' @param N maximum number of observed data points
#' @param estimate_mean Boolean indicating whether the pre-change mean should
#' be estimated, or assumed to be zero
#' @param seed seed for the RNG
#' @return Vector of length 2 with leading constants for the penalty function
#' for the meanDetector
#' @export
MC_mean <- function(p, false_alarm_prob, constant_penalty,
                    min_prechange_obs = 1, MC_reps, N,
                    estimate_mean = TRUE, seed = 123) {
  set.seed(seed)
  max_statistics <- matrix(NA, nrow = 2, ncol = MC_reps)
  cat("Running MC simulation for meanDetector\n")
  for (j in 1:MC_reps) {
    if (j %% 100 == 0) {
      cat("Iteration: ", j, "\n")
    }
    detector <- CHAD(p,
      method = "mean", leading_constant = c(1000, 1000),
      constant_penalty = constant_penalty,
      estimate_mean = estimate_mean,
      min_prechange_obs = min_prechange_obs
    )
    ys <- matrix(rnorm(N * p), nrow = p, ncol = N)
    for (i in 1:N) {
      detector <- getData(detector, ys[, i])
    }
    max_statistics[, j] <- attr(detector, "overall_max_statistics")
  }

  qq <- 1 - false_alarm_prob / 2 # Bonferroni correction

  return(c(
    quantile(max_statistics[1, ], qq),
    quantile(max_statistics[2, ], qq)
  ))
}

#' @title MC simulation for covarianceDetector
#' @description Chooses the leading constant for the penalty function for
#' the 'covarianceDetector' subclass of 'CHAD'.
#' @param p data dimension
#' @param false_alarm_prob desired false alarm probability
#' @param constant_penalty If \code{true}, the penalty function does not grow
#' with the sample size
#' @param MC_reps number of MC iterations
#' @param N maximum number of observed data points
#' @param estimate_mean Boolean indicating whether mean-centering should be
#' performed or if the mean of the data is assumed to be zero
#' @param min_prechange_obs Minimum number of observations before a candidate
#' changepoint location is considered
#' @param baseline_operatornorm A numeric value
#' specifying the baseline operator norm of the covariance matrix
#' for each coordinate of the data. If \code{baseline_operatornorm = NA}, the
#' operator norm is estimated in real-time
#' @param seed seed for the RNG
#' @return Numerical value of the leading constant for the penalty function
#' @export
MC_covariance <- function(p, false_alarm_prob, constant_penalty, MC_reps, N,
                          estimate_mean = TRUE,
                          min_prechange_obs = 1,
                          baseline_operatornorm = NA,
                          seed = 123) {
  set.seed(seed)
  max_statistics <- rep(NA, MC_reps)
  cat("Running MC simulation for covarianceDetector\n")
  for (j in 1:MC_reps) {
    if (j %% 100 == 0) {
      cat("Iteration: ", j, "\n")
    }
    detector <- CHAD(p,
      method = "covariance", leading_constant = c(100000),
      constant_penalty = constant_penalty,
      estimate_mean = estimate_mean,
      baseline_operatornorm = baseline_operatornorm,
      min_prechange_obs = min_prechange_obs
    )
    ys <- matrix(rnorm(N * p), nrow = p, ncol = N)
    for (i in 1:N) {
      detector <- getData(detector, ys[, i])
    }
    max_statistics[j] <- attr(detector, "overall_max_statistics")
  }

  qq <- 1 - false_alarm_prob
  return(quantile(max_statistics, qq))
}

#' Printing methods for the 'CHAD' class
#' @param x object of the 'CHAD' class
#' @param ... other arguments used in \code{print}
#' @export
print.CHAD <- function(x, ...) {
  detector <- x
  cat("Method:", method(detector), "\n\n")
  cat("Time =", n_obs(detector), "\n\n")

  mx <- rbind(statistics(detector), leading_constant(detector))
  mx <- round(mx, 3)
  row.names(mx) <- c("statistics", "critical values")
  print(mx)

  if (is.numeric(status(detector))) {
    cat("\nChangepoint declared at time =", status(detector), "\n\n")
  }
}
