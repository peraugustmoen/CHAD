# This is just FocusCH_HighDim (the R version of MdFOCuS, from
# MdFocus_MeanGaussian_md.R in the github repo github.com/guillemr/Focused)
# with timing added. The main loop is already sequential. We time in blocks
# in the same way as in simulation_study_runtime.R

FocusCH_HighDim_timed <- function(data,
                                  get_opt_cost = get_glo_opt,
                                  common_difference_step = 1,
                                  common_ratio_step = 1.5,
                                  first_step_qhull = ncol(data) + 5,
                                  threshold = Inf,
                                  dim_indexes = map2(seq(0, ncol(data) - 1, by = 2), seq(1, ncol(data), by = 2), \(f, s) c(f, s %% ncol(data)) + 1),
                                  seg_bounds = c(0, nrow(data))) {
  # Initialization-------------
  ## get data parameters
  n <- nrow(data)
  p <- ncol(data)

  ## get cnsts for costs
  data_cumsum_square <- cumsum(rowSums(data^2))

  ## get points P(i), i in {1,n}
  data_left_cumsum <- apply(cbind(1, data), 2, cumsum)

  ## pruning step
  next_step_qhull <- first_step_qhull

  ## Output candidate list
  list_cand <- list(
    index = integer(n), # set of candidates
    nb = 0, # number of candidates
    nb_at_step = integer(n), # vector of the number of candidates at each iteration
    opt.cost = vector(mode = "list", length = n), # vector of the functional cost at each iteration
    opt.change = vector(mode = "list", length = n)
  ) # vector of the optimal changepoint candidate at each iteration

  ## First iteration
  list_cand$nb <- list_cand$nb + 1
  list_cand$index[list_cand$nb] <- 1
  list_cand$nb_at_step[1] <- 1

  ## per-SEGMENT timing (the instrumentation): the loop iterations in each
  ## measurement segment are timed together with ONE Sys.time() pair, and we
  ## also record the number of stored candidates at the end of each segment.
  num_segments <- length(seg_bounds) - 1
  segment_times <- rep(NA_real_, num_segments) # elapsed seconds for the block
  nb_candidates <- rep(NA_integer_, num_segments) # candidates at end of segment

  #------------------------------------
  # For loop (over measurement segments)-
  broke <- FALSE
  for (b in seq_len(num_segments)) {
    lo <- max(seg_bounds[b] + 1, 2) # the loop starts at 2
    hi <- seg_bounds[b + 1]
    startt <- as.numeric(Sys.time()) # start of the timed block for this segment
    if (lo <= hi) {
      for (i in lo:hi) {
        ## add a candidate
        list_cand$nb <- list_cand$nb + 1
        list_cand$index[list_cand$nb] <- i
        list_cand$nb_at_step[i] <- list_cand$nb

        ## First step: search of optimal changepoint candidate (index and maximization)
        index_cand <- list_cand$index[1:list_cand$nb]
        left_mean <- data_left_cumsum[index_cand, ] # get cumsum for cost at time i

        out <- get_opt_cost(left_mean, data_cumsum_square[i]) # get optimal cost
        list_cand$opt.change[[i]] <- out$opt.change
        list_cand$opt.cost[[i]] <- out$opt.cost


        ## Second step: Pruning by the Quickhull algorithm from the package "geometry"
        if (list_cand$nb >= next_step_qhull | i == n) {
          on_the_hull <-
            map(dim_indexes, \(i) convhulln(left_mean[, c(1, i + 1)]) |>
              as.vector() |>
              unique()) |>
            unlist() |>
            unique() |>
            sort()

          list_cand$nb <- length(on_the_hull) # get number of vertices
          list_cand$index[1:list_cand$nb] <- index_cand[on_the_hull] # update the set of candidates
          ## update next_step_qhull
          next_step_qhull <- common_ratio_step * list_cand$nb + common_difference_step
        }

        ## if the (minus) optimal cost is over the threshold end the algorithm
        if (length(list_cand$opt.cost[[i]]) == length(threshold)) {
          change_detected <- sum(-list_cand$opt.cost[[i]] >= threshold)
        } else {
          warning("The number of statistics to test does not match the length of the threshold")
        }

        if (change_detected) {
          broke <- TRUE
          break
        }
      }
    }
    endd <- as.numeric(Sys.time()) # end of the timed block for this segment
    segment_times[b] <- endd - startt
    nb_candidates[b] <- list_cand$nb
    if (broke) break
  }
  return(list(
    list_cand = list_cand,
    segment_times = segment_times,
    nb_candidates = nb_candidates
  ))
}
