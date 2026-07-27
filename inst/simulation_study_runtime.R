####### Simulation study for run times #######
## ...for the paper 'A grid-based framework for fast online changepoint detection',
## Per August Jarval Moen, 2026.

## The run time study compares:
##   1. CHAD (Rcpp)     - the CHAD package with its C routine
##   2. ocd             - R package ocd (R)
##   3. Mei             - R package ocd (R)
##   4. XS              - R package ocd (R)
##   5. Chan            - R package ocd (R)
##   6. MdFOCuS (Rcpp)  - sequential Rcpp implementation in the R package
##                        'focus' (CRAN)
##   7. MdFOCuS (R)     - the sparsity-adaptive variant in R
##                        (GitHub repo "Focused"), as used in the statistical
##                        simulation studies. Timed per iteration of its main
##                        loop, see mdfocus_r_timed.R
##   8. CHAD (R)        - the CHAD method run in pure R


## To run the simulation, install the CHAD package from GitHub:
# devtools::install_github("peraugustmoen/CHAD").

# TO install focus and ocd from CRAN:
# install.packages("ocd")
# install.packages("focus")

## Imports
library(CHAD)
library(ocd)
library(focus)
library(foreach)
library(doSNOW)
library(abind)
library(ggplot2)
library(patchwork)

## Saving options
save <- FALSE # if results should be saved

## IMPORTANT! Specify the directory in which results should be saved:
## (in the maindir variable)
maindir <- ""
dateandtime <- gsub(" ", "--", as.character(Sys.time()))
dateandtime <- gsub(":", ".", dateandtime)
savedir <- file.path(maindir, dateandtime)


# Creating subfolder with current time as name:
if (save) {
  dir.create(savedir, showWarnings = FALSE)
  savedir <- file.path(
    maindir,
    sprintf("%s/computational_performance", dateandtime)
  )
  dir.create(savedir, showWarnings = FALSE)
  plotdir <- file.path(
    maindir,
    sprintf("%s/computational_performance/plots", dateandtime)
  )
  dir.create(plotdir, showWarnings = FALSE)
  datadir <- file.path(
    maindir,
    sprintf("%s/computational_performance/data", dateandtime)
  )
  dir.create(datadir, showWarnings = FALSE)
}


source(system.file("tuning_competing_methods.R", package = "CHAD"))
source(system.file("mdfocus_r_timed.R", package = "CHAD"))
source(system.file("plot_style.R", package = "CHAD")) # shared colors/line types for all plots

# Helper function for mdfocus:
mdfocus_cran_projections <- function(dim) {
  if (dim < 5) {
    return(NULL)
  }
  focus::generate_projection_indexes(dim, 2)
}


## Global params
num_sim_n <- 10 # number of iterations for the simulation varying n
num_sim_p <- 10 # number of iterations for the simulation varying p
num_methods <- 8
estimate_mean <- FALSE
estimate_mean_until <- 0
constant_penalty <- TRUE

if (!estimate_mean) {
  estimate_mean_until <- 0
}

## Measurement grids, geometrically spaced:
num_t_points <- 25 # number of measurement points over t (max t is N below)
t_min <- 100 # smallest measurement point over t
num_p_points <- 15 # number of values of p
p_min <- 8 # smallest value of p
p_max <- 200 # largest value of p

## Common threshold to have no detections:
runtime_thresh <- 80000


## Labels for the method variants. colors and line types are defined in
## plot_style.R
method_labels <- c(
  "CHAD (Rcpp)", "ocd", "Mei", "XS", "Chan",
  "MdFOCuS (Rcpp)", "MdFOCuS (R)", "CHAD (R)"
)

## The two plot sets are defined below.
## "main" shows only the pure-R variants.
## "all" additionally includes the C-accelerated
## CHAD (Rcpp) and MdFOCuS (Rcpp).
plot_subsets <- list(
  main = c(8, 2, 3, 4, 5, 7),
  all = c(1, 8, 2, 3, 4, 5, 6, 7)
)


### Simulation for N
{
  N <- 10000 # max number of data samples considered
  p_const <- 100


  t_points <- unique(round(exp(
    seq(log(t_min), log(N), length.out = num_t_points)
  )))
  seg_bounds <- c(0, t_points)
  num_segments <- length(t_points)


  runtimes_n_all <- array(NA_real_,
    dim = c(num_methods, num_segments, num_sim_n)
  )
  memory_n_all <- array(NA_real_,
    dim = c(num_methods, num_segments, num_sim_n)
  )

  for (v in 1:num_sim_n) {
    cat("v = ", v, "\n")
    ys <- matrix(rnorm(N * p_const), nrow = p_const, ncol = N)
    detector <- CHAD(p_const,
      method = "mean", leading_constant = rep(runtime_thresh, 2),
      constant_penalty = constant_penalty, estimate_mean = estimate_mean
    )
    detector_ocd <- ocd::ChangepointDetector(
      dim = p_const, method = "ocd", thresh = rep(runtime_thresh, 3)
    )
    detector_mei <- ocd::ChangepointDetector(
      dim = p_const, method = "Mei", thresh = rep(runtime_thresh, 2)
    )
    detector_xs <- ocd::ChangepointDetector(
      dim = p_const, method = "XS", thresh = runtime_thresh
    )
    detector_chan <- ocd::ChangepointDetector(
      dim = p_const, method = "Chan", thresh = runtime_thresh
    )
    detector_chad_r <- CHAD(p_const,
      method = "mean", leading_constant = rep(runtime_thresh, 2),
      constant_penalty = constant_penalty, estimate_mean = estimate_mean,
      pure_R = TRUE
    )

    detectors <- list()
    detectors[[1]] <- detector
    detectors[[2]] <- detector_ocd
    detectors[[3]] <- detector_mei
    detectors[[4]] <- detector_xs
    detectors[[5]] <- detector_chan
    detectors[[8]] <- detector_chad_r

    for (m in c(1:5, 8)) {
      for (j in 1:num_segments) {
        seg_inds <- (seg_bounds[j] + 1):seg_bounds[j + 1]
        ## start timing
        startt <- Sys.time()
        for (i in seg_inds) {
          if (m == 1 || m == 8) {
            detectors[[m]] <- CHAD::getData(detectors[[m]], ys[, i])
          } else {
            detectors[[m]] <- ocd::getData(detectors[[m]], ys[, i])
          }
        }
        endd <- Sys.time()
        runtimes_n_all[m, j, v] <-
          as.numeric(difftime(endd, startt, units = "secs")) /
            length(seg_inds)
        memory_n_all[m, j, v] <- as.numeric(object.size(detectors[[m]]))
      }
    }

    # run Rcpp version of mdfocus:
    det_mdfocus <- detector_create(
      type = "multivariate",
      dim_indexes = mdfocus_cran_projections(p_const)
    )
    theta0_mdfocus <- rep(0, p_const)
    for (j in 1:num_segments) {
      seg_inds <- (seg_bounds[j] + 1):seg_bounds[j + 1]
      startt <- Sys.time()
      for (i in seg_inds) {
        detector_update(det_mdfocus, ys[, i])
        r <- get_statistics(det_mdfocus,
          family = "gaussian",
          theta0 = theta0_mdfocus
        )
        if (!is.null(r$stat) && r$stat > runtime_thresh) {
          cat("mdfocus declared a change (should not happen)\n")
        }
      }
      endd <- Sys.time()
      runtimes_n_all[6, j, v] <-
        as.numeric(difftime(endd, startt, units = "secs")) /
          length(seg_inds)

      memory_n_all[6, j, v] <-
        detector_pieces_len(det_mdfocus) * (p_const + 1) * 8
    }

    # run R version of mdfocus:
    sparsity_levels <- 2^seq_len(floor(log2(p_const)))
    dat <- data.frame(t(ys))
    res <- FocusCH_HighDim_timed(dat,
      get_opt_cost = \(...) get_partial_opt(...,
        cost = cost_lr_partial0, which_par = sparsity_levels
      ),
      threshold = rep(Inf, 2 + length(sparsity_levels)),
      seg_bounds = seg_bounds
    )
    for (j in 1:num_segments) {
      seg_inds <- (seg_bounds[j] + 1):seg_bounds[j + 1]
      ## block timing, divided by the block length, exactly as for the other
      ## methods above:
      runtimes_n_all[7, j, v] <-
        res$segment_times[j] / length(seg_inds)
      memory_n_all[7, j, v] <-
        res$nb_candidates[j] * (p_const + 1) * 8
    }
  } # end for (v in 1:num_sim_n)

  ## Compute median:
  runtimes_n <- apply(runtimes_n_all, c(1, 2), median, na.rm = TRUE)
  memory_n <- apply(memory_n_all, c(1, 2), median, na.rm = TRUE)

  # convert to milliseconds:
  runtimes_n <- runtimes_n * 1000

  # convert to kB:
  memory_n <- memory_n / 1000

  num_obs_vector <- t_points

  if (save) {
    saveRDS(runtimes_n, file = sprintf("%s/runtimes_n.RDA", datadir))
    saveRDS(memory_n, file = sprintf("%s/memory_n.RDA", datadir))
    saveRDS(num_obs_vector, file = sprintf("%s/num_obs_vector.RDA", datadir))
  }
}

matplot(t(runtimes_n), type = "l", col = 1:num_methods, lty = 1:num_methods)
legend("topleft", legend = method_labels, col = 1:num_methods, lty = 1:num_methods)

matplot(t(memory_n), type = "l", col = 1:num_methods, lty = 1:num_methods)
legend("topleft", legend = method_labels, col = 1:num_methods, lty = 1:num_methods)

## Simulation for p
{
  n_const <- 500 # number of data samples considered

  ## Geometrically spaced values of p:
  p_values <- unique(round(exp(
    seq(log(p_min), log(p_max), length.out = num_p_points)
  )))

  runtimes_p_all <- array(NA_real_,
    dim = c(num_methods, length(p_values), num_sim_p)
  )
  memory_p_all <- array(NA_real_,
    dim = c(num_methods, length(p_values), num_sim_p)
  )

  for (v in 1:num_sim_p) {
    cat("v = ", v, "\n")
    for (j in seq_along(p_values)) {
      p <- p_values[j]
      cat("p = ", p, "\n")
      ys <- matrix(rnorm(n_const * p), nrow = p, ncol = n_const)
      detector <- CHAD(p,
        method = "mean", leading_constant = rep(runtime_thresh, 2),
        constant_penalty = constant_penalty, estimate_mean = estimate_mean
      )
      detector_ocd <- ocd::ChangepointDetector(
        dim = p, method = "ocd", thresh = rep(runtime_thresh, 3)
      )
      detector_mei <- ocd::ChangepointDetector(
        dim = p, method = "Mei", thresh = rep(runtime_thresh, 2)
      )
      detector_xs <- ocd::ChangepointDetector(
        dim = p, method = "XS", thresh = runtime_thresh
      )
      detector_chan <- ocd::ChangepointDetector(
        dim = p, method = "Chan", thresh = runtime_thresh
      )
      detector_chad_r <- CHAD(p,
        method = "mean", leading_constant = rep(runtime_thresh, 2),
        constant_penalty = constant_penalty, estimate_mean = estimate_mean,
        pure_R = TRUE
      )

      detectors <- list()
      detectors[[1]] <- detector
      detectors[[2]] <- detector_ocd
      detectors[[3]] <- detector_mei
      detectors[[4]] <- detector_xs
      detectors[[5]] <- detector_chan
      detectors[[8]] <- detector_chad_r


      ## process data point and take average over the time

      for (m in c(1:5, 8)) {
        startt <- Sys.time()
        for (i in 1:n_const) {
          if (m == 1 || m == 8) {
            detectors[[m]] <- CHAD::getData(detectors[[m]], ys[, i])
          } else {
            detectors[[m]] <- ocd::getData(detectors[[m]], ys[, i])
          }
        }
        endd <- Sys.time()
        runtimes_p_all[m, j, v] <-
          as.numeric(difftime(endd, startt, units = "secs")) / n_const
        memory_p_all[m, j, v] <- as.numeric(object.size(detectors[[m]]))
      }

      ## run Rccp version of mdfocus:
      det_mdfocus <- detector_create(
        type = "multivariate",
        dim_indexes = mdfocus_cran_projections(p)
      )
      theta0_mdfocus <- rep(0, p)
      startt <- Sys.time()
      for (i in 1:n_const) {
        detector_update(det_mdfocus, ys[, i])
        r <- get_statistics(det_mdfocus,
          family = "gaussian",
          theta0 = theta0_mdfocus
        )
        if (!is.null(r$stat) && r$stat > runtime_thresh) {
          cat("mdfocus declared a change (should not happen)\n")
        }
      }
      endd <- Sys.time()
      runtimes_p_all[6, j, v] <-
        as.numeric(difftime(endd, startt, units = "secs")) / n_const
      memory_p_all[6, j, v] <-
        detector_pieces_len(det_mdfocus) * (p + 1) * 8

      ## The R version of mdfocus:
      sparsity_levels <- 2^seq_len(floor(log2(p)))
      dat <- data.frame(t(ys))
      startt <- Sys.time()
      res <- FocusCH_HighDim(dat,
        get_opt_cost = \(...) get_partial_opt(...,
          cost = cost_lr_partial0, which_par = sparsity_levels
        ),
        threshold = rep(Inf, 2 + length(sparsity_levels))
      )
      endd <- Sys.time()
      runtimes_p_all[7, j, v] <-
        as.numeric(difftime(endd, startt, units = "secs")) / n_const
      memory_p_all[7, j, v] <-
        res$nb * (p + 1) * 8
    }
  }


  ## compute median:
  runtimes_p <- apply(runtimes_p_all, c(1, 2), median, na.rm = TRUE)
  memory_p <- apply(memory_p_all, c(1, 2), median, na.rm = TRUE)

  # convert to milliseconds:
  runtimes_p <- runtimes_p * 1000

  # convert to kB:
  memory_p <- memory_p / 1000

  p_vector <- p_values

  if (save) {
    saveRDS(runtimes_p, file = sprintf("%s/runtimes_p.RDA", datadir))
    saveRDS(memory_p, file = sprintf("%s/memory_p.RDA", datadir))
    saveRDS(p_vector, file = sprintf("%s/p_vector.RDA", datadir))
  }
}


## Plotting

## Plot helper:
make_runtime_panel <- function(mat, xvec, subset, ylab_text, xlab_expr,
                               loglog = TRUE,
                               legend_labels = method_labels[subset]) {
  nbin <- length(xvec)
  plotdata <- data.frame(
    x = rep(xvec, times = length(subset)),
    y = c(t(mat[subset, , drop = FALSE])),
    Method = factor(rep(method_labels[subset], each = nbin),
      levels = legend_labels
    )
  )
  gg <- ggplot(
    data = plotdata,
    aes(x = x, y = y, color = Method, linetype = Method)
  ) +
    geom_line(show.legend = TRUE) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_linetype_manual(values = method_linetypes, drop = FALSE) +
    theme_bw() +
    theme(legend.position = "right") +
    ylab(ylab_text) +
    xlab(xlab_expr) +
    theme(legend.title = element_blank())
  if (loglog) {
    ## Add a line with slope 1:
    ytop <- max(plotdata$y, na.rm = TRUE)
    x1 <- max(xvec)
    x0 <- x1 / 3
    y0 <- ytop * 1.2
    y1 <- y0 * (x1 / x0) # slope 1 on the log-log scale
    gg <- gg +
      scale_x_log10() +
      scale_y_log10(
        limits = c(NA, ytop * 5),
        labels = function(b) format(b, scientific = FALSE, drop0trailing = TRUE)
      ) +
      annotate("segment",
        x = x0, xend = x1, y = y0, yend = y1,
        color = "gray40", linewidth = 0.35
      ) +
      annotate("text",
        x = sqrt(x0 * x1) / 1.35, y = sqrt(y0 * y1) * 1.35,
        label = "slope 1", color = "gray40", size = 2.6
      )
  }
  return(gg)
}


for (sname in names(plot_subsets)) {
  sub <- plot_subsets[[sname]]
  ## show all methods in the memory panels. MdFOCuS memory (both variants) is
  ## the analytic active-candidate storage nb*(p+1)*8 bytes; object.size does
  ## not work for it..
  sub_mem <- sub

  plot1 <- make_runtime_panel(
    runtimes_n, num_obs_vector, sub,
    "Update time (ms)", bquote(t)
  )
  plot2 <- make_runtime_panel(
    memory_n, num_obs_vector, sub_mem,
    "Memory use (Kb)", bquote(t),
    legend_labels = method_labels[sub]
  )
  plot3 <- make_runtime_panel(
    runtimes_p, p_vector, sub,
    "Update time (ms)", bquote(p)
  ) +
    theme(axis.title.y = element_blank())
  plot4 <- make_runtime_panel(
    memory_p, p_vector, sub_mem,
    "Memory use (Kb)", bquote(p),
    legend_labels = method_labels[sub]
  ) +
    theme(axis.title.y = element_blank())

  combined_plot <- (plot1 + plot3) / (plot2 + plot4) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  print(combined_plot)

  if (save) {
    ggsave(
      filename = sprintf("%s/plot1_loglog_%s.pdf", plotdir, sname),
      plot = plot1, device = "pdf", width = 7, height = 5
    )
    ggsave(
      filename = sprintf("%s/plot2_loglog_%s.pdf", plotdir, sname),
      plot = plot2, device = "pdf", width = 7, height = 5
    )
    ggsave(
      filename = sprintf("%s/plot3_loglog_%s.pdf", plotdir, sname),
      plot = plot3, device = "pdf", width = 7, height = 5
    )
    ggsave(
      filename = sprintf("%s/plot4_loglog_%s.pdf", plotdir, sname),
      plot = plot4, device = "pdf", width = 7, height = 5
    )
    ggsave(
      filename = sprintf("%s/combined_runtime_loglog_%s.eps", plotdir, sname),
      plot = combined_plot, device = "eps", width = 7, height = 6
    )
    ggsave(
      filename = sprintf("%s/combined_runtime_loglog_%s.pdf", plotdir, sname),
      plot = combined_plot, device = "pdf", width = 7, height = 6
    )
  }
}

## Absolute plots on linear axes for all eight variants (for reference).
sub_all <- plot_subsets$all
sub_all_mem <- sub_all
plot1 <- make_runtime_panel(runtimes_n, num_obs_vector, sub_all,
  "Update time (ms)", bquote(t),
  loglog = FALSE
)
plot2 <- make_runtime_panel(memory_n, num_obs_vector, sub_all_mem,
  "Memory use (Kb)", bquote(t),
  loglog = FALSE, legend_labels = method_labels[sub_all]
)
plot3 <- make_runtime_panel(runtimes_p, p_vector, sub_all,
  "Update time (ms)", bquote(p),
  loglog = FALSE
) +
  theme(axis.title.y = element_blank())
plot4 <- make_runtime_panel(memory_p, p_vector, sub_all_mem,
  "Memory use (Kb)", bquote(p),
  loglog = FALSE, legend_labels = method_labels[sub_all]
) +
  theme(axis.title.y = element_blank())

combined_plot <- (plot1 + plot3) / (plot2 + plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot

if (save) {
  ggsave(
    filename = sprintf("%s/plot1.pdf", plotdir),
    plot = plot1, device = "pdf", width = 7, height = 5
  )
  ggsave(
    filename = sprintf("%s/plot2.pdf", plotdir),
    plot = plot2, device = "pdf", width = 7, height = 5
  )
  ggsave(
    filename = sprintf("%s/plot3.pdf", plotdir),
    plot = plot3, device = "pdf", width = 7, height = 5
  )
  ggsave(
    filename = sprintf("%s/plot4.pdf", plotdir),
    plot = plot4, device = "pdf", width = 7, height = 5
  )
  ggsave(
    filename = sprintf("%s/combined_runtime.eps", plotdir),
    plot = combined_plot, device = "eps", width = 7, height = 6
  )
  ggsave(
    filename = sprintf("%s/combined_runtime.pdf", plotdir),
    plot = combined_plot, device = "pdf", width = 7, height = 6
  )
}
