####### Simulation study for run times #######
## ...for the paper 'A grid-based methodology for fast online changepoint detection',
## Per August Jarval Moen, 2026.

## The run time study considers only CHAD


## To run the simulation, install the CHAD package from GitHub:
# devtools::install_github("peraugustmoen/CHAD").


## Imports
library(CHAD)


library(abind)
library(ggplot2)
library(patchwork)

## Seed
set.seed(123)
## Saving options
save <- TRUE # if results should be saved

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
    sprintf("%s/computational_performance_ONLY_CHAD", dateandtime)
  )
  dir.create(savedir, showWarnings = FALSE)
  plotdir <- file.path(
    maindir,
    sprintf("%s/computational_performance_ONLY_CHAD/plots", dateandtime)
  )
  dir.create(plotdir, showWarnings = FALSE)
  datadir <- file.path(
    maindir,
    sprintf("%s/computational_performance_ONLY_CHAD/data", dateandtime)
  )
  dir.create(datadir, showWarnings = FALSE)
}



source(system.file("plot_style.R", package = "CHAD")) # shared colors/line types for all plots

## Global params
num_sim_n <- 20 # number of iterations for the simulation varying n
num_sim_p <- 20 # number of iterations for the simulation varying p
num_methods <- 1
estimate_mean <- FALSE
estimate_mean_until <- 0
constant_penalty <- TRUE

if (!estimate_mean) {
  estimate_mean_until <- 0
}

## Measurement grids, geometrically spaced:
num_t_points <- 15 # number of measurement points over t (max t is N below)
t_min <- 100 # smallest measurement point over t
num_p_points <- 15 # number of values of p
p_min <- 10 # smallest value of p
p_max <- 20 # largest value of p

## Common threshold to have no detections:
runtime_thresh <- 80000


## Labels for the method variants. colors and line types are defined in
## plot_style.R
method_labels <- c(
  "CHAD"
)

## The two plot sets are defined below.
## "main" shows only the pure-R variants.
## "all" additionally includes the C-accelerated
## CHAD (Rcpp) and MdFOCuS (Rcpp).
plot_subsets <- list(
  main = c(1),
  all = c(1)
)


### Simulation for N
{
  N <- 1000000 # max number of data samples considered
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

    detector_chad_r <- CHAD(p_const,
      method = "mean", leading_constant = rep(runtime_thresh, 2),
      constant_penalty = constant_penalty, estimate_mean = estimate_mean,
      pure_R = TRUE
    )

    detectors <- list()
    detectors[[1]] <- detector_chad_r

    for (m in c(1)) {
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
    ## slope-1 reference segment in the top-right corner:
    yv <- plotdata$y[is.finite(plotdata$y) & plotdata$y > 0]
    ytop <- max(yv)
    lx0 <- log10(min(xvec))
    lx1 <- log10(max(xvec))
    xb <- lx0 + 0.94 * (lx1 - lx0)
    dlx <- min(0.22 * (lx1 - lx0), 0.6)
    xa <- xb - dlx # left end
    ya <- log10(ytop) + 0.08 # lower end, just above the tallest curve
    yb <- ya + dlx # slope 1: same number of decades in x and y
    gg <- gg +
      scale_x_log10() +
      scale_y_log10(
        limits = c(NA, ytop * 5),
        labels = function(b) format(b, scientific = FALSE, drop0trailing = TRUE)
      ) +
      annotate("segment",
        x = 10^xa, xend = 10^xb, y = 10^ya, yend = 10^yb,
        color = "gray40", linewidth = 0.35
      ) +
      annotate("text",
        x = 10^((xa + xb) / 2), y = 10^((ya + yb) / 2),
        label = "slope 1", color = "gray40", size = 2.6, vjust = -0.8
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
