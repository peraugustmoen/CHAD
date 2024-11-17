####### Simulation study for run times #######
## ...for the paper 'A general framework for fast online changepoint detection',
## Per August Jarval Moen, 2024.

## Install the CHAD package from GitHub:
# devtools::install_github("peraugustmoen/CHAD")

## Imports
library(CHAD)
library(ocd)
library(foreach)
library(doSNOW)
library(abind)
library(ggplot2)
library(patchwork)

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


## NOTE:
#  The package 'ocd' must be installed. Uncomment below line if needed:
# install.packages('ocd')





##  Methods included in the simulation study are:
#   1. mean detector in the paper
#   2. OCD
#   3. Mei
#   4. XS
#   5. Chan

## Global params
num_sim_n <- 200 # number of iterations for the simulation varying n
num_sim_p <- 5 # number of iterations for the simulation varying n
num_methods <- 5
estimate_mean <- FALSE
estimate_mean_until <- 0
constant_penalty <- TRUE

num_bins <- 25
if (!estimate_mean) {
  estimate_mean_until <- 0
}


### Simulation for N
{
  N <- 5000 # max number of data samples considered
  p_const <- 8
  binlength_N <- N / num_bins

  runtimes_n <- matrix(0, nrow = num_methods, ncol = N / binlength_N)
  memory_n <- matrix(0, nrow = num_methods, ncol = N / binlength_N)

  for (v in 1:num_sim_n) {
    cat("v = ", v, "\n")
    ys <- matrix(rnorm(N * p_const), nrow = p_const, ncol = N)
    detector <- CHAD(p_const,
      method = "mean", leading_constant = rep(80000, 2),
      constant_penalty = constant_penalty, estimate_mean = estimate_mean
    )
    detector_ocd <- ocd::ChangepointDetector(
      dim = p_const, method = "ocd", thresh = rep(80000, 3)
    )
    detector_mei <- ocd::ChangepointDetector(
      dim = p_const, method = "Mei", thresh = rep(80000, 2)
    )
    detector_xs <- ocd::ChangepointDetector(
      dim = p_const, method = "XS", thresh = 80000
    )
    detector_chan <- ocd::ChangepointDetector(
      dim = p_const, method = "Chan", thresh = 80000
    )

    detectors <- list()
    detectors[[1]] <- detector
    detectors[[2]] <- detector_ocd
    detectors[[3]] <- detector_mei
    detectors[[4]] <- detector_xs
    detectors[[5]] <- detector_chan

    for (m in 1:num_methods) {
      # cat("m = ", m, "\n")
      for (j in 1:dim(runtimes_n)[2]) {
        ## start timing
        startt <- proc.time()
        for (i in ((j - 1) * binlength_N + 1):(j * binlength_N)) {
          if (m == 1) {
            detectors[[m]] <- CHAD::getData(detectors[[m]], ys[, i])
          } else {
            detectors[[m]] <- ocd::getData(detectors[[m]], ys[, i])
          }
        }
        endd <- proc.time()
        runtimes_n[m, j] <- runtimes_n[m, j] +
          (endd - startt)[3] / num_sim_n / binlength_N
        memory_n[m, j] <- memory_n[m, j] + object.size(detectors[[m]]) /
          num_sim_n
      }
    }
  }

  # convert to milliseconds:
  runtimes_n <- runtimes_n * 800

  # convert to kB:
  memory_n <- memory_n / 824

  num_obs_vector <- 1:(dim(runtimes_n)[2]) * binlength_N

  if (save) {
    saveRDS(runtimes_n, file = sprintf("%s/runtimes_n.RDA", datadir))
    saveRDS(memory_n, file = sprintf("%s/memory_n.RDA", datadir))
    saveRDS(num_obs_vector, file = sprintf("%s/num_obs_vector.RDA", datadir))
  }
}

## Simulation for p
{
  n_const <- 500 # number of data samples considered
  P <- 200 # max value of p
  binlength_P <- P / num_bins

  runtimes_p <- matrix(0, nrow = num_methods, ncol = P / binlength_P)
  memory_p <- matrix(0, nrow = num_methods, ncol = P / binlength_P)

  for (v in 1:num_sim_p) {
    cat("v = ", v, "\n")
    for (j in 1:dim(runtimes_p)[2]) {
      p <- binlength_P * j
      cat("p = ", p, "\n")
      ys <- matrix(rnorm(2 * n_const * p), nrow = p, ncol = 2 * n_const)
      detector <- CHAD(p,
        method = "mean", leading_constant = rep(80000, 2),
        constant_penalty = constant_penalty, estimate_mean = estimate_mean
      )
      detector_ocd <- ocd::ChangepointDetector(
        dim = p, method = "ocd", thresh = rep(80000, 3)
      )
      detector_mei <- ocd::ChangepointDetector(
        dim = p, method = "Mei", thresh = rep(80000, 2)
      )
      detector_xs <- ocd::ChangepointDetector(
        dim = p, method = "XS", thresh = 80000
      )
      detector_chan <- ocd::ChangepointDetector(
        dim = p, method = "Chan", thresh = 80000
      )

      detectors <- list()
      detectors[[1]] <- detector
      detectors[[2]] <- detector_ocd
      detectors[[3]] <- detector_mei
      detectors[[4]] <- detector_xs
      detectors[[5]] <- detector_chan

      stopn <- n_const

      ## process the first n_const - 1 data points
      for (m in 1:num_methods) {
        for (i in 1:stopn) {
          if (m == 1) {
            detectors[[m]] <- CHAD::getData(detectors[[m]], ys[, i])
          } else {
            detectors[[m]] <- ocd::getData(detectors[[m]], ys[, i])
          }
        }
      }

      ## process data point and take average over the time

      for (m in 1:num_methods) {
        startt <- proc.time()
        for (i in (stopn + 1):(2 * n_const)) {
          if (m == 1) {
            detectors[[m]] <- CHAD::getData(detectors[[m]], ys[, i])
          } else {
            detectors[[m]] <- ocd::getData(detectors[[m]], ys[, i])
          }
        }
        endd <- proc.time()
        runtimes_p[m, j] <-
          runtimes_p[m, j] + (endd - startt)[3] / num_sim_p / n_const
        memory_p[m, j] <-
          memory_p[m, j] + object.size(detectors[[m]]) / num_sim_p
      }
    }
  }

  # convert to milliseconds:
  runtimes_p <- runtimes_p * 800

  # convert to Mb:
  memory_p <- memory_p / 824

  p_vector <- 1:(dim(runtimes_p)[2]) * binlength_P

  if (save) {
    saveRDS(runtimes_p, file = sprintf("%s/runtimes_p.RDA", datadir))
    saveRDS(memory_p, file = sprintf("%s/memory_p.RDA", datadir))
    saveRDS(p_vector, file = sprintf("%s/p_vector.RDA", datadir))
  }
}




## Plotting

# first, run time dependence on sample size:
plotdata1 <- data.frame(
  x = num_obs_vector,
  y = c(t(runtimes_n)),
  Method = factor(c(
    rep("CHAD", N / binlength_N),
    rep("ocd", N / binlength_N),
    rep("Mei", N / binlength_N),
    rep("XS", N / binlength_N),
    rep("Chan", N / binlength_N)
  ))
)


plot1 <- ggplot(
  data = plotdata1,
  aes(x = x, y = y, color = Method, linetype = Method)
) +
  geom_line() + # Plot lines
  scale_color_manual(values = c(
    "red", "blue",
    "green", "purple", "orange"
  )) + # Custom colors
  scale_linetype_manual(values = c(
    "solid",
    "dashed", "longdash", "dotdash", "twodash"
  )) + # Custom line types
  theme_bw() + # Add theme_bw()
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.32)) +
  # ggtitle(bquote(p == .(p_const))) +
  # theme(plot.title = element_text(hjust = 0.5))+
  ylab("Update time (ms)") +
  theme(legend.title = element_blank()) +
  xlab(bquote(t)) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

plot1

plotdata2 <- data.frame(
  x = num_obs_vector,
  y = c(t(memory_n)),
  Method = factor(c(
    rep("CHAD", N / binlength_N),
    rep("ocd", N / binlength_N),
    rep("Mei", N / binlength_N),
    rep("XS", N / binlength_N),
    rep("Chan", N / binlength_N)
  ))
)


plot2 <- ggplot(
  data = plotdata2,
  aes(x = x, y = y, color = Method, linetype = Method)
) +
  geom_line() + # Plot lines
  scale_color_manual(values = c(
    "red", "blue",
    "green", "purple", "orange"
  )) + # Custom colors
  scale_linetype_manual(values = c(
    "solid", "dashed",
    "longdash", "dotdash", "twodash"
  )) + # Custom line types
  theme_bw() + # Add theme_bw()
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 41)) +
  # theme(plot.title = element_text(hjust = 0.5))+
  # ggtitle(bquote(t == .(n_const))) +
  ylab("Memory use (Kb)") +
  xlab(bquote(t)) +
  theme(legend.title = element_blank())
plot2



# then, run time dependence on p:
plotdata3 <- data.frame(
  x = p_vector,
  y = c(t(runtimes_p)),
  Method = factor(c(
    rep("CHAD", P / binlength_P),
    rep("ocd", P / binlength_P),
    rep("Mei", P / binlength_P),
    rep("XS", P / binlength_P),
    rep("Chan", P / binlength_P)
  ))
)


plot3 <- ggplot(
  data = plotdata3,
  aes(x = x, y = y, color = Method, linetype = Method)
) +
  geom_line() + # Plot lines
  scale_color_manual(values = c(
    "red", "blue",
    "green", "purple", "orange"
  )) + # Custom colors
  scale_linetype_manual(values = c(
    "solid", "dashed",
    "longdash", "dotdash", "twodash"
  )) + # Custom line types
  theme_bw() + # Add theme_bw()
  theme(legend.position = "right") +
  # ggtitle(bquote(t == .(n_const))) +
  scale_y_continuous(limits = c(0, 3)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Update time (ms)") +
  xlab(bquote(p)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank())

plot3

plotdata4 <- data.frame(
  x = p_vector,
  y = c(t(memory_p)),
  Method = factor(c(
    rep("CHAD", P / binlength_P),
    rep("ocd", P / binlength_P),
    rep("Mei", P / binlength_P),
    rep("XS", P / binlength_P),
    rep("Chan", P / binlength_P)
  ))
)


plot4 <- ggplot(
  data = plotdata4,
  aes(x = x, y = y, color = Method, linetype = Method)
) +
  geom_line() + # Plot lines
  scale_color_manual(values = c(
    "red", "blue",
    "green", "purple", "orange"
  )) + # Custom colors
  scale_linetype_manual(values = c(
    "solid", "dashed",
    "longdash", "dotdash", "twodash"
  )) + # Custom line types
  theme_bw() + # Add theme_bw()
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 600)) +
  ylab("Memory use (Kb)") +
  # scale_y_continuous(limits = c(0,200)) +
  xlab(bquote(p)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank())
plot4



## combine plots
combined_plot <- (plot1 + plot3) / (plot2 + plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot

if (save) {
  ggsave(
    filename = sprintf("%s/plot1.eps", plotdir),
    plot = plot1,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/plot1.pdf", plotdir),
    plot = plot1,
    device = "pdf",
    width = 8,
    height = 8
  )

  ggsave(
    filename = sprintf("%s/plot2.eps", plotdir),
    plot = plot2,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/plot2.pdf", plotdir),
    plot = plot2,
    device = "pdf",
    width = 8,
    height = 8
  )

  ggsave(
    filename = sprintf("%s/plot3.eps", plotdir),
    plot = plot3,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/plot3.pdf", plotdir),
    plot = plot3,
    device = "pdf",
    width = 8,
    height = 8
  )

  ggsave(
    filename = sprintf("%s/plot4.eps", plotdir),
    plot = plot4,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/plot4.pdf", plotdir),
    plot = plot4,
    device = "pdf",
    width = 8,
    height = 8
  )

  ggsave(
    filename = sprintf("%s/combined_runtime.eps", plotdir),
    plot = combined_plot,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/combined_runtime.pdf", plotdir),
    plot = combined_plot,
    device = "pdf",
    width = 8,
    height = 8
  )
}



## Full plots, for supplementary material, where the range of the y axis
## is larger

plotdata1 <- data.frame(
  x = num_obs_vector,
  y = c(t(runtimes_n)),
  Method = factor(c(
    rep("CHAD", N / binlength_N),
    rep("ocd", N / binlength_N),
    rep("Mei", N / binlength_N),
    rep("XS", N / binlength_N),
    rep("Chan", N / binlength_N)
  ))
)


plot1 <- ggplot(
  data = plotdata1,
  aes(x = x, y = y, color = Method, linetype = Method)
) +
  geom_line() + # Plot lines
  scale_color_manual(values = c(
    "red", "blue",
    "green", "purple", "orange"
  )) + # Custom colors
  scale_linetype_manual(values = c(
    "solid", "dashed",
    "longdash", "dotdash", "twodash"
  )) + # Custom line types
  theme_bw() + # Add theme_bw()
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.32)) +
  # ggtitle(bquote(p == .(p_const))) +
  # theme(plot.title = element_text(hjust = 0.5))+
  ylab("Update time (ms)") +
  theme(legend.title = element_blank()) +
  xlab(bquote(t))

plot1

plotdata2 <- data.frame(
  x = num_obs_vector,
  y = c(t(memory_n)),
  Method = factor(c(
    rep("CHAD", N / binlength_N),
    rep("ocd", N / binlength_N),
    rep("Mei", N / binlength_N),
    rep("XS", N / binlength_N),
    rep("Chan", N / binlength_N)
  ))
)


plot2 <- ggplot(
  data = plotdata2,
  aes(x = x, y = y, color = Method, linetype = Method)
) +
  geom_line() + # Plot lines
  scale_color_manual(values = c(
    "red", "blue",
    "green", "purple", "orange"
  )) + # Custom colors
  scale_linetype_manual(values = c(
    "solid", "dashed",
    "longdash", "dotdash", "twodash"
  )) + # Custom line types
  theme_bw() + # Add theme_bw()
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 41)) +
  # theme(plot.title = element_text(hjust = 0.5))+
  # ggtitle(bquote(t == .(n_const))) +
  ylab("Memory use (Kb)") +
  xlab(bquote(t)) +
  theme(legend.title = element_blank())
plot2



# then, run time dependence on p:
plotdata3 <- data.frame(
  x = p_vector,
  y = c(t(runtimes_p)),
  Method = factor(c(
    rep("CHAD", P / binlength_P),
    rep("ocd", P / binlength_P),
    rep("Mei", P / binlength_P),
    rep("XS", P / binlength_P),
    rep("Chan", P / binlength_P)
  ))
)


plot3 <- ggplot(
  data = plotdata3,
  aes(x = x, y = y, color = Method, linetype = Method)
) +
  geom_line() + # Plot lines
  scale_color_manual(values = c(
    "red", "blue",
    "green", "purple", "orange"
  )) + # Custom colors
  scale_linetype_manual(values = c(
    "solid", "dashed",
    "longdash", "dotdash", "twodash"
  )) + # Custom line types
  theme_bw() + # Add theme_bw()
  theme(legend.position = "right") +
  # ggtitle(bquote(t == .(n_const))) +
  scale_y_continuous(limits = c(0, max(plotdata3$y))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Update time (ms)") +
  xlab(bquote(p)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank())

plot3

plotdata4 <- data.frame(
  x = p_vector,
  y = c(t(memory_p)),
  Method = factor(c(
    rep("CHAD", P / binlength_P),
    rep("ocd", P / binlength_P),
    rep("Mei", P / binlength_P),
    rep("XS", P / binlength_P),
    rep("Chan", P / binlength_P)
  ))
)


plot4 <- ggplot(
  data = plotdata4,
  aes(x = x, y = y, color = Method, linetype = Method)
) +
  geom_line() + # Plot lines
  scale_color_manual(values = c(
    "red", "blue",
    "green", "purple", "orange"
  )) + # Custom colors
  scale_linetype_manual(values = c(
    "solid", "dashed",
    "longdash", "dotdash", "twodash"
  )) + # Custom line types
  theme_bw() + # Add theme_bw()
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, max(plotdata4$y))) +
  ylab("Memory use (Kb)") +
  # scale_y_continuous(limits = c(0,200)) +
  xlab(bquote(p)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank())
plot4



## combine plots
combined_plot <- (plot1 + plot3) / (plot2 + plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
# plot_annotation(
#  title = "Computational costs",
#  theme = theme(plot.title = element_text(hjust = 0.5))
# )
combined_plot

if (save) {
  ggsave(
    filename = sprintf("%s/plot1_extended.eps", plotdir),
    plot = plot1,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/plot1_extended.pdf", plotdir),
    plot = plot1,
    device = "pdf",
    width = 8,
    height = 8
  )

  ggsave(
    filename = sprintf("%s/plot2_extended.eps", plotdir),
    plot = plot2,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/plot2_extended.pdf", plotdir),
    plot = plot2,
    device = "pdf",
    width = 8,
    height = 8
  )

  ggsave(
    filename = sprintf("%s/plot3_extended.eps", plotdir),
    plot = plot3,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/plot3_extended.pdf", plotdir),
    plot = plot3,
    device = "pdf",
    width = 8,
    height = 8
  )

  ggsave(
    filename = sprintf("%s/plot4_extended.eps", plotdir),
    plot = plot4,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/plot4_extended.pdf", plotdir),
    plot = plot4,
    device = "pdf",
    width = 8,
    height = 8
  )

  ggsave(
    filename = sprintf("%s/combined_runtime_extended.eps", plotdir),
    plot = combined_plot,
    device = "eps",
    width = 8,
    height = 8
  )
  ggsave(
    filename = sprintf("%s/combined_runtime_extended.pdf", plotdir),
    plot = combined_plot,
    device = "pdf",
    width = 8,
    height = 8
  )
}
