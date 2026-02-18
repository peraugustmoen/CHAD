####### Real data example #######
## for the paper 'A general methodology for fast online changepoint detection',
## Per August Jarval Moen, 2026.

## Install the CHAD package from GitHub:
# devtools::install_github("peraugustmoen/CHAD")


## Saving options
save <- TRUE # if results should be saved

## IMPORTANT! Specify the directory in which results should be saved:
## (in the maindir variable)
maindir <- "/Users/peraugustmoen/Library/CloudStorage/OneDrive-UniversitetetiOslo/project3/jrss-b_revision/final_plots_and_code_for_paper/dataexample"
dateandtime <- gsub(" ", "--", as.character(Sys.time()))
dateandtime <- gsub(":", ".", dateandtime)
savedir <- file.path(maindir, dateandtime)
if (save) {
  dir.create(savedir, showWarnings = FALSE)
}


## Imports
library(CHAD)
library(plotly)
library(tidyr)
library(dplyr)


## Load data (located in the inst/data folder in the package)
dat <- load_exchange_rates()

# Including the ten most traded currencies in 2022 except for the US Dollar.
# These are:
# EUR - Euro
# JPY - Japanese Yen
# GBP- British Pound
# CNY - Chinese Yuan
# AUD - Australian Dollar
# CAD - Canadian Dollar
# CHF - Swiss Franc
# HKD - Hong Kong Dollar
# SGD - Singapore Dollar
# SEK - Swedish Krona

keep_inds <- 1:11
dat <- dat[, keep_inds]
colnames(dat) <- c(
  "date", dat[3,2:11]
)
dat <- dat[-(1:5), ] # removing info columns
dat[dat == "ND"] <- NA
dat <- na.omit(dat) # remove rows with missing values

dat$date <- as.POSIXct(dat$date, format = "%Y-%m-%d", tz = "UTC")

# transform all columns except the first to numeric:
dat[, 2:ncol(dat)] <- sapply(dat[, 2:ncol(dat)], as.numeric)

# divide each row of dat by the first row:
dat[, 2:ncol(dat)] <-
  sweep(dat[, 2:ncol(dat)], 2, as.numeric(dat[1, 2:ncol(dat)]), `/`)

# apply diff to each column except the date column:
dat_diff <- dat[, ]
dat_diff[1:(nrow(dat) - 1), 2:ncol(dat)] <- apply(dat[, 2:ncol(dat)], 2, diff)
dat_diff[(nrow(dat_diff)), 2:ncol(dat_diff)] <- 0


## Plotting before changepoint analysis

# plot the raw data
matplot(dat$date, dat[, 2:ncol(dat)],
  type = "l", xlab = "Time",
  ylab = "Exchange rate", main = "Exchange rates"
)


# plot the raw data using plotly
# Convert to long format
df_long <- dat %>%
  pivot_longer(cols = -date, names_to = "series", values_to = "value")
plot <- plot_ly(df_long,
  x = ~date, y = ~value,
  color = ~series, type = "scatter", mode = "lines"
) %>%
  layout(
    xaxis = list(title = "Date"),
    yaxis = list(title = "Value")
  )
plot

## Online changepoint detection using the online covariance changepoint
## detector

p <- ncol(dat_diff) - 1
N <- nrow(dat_diff)
N_obs = N

# estimate baseline noise level from first year of data:
trainingdata <- dat_diff[
  dat_diff$date < as.POSIXct("2001-01-01", format = "%Y-%m-%d", tz = "UTC"),
]
baseline_op <- norm(cov(trainingdata[, 2:ncol(trainingdata)]), "2")


## Gaussian distribution

# Choose critical value via MC simulation with isotropic Gaussian noise,
# with 5% false alarm probability tolerance

critical_value <- MC_covariance(
  p = p, false_alarm_prob = 0.05,
  constant_penalty = FALSE, MC_reps = 1000, N = 1000,
  estimate_mean = FALSE,
  baseline_operatornorm = 1,
  seed = 102102
)
#critical_value = 6.039568


detector <- CHAD(p,
  method = "covariance", leading_constant = critical_value,
  estimate_mean = FALSE, baseline_operatornorm = baseline_op
)

teststats <- rep(NA, N)
changepoints <- c()
argmax <- c()
counter <- 0
for (i in 1:N) {
  counter <- counter + 1
  y_new <- as.numeric(dat_diff[i, 2:ncol(dat_diff)])
  detector <- CHAD::getData(detector, y_new)
  if (counter > 1) {
    teststats[i] <- attr(detector, "statistics")
  }
  if (!identical(status(detector), "monitoring")) {
    changepoints <- c(changepoints, i)
    argmax <- c(argmax, attr(detector, "allstats")$argmax)
    detector <- CHAD::reset(detector)
    counter <- 0
  }
}

dates <- dat_diff$date

cat("Changpoint detected at:\n")
dates[changepoints] +1 # plus one since we took a diff
cat("Candidate changepoint location where the change was flagged:\n")
dates[changepoints - argmax] +1 # plus one since we took a diff

# Interactive plot of exchange rates with changepoint declarations:
plot <- plot_ly(df_long,
  x = ~date, y = ~value, color = ~series,
  type = "scatter", mode = "lines"
) %>%
  layout(
    xaxis = list(title = "Date"),
    yaxis = list(title = "Value")
  )



shapes <- list()
for (i in 1:length(changepoints)) {
  cp <- dates[changepoints[i]]+1 #plus one since we took a diff
  shapes <- append(shapes, list(
    list(
      type = "line",
      x0 = cp, x1 = cp,
      y0 = 0, y1 = 1, # Using relative coordinates for y-axis
      xref = "x", # Reference x and y axes
      yref = "paper", # yref = 'paper' uses the whole plot height
      line = list(color = "red", dash = "dash")
    )
  ))
}
plot <- plot %>% layout(shapes = shapes)
plot

# Interactive plot with differentiated rates + changepoints
df_diff_long <- dat_diff %>%
  pivot_longer(cols = -date, names_to = "series", values_to = "value")
# Plotly plot with vertical lines at changepoints
plot <- plot_ly(df_diff_long,
                x = ~date, y = ~value, color = ~series,
                type = "scatter", mode = "lines"
) %>%
  layout(
    xaxis = list(title = "Date"),
    yaxis = list(title = "Value")
  )



shapes <- list()
for (i in 1:length(changepoints)) {
  cp <- dates[changepoints[i]]
  # print(cp)
  shapes <- append(shapes, list(
    list(
      type = "line",
      x0 = cp, x1 = cp,
      y0 = 0, y1 = 1, # Using relative coordinates for y-axis
      xref = "x", # Reference x and y axes
      yref = "paper", # yref = 'paper' uses the whole plot height
      line = list(color = "red", dash = "dash")
    )
  ))
}
plot <- plot %>% layout(shapes = shapes)
plot





# ggplot for nice pdf and eps output
line_types <- c(
  "solid", "dashed", "dotted", "dotdash", "longdash",
  "twodash", "solid", "dashed", "dotted", "dotdash"
)

plot <- ggplot(df_long, aes(
  x = date, y = value, color = series,
  linetype = series
)) +
  geom_line() +
  labs(x = "Date", y = "Exchange rate to USD") +
  theme_bw()
plot <- plot + theme(
  legend.position = "bottom", # Position the legend at the bottom
  legend.title = element_blank() # Remove the legend title
)
plot

cp_dates <- dates[changepoints]

for (i in 1:length(cp_dates)) {
  plot <- plot + geom_vline(
    xintercept = as.numeric(cp_dates[i]+1), # plus one since we took a diff
    color = "red", linetype = "dashed"
  )
}

print(plot)

if (save) {
  ggsave(
    filename = sprintf("%s/currencies.eps", savedir),
    plot = plot,
    device = "eps",
    width = 8,
    height = 6
  )
  ggsave(
    filename = sprintf("%s/currencies.pdf", savedir),
    plot = plot,
    device = "pdf",
    width = 8,
    height = 6
  )
}








## Gaussian distribution with long training period

# Choose critical value via MC simulation with isotropic Gaussian noise,
# with 5% false alarm probability tolerance, with sequences of length
# N = 10,000

critical_value_long <- MC_covariance(
  p = p, false_alarm_prob = 0.05,
  constant_penalty = FALSE, MC_reps = 1000, N = 10000,
  estimate_mean = FALSE,
  baseline_operatornorm = 1,
  seed = 102102
)
# critical_value_long = 6.146377

detector <- CHAD(p,
                 method = "covariance", leading_constant = critical_value_long,
                 estimate_mean = FALSE, baseline_operatornorm = baseline_op
)

teststats <- rep(NA, N)
changepoints <- c()
argmax <- c()
counter <- 0
for (i in 1:N) {
  counter <- counter + 1
  y_new <- as.numeric(dat_diff[i, 2:ncol(dat_diff)])
  detector <- CHAD::getData(detector, y_new)
  if (counter > 1) {
    teststats[i] <- attr(detector, "statistics")
  }
  if (!identical(status(detector), "monitoring")) {
    changepoints <- c(changepoints, i)
    argmax <- c(argmax, attr(detector, "allstats")$argmax)
    detector <- CHAD::reset(detector)
    counter <- 0
  }
}

dates <- dat_diff$date

cat("Changpoint detected at:\n")
dates[changepoints] +1 # plus one since we took a diff
cat("Candidate changepoint location where the change was flagged:\n")
dates[changepoints - argmax] +1 # plus one since we took a diff

# Plotly plot with vertical lines at changepoints
plot <- plot_ly(df_long,
                x = ~date, y = ~value, color = ~series,
                type = "scatter", mode = "lines"
) %>%
  layout(
    xaxis = list(title = "Date"),
    yaxis = list(title = "Value")
  )

shapes <- list()
for (i in 1:length(changepoints)) {
  cp <- dates[changepoints[i]]+1 #plus one since we took a diff
  # print(cp)
  shapes <- append(shapes, list(
    list(
      type = "line",
      x0 = cp, x1 = cp,
      y0 = 0, y1 = 1, # Using relative coordinates for y-axis
      xref = "x", # Reference x and y axes
      yref = "paper", # yref = 'paper' uses the whole plot height
      line = list(color = "red", dash = "dash")
    )
  ))
}
plot <- plot %>% layout(shapes = shapes)
plot


# ggplot for nice pdf and eps output
line_types <- c(
  "solid", "dashed", "dotted", "dotdash", "longdash",
  "twodash", "solid", "dashed", "dotted", "dotdash"
)

plot <- ggplot(df_long, aes(
  x = date, y = value, color = series,
  linetype = series
)) +
  geom_line() +
  labs(x = "Date", y = "Exchange rate to USD") +
  theme_bw()
plot <- plot + theme(
  legend.position = "bottom", # Position the legend at the bottom
  legend.title = element_blank() # Remove the legend title
)
plot

cp_dates <- dates[changepoints]

for (i in 1:length(cp_dates)) {
  plot <- plot + geom_vline(
    xintercept = as.numeric(cp_dates[i]+1), # plus one since we took a diff
    color = "red", linetype = "dashed"
  )
}

print(plot)

if (save) {
  ggsave(
    filename = sprintf("%s/currencies_long_training.eps", savedir),
    plot = plot,
    device = "eps",
    width = 8,
    height = 6
  )
  ggsave(
    filename = sprintf("%s/currencies_long_training.pdf", savedir),
    plot = plot,
    device = "pdf",
    width = 8,
    height = 6
  )
}










## calibration using t_5 distribution under the null

df = 5

critical_value_heavytail <- MC_covariance_heavytail(
  p = p, false_alarm_prob = 0.05,
  constant_penalty = FALSE, MC_reps = 1000, N = 1000,
  estimate_mean = FALSE,
  baseline_operatornorm = 1,
  seed = 102102,
  df = df
)
# critical_value_heavytail = 19.09366

detector <- CHAD(p,
                 method = "covariance", leading_constant = critical_value_heavytail,
                 estimate_mean = FALSE, baseline_operatornorm = baseline_op
)

teststats <- rep(NA, N)
changepoints <- c()
argmax <- c()
counter <- 0
for (i in 1:N) {
  counter <- counter + 1
  y_new <- as.numeric(dat_diff[i, 2:ncol(dat_diff)])
  detector <- CHAD::getData(detector, y_new)
  if (counter > 1) {
    teststats[i] <- attr(detector, "statistics")
  }
  if (!identical(status(detector), "monitoring")) {
    changepoints <- c(changepoints, i)
    argmax <- c(argmax, attr(detector, "allstats")$argmax)
    detector <- CHAD::reset(detector)
    counter <- 0
  }
}

dates <- dat_diff$date

cat("Changpoint detected at:\n")
dates[changepoints] +1 # plus one since we took a diff
cat("Candidate changepoint location where the change was flagged:\n")
dates[changepoints - argmax] +1 # plus one since we took a diff

# Plotly plot with vertical lines at changepoints
plot <- plot_ly(df_long,
                x = ~date, y = ~value, color = ~series,
                type = "scatter", mode = "lines"
) %>%
  layout(
    xaxis = list(title = "Date"),
    yaxis = list(title = "Value")
  )

shapes <- list()
for (i in 1:length(changepoints)) {
  cp <- dates[changepoints[i]]+1 #plus one since we took a diff
  # print(cp)
  shapes <- append(shapes, list(
    list(
      type = "line",
      x0 = cp, x1 = cp,
      y0 = 0, y1 = 1, # Using relative coordinates for y-axis
      xref = "x", # Reference x and y axes
      yref = "paper", # yref = 'paper' uses the whole plot height
      line = list(color = "red", dash = "dash")
    )
  ))
}
plot <- plot %>% layout(shapes = shapes)
plot


# ggplot for nice pdf and eps output
line_types <- c(
  "solid", "dashed", "dotted", "dotdash", "longdash",
  "twodash", "solid", "dashed", "dotted", "dotdash"
)

plot <- ggplot(df_long, aes(
  x = date, y = value, color = series,
  linetype = series
)) +
  geom_line() +
  labs(x = "Date", y = "Exchange rate to USD") +
  theme_bw()
plot <- plot + theme(
  legend.position = "bottom", # Position the legend at the bottom
  legend.title = element_blank() # Remove the legend title
)
plot

cp_dates <- dates[changepoints]

for (i in 1:length(cp_dates)) {
  plot <- plot + geom_vline(
    xintercept = as.numeric(cp_dates[i]+1), # plus one since we took a diff
    color = "red", linetype = "dashed"
  )
}

print(plot)

if (save) {
  ggsave(
    filename = sprintf("%s/currencies_heavytail.eps", savedir),
    plot = plot,
    device = "eps",
    width = 8,
    height = 6
  )
  ggsave(
    filename = sprintf("%s/currencies_heavytail.pdf", savedir),
    plot = plot,
    device = "pdf",
    width = 8,
    height = 6
  )
}
