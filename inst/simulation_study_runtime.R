####### Simulation study #######
## ...for the paper 'A general framework for fast online changepoint detection',
## Per August Jarval Moen, 2024.

## Imports
library(CHAD)

#install.packages("foreach")
#install.packages("doParallel")
#install.packages("abind")
library(foreach)
library(doSNOW)
library(abind)
library(ggplot2)
library(patchwork)

## Saving options
save = TRUE # if results should be saved
# Specify directory in which results should be saved:
maindir = "... fill in ... "
maindir = "/Users/peraugustmoen/Library/CloudStorage/OneDrive-UniversitetetiOslo/project3/simulations"
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)

load_threshes_dir = ""
load_results_dir = ""

# Creating subfolder with current time as name:
if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/computational_performance",dateandtime))
  dir.create(savedir, showWarnings =FALSE )
  plotdir = file.path(maindir, sprintf("%s/computational_performance/plots",dateandtime))
  dir.create(plotdir, showWarnings =FALSE )
  datadir = file.path(maindir, sprintf("%s/computational_performance/data",dateandtime))
  dir.create(datadir, showWarnings =FALSE )

}

## NOTE:
#  The working directory (as specified by setwd()) should be a parent
#  directory of the inst/ folder in which tuning_competing_methods.R
#  can be found, for instance in the source file directory of the
#  CHAD package

source("inst/tuning_competing_methods.R")

## NOTE:
#  The package 'ocd' must be installed. Uncomment below line if needed:
#install.packages('ocd')


## Global parameters
# N = 1500 # number of data samples considered
# chgptloc = round(N/3)
# num_sim = 1000 # number of iterations in the simulation
# ps = c(100,1000)
# sparsities100 = c(1,5,10,100)
# sparsities1000 = c(1,5,30,1000)
# thetas = seq(0.0, 4.0, by=0.2)
# num_methods = 5
# num_cores = 6
# MC_reps = 1000
# false_alarm_prob = 0.05
# estimate_mean = FALSE
# estimate_mean_until = 0

# testing

num_sim = 5 # number of iterations in the simulation
num_methods = 5
estimate_mean = FALSE
estimate_mean_until = round(chgptloc/2)
constant_penalty = TRUE

if(!estimate_mean){
  estimate_mean_until=0
}



##  Methods included in the simulation study are:
#   1. mean detector in the paper
#   2. OCD
#   3. Mei
#   4. XS
#   5. Chan



### Simulation for N
N = 500 # number of data samples considered
p = 100
binlength = N / 100

runtimes_n = matrix(NA, nrow = num_methods, ncol = N / binlength)

for (v in 1:num_sim) {
  ys = matrix(rnorm(N*p), nrow = p, ncol = N)
  detector = CHAD(p, method = "mean", leading_constant = rep(100000,2),
                  constant_penalty=constant_penalty,estimate_mean = estimate_mean)
  detector_ocd <- ocd::ChangepointDetector(dim=p, method='ocd',thresh=rep(100000,3))
  detector_mei <- ocd::ChangepointDetector(dim=p, method='Mei',thresh=rep(100000,2))
  detector_xs <- ocd::ChangepointDetector(dim=p, method='XS',thresh=100000)
  detector_chan <- ocd::ChangepointDetector(dim=p, method='Chan',thresh=100000)


  for (m in 1:num_methods) {
    for (j in 1:dim(runtimes_n)[2]) {
      ## start timing
      for (i in ((j-1)*binlength + 1):(j*binlength)){
        if (m==1){
          detectors[[v]] <-getData(detectors[[v]], ys[,i])
        }else{
          detectors[[v]] <-ocd::getData(detectors[[v]], ys[,i])
        }

      }
      ## end timing
    }

  }


}


results <- provideDimnames(results , sep = "_", base = list("p", "sparsity", "theta", "method", "iteration"))

results[1,3,4,,]-chgptloc

meanabove = function(v) mean(v[v>0])

s_ind = 4
p_ind = 1

plot(thetas,apply(results[p_ind,s_ind,,1,]-chgptloc, 1, meanabove),type="l"  )
lines(thetas,apply(results[p_ind,s_ind,,2,]-chgptloc, 1, meanabove),type="l" ,col=2 )
lines(thetas,apply(results[p_ind,s_ind,,3,]-chgptloc, 1, meanabove),type="l" ,col=3 )
lines(thetas,apply(results[p_ind,s_ind,,4,]-chgptloc, 1, meanabove),type="l" ,col=4 )
lines(thetas,apply(results[p_ind,s_ind,,5,]-chgptloc, 1, meanabove),type="l" ,col=5 )






## Plotting
p_ind = 1
lenn = length(apply(results[p_ind,s_ind,,1,]-chgptloc, 1, meanabove))
plots = list()
for (i in 1:length(sparsities100)) {
  s_ind = i
  plotdata <- data.frame(
    x = thetas,
    y = c(apply(results[p_ind,s_ind,,1,]-chgptloc, 1, meanabove),
          apply(results[p_ind,s_ind,,2,]-chgptloc, 1, meanabove),
          apply(results[p_ind,s_ind,,3,]-chgptloc, 1, meanabove),
          apply(results[p_ind,s_ind,,4,]-chgptloc, 1, meanabove),
          apply(results[p_ind,s_ind,,5,]-chgptloc, 1, meanabove)),
    Method = factor(c(rep("CHAD", lenn), rep("ocd", lenn), rep("Mei",lenn), rep("XS", lenn),
                      rep("Chan",lenn)))

  )
  if(p_ind == 1){
    ss = sparsities100[i]
  }else{ss = sparsities1000[i]}

  plot_base <- ggplot(data = plotdata, aes(x = x, y = y, color = Method, linetype = Method)) +
    geom_line() +              # Plot lines
    scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) + # Custom colors
    scale_linetype_manual(values = c("solid", "dashed", "longdash", "dotdash", "twodash")) + # Custom line types
    theme_bw() +               # Add theme_bw()
    theme(legend.position = "right") +
    scale_y_continuous(limits = c(0,N - chgptloc)) +
    ggtitle(bquote(s == .(ss))) +
    theme(plot.title = element_text(hjust = 0.5))+
    ylab("Detection delay") +
    xlab(bquote(phi)) +
    theme(legend.title = element_blank())


  if (i %in% c(2,4)) {
    plot_base <- plot_base + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  }
  if (i %in% c(1,2)) {
    plot_base <- plot_base + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  }

  plots[[i]] = plot_base
}
# Combine the plots using patchwork, and share the legend
combined_plot <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom") &
  plot_annotation(
    title = "Average detection delay",
    theme = theme(plot.title = element_text(hjust = 0.5))
  )
combined_plot












## Array of method outputs
## Dimensions: [p, sparsity, theta, method, iteration]
outputs = array(NA, dim = c(length(ps), length(sparsities100), length(thetas),
                            num_methods, num_sim))


#### Step 3: Simulation study investigating computational performance

# use object.size()

