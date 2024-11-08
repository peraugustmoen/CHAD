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
  savedir = file.path(maindir, sprintf("%s/statistical_performance",dateandtime))
  dir.create(savedir, showWarnings =FALSE )
  plotdir = file.path(maindir, sprintf("%s/statistical_performance/plots",dateandtime))
  dir.create(plotdir, showWarnings =FALSE )
  datadir = file.path(maindir, sprintf("%s/statistical_performance/data",dateandtime))
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
N = 30 # number of data samples considered
chgptloc = round(N/3)
num_sim = 3*6 # number of iterations in the simulation
#ps = c(100,1000)
ps = c(2000)
sparsities100 = c(1,5,10,100)
sparsities1000 = c(1,5,30,1000)
thetas = seq(0.0, 8.0, by=0.4)
num_methods = 5
num_cores = 6
MC_reps = 1000
false_alarm_prob = 0.05
estimate_mean = FALSE
estimate_mean_until = round(chgptloc/2)
constant_penalty = TRUE

if(!estimate_mean){
  estimate_mean_until=0
}


## print parameters to file
{
paramfile = sprintf("%s/parameters.txt", savedir)
cat("Parameters:\n", file = paramfile, append = TRUE)
cat("N = ", N, " \n", file = paramfile, append=TRUE)
cat("chgptloc = ", chgptloc, "\n",file = paramfile, append=TRUE)
cat("ps = ", ps, "\n",file = paramfile, append=TRUE)
cat("sparsities100 = ", sparsities100, "\n",file = paramfile, append=TRUE)
cat("sparsities1000 = ", sparsities1000, "\n",file = paramfile, append=TRUE)
cat("thetas = ", thetas, "\n",file = paramfile, append=TRUE)
cat("num_methods = ", num_methods, "\n",file = paramfile, append=TRUE)
cat("num_cores = ", num_cores, "\n",file = paramfile, append=TRUE)
cat("MC_reps = ", MC_reps, "\n",file = paramfile, append=TRUE)
cat("false_alarm_prob = ", false_alarm_prob, "\n",file = paramfile, append=TRUE)
cat("estimate_mean = ", estimate_mean, "\n",file = paramfile, append=TRUE)
cat("constant_penalty = ", constant_penalty, "\n",file = paramfile, append=TRUE)
}

##  Methods included in the simulation study are:
#   1. mean detector in the paper
#   2. OCD
#   3. Mei
#   4. XS
#   5. Chan


#### Step 1: Choosing thresholds for the methods via MC simulations ####
thresholds = NULL
if(!identical(load_threshes_dir,"")){
  thresholds = readRDS( file=sprintf("%s/thresholds.RDA", load_threshes_dir))
}else{
  ## List of MC thresholds for the methods:
  thresholds = list(mean = list(), ocd = list(), mei = list(), xs = list(), chan = list())
  ## .. so e.g. thresholds[[2]][[1]] is a vector of the MC simulated thresholds for
  ## the OCD method for the first value in the vector ps

  for (v in 1:length(ps)) {
    p = ps[v]

    thresholds[[1]][[v]] = MC_mean(p, false_alarm_prob, constant_penalty=constant_penalty,
                                   estimate_mean = estimate_mean,
                                   MC_reps=MC_reps, N=N,seed = 123)

    thresholds[[2]][[v]] = MC_ocd_FA(dim = p, false_alarm_prob = false_alarm_prob,
                                     MC_reps=MC_reps, N = N, est_length = estimate_mean_until, seed = 123)
    thresholds[[3]][[v]] = MC_Mei_FA(p, false_alarm_prob=false_alarm_prob, N=N, MC_reps = MC_reps,
                           seed = 123)
    thresholds[[4]][[v]] <- MC_XS_FA(p, false_alarm_prob = false_alarm_prob, N=N, MC_reps = MC_reps,
                          seed=123)

    thresholds[[5]][[v]] <- MC_Chan_FA(p, false_alarm_prob = false_alarm_prob, N = N,
                              MC_reps = MC_reps, seed = 123)
  }
  if(save){
    saveRDS(thresholds, file=sprintf("%s/thresholds.RDA", datadir))
  }
}


#### Step 2: Simulation study investigating statistical performance

## Performed in parallel for faster execution
if(!identical(load_results_dir, "")){
  results = readRDS( file=sprintf("%s/results.RDA", load_results_dir))
}else{
  cl <- makeCluster(num_cores, type = "SOCK")
  registerDoSNOW(cl)

  # Set up progress bar
  pb <- txtProgressBar(max = num_sim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  results <- foreach(z = 1:num_sim,
                   .combine = function(...) abind(..., along=5),
                    .multicombine = TRUE, .options.snow = opts) %dopar% {
    #for(i in 1:num_sim){

    library(CHAD)
    set.seed(z+1000)
    result_array <- array(NA, dim = c(length(ps), length(sparsities100), length(thetas),
                                      num_methods))
    {sink("/dev/null") ; # avoiding prints to terminal
    for (v in 1:length(ps)) {
      p = ps[v]
      for (j in 1:length(sparsities100)) {
        s = 1 #sparsity
        if(p <100){
          s = j
          if(j == length(sparsities100)){
            s = p
          }
        }else if(v==1){
          s = sparsities100[j]
        }else{
          s = sparsities1000[j]
        }
        for (t in 1:length(thetas)) {

          theta = thetas[t]
          ys = matrix(rnorm(N*p), nrow = p, ncol = N)
          ys[1:s,(chgptloc+1):N] = ys[1:s,(chgptloc+1):N] + theta/sqrt(s)

          detector = CHAD(p, method = "mean", leading_constant = thresholds[[1]][[v]],
                          constant_penalty=constant_penalty,estimate_mean = estimate_mean)
          detector_ocd <- ocd::ChangepointDetector(dim=p, method='ocd',thresh=thresholds[[2]][[v]])
          detector_mei <- ocd::ChangepointDetector(dim=p, method='Mei',thresh=thresholds[[3]][[v]])
          detector_xs <- ocd::ChangepointDetector(dim=p, method='XS',thresh=thresholds[[4]][[v]])
          detector_chan <- ocd::ChangepointDetector(dim=p, method='Chan',thresh=thresholds[[5]][[v]])

          if(estimate_mean){
            mean_est = rowSums(ys[,1:estimate_mean_until])/estimate_mean_until
          }
          for (i in 1:N) {

            detector <- getData(detector, ys[,i])


            if(i>estimate_mean_until){
              detector_ocd <- ocd::getData(detector_ocd, ys[,i])

              detector_mei <- ocd::getData(detector_mei, ys[,i])

              detector_xs <- ocd::getData(detector_xs, ys[,i])

              detector_chan <- ocd::getData(detector_chan, ys[,i])
            }

          }

          detectors = list()
          detectors[[1]] = detector
          detectors[[2]] = detector_ocd
          detectors[[3]] = detector_mei
          detectors[[4]] = detector_xs
          detectors[[5]] = detector_chan


          if(identical(status(detector), 'monitoring')){
            result_array[v,j,t,1] = N
          }else{
            result_array[v,j,t,1] = (status(detector))
          }

          for (i in 2:length(detectors)) {
            det = detectors[[i]]
            if(identical(status(det), 'monitoring')){
              result_array[v,j,t,i] = N
            }else{
              result_array[v,j,t,i] = (status(det)) + estimate_mean_until
            }
          }

        }

      }

    }
    ; sink();}

    result_array
  }
  close(pb)
  stopCluster(cl)
  if(save){
    saveRDS(results, file=sprintf("%s/results.RDA", datadir))
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

