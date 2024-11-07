
# install.packages("ocd)

## Ensure WD is the same as the source package location
source("inst/tuning_competing_methods.R")


## get thresholds
N = 500
p = 50
MC_reps = 500
false_alarm_prob = 0.05
ocd_mc_res = MC_ocd_FA(p, false_alarm_prob, MC_reps=MC_reps,
                       N=N,seed = 123)
mei_mc_res = MC_Mei_FA(p, false_alarm_prob=false_alarm_prob, N=N, MC_reps = MC_reps,
                      seed = 123)
xs_mc_res <- MC_XS_FA(p, false_alarm_prob = false_alarm_prob, N=N, MC_reps = MC_reps,
                  seed=123)

chan_mc_res <- MC_Chan_FA(p, false_alarm_prob = false_alarm_prob, N = N,
                      MC_reps = MC_reps, seed = 123)

mean_mc_res = MC_mean(p, false_alarm_prob, constant_penalty=TRUE, MC_reps=MC_reps, N=N,seed = 123)



K = 200
times = array(data = NA, dim = c(5, N))
dd = array(data = NA, dim = c(5,K))
chgptloc = N / 2
sparsity = 5
theta = 1
for (j in 1:K) {
  if(j %% 10==0){
    print(j)
  }
  detector = CHAD(p, method = "mean", leading_constant = mean_mc_res,
                  constant_penalty=TRUE)
  detector_ocd <- ocd::ChangepointDetector(dim=p, method='ocd',thresh=ocd_mc_res)
  detector_mei <- ocd::ChangepointDetector(dim=p, method='Mei',thresh=mei_mc_res)
  detector_xs <- ocd::ChangepointDetector(dim=p, method='XS',thresh=xs_mc_res)
  detector_chan <- ocd::ChangepointDetector(dim=p, method='Chan',thresh=chan_mc_res)


  ys = matrix(rnorm(N*p), nrow = p, ncol = N)
  ys[1:sparsity,(chgptloc+1):N] = ys[1:sparsity,(chgptloc+1):N] + theta/sqrt(sparsity)

  for (i in 1:N) {
    startt = proc.time()
    detector <- getData(detector, ys[,i])
    endd = proc.time()
    times[1,i] = times[1,i] + (endd - startt)[3]/K

    startt = proc.time()
    detector_ocd <- ocd::getData(detector_ocd, ys[,i])
    endd = proc.time()
    times[2,i] = times[2,i] +  (endd - startt)[3]/K

    startt = proc.time()
    detector_mei <- ocd::getData(detector_mei, ys[,i])
    endd = proc.time()
    times[3,i] = times[3,i] + (endd - startt)[3]/K

    startt = proc.time()
    detector_xs <- ocd::getData(detector_xs, ys[,i])
    endd = proc.time()
    times[4,i] = times[4,i] + (endd - startt)[3]/K

    startt = proc.time()
    detector_chan <- ocd::getData(detector_chan, ys[,i])
    endd = proc.time()
    times[5,i] = times[5,i] + (endd - startt)[3]/K


  }

  detectors = list()
  detectors[[1]] = detector
  detectors[[2]] = detector_ocd
  detectors[[3]] = detector_mei
  detectors[[4]] = detector_xs
  detectors[[5]] = detector_chan

  for (i in 1:length(detectors)) {
    det = detectors[[i]]
    if(identical(status(det), 'monitoring')){
      dd[i, j] = N
    }else{
      dd[i,j] = (status(det)-chgptloc)
    }
  }
}

dd = na.omit(dd)
if(identical(theta,0)){
  falsealarms = rowSums(dd<N) / K
  falsealarms
}else{

  avg_detection_delays = rep(NA, dim(dd)[1])
  for (i in 1:dim(dd)[1]) {
    ddd = dd[i,]
    avg_detection_delays[i] = mean(ddd[ddd>0])
  }
  print(avg_detection_delays)
}




dd_ocd = na.omit(dd_ocd)
sum(dd<N) / K
sum(dd_ocd<N)/K
mean(abs(dd[dd>0]))
mean(abs(dd_ocd[dd_ocd>0]))

hist(dd[dd>0])
hist(abs(dd_ocd[dd_ocd>0]))

plot(times) #helt rått
points(times_ocd,col=2)

plot(1:n, sqrt(2*log(exp(1)*p*log((1:n + 20))/4^2)))
