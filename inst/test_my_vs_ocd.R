MC_ocd_FA <- function(dim, false_alarm_prob, N, beta=1, sparsity='auto', MC_reps){
  peak_stat <- matrix(0, MC_reps, 3)
  colnames(peak_stat) <- c('diag','off_d','off_s')
  if (sparsity == 'sparse') peak_stat <- peak_stat[,-2]
  if (sparsity == 'dense') peak_stat <- peak_stat[,-3]

  # run MC_reps simulations for peak statistics of S_diag, S_{off,d} and S_{off,s}
  for (rep in 1:MC_reps){
    if(rep %%100==0){
      cat("Iteration: ", rep, "\n")
    }
    A <- matrix(0, dim, 1)
    tail <- matrix(0, dim, floor(log2(dim))*2+4)

    for (i in 1:N){
      x_new <- rnorm(dim)
      ret <- ocd::ocd_update(x_new, A, tail, beta, sparsity)
      A <- ret$A; tail <- ret$tail
      peak_stat[rep,] <- pmax(peak_stat[rep,], ret$stat)
    }
  }
  qq = 1-false_alarm_prob/ dim(peak_stat)[2]
  # compute the MC thresholds from the peak statistics
  thresh_est <- function(v) quantile(sort(v), qq)
  th_individual <- apply(peak_stat, 2, thresh_est)

  #th_multiplier <- thresh_est(apply(t(peak_stat)/th_individual, 2, max))
  #th <- th_individual * th_multiplier
  th <- th_individual
  names(th) <- colnames(peak_stat)
  return(th)
}

MC_Mei_FA <- function(dim, false_alarm_prob, N, b=beta/sqrt(dim), MC_reps){
  peak_stat <- matrix(0, MC_reps, 2)
  colnames(peak_stat) <- c('max','sum')

  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps){
    if(rep %%100==0){
      cat("Iteration: ", rep, "\n")
    }
    R <- matrix(0, dim, 2)

    for (i in 1:N){
      x_new <- rnorm(dim)
      ret <- ocd::Mei_update(x_new, R, b)
      R <- ret$R
      peak_stat[rep,] <- pmax(peak_stat[rep,], ret$stat)
    }
  }

  qq = 1-false_alarm_prob/ dim(peak_stat)[2]

  # compute the MC thresholds from the peak statistics
  thresh_est <- function(v) quantile(sort(v), qq)
  th_individual <- apply(peak_stat, 2, thresh_est)

  th <- th_individual
  names(th) <- colnames(peak_stat)
  return(th)

}

MC_XS <- function(dim, false_alarm_prob, N, p0=1/sqrt(dim), w=200, MC_reps){
  peak_stat <- rep(-Inf, MC_reps)

  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps){
    if(rep %%100==0){
      cat("Iteration: ", rep, "\n")
    }
    X_recent <- CUSUM <- matrix(0, dim, w)

    for (i in 1:N){
      x_new <- rnorm(dim)
      ret <- ocd::XS_update(x_new, X_recent, CUSUM, p0, w)
      X_recent <- ret$X_recent; CUSUM <- ret$CUSUM
      peak_stat[rep] <- pmax(peak_stat[rep], ret$stat)
    }
  }

  qq = 1-false_alarm_prob
  # compute the MC thresholds from the peak statistics
  th <- quantile(sort(peak_stat), qq)
  return(th)
}

MC_Chan <- function(dim, false_alarm_prob, N,  p0=1/sqrt(dim), w=200,
                    lambda=sqrt(8)-2, MC_reps){
  peak_stat <- rep(-Inf, MC_reps)

  # run MC_reps simulations for peak statistics
  for (rep in 1:MC_reps){
    if(rep %%100==0){
      cat("Iteration: ", rep, "\n")
    }
    X_recent <- CUSUM <- matrix(0, dim, w)

    for (i in 1:N){
      x_new <- rnorm(dim)
      ret <- ocd::Chan_update(x_new, X_recent, CUSUM, p0, w, lambda)
      X_recent <- ret$X_recent; CUSUM <- ret$CUSUM
      peak_stat[rep] <- pmax(peak_stat[rep], ret$stat)
    }
  }

  qq = 1-false_alarm_prob
  # compute the MC thresholds from the peak statistics
  th <- quantile(sort(peak_stat), qq)
  return(th)
}

# OCD penalty
p = 100
N = 100
ocd_mc_res = MC_ocd_FA(p, 0.05, beta=1, sparsity='auto', MC_reps=100, N = N)

# mean method penalty manual
MC_reps = 100
maxes = matrix(NA, nrow = 2, ncol = MC_reps)
argmax_n_obs = rep(NA, MC_reps)
argmax_s = rep(NA, MC_reps)
argmax_g = rep(NA, MC_reps)
set.seed(123)
for (j in 1:MC_reps) {
  if(j%% 10 == 0){
    print(j)
  }
  detector = CHAD(p, method = "mean", leading_constant = c(100,100),constant_penalty=TRUE)
  ys = matrix(rnorm(N*p), nrow = p, ncol = N)
  for (i in 1:N) {
    detector <- getData(detector, ys[,i])
  }
  mm = attr(detector, 'overall_max_statistics')
  maxes[,j] = mm
}

leading_const_manual = c(quantile(maxes[1,],0.95), quantile(maxes[2,],0.95))

leading_const = MC_mean(p, 0.05, constant_penalty=FALSE, MC_reps=N, N=N,seed = 123)

## check where the upper bound on the constant factor is tight
MC_reps = 1000
N = 1000
p=1000
maxes = matrix(NA, nrow = 2, ncol = MC_reps)
argmax_n_obs_sparse = rep(NA, MC_reps)
argmax_s_sparse = rep(NA, MC_reps)
argmax_g_sparse = rep(NA, MC_reps)
argmax_n_obs_dense = rep(NA, MC_reps)
argmax_g_dense = rep(NA, MC_reps)
set.seed(123)
for (j in 1:MC_reps) {
  if(j%% 10 == 0){
    print(j)
  }
  detector = CHAD(p, method = "mean", leading_constant = c(100,100),constant_penalty=TRUE)
  ys = matrix(rnorm(N*p), nrow = p, ncol = N)
  for (i in 1:N) {
    detector <- getData(detector, ys[,i])
  }
  mm = attr(detector, 'overall_max_statistics')
  maxes[,j] = mm
  argmax_n_obs_sparse[j] = attr(detector, 'allstats')$overall_argmax_n_obs_sparse
  argmax_s_sparse[j] = attr(detector, 'allstats')$overall_argmax_s_sparse
  argmax_g_sparse[j] = attr(detector, 'allstats')$overall_argmax_g_sparse

  argmax_n_obs_dense[j] = attr(detector, 'allstats')$overall_argmax_n_obs_dense
  argmax_g_dense[j] = attr(detector, 'allstats')$overall_argmax_g_dense
}

hist(argmax_n_obs_sparse)
hist(argmax_n_obs_dense)
hist(argmax_g_sparse)
hist(argmax_g_dense)
hist(argmax_s_sparse)

# argmax is attained at


## check false alarm rate approximately OK and detection delay
N = 1000
p = 5
N = 1000
MC_reps = 1000
false_alarm_prob = 0.05
ocd_mc_res = MC_ocd_FA(p, false_alarm_prob, beta=1, sparsity='auto', MC_reps=MC_reps,
                       N=N)
leading_const = MC_mean(p, false_alarm_prob, constant_penalty=TRUE, MC_reps=MC_reps, N=N,seed = 123)



K = 100
times = rep(NA, N)
times_ocd = rep(NA, N)
dd = rep(NA, K)
dd_ocd = rep(NA, K)
chgptloc = N / 2
sparsity = 1
theta = 4
for (j in 1:K) {
  if(j %% 10==0){
    print(j)
  }
  detector = CHAD(p, method = "mean", leading_constant = leading_const,
                  constant_penalty=TRUE)
  detector_ocd <- ocd::ChangepointDetector(dim=p, method='ocd',thresh=ocd_mc_res,
                                           sparsity= 'auto', beta=1)
  ys = matrix(rnorm(N*p), nrow = p, ncol = N)
  ys[1:sparsity,(chgptloc+1):N] = ys[1:sparsity,(chgptloc+1):N] + theta/sqrt(sparsity)
  #ys[1:3,20:n] = ys[1:3,20:n] + 0.5
  for (i in 1:N) {
    startt = proc.time()
    detector <- getData(detector, ys[,i])
    endd = proc.time()
    times[i] = (endd - startt)[3]

    startt = proc.time()
    detector_ocd <- ocd::getData(detector_ocd, ys[,i])
    endd = proc.time()
    times_ocd[i] = (endd - startt)[3]
  }
  if(identical(status(detector), 'monitoring')){
    dd[j] = N
  }else{
    dd[j] = (status(detector)-chgptloc)
  }

  if(identical(status(detector_ocd), 'monitoring')){
    dd_ocd[j] = N
  }else{
    dd_ocd[j] = (status(detector_ocd)-chgptloc)
  }
}

dd = na.omit(dd)
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
