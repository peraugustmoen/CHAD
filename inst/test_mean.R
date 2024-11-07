
p = 100

detector = CHAD(p, method = "mean", leading_constant = c(2,2))

n = 1000
ys = matrix(rnorm(n*p), nrow = p, ncol = n)

for (i in 1:n) {
  detector = getData(detector, ys[,i])
}

attr(detector, 'grid')

attr(detector, 'allstats')$cumsum



# check grid and cumulative sum saving
gg = get_grid(n)
sum(abs(gg-attr(detector, 'grid'))) #grid OK

cumsum =  t(apply(ys, 1, cumsum))
sum(abs(cumsum[,n-gg] - attr(detector, 'allstats')$cumsums)) #cumulative sums are OK

sum(abs(attr(detector, 'allstats')$cumsum - cumsum[,n])) # last cumulative sum is also OK


# check rest


p = 100

detector = CHAD(p, method = "mean", leading_constant = c(2,2))

n = 100
ys = matrix(rnorm(n*p), nrow = p, ncol = n)
ys[,2:100] = ys[,2:100] + 1

detector = getData(detector, ys[,1])
attr(detector, "allstats")

detector = getData(detector, ys[,2])
attr(detector, "allstats")
attr(detector, "baseline_sd")

detector = getData(detector, ys[,3])
detector = getData(detector, ys[,4])
attr(detector, "allstats")

detector = CHAD(p, method = "mean", leading_constant = c(2,2))
for (i in 1:n) {
  detector = getData(detector, ys[,i])
}
attr(detector, "allstats")$Amatrix



cc = HDCD::CUSUM(ys[,])
sum(cc[,5]^2)-p
attr(detector, "allstats")$Amatrix[1,95]


sum(cc[,2]^2)-p
attr(detector, "allstats")$Amatrix[1,2]

# seems OK
print(detector)

attr(detector, "allstats")




p = 100

detector = CHAD(p, method = "mean", leading_constant = c(2.4,0.8))

n = 1000
ys = matrix(rnorm(n*p), nrow = p, ncol = n)
ys[,500:n] = ys[,500:n] + 0.0
ys[1:2,500:n] = ys[1:2,500:n] + 1

times = rep(NA, n)
times_ocd = rep(NA, n)
detector_ocd <- ocd::ChangepointDetector(dim=100, method='ocd',thresh=c(11.6, 179.5, 54.9), beta=1)
for (i in 1:n) {
  startt = proc.time()
  detector <- getData(detector, ys[,i])
  endd = proc.time()
  times[i] = (endd - startt)[3]

  startt = proc.time()
  detector_ocd <- ocd::getData(detector_ocd, ys[,i])
  endd = proc.time()
  times_ocd[i] = (endd - startt)[3]
}
plot(times) #helt rått
points(times_ocd,col=2)


attr(detector, 'grid')

attr(detector, 'allstats')$cumsum



## check if the test statistic is correct
p = 100
detector = CHAD(p, method = "mean", leading_constant = c(2,2))

n = 1000
ys = matrix(rnorm(n*p), nrow = p, ncol = n)
ys[,500:n] = ys[,500:n] + 0.0
ys[1:2,500:n] = ys[1:2,500:n] + 1

times = rep(NA, n)
times_ocd = rep(NA, n)
for (i in 1:n) {
  startt = proc.time()
  detector <- getData(detector, ys[,i])
  endd = proc.time()
  times[i] = (endd - startt)[3]
}
gg = get_grid(n)
cumsum =  t(apply(ys, 1, cumsum))
g = gg[100]

cc = HDCD::CUSUM(ys)
cc = cc[,n-g]

mean_penalty_and_threshold = mean_penalty_and_threshold(n,p)
as = mean_penalty_and_threshold[[1]]
#cat("as : ")
#print(as)
nu_as = mean_penalty_and_threshold[[2]]
r_values = mean_penalty_and_threshold[[3]]
s_values = mean_penalty_and_threshold[[4]]

Arec = rep(NA, length(nu_as))
for (i in 1:length(nu_as)) {
  a = as[i]
  nu_a = nu_as[i]
  ccc = cc[]
  ccc = ccc[abs(ccc)>a]
  if(identical(ccc,numeric(0))){
    Arec[i] = 0
  }else{
    ccc = ccc^2- nu_a
    Arec[i] = sum(ccc)
  }
}
Arec
attr(detector,'allstats')$Amatrix[,100]

cc = HDCD::CUSUM(ys)
plot(rev(colSums(cc^2)-p))
plot(attr(detector,'allstats')$Amatrix[1,])

plot(rev(colSums(cc^2)-p)/ r_values[1])
plot(attr(detector,'allstats')$A_scaled[1,])


v=6
a = as[v]
nu_a = nu_as[v]
ccc = cc[]
for (i in 1:dim(ccc)[1]) {
  for(j in 1:dim(ccc)[2]){
    if(abs(ccc[i,j])>a){
      ccc[i,j] = ccc[i,j]^2 - nu_a
    }else{
      ccc[i,j] = 0
    }

  }
}

plot(rev(colSums(ccc)))
plot(attr(detector,'allstats')$Amatrix[3,])

plot(rev(colSums(ccc)) / r_values[3])
plot(attr(detector,'allstats')$A_scaled[3,])

max((rev(colSums(ccc)) / r_values[3]))
max((attr(detector,'allstats')$A_scaled[3,]))

attr(detector,'statistics')
max((colSums(ccc)) / r_values[v])
max((attr(detector,'allstats')$A_scaled[v,]))



## check why it sucks so much:
N = 150 # number of data samples considered
chgptloc = round(N/3)
num_sim = 2*6 # number of iterations in the simulation
ps = c(100)
p = ps[1]
sparsities100 = c(1,5,10,100)
sparsities1000 = c(1,5,30,1000)
thetas = seq(0.0, 8.0, by=0.4)
num_methods = 5
num_cores = 6
MC_reps = 1000
false_alarm_prob = 0.05

s = p
theta = 1.6
ys = matrix(rnorm(N*p), nrow = p, ncol = N)
ys[1:s,(chgptloc+1):N] = ys[1:s,(chgptloc+1):N] + theta/sqrt(s)

detector = CHAD(p, method = "mean", leading_constant = thresholds[[1]][[v]],
                constant_penalty=TRUE)
detector_ocd <- ocd::ChangepointDetector(dim=p, method='ocd',thresh=thresholds[[2]][[v]])


for (i in 1:N) {
  detector <- getData(detector, ys[,i])
  detector_ocd <- ocd::getData(detector_ocd, ys[,i])

  detectors = list()
  detectors[[1]] = detector
  detectors[[2]] = detector_ocd

  if(!identical(status(detector_ocd),'monitoring')){
    print(i)
    break
  }

}
status(detector_ocd)
status(detector)
attr(detector, 'allstats')$Amatrix
attr(detector, 'allstats')$A_scaled

attr(detector, 'grid')
attr(detector,'statistics')
leading_constant(detector)
44.8688832 / mean_penalty_and_threshold(2,p)[[3]][1]








### simu to find appropriate thresholds
p = 100
n = 5000
N = 100
maxes = matrix(NA, nrow = 2, ncol = N)
for (j in 1:N) {
  if(j%% 10 == 0){
    print(j)
  }
  detector = CHAD(p, method = "mean", leading_constant = c(100,100))
  ys = matrix(rnorm(n*p), nrow = p, ncol = n)
  for (i in 1:n) {
    detector <- getData(detector, ys[,i])
  }
  mm = attr(detector, 'overall_max_statistics')
  maxes[,j] = mm
}

max(maxes[1,])
max(maxes[2,])

plot(times) #helt rått

attr(detector, 'grid')

attr(detector, 'allstats')$cumsum





