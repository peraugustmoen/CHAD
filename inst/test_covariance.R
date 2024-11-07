
p = 10

detector = CHAD(p, method = "covariance", leading_constant = c(2))

n = 100
ys = matrix(rnorm(n*p), nrow = p, ncol = n)

cumsum = array(NA, dim = c(n, p,p))
cumsum[1,,] = ys[,1] %*% t(ys[,1])

for (i in 1:n) {
  detector = getData(detector, ys[,i])
  if(i>1){
    cumsum[i,,] = cumsum[i-1,,] + ys[,i] %*% t(ys[,i])
  }
}

attr(detector, 'grid')

attr(detector, 'allstats')$cumsum


# check grid and cumulative sum saving
gg = get_grid(n)
sum(abs(gg-attr(detector, 'grid'))) #grid OK

sum(abs(cumsum[n-gg,,] - attr(detector, 'allstats')$cumsums)) #cumulative sums are OK

sum(abs(attr(detector, 'allstats')$cumsum - cumsum[n,,])) # last cumulative sum is also OK

# check log2cumulativesum
dim(attr(detector, 'allstats')$log2cumsum)
2^(0:6)

sum(abs(attr(detector, 'allstats')$log2cumsum - cumsum[2^(0:6),,]))


## check the operator norms

p = 10

detector = CHAD(p, method = "covariance", leading_constant = c(2),mean_centering=TRUE)

n = 100
ys = matrix(rnorm(n*p), nrow = p, ncol = n)

cumsum = array(NA, dim = c(n, p,p))
cumsum[1,,] = ys[,1] %*% t(ys[,1])

teststats = rep(NA, n)
for (i in 1:n) {
  detector = getData(detector, ys[,i])
  if(i>1){
    cumsum[i,,] = cumsum[i-1,,] + ys[,i] %*% t(ys[,i])
    teststats[i] = attr(detector, 'statistics')
  }
}
teststats
max(na.omit(teststats))
attr(detector, 'allstats')$teststats

# check cumsum on means
attr(detector, 'allstats')$cumsum_mean
cumsum =  t(apply(ys, 1, cumsum))
gg = get_grid(n)
sum(abs(cumsum[,n-gg] - attr(detector, 'allstats')$cumsums_mean)) #cumulative sums are OK

sum(abs(attr(detector, 'allstats')$cumsum_mean - cumsum[,n])) # last cumulative sum is also OK

# check log2cumulativesum
dim(attr(detector, 'allstats')$log2cumsum)
2^(0:6)

sum(abs(attr(detector, 'allstats')$log2cumsum_mean - cumsum[,2^(0:6)]))







dat = read.csv("/Users/peraugustmoen/Library/CloudStorage/OneDrive-UniversitetetiOslo/project3/exchange_rates.csv")
dat[dat=="ND"] = NA
dat$time = as.POSIXct(dat$Series.Description, format="%Y-%m-%d", tz = "UTC")


plot(dat$time, dat$INDIA....SPOT.EXCHANGE.RATE..RUPEES.US..)

old_dat = dat[,]

dat = na.omit(dat)


dat = dat[dat$time >=as.POSIXct("2003-01-01", format="%Y-%m-%d", tz = "UTC"),]


#vars = c(2,3)
vars = c(2,3,5,6,8,9,12,13,17,23,25)
p = length(vars)
X = dat[,vars]
X = as.matrix(sapply(X, as.numeric))
X = t(X)

n = dim(X)[2]
detector = CHAD(p, method = "covariance", leading_constant = c(5),mean_centering=TRUE)

teststats = rep(NA, n)
changepoints = c()
argmax = c()
counter = 0
for (i in 1:n) {
  counter = counter+1
  detector = getData(detector, X[,i])
  if(counter>1){
    teststats[i] = attr(detector, 'statistics')
  }
  if(!identical(status(detector), 'monitoring')){
    changepoints = c(changepoints, i)
    argmax = c(argmax, attr(detector, "allstats")$argmax)
    detector = reset(detector)
    counter = 0
  }
}
teststats
status(detector)

attr(detector, 'allstats')$teststats
plot(teststats)
which.max(teststats)
dat$time[which.max(teststats)]
dat$time[changepoints]
dat$time[changepoints - argmax]
