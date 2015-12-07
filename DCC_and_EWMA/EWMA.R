# EWMA estimation and forecasting of the Dow-Jones covariance matrix.
#
# Author: Laurent Callot
# 17/05/2015
library(parallel)
library(magrittr)
source('../subs//ewma_subs.R')

ncores <- 2

# Loading the data
dj <- 100*read.table('../data/returns_dja.dat') %>% tail(.,-21)
nstock <- 30
diag.ind  <-rep(0,nstock*(nstock+1)/2)
diag.ind[cumsum(1:nstock)]	<-1
# dates
dates <- as.Date(read.table('..//data/dates')$V1)
# nbr forecasts
nfc <- 455

cat('\nEWMA estimation start.\n')
ts <- proc.time()
# covfc <-  mclapply(1:nfc,mc_ewma,dj,mc.cores = ncores)
# save(file='EWMA_forecasts_scaled',covfc)
te <- proc.time()

cat('\nEWMA estimation time: ',round(te-ts,2)[3],'sec.\n')


# WEEKLY and MONTHLY FORECASTS

# Dates in week and month format
week  <- strftime(dates,format='%Y-%W')
month <- strftime(dates,format='%Y-%m')

# daily data untrimmed
ddj <- read.table('../data/returns_dja.dat')
#		return aggregation
cumpct <- function(x) {
#  y <- cumprod(1+ x) - 1
  y <- mean( x) 
  return(as.numeric(tail(y,1)))
}

wdj <- 100*apply(ddj,2,by,week,cumpct)
mdj <- 100*apply(ddj,2,by,month,cumpct)


cat('\nAggregated EWMA estimation start.\n')
ts <- proc.time()
 wcovfc <-  mclapply(1:52,mc_ewma,wdj,agg='w',mc.cores = ncores)
 mcovfc <-  mclapply(1:12,mc_ewma,mdj,agg='m',mc.cores = ncores)

 save(file='W_EWMA_forecasts_scaled',wcovfc)
 save(file='M_EWMA_forecasts_scaled',mcovfc)
te <- proc.time()

cat('\nAggregated EWMA estimation time: ',round(te-ts,2)[3],'sec.\n')