# DCC estimation and forecasting of the Dow-Jones covariance matrix.
#
# Author: Laurent Callot
# 11/04/2015

library(plyr)
library(parallel)
library(rmgarch)
library(magrittr)

# Specs inspired from :
# http://unstarched.net/2013/01/03/the-garch-dcc-model-and-2-stage-dccmvt-estimation/

nstock <- 30
diag.ind  <-rep(0,nstock*(nstock+1)/2)
diag.ind[cumsum(1:nstock)]	<-1
# dates
dates <- as.Date(read.table('..//data/dates')$V1)

# Dates in week and month format
week  <- strftime(dates,format='%Y-%W')
month <- strftime(dates,format='%Y-%m')

# daily data untrimmed
ddj <- read.table('../data/returns_dja.dat')
#  	return aggregation
cumpct <- function(x) {
  y <- mean(x)
 # y <- cumprod(1+ x) - 1
  return(as.numeric(tail(y,1)))
}

#  	return aggregation
wdj <- 100*apply(ddj,2,by,week,cumpct) 
mdj <- 100*apply(ddj,2,by,month,cumpct) 

# Creating the cluster
ncores <- 16
cl <- makeForkCluster(ncores)

# timer init
ti <- proc.time()
# GARCH-spec
# univariate GARCH specification
#ugs <- ugarchspec() 
ugs <- ugarchspec(mean.model = list(include.mean=TRUE),
                  variance.model = list(garchOrder = c(1,1), model = 'sGARCH'),
                  distribution.model = 'norm')
mspec <- multispec(replicate(nstock, ugs))

# fitting the univariate GARCHS
#ufit  <- multifit(mspec,data = dj,cluster = cl)
wufit  <- multifit(mspec,data = wdj,cluster = cl)
#mufit  <- multifit(mspec,data = mdj,cluster = cl)
# Uni-garch time
tu <- proc.time()
cat('\nUnivarite GARCH estimation time: ',round(tu-ti,2)[3],'sec.')

# DCC spec
#dccs <- dccspec(uspec = mspec, dccOrder = c(1, 1), distribution = 'mvnorm')
wdccs <- dccspec(uspec = mspec, dccOrder = c(1, 1), distribution = 'mvnorm')
#mdccs <- dccspec(uspec = mspec, dccOrder = c(1, 1), distribution = 'mvnorm')

# DCC fit
#dccf <- dccfit(dccs, data = dj, fit.control = list(eval.se = FALSE), fit = ufit, cluster = cl)
wdccf <- dccfit(wdccs, data = wdj, fit.control = list(eval.se = FALSE), fit = wufit, cluster = cl)
#mdccf <- dccfit(mdccs, data = mdj, fit.control = list(eval.se = FALSE), fit = mufit, cluster = cl)

# DCC time
td <- proc.time()
cat('\nDCC estimation time: ',round(td-tu,2)[3],'sec.')

#rolling forecasts with full data
dj <- 100*read.table('../data/returns_dja.dat') %>% tail(.,-21)
dj <- dj[,1:nstock]


#dccr <- dccroll(dccs,dj,n.ahead=1,forecast.length = 455,refit.every = 1,window.size = 1000,cluster = cl,fit.control = list('eval.se'=FALSE))
wdccr <- dccroll(wdccs,wdj,n.ahead=1,forecast.length = 52,refit.every = 1,window.size = 1000,cluster = cl,fit.control = list('eval.se'=FALSE))
#mdccr <- dccroll(mdccs,mdj,n.ahead=1,forecast.length = 12,refit.every = 1,window.size = 1000,cluster = cl,fit.control = list('eval.se'=FALSE))
# DCC-roll time
tr <- proc.time()
cat('\nDCC estimation time: ',round(tr-td,2)[3],'sec.')

stopCluster(cl)

# Extract the forecasted covariance matrix
#covfc <- rcov(dccr)
wcovfc <- rcov(wdccr)
#mcovfc <- rcov(mdccr)

# Transforms a covariance matrix to a vector.
cov2vec<-function(x){return(x[upper.tri(x,diag=TRUE)])}


# vectorize
#covfc <- t(apply(covfc,3,cov2vec))
#save(file='DCC_forecasts_scaled',covfc)

wcovfc <- t(apply(wcovfc,3,cov2vec))
save(file='W_DCC_forecasts_scaled',wcovfc)

#mcovfc <- t(apply(mcovfc,3,cov2vec))
#save(file='M_DCC_forecasts_scaled',covfc)


cat('\nTotal estimation time: ',round(tr-ti,2)[3],'sec.')

