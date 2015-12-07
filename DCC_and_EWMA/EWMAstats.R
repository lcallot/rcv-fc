
library(magrittr)

# Loading the data
nstock <- 30
diag.ind  <-rep(0,nstock*(nstock+1)/2)
diag.ind[cumsum(1:nstock)]  <-1

# dates
dates <- as.Date(read.table('..//data/dates')$V1)
# Dates in week and month format
week  <- strftime(dates,format='%Y-%W')
month <- strftime(dates,format='%Y-%m')

load('EWMA_forecasts_scaled')
load('W_EWMA_forecasts_scaled')
load('M_EWMA_forecasts_scaled')

# Load the RCovs:
rcv <- get(load('..//data/CRK-har-dj-cens-lcov'))$d[-c(1:1019),]
wrcv <- get(load('..//data/CRK-dj-cens-lcov-W'))[-c(1:263),]
mrcv <- get(load('..//data/CRK-dj-cens-lcov-M'))[-c(1:60),]



rcv[,diag.ind==1] <- exp(rcv[,diag.ind==1])
wrcv[,diag.ind==1] <- exp(wrcv[,diag.ind==1])
mrcv[,diag.ind==1] <- exp(mrcv[,diag.ind==1])


#forecast error
dfcerr <- rcv-do.call(rbind,covfc)
wfcerr <- wrcv-do.call(rbind,wcovfc)
mfcerr <- mrcv-do.call(rbind,mcovfc)

#lazy hack...
for(fcerr in list(dfcerr,wfcerr,mfcerr)){
  # storage
  stats <- matrix(NA,nrow=3,ncol=3)
  rownames(stats) <- c('MedAFE','MaxAFE','Frobenuis')
  colnames(stats) <- c('A','D','O')
  
  for(di in list(0,1,c(0,1))){
    ci <- ifelse(prod(di==0),3,ifelse(prod(di==1),2,1))
    
    stats[1,ci] <- mean(apply(abs(fcerr[,diag.ind%in% di]),1,median)) #MedAFE
    stats[2,ci] <- mean(apply(abs(fcerr[,diag.ind%in% di]),1,max)) #MaxAFE
    stats[3,ci] <- mean(apply(fcerr[,diag.ind%in% di]^2,1,function(x)sqrt(sum(x)))) #Frobenius
  }
  print(stats)
  
  stats %>% t %>% c %>% round(.,2) %>% paste(.,collapse=' & ') %>% print
  
  
  
}
  
