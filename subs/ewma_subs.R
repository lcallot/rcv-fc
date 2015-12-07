
# EWMA function adapted from EWMAvol.R from Ruey Tsay
# source: http://faculty.chicagobooth.edu/ruey.tsay/teaching/mts/sp2013/EWMAvol.R

## Compute exponentially weighted moving average covariance matrix.
ewma <- function(rtn,lambda=0.96,agg='d'){
  if(!is.matrix(rtn))rtn=as.matrix(rtn)
  Mean=apply(rtn,2,mean)
  nT=dim(rtn)[1]; 
  k=dim(rtn)[2]
  x=rtn
  for (i in 1:k){
    x[,i]=rtn[,i]-Mean[i]
  }

  # init EWMA
  Sigt=cov(rtn)
  V1=cov2vec(Sigt)
  for (t in 2:nT){
    xl=x[t-1,]
    for (i in 1:k){
    Sigt[i,]= (1-lambda)*xl*xl[i]+lambda*Sigt[i,]
    }
  V1=rbind(V1,cov2vec(Sigt))
  }

return(V1)
}

mc_ewma <- function(i,dj,agg){
  
  if(agg=='d') trainobs <- 1000
  if(agg=='w') trainobs <- 263
  if(agg=='m') trainobs <- 60
  
  rtn <- dj[i+(0:(trainobs-1)),]
  fc <- ewma(rtn,lambda=0.96)[trainobs,]
  return(fc)
}


# Transforms a sym matrix into a vector of unique elements
cov2vec <- function(covmat){covmat[upper.tri(covmat,diag=TRUE)]}
