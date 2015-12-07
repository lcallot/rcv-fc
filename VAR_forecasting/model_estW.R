# Libraries
require('reshape2')
require('ggplot2')
require('Matrix')
require('SparseM')
require('glmnet')
require('xtable')
require('expm')
require('plyr')
require('lassovar')


# Sourcing the subs
source('../subs/plot_subs.R')
source('../subs/rv_subs.R')
source('../subs/ptf_subs.R')


# catching the nbr of cores from the command line
argv = commandArgs(trailingOnly=TRUE)
cat('argv',argv,'\n',sep='')

ncores <- 12
if(length(argv)) ncores <- argv[1]

cat('ncores: ',ncores,'\n',sep='')


# Estimation settings
horizon		<-1
fc.window	<-'fix'
crit		<-'BIC'


# Creating the list of diagonal indices. 
nstock     <- 30
diag.ind	<-rep(0,nstock*(nstock+1)/2)
diag.ind[cumsum(1:nstock)]	<-1

models <- matrix(rbind(
  c('var',1,'Lasso','none','dj.cens.lcov.M',263,'none'),
  c('var',1,'Lasso','none','dj.cens.lcov.M',263,'ar'),
  c('var',1,'Lasso','lasso','dj.cens.lcov.M',263,'none'),
  c('var',5,'Lasso','none','dj.cens.lcov.M',263,'none'),
  c('var',5,'Lasso','none','dj.cens.lcov.M',263,'ar'),
  c('var',5,'Lasso','lasso','dj.cens.lcov.M',263,'none')
),ncol=7,
dimnames=c(list('Model'=NULL,'spec'=c('Model','Lag','Estimator','Adaptive','Data','Est.smpl','Restrictions'))))



# Main estimation function, calls forecast.lassovar and sorts out the results. 
fm <- fc.mod(models,horizon=1,fc.window=fc.window,crit=crit,ncores=ncores,agg=TRUE,diag.ind=diag.ind,post=TRUE)


