
#	Computes the forecasts and save the R object returned for a set of models in mod.smpl.

fc.mod <- function(mod.smpl,horizon=1,fc.window='fix',crit='BIC',ncores=ncores,agg=FALSE,diag.ind,post=TRUE){
  
  #Loop across models
	for(m in 1:nrow(mod.smpl)){
		m.vec <- mod.smpl[m,]	# The model row
		cat('\n Forecasting with Model:\n')
		cat(paste(m.vec,collapse=', '))
		mn <- paste(mod.smpl[m,],collapse='.')	#model names
		cat('\n Model Name: \n')
		cat(mn)
		if(!agg)dn <- load(paste('../data/CRK-har',gsub('\\.','-',mod.smpl[m,'Data']),sep='-')) #data load + assign name
		if(agg)dn <- load(paste('../data/CRK',gsub('\\.','-',mod.smpl[m,'Data']),sep='-')) #data load + assign name
		
		cat('\n Data set name: ',dn)
		cat('\n')
		ada <- NULL
		if(m.vec['Adaptive']!='none')ada <- m.vec['Adaptive']

    trans=NULL
		if(length(grep('lmat',mn)))trans <- 'lmat'
		if(length(grep('lcov',mn)))trans <- 'lcov'
    
    
		if(!agg){
			CRK	 <- get(dn)[['d']]
			CRK.w	 <- get(dn)[['w']]
			CRK.m	 <- get(dn)[['m']]
		}
		if(agg){
			CRK <- get(dn)
			CRK.w <- NULL
			CRK.m <- NULL
			}

		nstock	 <- (-1+sqrt(1+8*ncol(CRK)))/2
		assign('nbr.stock',nstock,envir=.GlobalEnv)




		if(m.vec['Model']=='var'){
      
			if(m.vec['Estimator']=='Lasso'){					   
				assign(paste('time',mn,sep='.'),system.time(assign(paste('fc',mn,sep='.'),
					forecast.lassovar(CRK,ntrain=as.numeric(m.vec['Est.smpl']),
					horizon=horizon,fc.window=fc.window,ic=crit,
					lags=as.numeric(m.vec['Lag']),type='equation',adaptive=ada,
					mc=TRUE,mclas=FALSE,ncores=ncores,
					rest=m.vec['Restrictions'],post=post))))}
      
			if(m.vec['Estimator']=='ols'){	
				assign(paste('time',mn,sep='.'),system.time(assign(paste('fc',mn,sep='.'),
					forecast.olsvar(CRK,ntrain=as.numeric(m.vec['Est.smpl']),
					horizon=horizon,fc.window=fc.window,ic=crit,
					lags=as.numeric(m.vec['Lag']),type='equation',adaptive=ada,
					mc=TRUE,ncores=ncores,exo=NULL,
					rest=m.vec['Restrictions']))))}			   

		}
		if(m.vec['Model']=='har'){
      
			if(m.vec['Lag']=='d.w')xo <- CRK.w
			if(m.vec['Lag']=='d.w.m')xo <- cbind(CRK.w,CRK.m)
      
			assign(paste('time',mn,sep='.'),system.time(assign(paste('fc',mn,sep='.'),
					forecast.lassovar(CRK,ntrain=as.numeric(m.vec['Est.smpl']),
						horizon=horizon,fc.window=fc.window,ic=crit,
						lags=1,type='equation',
						adaptive=ada,mc=TRUE,mclas=FALSE,ncores=ncores,
						covrest=ifelse(m.vec['Restrictions']=='rest',TRUE,FALSE),
						exo=xo))))
		}
    
		if(m.vec['Model']=='NoChange'){
			assign(paste('time',mn,sep='.'),system.time(assign(paste('fc',mn,sep='.'),
					  fc.nochange(CRK,training.init=as.numeric(m.vec['Est.smpl']),horizon=horizon))))
		}
		
    fc.est <- get(paste('fc',mn,sep='.'))
		fc.est$CRK <- CRK 
		assign(paste('fc',mn,sep='.'),fc.est)

		save(list=paste('fc',mn,sep='.'),file=paste('fc_data/fc',mn,sep='.'))
		save(list=paste('time',mn,sep='.'),file=paste('time/time',mn,sep='.'))
		rm(list=c(paste('fc',mn,sep='.'),'CRK','fc.est','CRK.w','CRK.m',dn))
		
		gc()
	}
	return(TRUE)
}



#No change forecasts:
fc.nochange <- function(rv.data,training.init,horizon)
{
	
	rv.diff<-diff(rv.data[-c(1:(training.init-1)),],lag=horizon,diff=1)
	rv.ncfc<-rv.data[-c(1:(training.init-1),nrow(rv.data)),]

	return(list('err'=rv.diff,'pred'=rv.ncfc))
}


# fixed window average, internal function
wd.avg<-function(x,wd.length){
	ldat	<-length(x)
	xl	<-x

	for(w in 1:(wd.length-1)){
		xw.tmp		<-x[-c(1:w)]
		length(xw.tmp)	<-ldat
		xl		<-cbind(xl,xw.tmp)
	}
	wavg	<-rowMeans(xl)

	#The padding with NA is not done in the right place, gotta manually move them from tail to end:
	wavg	<-c(rep(NA,wd.length-1),wavg[-c((length(wavg)-wd.length+2):length(wavg))])

	return(wavg)
}

# Lag selection, internal function
lagseleq<-function(beq.grp,neq,nlag){
	lnz	<-(unique(ceiling(which(beq.grp!=0)/neq)))
	z	<-rep(0,nlag)
	z[lnz]<-1
	return(z)
}


# Takes a forecast.lassovar object, the transformation applied to the data, the index of diagonal elements,
# and the observed covariances.
# Returns the untransformed predicted and observed covariance matrices.
fc.covmat<- 
function(fc.lv,transfo=c(NULL,'lcov','lcor','cov','lmat'),diag.ind,obs.cov){
	
	fc.raw	<- fc.lv$pred

	if(transfo=='lcor'|transfo=='lcov'){
		#rescaling the log variances
		fc.raw[,diag.ind==1]	<-exp(fc.raw[,diag.ind==1])
		obs.cov[,diag.ind==1]	<-exp(obs.cov[,diag.ind==1])
	}
	
	if(transfo=='lcor'){
		#Transforming back the correlations into covariances. that's going to be messy...
		fc.raw	<-t(apply(fc.raw,1,varcor2cov,diag.ind))
		obs.cov	<-t(apply(obs.cov,1,varcor2cov,diag.ind))
	}

	
	if(transfo=='lcor'|transfo=='lcov'){
		fc.raw	<-array(apply(fc.raw,1,covvec2covmat,diag.ind),dim=c(sum(diag.ind==1),sum(diag.ind==1),nrow(fc.raw)))
		obs.cov	<-array(apply(obs.cov,1,covvec2covmat,diag.ind),dim=c(sum(diag.ind==1),sum(diag.ind==1),nrow(obs.cov)))
	}

	if(transfo=='lmat'){
		fc.raw	<-array(apply(fc.lv$pred,1,lmatvec2covmat,diag.ind),dim=c(sum(diag.ind==1),sum(diag.ind==1),nrow(fc.lv$pred)))
		obs.cov	<-array(apply(obs.cov,1,lmatvec2covmat,diag.ind),dim=c(sum(diag.ind==1),sum(diag.ind==1),nrow(fc.lv$pred)))
	}
	return(list('obs.cov'=obs.cov,'pred'=fc.raw))
}


# Takes a log-matrix transformed covariance matrix in vector form, returns its matrix exponential in matrix form
lmatvec2covmat <- function(x,diag.ind){

	mx <- matrix(0,ncol=sum(diag.ind==1),nrow=sum(diag.ind==1))
	mx[upper.tri(mx,diag=TRUE)]	<-x
	mx	<-mx+ t(mx)
	diag(mx) <- diag(mx)/2

	return(expm(mx))	

}

# Takes a covariance matrix in vector form, returns it in matrix form
covvec2covmat <- 
function(x,diag.ind){
	
	mx <- matrix(0,ncol=sum(diag.ind==1),nrow=sum(diag.ind==1))
	mx[upper.tri(mx,diag=TRUE)]	<-x
	mx	<-mx+ t(mx)
	diag(mx) <- diag(mx)/2

	return(mx)
	}




#Takes a varcor vector (1 day) and transforms it back to a covariance vector 
varcor2cov<-
function(x,diag.ind){
	
	#creating index vectors
	N	<-sum(diag.ind==1)	

	c1<-c2<-varcov<-rep(0,length(diag.ind))

	for(n in 1:N)	{c1[cumsum(1:(n-1))[n-1] + c(1:(n-1))]<-1:(n-1)}
	for(n in 2:N)		{c2[cumsum(1:(n-1))[n-1] + c(1:(n-1))]<-n}
	
	varcov[diag.ind==2] <- sqrt(x[diag.ind==1][c2]*x[diag.ind==1][c1])*x[diag.ind==2]
	varcov[diag.ind==1] <- x[diag.ind==1]

	return(varcov)
}

# Transforms a covariance matrix to a vector.
cov2vec<-function(x){return(x[upper.tri(x,diag=TRUE)])}



# Takes a forecast.lassovar object and some extra information.
# Returns a vector with the MSE medSE and MaxSE untransformed. 
fc.statvec <- function(fc,diag.ind,xtr.ind=NULL,trans='lcov',ind='dj')
{
	fc.stat <- dim(fc$err)
	if(is.null(xtr.ind))	err	<- .raw.err(fc,diag.ind[[ind]],fc$CRK,trans)
	if(!is.null(xtr.ind))	err	<- .raw.err(fc,diag.ind[[ind]],fc$CRK,trans)[-xtr.ind,]

	for(i in unique(diag.ind[[ind]])){
		fc.stat <- c(fc.stat,mean(colMeans(err^2)[diag.ind[[ind]]==i]))} #MSE diag and off
	
	for(i in unique(diag.ind[[ind]])){
		fc.stat <- c(fc.stat,median(apply(err[,diag.ind[[ind]]==i]^2,2,median)))} #  median of median squared error

	for(i in unique(diag.ind[[ind]])){
		fc.stat <- c(fc.stat,mean(apply(err[,diag.ind[[ind]]==i]^2,2,max)))} # Largest squared error

	return(fc.stat)
}


.raw.err <- function(fc,diag.ind,CRK,trans)
{
	fc.nfc<- nrow(fc$err)
	fc.cm <- fc.covmat(fc,transfo=trans,diag.ind,obs.cov=CRK[-c(1:(nrow(CRK)-fc.nfc)),])	
	rm(fc)
	gc()

		err <- t(apply(fc.cm$obs.cov-fc.cm$pred,3,cov2vec))
	return(err)
}



.xtrm.obs <- function(x,prob=0.99)
	{
	co <- quantile(x,probs=prob,na.rm=TRUE)
	xtr.ind<-which(x>=co,arr.ind=TRUE)
	return(xtr.ind)
	}

fc.xtrm.err <- function(fc,diag.ind,CRK,nc=FALSE,trans)
	{
	err	<- .raw.err(fc,diag.ind,CRK,nc,trans)
	if(!nc)err	<-rbind(matrix(NA,ncol=ncol(err),nrow=nrow(CRK)-nrow(err)),err)
	if(nc)err	<-rbind(matrix(NA,ncol=ncol(err),nrow=dim(CRK)[3]-nrow(err)),err)

	xtr.ind <- apply(err^2,2,.xtrm.obs,prob=0.99)
	
	return(xtr.ind)
	}



# Computes ols based forecasts.
forecast.olsvar<-function(dat,training.init,horizon=1
                          ,fc.window=c('expanding','rolling'),fc.type=c('direct','recursive'),mc=FALSE,... )
{
	#setting up stuff
	start.date	<-start(dat)
	end.date	<-end(dat)
	freq		<-frequency(dat)
	fc.lassovar	<-list(err=NULL)

	#Optional args
	argnames<- names(list(...)) 
	argList	<- list(...)
	ic	<-argList$ic
	lags	<-argList$lags
	type	<-argList$type
	exo	<-argList$exo
  rest	<-argList$rest
  trans <-argList$trans
  if(is.null(trans))trans<-'lcov'
  if(!(trans%in%c('lcov','lmat')))stop(cat('\n unknown transformation: ',trans))
  
	if(is.null(lags))lags<-1
	if(is.null(type))type<-'equation'
	if(is.null(rest))rest<-'none'

	fc.lv	<-list('call'=match.call(),'err'=NULL,'pred'=NULL,'coefficients'=list())	

	if(is.null(argList$ncores)){ncores=1}
	else ncores=argList$ncores

	fc.init	<- window(dat,start=start.date,end=training.init)
	fc.date	<-end(fc.init)

	nbr.fc	<-nrow(dat)-nrow(fc.init)

	cat('\n		-----------------------------		\n')	
	cat(fc.window,' window forecasts\n',sep='')
	cat(horizon,'-steps ahead ',fc.type,' forecasts. Initial training sample: ',fc.date[1],'-',fc.date[2],'\n','Number of forecasts: ',nbr.fc,'\n',sep='')
	cat('Full sample start:',start.date[1],'-',start.date[2],' end:',end.date[1],'-',end.date[2],'\n',sep='')



	if(!mc)fc.olsvar<-lapply(0:(nbr.fc-horizon),.fc.loop.olsvar
                           ,dat,fc.init,start.date,freq,fc.window,lags,horizon
                           ,exo=exo,rest=rest)
	if(mc)fc.olsvar<-mclapply(0:(nbr.fc-horizon),.fc.loop.olsvar,dat,fc.init
                            ,start.date,freq,fc.window,lags,horizon
                            ,exo=exo,rest=rest
                            ,mc.cores=ncores)
	


	fc.lv$spectest <- array(NA,dim=c(3,ncol(dat),nbr.fc),dimnames=list('Test'=c('Autocorr','Normality','Rsquared'),'Equation'=1:ncol(dat),'Forecast'=1:nbr.fc))

	fc.i <- 1
	for(lv.tmp in fc.olsvar)
		{
			fc.lv$err<-rbind(fc.lv$err,lv.tmp$err)
			fc.lv$pred<-rbind(fc.lv$pred,lv.tmp$pred)
			fc.lv$coefficients<-c(fc.lv$coefficients,Matrix(lv.tmp$coefficients,sparse=TRUE))
		
			#fc.lv$spectest[,,fc.i] <- lv.tmp$spectest
			fc.i <- fc.i+1
		}


return(fc.lv)
}



# internal, produces direct ols forecasts for benchmarking.
.fc.loop.olsvar <-
function(fc,dat,fc.init,start.date,freq
         ,fc.window,lags,horizon,exo
         ,rest)
{
	if(fc.window=='fix'){start.fc	<- incdate.ols(start.date,fc,freq)}
	else{start.fc	<-start.date}

	end.fc	<-incdate.ols(start.date,fc+nrow(fc.init)-1,freq)
	fc.date	<-incdate.ols(end.fc,1,freq)

	train.dat	<-window(dat,start=start.fc,end=end.fc)
	fc.dat		<-window(dat,start=end.fc,end=fc.date)[2,]

	y.var		<-mkvar.ols(train.dat,lags,horizon,exo)

	if(rest=='none'){
    dat     <- data.frame(y=y.var$y,x=y.var$x)
    fmla    <- as.formula(paste('cbind(',paste(colnames(dat)[grep('y',colnames(dat))],collapse=',')
                                ,') ~ ', paste(colnames(dat)[grep('x',colnames(dat))], collapse= "+")))
    ols.mod	<-lm(fmla,data=dat)
    
		ols.pred<-predict(ols.mod,tail(dat,1))
    ols.err <- fc.dat - ols.pred
		ols.coeff <- coef(ols.mod)
	}
  if(rest=='ar'){
    lind <- ncol(y.var$y) * (0 : (ncol(y.var$x)/ncol(y.var$y)-1)) 
		ols.err<-matrix(0,nrow=1,ncol=ncol(y.var$y))
		ols.pred<-matrix(0,nrow=1,ncol=ncol(y.var$y))
		ols.coeff<-matrix(0,nrow=(ncol(y.var$x)+1),ncol=ncol(y.var$y))
    for(i in 1: ncol(y.var$y)){
      # estimation using biglm
      dat     <- data.frame(y=y.var$y[,i],x=y.var$x[,lind+i])
      fmla    <- as.formula(paste("y ~ ", paste(colnames(dat)[-1], collapse= "+")))
      ols.mod	<-biglm(fmla,data=dat)
           
  		ols.pred[i]<-predict(ols.mod,tail(dat,1))
      ols.err[i] <-(fc.dat[i] - ols.pred[i])
  		ols.coeff[c(1,lind + i + 1),i] <- coef(ols.mod)
    }
	}
    
	if(rest=='covrest'){
		ols.err<-matrix(0,nrow=1,ncol=ncol(y.var$y))
		ols.pred<-matrix(0,nrow=1,ncol=ncol(y.var$y))
		ols.coeff<-matrix(0,nrow=(ncol(y.var$x)+1),ncol=ncol(y.var$y))
		for(i in 1:ncol(y.var$y)){
			cov.keep <-cov.keep2(i,nbr.stock)	
			cov.excl <-setdiff(1:ncol(y.var$y),cov.keep)
			all.excl <-cov.excl
			if(ncol(y.var$x)>ncol(y.var$y))for(l in 2:ncol(y.var$x)/ncol(y.var$y))cov.excl<-c(all.excl,cov.excl+(l-1)*ncol(y.var$y))
			ols.mod	<-lm(y.var$y[,i]~y.var$x[,-c(all.excl)])
			ols.pred[i]<- predict(ols.mod,y.var$x)[nrow(y.var$x)]
			ols.err[i]<-(fc.dat-predict(ols.mod,y.var$x))[nrow(y.var$x)]
			ols.coeff[c(1,1+cov.keep),i]<-ols.mod$coefficients
		}
	}
  sptest<-NULL #OBSOLETE
	return(list('err'=ols.err,'pred'=ols.pred,'coefficients'=ols.coeff,'spectest'=sptest))
}



# BEURK!! This could be done in 5 lines spitting out dataframe. also change in lassovar.  
mkvar.ols <-
function(data.var,lags,horizon,exo=NULL)
{
	nbrser	<-ncol(data.var)
	y	<-data.var
	if(!is.null(colnames(data.var)))var.names<-colnames(data.var)
	if(is.null(colnames(data.var)))	var.names<-paste('Eq_',1:nbrser,sep='')

	if(!('ts'%in%class(y))) {y<-ts(y,start=1,frequency=1);data.var<-ts(data.var,start=1,frequency=1)}

	for(l in 1:lags)
		{y	<-ts.intersect(y,lag(data.var,-l-horizon+1))
		var.names<-c(var.names,paste(l,'L_',var.names[1:nbrser],sep=''))
		}
	colnames(y)<-var.names
	
	if(is.null(exo))x<-cbind(y[,-c(1:nbrser)])
	if(!is.null(exo)){
		if(!('ts'%in%class(exo))){exo<-ts(exo,start=1,frequency=1)}
		x<-ts.intersect(y[,-c(1:nbrser)],lag(exo,-1))
	}


	y.var<-list('y'=y[,c(1:nbrser)],'x'=x,'lags'=lags,'horizon'=horizon,'nbrser'=nbrser)
	return(y.var)
}

cov.keep2<- function(i,N)
{
	ind.mat	<-matrix(NA,N,N)
	ind.mat[upper.tri(ind.mat,diag=TRUE)]<-1:(N*(N+1)/2)
	x	<-ind.mat[upper.tri(ind.mat)]
	ind.mat	<-t(ind.mat)
	ind.mat[upper.tri(ind.mat,diag=FALSE)]<-x

	if(i%in%diag(ind.mat))incl<-unique(c(ind.mat[which(ind.mat==i,arr.ind=TRUE)[,1],],diag(ind.mat)))
	else incl<-unique(c(ind.mat[which(ind.mat==i,arr.ind=TRUE)[,1],],diag(ind.mat)))

	return(incl)
}



# A little function to increment a 2-element date by inc period given a frequency 
incdate.ols<-
function(date,inc,freq)
{
	m	<-date[2]
	y	<-date[1]

	nd	<-freq*y + m-1
	nnd	<-nd+inc

	ny	<-nnd%/%freq
	nm	<-nnd%%freq +1
	dateinc	<-c(ny,nm)
	return(dateinc)
}


# de-transform the  predictions     
# OBSOLETE??
.rawpred <- function(pred,trans,diag.ind)
{
  hmax <- ifelse(length(dim(pred))==3,dim(pred)[3],1)
  
  if(hmax>1)rawpred <- array(NA,dim=dim(pred))
  if(hmax==1)rawpred <- array(NA,dim=c(dim(pred),1))
  
  
  
  for(h in 1:hmax){
    for(r in 1:dim(pred)[1]){
      if(trans=='lmat')
      {
        if(hmax>1)rawpredtmp <- lmatvec2covmat(pred[r,,h],diag.ind)
        if(hmax==1)rawpredtmp <- lmatvec2covmat(pred[r,],diag.ind)
        rawpred[r,,h]  <- cov2vec(rawpredtmp)
      }
      if(trans=='lcov')
      {
        if(hmax>1)rawpred[r,,h] <- pred[r,,h]
        if(hmax==1)rawpred[r,,h] <- pred[r,]
        rawpred[r,diag.ind==1,h] <- exp(rawpred[r,diag.ind==1,h])
      }
      if(trans=='lcpd')
      {
        if(hmax>1) prd <- pred[r,,h]
        if(hmax==1) prd <- pred[r,]
        prd[,diag.ind==1] <- exp(prd[,diag.ind==1])
        # transform rawpred to pd.
        prd <- vecreg(prd,diag.ind=diag.ind)
        rawpred[r,,h] <- prd #storing
      }    }
  }
  return(rawpred)
}


# de-transform and compute raw predictions and errors     
.raw <- function(pred,obs,trans,diag.ind)
  {
  if(trans=='lmat')
    {
      rawpred <- lmatvec2covmat(pred,diag.ind)
      rawobs  <- lmatvec2covmat(obs,diag.ind)
      rawerr  <- rawobs-rawpred 
      rawerr  <- cov2vec(rawerr)
    }
  if(trans=='lcov')
    {
      rawpred <- pred
      rawpred[diag.ind==1] <- exp(rawpred[diag.ind==1])
      rawobs <- obs
      rawobs[diag.ind==1]  <- exp(rawobs[diag.ind==1])
      rawerr  <- rawobs-rawpred  
      }
  if(trans=='lcpd')
  {
    rawpred <- pred
    rawpred[diag.ind==1] <- exp(rawpred[diag.ind==1])
    rawobs <- obs
    rawobs[diag.ind==1]  <- exp(rawobs[diag.ind==1])
    # transform rawpred to pd.
    rawpred <- vecreg(rawpred,diag.ind=diag.ind)
    
    # compute errors
    rawerr  <- rawobs-rawpred  
  }
        return(rawerr)
  }


# Test if the input matrix is ill conditionned or non positive definite
# If one is true, regularize as in Hautsch, Kyj, and Oomen
vecreg <- function(pred,diag.ind){
  
  nstk <- sum(diag.ind)
  cx <- covvec2covmat(pred,diag.ind)
  
  # Check ill conditionned
  ill <- kappa(cx) > 10*nstk
  
  # Getting the eigenvalues and vectors of cx
  ev   <- eigen(cx)
  eval <- ev$value
  evec <- ev$vector
  
  # Check !PD
  NPD <- !(all(eval > 0))         
           
  # Regularization         
  if(ill | NPD ){

    # The other way
    # Smallest positive eigenvalue
    lmp  <- min(eval[eval>0])
    evreg <- eval*(eval>lmp) + lmp*(eval<=lmp)
    # PD covariance matrix
    cxpd <- evec %*% diag(evreg) %*% t(evec)
    # cat('\n Smallest positive EV:',lmp,'nbr regularized eigenvalues: ',sum(eval<lmp),'\n',sep=' ')
    # cat('\n Condition number of the original matrix:',kappa(cx),' reg matrix: ',kappa(cxpd),'\n',sep=' ')
    
  }
  else {cxpd <- cx}        
  
    return(cov2vec.none(cxpd))
}


# Using the forecast object to construt hmax-step ahead forecasts
fc.roll <- function(stat.smpl,hmax=2,diag.ind){
	hfc <- list()
  
  # Looping over samples
	for(r in 1:nrow(stat.smpl)){
		#Loading the forecasts and setting some variables
		mn	 <- paste(stat.smpl[r,],collapse='.')
    
    # Is it a no-change forecast?
    nc <- (stat.smpl[r,'Model']=='NoChange')
    
    post <- FALSE
    
    # Checking if the roll file already exists, if it does not, compute the stuff 
    if(!file.exists(paste('fc_roll/roll',mn,sep='.')))
    {
      datpath <- paste('fc_data/fc',mn,sep='.')
      if(length(grep('lcpd',mn))) datpath <- gsub('lcpd','lcov',datpath)
  		fcn	 <- load(datpath)
  		fc	 <- get(fcn)
      rm(list=fcn)
      gc()
  		nfc	 <- nrow(fc$err)
  
      # Is there a post estimator in the forecast model?
  		if( (!is.null(fc$post)) & length(fc$post)) post <- TRUE
      
      
      # Getting the transformation:
  		if(length(grep('lmat',mn)))trans <- 'lmat'
  		if(length(grep('lcov',mn)))trans <- 'lcov'
  		if(length(grep('lcpd',mn)))trans <- 'lcpd'
      
  		# Building the forecast array
  		fc.storage <- array(NA,dim=c(dim(fc$pred)-c(hmax-1,0),hmax)
                    ,dimnames=list('fc'=1:(nfc-hmax+1),'equation'=dimnames(fc$pred)[[2]],'h'=1:hmax))
      
      fc.h    <- list('err'=fc.storage,'pred'=fc.storage,'obs'=fc.storage)
      fc.post <- list('err'=fc.storage,'pred'=fc.storage,'obs'=fc.storage)
      rm(fc.storage)
      gc()
      
      # Initializing the forecats
      # everything in fc is untransformed
      # getting the raw errors:
      if(hmax==1)predh1 <- fc$pred
      else predh1 <- head(fc$pred,-(hmax-1))
      
      if(hmax==1)obsh1  <- tail(fc$CRK,nfc)
      else obsh1  <- head(tail(fc$CRK,nfc),-(hmax-1)) 
      
      rawerrh1 <- array(NA,dim=dim(predh1)) 
      for(p in 1: nrow(predh1)){rawerrh1[p,]<-.raw(predh1[p,],obsh1[p,],trans=trans,diag.ind=diag.ind)}
      
      # Storig the prediction and raw errors for horizon 1.
      fc.h$pred[,,1]  <- predh1
    	fc.h$err[,,1]   <- rawerrh1
      fc.h$obs[,,1]   <- obsh1  
      

      # Looping over forecasts
  		for(i in 1:(nfc-hmax+1))
      {      
        refday <- nrow(fc$CRK) - nfc + i - 1 
        
        # No change forecast at every horizon, equal to first forecast.
        if(nc&&(hmax>1)){
          fc.h$pred[i,,2:hmax] <- matrix(fc.h$pred[i,,1],nrow=ncol(fc$CRK),ncol=hmax-1,byrow=FALSE)
        }
      
        # Storing the parameters:
        if(!nc)
          {
          # Using coeff if no post estimator
          fc.par <- as.matrix(fc$coefficients[[i]])
          lags <- (nrow(fc.par)-1)/ncol(fc$CRK)
          
          if(post) 
            {
            # Storing the post lasso parameters
            post.par <- as.matrix(fc$post[[i]])
            post.par[which(is.na(post.par))] <- 0  
            # Constructing the past observation vector
            post.rec <- 1
            for(lh in 1:lags){post.rec<-c(post.rec,fc$CRK[refday-lh+1,])}
            
            # post Forecast at horizon 1
        		post.tmp <- post.rec%*%post.par       
      			fc.post$pred[i,,1] <- post.tmp 
      			fc.post$err[i,,1] <- .raw(pred=post.tmp,fc$CRK[refday+1,],trans=trans,diag.ind=diag.ind)  
            fc.post$obs[i,,1]   <- fc$CRK[refday+1,]
            }
        }
        
      if(hmax>1){
        # Rolling over horizon
  			for(h in 2:hmax)
          {
  				fcday <- refday+h
          
  			 # No change forecasts errors
  			  if(nc){
  			    fc.h$err[i,,h] <-.raw(
              pred=fc.h$pred[i,,h]
              ,obs=fc$CRK[refday + h, ]
              ,trans=trans,diag.ind=diag.ind)
  			  }
        
      # 'model' forecast errors
        if(!nc)
          {
            
            # construct the recursive forecast vector for prediction
            fcy.rec <- 1
            for(l in 1:min(h-1,lags)){fcy.rec<-c(fcy.rec,fc.h$pred[i,,h-l])}     
            
            # if some lags are not available from previous forecasts take them from the observed data:
            if(h<=lags)for(lh in 1:(lags-h+1)){fcy.rec<-c(fcy.rec,fc$CRK[refday-lh+1,])}
            
            # Forecast
    				fcy.tmp <- fcy.rec%*%fc.par        
    				fc.h$pred[i,,h] <- fcy.tmp 
    				fc.h$err[i,,h] <- .raw(pred=fcy.tmp,fc$CRK[fcday,],trans=trans,diag.ind=diag.ind)
            
            if(post)
            {
              # Forecast with post lasso if available
              post.rec <- 1
              for(l in 1:min(h-1,lags)){post.rec<-c(post.rec,fc.post$pred[i,,h-l])}            
              # if some lags are not available from previous forecasts take them from the observed data:
              if(h<=lags)for(lh in 1:(lags-h+1)){post.rec<-c(post.rec,fc$CRK[refday-lh+1,])}
              
              # post Forecast
      				post.tmp <- post.rec%*%post.par       
      				fc.post$pred[i,,h] <- post.tmp 
      				fc.post$err[i,,h] <- .raw(pred=post.tmp,fc$CRK[fcday,],trans=trans,diag.ind=diag.ind)
            }
            
            
    			}
          fc.h$obs[i,,h] <- fc$CRK[fcday,]
          fc.post$obs[i,,h] <- fc$CRK[fcday,]
        }
      }
  		}
  		rm(fc)
  		gc()

      # adding the parameter list to the forecast object if not a no change
      if(!nc) fc.h$coefficients <- fc.par
      save(fc.h,file=paste('fc_roll/roll',mn,sep='.'))
      
      # Saving a post lasso forecast object if post exists
      if(post)
      {
        fc.post$coefficients <- post
        save(fc.post,file=paste('fc_roll/post',mn,sep='.'))
      }
    rm(fc.h)
    rm(fc.post)
    gc()
    }

	}
	
	return(TRUE)
}




roll.stattab <- function(stat.smpl,diag.ind,hsel=1){

  # Create a full model list with the post models
  full.smpl <- NULL
  
  for(m in 1:nrow(stat.smpl))
    {
    full.smpl <- rbind(full.smpl,paste('roll',paste(stat.smpl[m,],collapse='.'),sep='.'))
    if(file.exists(paste('fc_roll/post',paste(stat.smpl[m,],collapse='.'),sep='.')))
      { 
      full.smpl <- rbind(full.smpl,paste('post',paste(stat.smpl[m,],collapse='.'),sep='.'))
      }
    }
  
	# forecast dimensions
	nmod <- nrow(full.smpl)
  
  hfc <- get(load(paste('fc_roll/',full.smpl[1,],sep='')))
  
	nfc <- dim(hfc[[1]])[1]
	hmax <- dim(hfc[[1]])[3]

	# The stats:
	statn <- c('beat','RMSE','MdSE','MxSE','frob')
	nstat <- length(statn)

	# The stat array
	fctab <- array(NA,dim=c(nmod,nstat,3,length(hsel)),dimnames=list('model'=full.smpl,'stat'=statn,'diag'=c('A','D','O'),'h'=hsel))


	#Time to fill up the table:
	for(m in 1:nmod){
    
	  mn	 <- full.smpl[m,]
    hfc <- get(load(paste('fc_roll/',mn,sep='')))
    
    
	  if(length(grep('lmat',mn)))trans <- 'lmat'
	  if(length(grep('lcov',mn)))trans <- 'lcov'
	  
    # Storing the benchmark
		if(m==1)fc.bmk <- hfc
    # Loop over the element subsets
		for(d in c('A','D','O')){
      ind <- ifelse(d=='D',1,0)
      if(d=='A') ind<-c(0,1)
			if(m>1){
			fctab[m,'beat',d,] <- apply(abs(fc.bmk$err[
        ,diag.ind %in% ind
        ,hsel])>abs(hfc$err[,diag.ind%in%ind,hsel]),3,
        mean)
			}
			fctab[m,'RMSE',d,] <- apply(hfc$err[
        ,diag.ind%in%ind,hsel],3
        ,function(x)mean(sqrt(rowMeans(x^2))))
      
			fctab[m,'MdSE',d,] <- apply(hfc$err[
        ,diag.ind%in%ind,hsel],3
        ,function(x)mean(apply(abs(x),1,median)))
      
			fctab[m,'MxSE',d,] <- colMeans(apply(hfc$err[
        ,diag.ind%in%ind,hsel],c(1,3)
        ,function(x)max(abs(x))))
		}

    frobenius<-frobnorm(hfc$err[,,],diag.ind)
    
    fctab[m,'frob',,] <-frobenius[,hsel]
	}
	return(fctab)

}


frobnorm<-function(err,diag.ind)
{
  frobenius<-rbind(aaply(err,c(3),function(x)mean(sqrt(rowSums(x^2))))
                   ,aaply(err[,diag.ind==1,],c(3),function(x)mean(sqrt(rowSums(x^2))))
                   ,aaply(err[,diag.ind==0,],c(3),function(x)mean(sqrt(rowSums(x^2))))
  )
  rownames(frobenius) <- c('all','diag','off-diag')      
  
  return(frobenius)
}

