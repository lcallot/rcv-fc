# The main function to construct tge portfolio selection and risk tables.


# the main function to generate portfolio risks.
# returns an object that can be called by the plotting and table generating functions.
portfolio <- function(stat.smpl,diag.ind,trans,wtypes)
{
  
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
  nstk <- sum(diag.ind) 

	# The risk array
	rsk <- array(NA,
               dim=c(nmod,length(wtypes),2,nfc,hmax),
               dimnames=list('model'=full.smpl,'wtype'=wtypes,'type'=c('pred','obs'),'fc'=NULL,'h'=NULL))
  # The ptf weights array
	wall  <- array(NA,
               dim=c(nmod,length(wtypes),nfc,hmax,nstk),
               dimnames=list('model'=full.smpl,'wtype'=wtypes,'fc'=NULL,'h'=NULL,'stk'=NULL))

	#Time to fill up the table:
	for(m in 1:nmod){
    
	  mn	 <- full.smpl[m,]
    hfc <- get(load(paste('fc_roll/',mn,sep='')))
    
    # Getting the transformation:
    trans <- 'none'
  	if(length(grep('lmat',mn)))trans <- 'lmat'
  	if(length(grep('lcov',mn)))trans <- 'lcov'
	  
    # detransforming the predictions and observations
	  if(trans!='none'){
	    pred <- .rawpred(hfc$pred[,,],diag.ind=diag.ind,trans=trans)
	    obs <- .rawpred(hfc$obs[,,],diag.ind=diag.ind,trans=trans)
	  }
	  else {
	    pred  <-  hfc$pred[,,]
	    obs   <-  hfc$obs[,,]
	  }
    # cleaning	  
	  rm(hfc)
	  gc()
	  
	  # Cumulated sum over horizons of predicted and observed covariances
	  cpred <- aaply(pred,c(1,2),cumsum,.parallel=TRUE)
	  cobs  <- aaply(obs,c(1,2),cumsum,.parallel=TRUE)
	  
	  # array of cumulated predicted covariances
	  cvpred <- aaply(cpred,c(1,3),covvec2covmat,diag.ind=diag.ind,.parallel=TRUE)
	  # array of cumulated observed covariances
	  cvobs <- aaply(cobs,c(1,3),covvec2covmat,diag.ind=diag.ind,.parallel=TRUE)
	  # Regularization of the predicted covariance matrices.  
	  cvpred <- aaply(cvpred,c(1,2),cvregularize,nstk=nstk,.parallel=FALSE)
	   
    # Looping over ptf weight types
    for(wn in wtypes){
      prsk <- pftrsk(cvpred,cvobs,nstk,wn)
	    rsk[m,wn,,,] <- prsk$rsk 
	    wall[m,wn,,,] <- prsk$w     
    }  
  }
  
  return(list('rsk'=rsk,'wall'=wall))
  
}

# computes the weights and the associated predicted and observed risks.
pftrsk<-function(cvpred,cvobs,nstk,wn)
{ 
  # getting ptf weights
  if(wn!='true-min-risk') wc <- aaply(cvpred,c(1,2),ptfweights,wn=wn,nstk=nstk,.parallel=TRUE)
  if(wn=='true-min-risk') wc <- aaply(cvobs,c(1,2),ptfweights,wn=wn,nstk=nstk,.parallel=TRUE)
  
  rsk <- array(NA,dim=c(2,dim(wc)[-3]),dimnames = list('type'=c('pred','obs'),'fc'=NULL,'h'=NULL))
 
  for(h in 1:dim(wc)[2])
  {
    rskop <- foreach(fc=1:dim(wc)[1],.combine='rbind') %dopar%
    {
      cat(h)
      # Cumulated
      wchf <-t(wc[fc,h,])
      cbind(wchf%*%cvobs[fc,h,,]%*%t(wchf),wchf%*%cvpred[fc,h,,]%*%t(wchf))
    }   
      rsk['obs',,h] <-rskop[,1]
      rsk['pred',,h] <- rskop[,2]
  }

return(list('rsk'=rsk,'w'=wc))
}

# Test if the inpute matrix is ill conditionned or non positive definite
# If one is true, regularize as in Hautsch, Kyj, and Oomen
cvregularize <- function(cx,nstk){
  
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
  
    return(cxpd)
}


# Calculate the portfolio weights for a given covariance matrix 
ptfweights <- function(cx,wn,nstk){
  
  nstk <- dim(cx)[1]
  
  if(wn%in%c('min-risk','true-min-risk'))icx <- try(solve(cx))
  else icx <- NULL
  
  # Computing whichever fancy weight is demanded, falling back to equal weights in case of error.
    if((!("try-error" %in% class(icx))) & (wn%in%c('min-risk','true-min-risk')))
      { # minrisk weights only if we can invert the covariance
      w <- matrix(rowSums(icx)/sum(icx),ncol=nstk,nrow=1)}
    else {w <- matrix(1/nstk,ncol=nstk,nrow=1)}
      
  if(rcond(cx)<0.0001) cat('Exposure: ',sum(abs(w)),' kappa: ',kappa(cx),' rcond ', rcond(cx), '\n',sep='')
  
  return(w)
}




cat.ptftab<- function(ptf,hsel){
  
  # weight names
  wn <- dimnames(ptf[[1]])[[2]]
  # model names
  mn <- dimnames(ptf[[1]])[[1]]

  # Constructing the avg risk
  rsk <- aaply(ptf$rsk,c(1,2,3,5),mean,na.rm=TRUE)
  
  # Constructing the exposure
  expo <- aaply(ptf$wall,c(1,2,4),function(x){mean(rowSums(abs(x),na.rm=TRUE))})
	cat('\\begin{sidewaystable}\n')
  cat('\\tiny\n')
	cat('\\begin{tabular}{',rep('l',1),rep('r',1+3*length(wn)), '}\n',sep=' ')
	cat('\\toprule\n')

  cat( ' & & \\multicolumn{',2*length(wn),'}{c}{Portfolio Risk} & \\multicolumn{',length(wn),'}{c}{Exposure} \\\\ \n')

  cat('\\cmidrule(r){',paste(3,2*length(wn)+2,sep='-'),'}')
  cat('\\cmidrule(r){',paste(2*length(wn)+3,2+3*length(wn),sep='-'),'}')
  
  cat( '& weights ')
	for(n in wn){cat('& \\multicolumn{2}{c}{',n,'}')}
	for(n in wn){cat('& \\multicolumn{1}{c}{',n,'}')}
	cat('\\\\\n')
  for(i in 1:length(wn))cat('\\cmidrule(r){',paste(2*(i-1)+3,2*i+2,sep='-'),'}  ')
  
	cat(c('Model','Horizon'),rep(c('Predicted','Observed'),length(wn)),rep('',length(wn)),sep=' & ')
	cat('\\\\\n')
	cat('\\midrule\n')

  for(m in mn){
    for(h in hsel){
      
      if(length(grep(m,pattern='\\.ar'))>0)cat('\\rowcolor{LightRed}\n')
      if(length(grep(m,pattern='cens'))>0)cat('\\rowcolor{LightBlue}\n')
      if(length(grep(m,pattern='NoChange'))>0)cat('\\rowcolor{LightGreen}\n')
      
      cat(m,h,'',sep='&')
      
      rrsk <- NULL
      for(w in wn) rrsk <- cbind(rrsk,t(rsk[m,w,,h]))
      rrsk <- round(rrsk,3)
      
      rexpo <-round(expo[m,,h],3)
      
      cat(rrsk,sep=' & ')
      cat(' & ')
      cat(rexpo,sep=' & ')
      cat(' \\\\\n ')
      
      if(h==hsel[length(hsel)])cat('\\midrule\n')
    }
  }
  
  
  
	cat('\n')
	cat('\\bottomrule\n')			
	cat('\\end{tabular}')
	cat('\\label{ptf.tab}')
	cat('\\caption{Summary statistics portfolio selection.}')
	cat('\\end{sidewaystable}')
}



