

# Applies the desired transformation (NULL,lcov,lmat) on a CRK data set. returns the transformed set, nothing fancy. 
transCRK <- function(CRK,transfo='none'){

	
	CRK.trans <- as.matrix(t(apply(CRK,3,get(paste('cov2vec',transfo,sep='.')))))

	return(CRK.trans)
}



# A function to aggregate and save a set of transformed subsets of the raw data in vectorized form. Can get heavy if the number of stock is large. 
agg.CRK <- function(data.smpl,bad.scale,bad.share,dat.path){
	#	Loading the file containing the stock indices and names:
	for(ind in unique(data.smpl[,'Stocks'])){
		assign(paste(ind,'.ind',sep=''),read.table(paste(ind,'-ind',sep=''))[,2])
		assign(paste(ind,'.names',sep=''),read.table(paste(ind,'-ind',sep=''))[,1])
		assign(paste(ind,'.nstk',sep=''),length(get(paste(ind,'ind',sep='.'))))
	#	Load a subset of the raw data	
		CRK.tmp <- h5read(dat.path,name='CRK_corrected'
                   ,index=list(get(paste(ind,'ind',sep='.'))
                               ,get(paste(ind,'ind',sep='.')),NULL))
		       
  
		bad.nbr		<-floor(bad.share*prod(dim(CRK.tmp)[1]
                                    ,dim(CRK.tmp)[2]+1)/2)
  
  	#	Performing the censoring on the loaded subset if required.
  	if(sum(data.smpl[,'Censoring']=='cens')!=0)
  	{
  	  assign(paste('baddays',ind,sep='.'),bad.days(CRK.tmp,bad.scale=bad.scale,bad.nbr=bad.nbr) )
  	  CRK.tmp <- censor.data(CRK.tmp,get(paste('baddays',ind,sep='.')),avg.wd=5)
  	  cat('\n\t\t.....................................................\n')
  	  cat('\t\t		 Data Censoring: ')
  	  cat('\n\t ', length(get(paste('baddays.',ind,sep=''))),' days flagged for censoring. \n',sep='')
  	  cat('Criterion: ',bad.nbr,' variables (',round(100*bad.share,0),'%) with values more than\n ', bad.scale,' std.dev. away from the sample avg. of the variable on a day.\n',sep='')
  	  cat('Days excluded: ',as.character(read.table('dates')[get(paste('baddays.',ind,sep='')),]),sep='\n\t\t')
	}
	
  
	#	Get the dates	
		dates <- read.table('dates')$V1
		week  <- strftime(dates,format='%Y-%W')
		month <- strftime(dates,format='%Y-%m')

	#		
		
		CRK.W <- aaply(CRK.tmp,.margins=c(1,2),.fun=by,week,mean)
		CRK.M <- aaply(CRK.tmp,.margins=c(1,2),.fun=by,month,mean)

		assign(paste('CRK',ind,'W',sep='.'),CRK.W)
		assign(paste('CRK',ind,'M',sep='.'),CRK.M)

		cat('...................................................................\n')
		cat(' Extracting data for stock set:',ind,'. Contains ',get(paste(ind,'.nstk',sep='')),'\n stocks.\n')
		cat('Stock Names and Index:')
		for(i in 1:length(get(paste(ind,'names',sep='.')))){
			cat('\n			',paste(as.character(read.table(paste(ind,'-ind',sep=''))[i,1]),read.table(paste(ind,'-ind',sep=''))[i,2],sep=' '))
		}
		cat('\n Aggregating by weeks and month. \n')
	}
	#	Performing required transfo for each dataset in data.smpl. aggregates the data , format the names, and save the whole mess. 
	cat('\n')
	cat('\n')
	cat('\n\t\t.....................................................\n')
	cat('\t\t		 Data Saving: \n')
	cat('\n Total number of specifications: ',nrow(data.smpl),'\n',sep='')
	for(r in 1:nrow(data.smpl)){
		sp		 <- data.smpl[r,]
		# Transforming
		CRK.tr.W	<- transCRK(get(paste('CRK',sp[1],'W',sep='.')),sp[3]) 
		assign(paste('CRK',sp[1],sp[2],sp[3],'W',sep='.'),CRK.tr.W)
		CRK.tr.M	<- transCRK(get(paste('CRK',sp[1],'M',sep='.')),sp[3]) 
		assign(paste('CRK',sp[1],sp[2],sp[3],'M',sep='.'),CRK.tr.M)

		
		
		# Saving
		save(list=paste('CRK',sp[1],sp[2],sp[3],'W',sep='.'),file=paste('CRK',sp[1],sp[2],sp[3],'W',sep='-'))
		save(list=paste('CRK',sp[1],sp[2],sp[3],'M',sep='.'),file=paste('CRK',sp[1],sp[2],sp[3],'M',sep='-'))

		cat('\n Stock Index: ',sp[1],' censoring: ',sp[2],', Transformation: ',sp[3],'.',sep='') 
		cat('\n File names: ',paste('CRK',sp[1],sp[2],sp[3],sep='.'),'\n',sep='') 
	}

	cat('\n\nAnd we\'re done\n')

}


# A function to save a set of transformed subsets of the raw data in vectorized form. Can get heavy if the number of stock is large. 
save.CRK <- function(data.smpl,bad.scale,bad.share,dat.path){
	#	Loading the file containing the stock indices and names:
	for(ind in unique(data.smpl[,'Stocks'])){
		assign(paste(ind,'.ind',sep=''),read.table(paste(ind,'-ind',sep=''))[,2])
		assign(paste(ind,'.names',sep=''),read.table(paste(ind,'-ind',sep=''))[,1])
		assign(paste(ind,'.nstk',sep=''),length(get(paste(ind,'ind',sep='.'))))
		assign(paste('CRK',ind,'none',sep='.'),h5read(dat.path,name='CRK_corrected',index=list(get(paste(ind,'ind',sep='.')),get(paste(ind,'ind',sep='.')),NULL)))
		bad.nbr		<-floor(bad.share*prod(dim(get(paste('CRK',ind,'none',sep='.')))[1],dim(get(paste('CRK',ind,'none',sep='.')))[2]+1)/2)
		       

		cat('...................................................................\n')
		cat(' Extracting data for stock set:',ind,'. Contains ',get(paste(ind,'.nstk',sep='')),'\n stocks.\n')
		cat('Stock Names and Index:')
		for(i in 1:length(get(paste(ind,'names',sep='.')))){
			cat('\n			',paste(as.character(read.table(paste(ind,'-ind',sep=''))[i,1]),read.table(paste(ind,'-ind',sep=''))[i,2],sep=' '))
		}

	#	Performing the censoring on the loaded subset if required.
		if(sum(data.smpl[,'Censoring']=='cens')!=0)
			{
			assign(paste('baddays',ind,sep='.'),bad.days(get(paste('CRK',ind,'none',sep='.')),bad.scale=bad.scale,bad.nbr=bad.nbr) )
			assign(paste('CRK',ind,'cens',sep='.'),censor.data(get(paste('CRK',ind,'none',sep='.')),get(paste('baddays',ind,sep='.')),avg.wd=5))
			cat('\n\t\t.....................................................\n')
			cat('\t\t		 Data Censoring: ')
			cat('\n\t ', length(get(paste('baddays.',ind,sep=''))),' days flagged for censoring. \n',sep='')
			cat('Criterion: ',bad.nbr,' variables (',round(100*bad.share,0),'%) with values more than\n ', bad.scale,' std.dev. away from the sample avg. of the variable on a day.\n',sep='')
			cat('Days excluded: ',as.character(read.table('dates')[get(paste('baddays.',ind,sep='')),]),sep='\n\t\t')
			}
	}
	#	Performing required transfo for each dataset in data.smpl. Creates the HAR, format the names, and save the whole mess. 
	cat('\n')
	cat('\n')
	cat('\n\t\t.....................................................\n')
	cat('\t\t		 Data Saving: \n')
	cat('\n Total number of specifications: ',nrow(data.smpl),'\n',sep='')
	for(r in 1:nrow(data.smpl)){
		sp		 <- data.smpl[r,]
		CRK.trans	 <- transCRK(get(paste('CRK',sp[1],sp[2],sep='.')),sp[3]) 
		assign(paste('CRK',sp[1],sp[2],sp[3],sep='.'),CRK.trans)
		

		# creating weekly and monthly averages
		assign(paste('CRK.5',sp[1],sp[2],sp[3],sep='.'),apply(get(paste('CRK',sp[1],sp[2],sp[3],sep='.')),2,wd.avg,5))
		assign(paste('CRK.22',sp[1],sp[2],sp[3],sep='.'),apply(get(paste('CRK',sp[1],sp[2],sp[3],sep='.')),2,wd.avg,22))

		# trimming and listing
		assign(paste('CRK.har',sp[1],sp[2],sp[3],'d',sep='.'),get(paste('CRK',sp[1],sp[2],sp[3],sep='.'))[-c(1:21),])
		assign(paste('CRK.har',sp[1],sp[2],sp[3],'w',sep='.'),get(paste('CRK.5',sp[1],sp[2],sp[3],sep='.'))[-c(1:21),])
		assign(paste('CRK.har',sp[1],sp[2],sp[3],'m',sep='.'),get(paste('CRK.22',sp[1],sp[2],sp[3],sep='.'))[-c(1:21),])
		# Saving
		x <- list()
		for(p in c('d','w','m')){x[[p]] <- get(paste('CRK.har',sp[1],sp[2],sp[3],p,sep='.'))}
		assign(paste('CRK.har',sp[1],sp[2],sp[3],sep='.'),x)
		save(list=paste('CRK.har',sp[1],sp[2],sp[3],sep='.'),file=paste('CRK-har',sp[1],sp[2],sp[3],sep='-'))

		cat('\n Stock Index: ',sp[1],', Censoring: ',sp[2],', Transformation: ',sp[3],'.',sep='') 
		cat('\n File names: ',paste('CRK',sp[1],sp[2],sp[3],sep='.'),'\n',sep='') 
		cat('\n		    ',paste('CRK.har',sp[1],sp[2],sp[3],sep='.'),'\n',sep='') 
	}

	cat('\n\nAnd we\'re done\n')

}






#Computes the number of series for which an observation is larger (in absolute value) than
bad.days	<-function(CRK,bad.scale=2,bad.nbr=100){

	CRK.var	<-sqrt((apply(CRK,c(1:2),var,na.rm=TRUE)))
	CRK.avg	<- (apply(CRK,c(1:2),mean,na.rm=TRUE))

	CRK.big <- ((CRK>(array(CRK.avg,dim=dim(CRK))+bad.scale*array(CRK.var,dim=dim(CRK))))|(CRK<array(CRK.avg,dim=dim(CRK))-bad.scale*array(CRK.var,dim=dim(CRK))))

	CRK.badday <- apply(CRK.big,3,sum)
	CRK.badday <- which(CRK.badday>min(dim(CRK)[3],bad.nbr))
	return(CRK.badday)
	}


# Covariance to transformed vector functions
cov2vec.none<-function(x){x[upper.tri(x,diag=TRUE)]}
cov2vec.lmat<-function(x){x <- logm(x);x[upper.tri(x,diag=TRUE)]}
cov2vec.lcov<-function(x){diag(x)<-log(diag(x));return(x[upper.tri(x,diag=TRUE)])}




censor.data <- function(CRK,baddays,avg.wd=5){
	
	avg.block  	<- get.avgblock(baddays)
	CRK.censored	<-CRK

	#Constructing the matrix of index for the variables used to construct averages. 
	#	The averages are identical within blocks.
	for(b in 1:nrow(avg.block)){
		block	<-avg.block[b,]
		block[which(block>dim(CRK[3]))] <- NA
		cens.day <- apply(CRK[,,block],c(1,2),mean,na.rm=TRUE)
		CRK.censored[,,baddays[b]]	<- cens.day 
	}
	return(CRK.censored)
}


#Computes the vector of indicies to use for each average..
get.avgblock <- function(baddays){	
	bad.block	<- NULL
	for(b in baddays){
		bad.left	 <- NULL
		count		 <- 1
		while(length(bad.left)<5){x<-ifelse((b-count) %in% baddays,NA,b-count);if(!is.na(x)){bad.left <- c(bad.left,x)};count<-count+1;} 
		bad.left[which(bad.left<1)] <- NA
		bad.right	 <- NULL
		count		 <- 1
		while(length(bad.right)<5){x<-ifelse((b+count) %in% baddays,NA,b+count);if(!is.na(x)){bad.right <- c(bad.right,x)};count<-count+1;} 
		bad.block <- rbind(bad.block,c(bad.left,bad.right))
		}
	return(bad.block)
	}

		
# Extracts the forecasts.
# Save them as csv files.
fc.xtract <- function(stat.smpl,h,diag.ind,dates.all)
{
  
  # Create a full model list with the post models
  full.smpl <- NULL
  
  # Checking if a post lasso has also been estimated
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

  #Looping over the models (incl post estimators):
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
    
  # obs and pred dim  
  # 1 nbr fc
  # 2 series
  # 3 horizons
    write.csv(pred[,,h],file=paste('fc_marcelo/forecast_h',h,'_',mn,'.csv',sep=''),row.names=as.character(dates.all))
    write.csv(obs[,,h],file=paste('fc_marcelo/observed_h',h,'_',mn,'.csv',sep=''),row.names=as.character(dates.all))
    
  }  
}




# Extracts the paramters.
fc.xtpar <- function(stat.smpl,dates.all)
{  
  # storage list
  parmat <- list()
  # Create a full model list with the post models
  full.smpl <- NULL
  # formatting the names
  for(m in 1:nrow(stat.smpl)) { full.smpl <- rbind(full.smpl,paste('fc',paste(stat.smpl[m,],collapse='.'),sep='.'))}
  # nbr forecast models
  nmod <- nrow(full.smpl)
  #Looping over the models:
  for(m in 1:nmod){
    mn	 <- full.smpl[m,]
    hfc <- get(load(paste('fc_data/',mn,sep='')))
    parmat[[mn]] <- hfc$coefficients
    names(parmat[[mn]]) <- dates.all
  }  
  
  return(parmat)
}




# Extracts the forecasts.
# Save them as csv files.
fc.xtlambda <- function(stat.smpl,dates.all)
{  
  # storage list
  lambda <- list()
  
  # Create a full model list with the post models
  full.smpl <- NULL
  
  # formatting the names
  for(m in 1:nrow(stat.smpl)) { full.smpl <- rbind(full.smpl,paste('fc',paste(stat.smpl[m,],collapse='.'),sep='.'))}
  
  # nbr forecast models
  nmod <- nrow(full.smpl)
  
  #Looping over the models:
  for(m in 1:nmod){
    mn   <- full.smpl[m,]
    hfc <- get(load(paste('fc_data/',mn,sep='')))
        
    lambda[[mn]] <- hfc$lambda
    names(lambda[[mn]]) <- dates.all
  }   
  
  return(lambda)
}



# Extracts the forecasts.
# Save them as csv files.
fc.xterr <- function(stat.smpl,diag.ind,dates.all)
{
  # Create a full model list with the post models
  full.smpl <- NULL
  
  errlst <- list()
  
  # Checking if a post lasso has also been estimated
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
  
  #Looping over the models (incl post estimators):
  for(m in 1:nmod){
    
    
    mn	 <- full.smpl[m,]
    hfc <- get(load(paste('fc_roll/',mn,sep='')))
    
    # Getting the transformation:
    trans <- 'none'
    if(length(grep('lmat',mn)))trans <- 'lmat'
    if(length(grep('lcov',mn)))trans <- 'lcov'
    
    # detransforming the predictions and observations
    if(trans!='none'){err <- .rawpred(hfc$obs[,,],diag.ind=diag.ind,trans=trans) - .rawpred(hfc$pred[,,],diag.ind=diag.ind,trans=trans) }
    else {err  <-   hfc$obs[,,] - hfc$pred[,,]}
    # cleaning	  
    rm(hfc)
    gc()
    
    errlst[[mn]] <- err
  }  
 
  return(errlst)
}

