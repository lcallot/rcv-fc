
# Constructs the  stat table 
mk.stattab <- function(stat.smpl,diag.ind){
  # Constructing the storage matrix
  stat.names <- c('Nbr.param','Nbr.non.0' ,'beat.bnchmk','Rsquared','Autocorr','Normality', # Split in O and D
                  'MSE','Median','Max.err') # As before.
  
  OD.n <- NULL
  for(n in stat.names){OD.n <- c(OD.n,paste(c('D','O'),n,sep='.'))}
  stat.tab <- matrix(NA,nrow=nrow(stat.smpl),ncol=(ncol(stat.smpl)+3+2*length(stat.names)),
                     dimnames=c(list('Models'=NULL,'stat'=c(colnames(stat.smpl),'Time','Nbr.fc','Nbr.eq',OD.n))))
  
  #Constructing the plot
  for(s in 1:nrow(stat.smpl)){
    sn <- paste(stat.smpl[s,],collapse='.')	#model names
    
    stat.tab[s,colnames(stat.smpl)] <- stat.smpl[s,] #begining to fill up the matrix
    fc.sn	 <- load(paste('fc_data/fc',sn,sep='.'))
    fc.s	 <- get(fc.sn)
    if(s==1)fc.bchmk <- fc.s
    transfo	 <- 'unknown'
    if(length(grep('lmat',sn)))transfo <- 'lmat'
    if(length(grep('lcov',sn)))transfo <- 'lcov'
    ind	 <- 'unknown'
    if(length(grep('dj',sn)))ind<- 'dj'
    if(length(grep('test',sn)))ind <- 'test'
    fc.sv	 <- fc.statvec(fc.s,diag.ind,xtr.ind=NULL,trans=transfo,ind=ind)
    n0 <- NULL
    if(stat.smpl[s,'Model']!='NoChange')for(i in 1:length(fc.s$coefficients)){n0 <- rbind(n0,colSums(fc.s$coefficients[[i]]!=0))}
    
    nstock <- sum(diag.ind[[ind]]==1)
    
    for(d in c('D','O')){
      di <- ifelse(d=='D',1,2)
      # spec tests
      #for(st in c('Autocorr','Normality','Rsquared')){if(stat.smpl[s,'Model']!='NoChange')stat.tab[s,paste(d,st,sep='.')] <- mean(fc.s$spectest[st,diag.ind[[ind]]==di,])	}
      # Number of parameters, and non zero
      if(stat.smpl[s,'Model']!='NoChange'){
        if(stat.smpl[s,'Restrictions']=='none')stat.tab[s,paste(d,'Nbr.param',sep='.')] <- dim(fc.s$coefficients[[1]])[1]	
        if(stat.smpl[s,'Restrictions']=='rest')stat.tab[s,paste(d,'Nbr.param',sep='.')] <- 1+((dim(fc.s$coefficients[[1]])[1]-1)/dim(fc.s$coefficients[[1]])[2])*ifelse(d=='D',nstock+(nstock-1),nstock+(nstock-1)+(nstock-2))	
      }
      if(stat.smpl[s,'Model']!='NoChange')stat.tab[s,paste(d,'Nbr.non.0',sep='.')] <- mean(n0[,diag.ind[[ind]]==di])
      if(s>1)stat.tab[s,paste(d,'beat.bnchmk',sep='.')] <- mean((abs(fc.s$err)<abs(fc.bchmk$err))[,diag.ind[[ind]]==di])
    }
    
    #Getting the stats from the statvec function
    n.sv <- c('Nbr.fc','Nbr.eq')
    for(n in c('MSE','Median','Max.err')){n.sv <- c(n.sv,paste(c('D','O'),n,sep='.'))}
    
    names(fc.sv) <- n.sv
    for(n in names(fc.sv))stat.tab[s,n] <- fc.sv[n]
    
    #Saving the time
    
    tn <- load(paste('time/time',sn,sep='.'))
    stime <- (summary(get(tn))[[3]]*min(fc.sv[1],ncores))/(60*fc.sv[1])
    stat.tab[s,'Time'] <- stime
    rm(fc.s)
    gc()
  }
  return(stat.tab)
}


# Prints the stat table
cat.stattab <- function(stat.tab,stat.smpl){
  stat.n <- c('Par per Eq.','Non-0 per Eq.' ,'Beat Bchmk.','R-squared','Autocorr.','Normality.', 'Mean SFE','Med SFE','Max AFE')
  
  cat('\\begin{sidewaystable}\n')
  cat('\\tiny\n')
  cat('\\begin{tabular}{',rep('l',1),rep('r',3+2*length(stat.n)), '}\n',sep=' ')
  cat('\\toprule\n')
  
  #column label 
  
  smpl.n <- c('Sample','Time','$\\#$fc.','$\\#$Eq.')
  
  cat(paste(smpl.n,collapse=' &'))
  for(n in stat.n){cat('& \\multicolumn{2}{c}{',n,'}')}
  cat('\\\\\n')
  cat(c(rep('',1),'$\\left(\\frac{min}{fc/cpu}\\right)$',rep('',2),rep(c('D','O'),length(stat.n))),sep=' & ')
  cat('\\\\\n')
  cat('\\midrule\n')
  
  rest.ind <- which(stat.smpl[,'Restrictions']=='rest')
  lmat.ind <- grep('lmat',stat.smpl[,'Data'])
  ncens.ind <- grep('none',stat.smpl[,'Data'])
  
  
  #Some cols of the table have to be transformed back to numeric.
  stat.tab <- data.frame(stat.tab)
  #stat.tab[,-c(1:(col(stat.smpl)+3))] <- as.numeric(stat.tab[,-c(1:(col(stat.smpl)+3))])
  for(i in 1:nrow(stat.tab)){
    stat.vec <- stat.tab[i,]
    
    str.ind <- c('s',rep('f',3+2*length(stat.n)))
    rnd.ind <- c(0,2,0,0,rep(2,2*length(stat.n)))
    
    if(i%in%rest.ind)cat('\\rowcolor{LightRed}\n')
    if(i%in%lmat.ind)cat('\\rowcolor{LightGreen}\n')
    if(i==1)cat('\\rowcolor{LightBlue}\n')
    if(i%in%ncens.ind)cat('\\rowstyle{\\bfseries}\n')
    
    cat(paste(as.matrix(stat.vec)[1:ncol(stat.smpl)],collapse='.'))
    for(r in 2:length(str.ind)){
      cat(' & ')
      if(!is.na(stat.tab[[r+ncol(stat.smpl)-1]][i])){
        if(str.ind[r]=='s') cat(as.character(stat.tab[[r+ncol(stat.smpl)-1]][i]),sep='')
        if(str.ind[r]=='f') cat(round(as.numeric(as.character(stat.tab[[r+ncol(stat.smpl)-1]][i])),rnd.ind[r]),sep='')
      }
    }
    cat('\\\\\n')
    if(i%in%mdru.ind)cat('\\midrule\n')
  }
  cat('\n')
  cat('\\bottomrule\n')			
  cat('\\end{tabular}')
  cat('\\label{stat.tab}')
  cat('\\caption{Summary statistics, specification tests, and forecast statistics for every model.}')
  cat('\\end{sidewaystable}')
  
  return(TRUE)
}



#Takes the matrix of residuals, apply some specification tests (autocorr, heteroskedasticity, normality, ...)
spectest<-function(mres,my,lags=min(floor(3*fitdf),nrow(mres)/2),fitdf=1)
{
  sptest <- matrix(NA,nrow=3,ncol=ncol(my))
  
  
  BP <- apply(mres,2,function(x){Box.test(x,lag=lags,type="L",fitdf=floor(fitdf))->bt;return( bt$p.val)})
  SW <- apply(mres,2,function(x){shapiro.test(x)$p.value})
  R2 <- 1-((nrow(mres)*apply(mres,2,var,na.rm=TRUE))/(nrow(my)*apply(my,2,var)))
  
  sptest[1,]	 <- BP
  sptest[2,]	 <- SW
  sptest[3,]	 <- R2
  
  rownames(sptest) <- c('Autocorr','Normality','Rsquared')
  
  return(sptest)
}




catstattab <- function(astat,coltype,statn,hsel)
{
  cat('\\begin{sidewaystable}\n')
  cat('\\tiny\n')
  cat('\\begin{tabular}{',coltype, '}\n',sep=' ')
  cat('\\toprule\n')
  
  cat(' & ')
  #column label 
  for(n in c(stat.n)){cat('& \\multicolumn{3}{c}{',n,'}')}
  cat('\\\\\n')
  cat('Model&h&')
  cat(paste(rep('A & D & O',length(stat.n)),collapse='&'))
  cat('\\\\\n')
  cat('\\midrule\n')
  
  for(m in 1:nrow(astat)){
    if(length(grep(astat[m,1],pattern='post'))>0)cat('\\rowcolor{LightBlue}\n')
    if(length(grep(astat[m,1],pattern='lcov'))>0)cat('\\rowcolor{LightRed}\n')
    if(length(grep(astat[m,1],pattern='lasso'))>0)cat('\\rowcolor{LightGreen}\n')
    
    cat(gsub(pattern='\\.',x=astat[m,1],replacement=' '),astat[m,2],sep='&')
    cat('&')
    cat(round(as.numeric(astat[m,-c(1,2)]),2),sep='&')
    cat('\\\\\n')
    if(((m)%%length(hsel))==0)cat('\\midrule\n')
  }
  
  cat('\n')
  cat('\\bottomrule\n')  		
  cat('\\end{tabular}')
  cat('\\caption{Summary statistics, h-step ahead recursive forecasts, all statistics averaged across forecast iterations.}')
  cat('\\end{sidewaystable}')
}


selmat <- function(pm,djcat,nstock=30){
  # Stock categories
  # stock indices
  dj <- cbind(which(diag.ind==1),djcat)
  
  # constructing that category vectors
  catm1 <- matrix(djcat[,3],nstock,nstock)
  catv1 <- matrix(catm1[upper.tri(catm1,diag=TRUE)],ncol=1)
  catm2 <- matrix(djcat[,3],nstock,nstock,byrow = TRUE)
  catv2 <- matrix(catm2[upper.tri(catm2,diag=TRUE)],ncol=1)
    
  # Count diag and offdiaf per cat
  catcd  <- table(catv1[diag.ind==1])
  catco1 <- table(catv1[diag.ind==0])
  catco2 <- table(catv2[diag.ind==0]) 
  catn <- names(catcd)
  
  if(length(catco1)<length(catco2)){catco1 <- c(catco1,0);names(catco1) <- names(catco2)}
  # Storage matrices
  ncat <- length(unique(djcat[,3]))
  oeq <- deq <- matrix(0,2*ncat,ncat,dimnames = list(NULL,catn))  
  
  # parameter selection frequency
  sel <- matrix(0,length(diag.ind),length(diag.ind))
  for(fc in pm){ sel <- sel + 1*(as.matrix(fc[-1,])!=0) }
  
  sel <- sel/length(pm)
  
  # Loop over rows
  for(r in 1:ncol(sel)){
    rsel <- sel[,r]
    dv   <- ov <- NULL
    
    
    for(ct in 1:ncat){
      #diag cound
      dc <- sum(rsel[ (catv1==catn[ct]) & (diag.ind==1)])
      #dc <- dc + sum(rsel[ (catv2==catn[ct]) & (diag.ind==1)])
      dv <- c(dv, dc/catcd[ct])
      #off-diag cound
      oc <- sum(rsel[ (catv1==catn[ct]) & (diag.ind==0)])
      oc <- oc + sum(rsel[ (catv2==catn[ct]) & (diag.ind==0)])
      ov <- c(ov, oc/(catco1+catco2)[ct])
    }
    # Aggregation by categories
    if(diag.ind[r]==1) {
      rcat <- catv1[r]
      deq[1:ncat,rcat] <- deq[1:ncat,rcat] + dv/catcd[rcat]
      deq[ncat + (1:ncat),rcat] <- deq[ncat + (1:ncat),rcat] + ov/catcd[rcat]
    }
    if(diag.ind[r]==0) {
      rcat1 <- catv1[r]
      rcat2 <- catv2[r]
      oeq[1:ncat,rcat1] <- oeq[1:ncat,rcat1] + dv/catco1[rcat1]
      oeq[1:ncat,rcat2] <- oeq[1:ncat,rcat2] + dv/catco2[rcat2]
      oeq[ncat + (1:ncat),rcat1] <- oeq[ncat + (1:ncat),rcat1] + ov/catco1[rcat1]
      oeq[ncat + (1:ncat),rcat2] <- oeq[ncat + (1:ncat),rcat2] + ov/catco2[rcat2]
    } 
  }
  return(list('deq'=deq,'oeq'=oeq))
}

