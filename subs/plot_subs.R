
# trying to get fancy facet labels with facet_wrap
facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  g <- ggplotGrob(gg.plot)
	gg <- g$grobs      
	strips <- grep("strip_t", names(gg))
	
	for(ii in seq_along(labels))  {
		modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
						   grep=TRUE, global=TRUE)
		gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[[ii]])
	}
	
	g$grobs <- gg
	class(g) = c("arrange", "ggplot",class(g)) 
	g
}

#	Plotting from a set of models:
plot.mod.agg <- function(plot.smpl,diag.names,diag.ind,dates.all)
  {
	cat('\n\n')
	cat('\nCreating plots for the models in plot.smpl:\n')
	cat(apply(plot.smpl,1,paste,collapse=', '),sep='\n')
	cat('\n\nPlot types: Parameter Stability (ps) and Model Dimension (md).\n')
  
  
	for(p in 1:nrow(plot.smpl)){
		pn <- paste(plot.smpl[p,],collapse='.')	#model names
		dates	<-dates.all [-c(1:(as.numeric(plot.smpl[p,'Est.smpl'])))]
		fc.pn <- load(paste('fc_data/fc',pn,sep='.'))
		fc.p  <- get(fc.pn)

    
		mod.strg <-  paste(plot.smpl[p,],collapse=', ')
		assign(paste('gg.md',pn,sep='.'),gg.moddim(fc.p,est.name=paste(plot.smpl[p,'Estimator'],plot.smpl[p,'Adaptive'],sep=' '),grp.ind=diag.ind,grp.names=diag.names,dates=dates,model.string=mod.strg))
		assign(paste('gg.ps',pn,sep='.'),gg.parstab(fc.p,diag.ind,diag.names,model.string=mod.strg))


		for(ptype in c('ps','md')){
			ggsave(get(paste('gg',ptype,pn,sep='.')),file=paste('plots/',gsub('\\.','-',paste(ptype,pn,sep='.')),'.pdf',sep=''))
		}
	
	}
return(TRUE)
}



#	Plotting from a set of models:
plot.mod <- function(plot.smpl,diag.names,diag.ind,dates.all){
	cat('\n\n')
	cat('\nCreating plots for the models in plot.smpl:\n')
	cat(apply(plot.smpl,1,paste,collapse=', '),sep='\n')
	cat('\n\nPlot types: Parameter Stability (ps) and Model Dimension (md).\n')
  
	for(p in 1:nrow(plot.smpl)){
		pn <- paste(plot.smpl[p,],collapse='.')	#model names
		dates	<-dates.all[-c(1:(21+as.numeric(plot.smpl[p,'Est.smpl'])))] # Forecasts dates accounting for the HAR truncation. 
		fc.pn <- load(paste('fc_data/fc',pn,sep='.'))
		fc.p  <- get(fc.pn)


		mod.strg <-  paste(plot.smpl[p,],collapse=', ')
		assign(paste('gg.md',pn,sep='.'),gg.moddim(fc.p,est.name=paste(plot.smpl[p,'Estimator'],plot.smpl[p,'Adaptive'],sep=' '),grp.ind=diag.ind,grp.names=diag.names,dates=dates,model.string=mod.strg))
		assign(paste('gg.gw',pn,sep='.'),gg.gw(fc.p,est.name=paste(plot.smpl[p,'Estimator'],plot.smpl[p,'Adaptive'],sep=' '),grp.ind=diag.ind,grp.names=diag.names,mse.scale=NULL,dates=dates))



		for(ptype in c('md','gw')){
			ggsave(get(paste('gg',ptype,pn,sep='.')),file=paste('plots/',gsub('\\.','-',paste(ptype,pn,sep='.')),'.pdf',sep=''))
			ggsave(get(paste('gg',ptype,pn,sep='.')),file=paste('plots/',gsub('\\.','-',paste(ptype,pn,sep='.')),'.pdf',sep=''))
		}
	
	}
return(TRUE)
}



plot.print <- function(plot.smpl){
	for(p in 1:nrow(plot.smpl)){
		cat('\n\\clearpage\n\\newpage\n')
		pn <- paste(plot.smpl[p,],collapse='.')	#model names
		cat('\\subsection*{Plots for model: ',pn ,'.}\n',sep='')
		for(ptype in c('ps','md')){
			pdf.n <- (paste('plots/',gsub('\\.','-',paste(ptype,pn,sep='.')),'.pdf',sep=''))
			cat('\n\\begin{figure}[h!]')
			cat('\n\\begin{center}')
			cat('\n\\includegraphics[width=\\textwidth]{',pdf.n,'}',sep='')
			cat('\n\\end{center}')
			cat('\n\\caption{',paste(ptype,pn,sep='.'),'}',sep='')
			cat('\n\\label{fig:',paste(ptype,pn,sep='.'),'}',sep='')
			cat('\n\\end{figure}\n')
			}
		}
}



# parameter stability plot, main function
gg.parstab<-function(fc.lv,grp.ind=NULL,grp.names=NULL,discrete=FALSE,model.string=NULL)
{
	
	dmmod<-matrix(0,nrow=nrow(fc.lv$coefficients[[1]])-1,ncol=ncol(fc.lv$coefficients[[1]]))
	for(fccoef in fc.lv$coefficients)dmmod<-dmmod+1*(matrix(fccoef,nrow=nrow(fccoef),ncol=ncol(fccoef))[-1,]!=0)
	par.nz <- dmmod/length(fc.lv$coefficients) 

	par.nz[which(par.nz==0)]<-NA

	if(!is.null(grp.ind))	gg.pstab <- gg.par.grp(par.nz,grp.ind,grp.names,discrete,model.string)
	if(is.null(grp.ind))	gg.pstab <- gg.par(par.nz,model.string)
	return(gg.pstab)
}

gg.par<-function(par.nz,model.string)
{
	mpar<-melt(par.nz)
	gg0 <- ggplot(mpar,aes(Var2,Var1))+ geom_tile(aes(fill = value),colour = "white")
	gg0 <- gg0 + theme_minimal()
	gg0 <- gg0 + scale_fill_gradient(low = "white", high = "midnightblue",na.value='white',limits=c(0,1),name='Variable selection frequency \n')
	gg0 <- gg0 + xlab('Equation')+ylab('Covariate')+ggtitle('Parameter Stability.')
	#gg0 <- gg0 + theme(axis.text.y=element_blank())

	return(gg0)
}

gg.par.grp<-function(par.nz,grp.ind,grp.names,discrete,model.string)		 
{
  
  grp.ind <- 1-grp.ind
	lags	<-nrow(par.nz)/ncol(par.nz)

	rownames(par.nz)<-(1:nrow(par.nz))
	tB	<-t(par.nz)
	dfB	<-data.frame(tB)
	dfB$grp	<-grp.ind
	dfB$eq.id <-1:nrow(tB)
	dfB$pl.id <-order(1:nrow(tB))[order(grp.ind)]

	bp	<-c(0.5,0.5+which((sort(grp.ind)[-length(grp.ind)]-sort(grp.ind)[-1])!=0),length(grp.ind)+0.5)  #break points for index
	glp	<-(bp[-1]-bp[-length(bp)])/2+bp[-length(bp)]

	#rownames(par.nz)<-1:nrow(par.nz)
	mB0 <- melt(dfB,id.vars=c('grp','eq.id','pl.id'))

	grp.ind2<-grp.ind
	bp2	<-bp	
	glp2	<-glp
	grp.names2<-grp.names

	if(lags>1)for(i in 2:lags){
		grp.ind2<-c(grp.ind2,max(grp.ind2)+grp.ind)
		bp2	<-c(bp2,bp+(i-1)*nrow(tB))
		glp2	<-c(glp2,glp+(i-1)*nrow(tB))
		grp.names2	<-c(grp.names2,grp.names)
	}
	mB0$grp2	<-rep(grp.ind2,each=ncol(par.nz))
	mB0$variable	<-as.double(sub('X','',mB0$variable))

	mB0$var.pl	<-rep(order((1:(lags*nrow(tB)))[order(grp.ind2)]),each=ncol(par.nz))

	if(discrete)mB0$val.dis	<- cut(mB0$value,breaks=c(0,0.2,0.8,1),right=TRUE)

	gg0 <- ggplot(mB0,aes(x=pl.id,y=var.pl))
	gg0 <- gg0 + theme_minimal()
	if(!discrete){
		gg0 <- gg0 + geom_tile(aes(fill = value),colour = "white")
		gg0 <- gg0 + scale_fill_gradient(low = "white", high = "midnightblue",na.value='white',limits=c(0,1),name='Variable selection frequency \n')
		#legend
		gg0 <- gg0 + theme(legend.position='bottom')+ggtitle(paste('Parameter stability: ',model.string,sep=''))
		gg0 <- gg0 + guides(fill=guide_colorbar(barwidth=10,barheight=0.5))
	}
	if(discrete){
		gg0 <- gg0 + geom_tile(aes(fill = val.dis),colour = "white")
		#gg0 <- gg0 + scale_fill_brewer(palette='RdPu',na.value='white',name='sel freq')
		gg0 <- gg0 + scale_fill_manual(values=c('orange','red','purple'),na.value='white',name='Variable selection frequency \n')
		#legend
		gg0 <- gg0 + theme(legend.position='bottom')+ggtitle('Parameter Stability.')
		#gg0 <- gg0 + guides(fill=guide_colorbar(barwidth=10,barheight=0.5))
	}

	#grp marks X
	gg0 <- gg0 + geom_segment(data=data.frame(grp.ind),aes(x=0,xend=length(grp.ind)+.5,y=0.5,yend=0.5))	
	gg0 <- gg0 + geom_segment(data=data.frame(bp),aes(x=bp,xend=bp,y=-1,yend=1))	

	#grp marks Y
	gg0 <- gg0 + geom_segment(data=data.frame(grp.ind2),aes(y=0,yend=length(grp.ind2)+.5,x=0.5,xend=0.5))	
	gg0 <- gg0 + geom_segment(data=data.frame(bp2),aes(y=bp2,yend=bp2,x=-1,xend=1))	

	#gg0 <- gg0 + scale_y_discrete('Covariate',breaks=NULL,labels=NULL)
	gg0 <- gg0 + scale_y_continuous('Covariate',breaks=floor(glp2),labels=grp.names2)
	gg0 <- gg0 + theme(axis.text.y = element_text(angle = 45, hjust = 1))
	gg0 <- gg0 + scale_x_continuous('Equation',breaks=floor(glp),labels=grp.names)
	gg0 <- gg0 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	gg0 <- gg0 + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

	#fixing the ratio to y:lag, x:1
	#gg0 <- gg0 + coord_fixed(ratio=lags)
	# For fixed ratio set manual axis in ggsave...

	return(gg0)
}





gg.gw<-function(l.fc.lv,est.name,grp.ind=NULL,grp.names=NULL,mse.scale=NULL,dates)
{
	#need to allow for several estimators passed as list or only one passed raw

	if(is.list(l.fc.lv[[1]])){
		grp.mse<-data.frame()
		for(i in 1:length(l.fc.lv))
		{
		fc.lv	<-l.fc.lv[[i]]
		e.name	<-est.names[i]
		ggmse<-.gg.gw.dataprep(fc.lv,e.name,grp.ind,grp.names,mse.scale,dates=dates)
		grp.mse<-rbind(grp.mse,ggmse)
		}
	}
	else grp.mse<-.gg.gw.dataprep(l.fc.lv,est.name,grp.ind,grp.names,mse.scale,dates=dates)


	m.mse	<-melt(grp.mse,id.vars=c('Estimator','Date'))
	gw<-ggplot(m.mse,aes(x=Date,y=value))+geom_line(data=m.mse,aes(colour=Estimator))+facet_wrap(~variable,ncol=1,scale='free_y')
	gw<-gw+theme_minimal()+ylab('Cumulative root squared forecast error') + theme(legend.position='bottom')
	return(gw)
}

.gg.gw.dataprep<-function(fc.lv,e.name,grp.ind=NULL,grp.names=NULL,mse.scale=NULL,dates)
	{
	#cumulating
	cum.mse <- apply(fc.lv$err^2,2,cumsum)
  #cum.mse <- sqrt(cum.mse)
	#scaling if scale provided
	if(!is.null(mse.scale)){sca.mse<-t(apply(cum.mse,1,function(x){return(x/mse.scale)}))}
	else	sca.mse<-cum.mse
	#aggregation
	grp.mse	<-NULL
	if(!is.null(grp.ind)){
		for(i in unique(grp.ind)){grp.mse<-cbind(grp.mse,rowMeans(sca.mse[,i==grp.ind]))}
		}
	else grp.mse<-sca.mse
  
  # square root
  grp.mse <- sqrt(grp.mse)
  
	#creating some names if not provided
	if(is.null(grp.names))grp.names<-paste('grp_',1:ncol(grp.mse),sep='')
	colnames(grp.mse)<-grp.names
	grp.mse<-data.frame(grp.mse)

	#adding estimator name
	grp.mse$Estimator<-e.name

	#adding date vector
	datevec<-(as.Date(dates,format='%Y-%m-%d'))
	grp.mse$Date<-datevec

	return(grp.mse)}



gg.moddim<-function(l.fc.lv,est.name,grp.ind=NULL,grp.names=NULL,mse.scale=NULL,dates=dates,model.string=NULL)
{
	grp.md<-.gg.moddim.dataprep(l.fc.lv,est.name,grp.ind,grp.names,mse.scale,dates=dates)

  #colnames(grp.md)[1:6] <- c('min','Q10','median','Q90','max','mean')
  
	m.md	<-melt(grp.md,id.vars=c('Estimator','Date'))
	#m.md	<-melt(grp.md,id.vars=c('Estimator','Date','grp'))
  
	gw<-ggplot(m.md,aes(x=Date,y=value))
	#gw<-ggplot(m.md,aes(x=Date,y=value,colour=variable) )
  gw <- gw + geom_point()
  #gw <- gw + geom_line()
  #gw <- gw + facet_wrap(~grp,ncol=1,scale='free_y')
  gw <- gw + facet_wrap(~variable,ncol=1,scale='free_y')
  gw <- gw + stat_smooth(data=m.md,se=TRUE,colour='blue',level=0.99,method='loess') 
  gw <- gw + ggtitle('Model Dimensions')
	gw <- gw + theme_minimal()+ylab(' Model size (average across equations)') + theme(legend.position='bottom')
	return(gw)


}

.gg.moddim.dataprep<-function(fc.lv,e.name,grp.ind=NULL,grp.names=NULL,mse.scale=NULL,dates)
{
	dmmod<-NULL
	for(fccoef in fc.lv$coefficients)dmmod<-rbind(dmmod,colSums((1*(matrix(fccoef,nrow=nrow(fccoef),ncol=ncol(fccoef))[-1,])!=0)))
	#aggregation + desparsifications
	grp.md	<-NULL

  if(!is.null(grp.ind)){
		for(i in unique(grp.ind)){grp.md<-cbind(grp.md,rowMeans(dmmod[,i==grp.ind]))}
		}
	else grp.md<-dmmod
	
  #creating some names if not provided
	if(is.null(grp.names))grp.names<-paste('grp_',1:ncol(grp.md),sep='')
	colnames(grp.md)<-grp.names
  grp.md<-data.frame(grp.md)
	#adding estimator name
	grp.md$Estimator<-e.name
	#adding date vector
	grp.md$Date <- dates 

	return(grp.md)

}



parsel<-function(lv.fit.full,ic)
{
  if(is.null(ic))lv.fit<-lv.fit.full
  else lv.fit<-lv.fit.full[[ic]]
  
  rownames(lv.fit$coefficients)<-c('cste',(1:nrow(lv.fit$coefficients[-1,])))
  mcoef	<-melt((lv.fit$coefficients[-1,]!=0))
  gglvsel	<-ggplot(mcoef,aes(x=Var2,y=Var1))+theme_minimal()+xlab('Equation')+ylab('Covariate')+ geom_tile(aes(fill = value))+scale_fill_manual(values=c('white','navyblue'))
  
  return(gglvsel)
}






gg.fc.dens <- function(fc.lv,diag.ind,diag.names=NULL)
{
  if(is.null(diag.names))diag.names<-c('Diag','Off-diag')
  
  m.dens	<-NULL
  for(type in c('err','pred')){
    dens	<-	fc.lv[[type]]
    colnames(dens)<-NULL
    
    m.fc.dens	<-melt(dens)
    Diag.ind	<-(m.fc.dens$Var1==m.fc.dens$Var2)
    m.fc.dens$diag	<-rev(diag.names)[as.numeric(Diag.ind)+1]
    m.tot		<-m.fc.dens
    m.tot$diag	<-'All'
    m.t.dens	<-rbind(m.fc.dens,m.tot)
    m.t.dens$type	<-type
    m.dens		<-rbind(m.dens,m.t.dens)
  }
  m.dens$type	<-	factor(m.dens$type,levels=c('err','pred'))	
  
  gg.dens<-ggplot(m.dens,aes(x=value,colour=diag))+facet_wrap(~type)+geom_line(stat='density',size=1)+theme_minimal()+scale_colour_manual(values=c('lightblue','black','steelblue'))
  
  return(gg.dens)
}



