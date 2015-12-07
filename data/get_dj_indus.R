
dj.ind <- read.table('dj-ind')
sp.cat <- read.csv('sp_indus.csv')

dj.cat <- NULL

for(dj in 1:nrow(dj.ind)){
	spi <- which(sp.cat[,1]%in%paste(dj.ind[dj,1],' ',sep=''))
	dj.cat <- rbind(dj.cat,data.frame('rind'=spi,(sp.cat[spi,])))
}

write.table(file='dj-cat',dj.cat)

