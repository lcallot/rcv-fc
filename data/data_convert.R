# Convert the R data files to csv with dates and 


# get list of the R data files. 
CRKnames <- list.files('.',pattern = 'CRK')

# get the dates 
dates <- read.table('dates')$V1 # raw dates
day   <- tail(dates,-21) # legacy oddity 
week  <- unique(strftime(dates,format='%Y-%W')) # week dates
month <- unique(strftime(dates,format='%Y-%m')) # month dates
 

# load R data file and save as csv 
for (fn in CRKnames){
  # loading the data
  lfn <- load(fn)
  rcv <- get(lfn)
  if(length(rcv)==3) rcv <- rcv[[1]]
  
  # Assigning dates
  wk <- grep('W',fn)
  mt <- grep('M',fn)
  if(length(wk)) rownames(rcv) <- week
  if(length(mt)) rownames(rcv) <- month
  if(!(length(wk)|length(mt))) rownames(rcv) <- day
  
  # saving  
  write.csv(rcv,paste0(fn,'.csv'))
}

