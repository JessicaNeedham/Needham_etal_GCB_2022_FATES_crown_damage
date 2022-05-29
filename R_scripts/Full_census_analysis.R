##########################
### Load full data set
# Get size and mortality 
################################################
################################################
rm(list = ls())
# LIBRARIES #
library(stringr)

##############################################################################
# User defined variables      
site.name <- 'bci'  # CTFS site name used in r table 
censuses <- c(7,8)      # which censuses to use - just use the most recent
n.census <- length(censuses) # number of census to include in analysis

##############################################################################
set.seed(1)
# load in all censuses
for (i in censuses){ 
  load(paste("Data/bci.tree",i,".rdata", sep=""))
}

# put them in a list
cns <- vector('list', n.census)

for(i in 1:n.census){
  cns[[i]] = eval(as.name(paste( site.name, ".tree",censuses[i], sep="")))
}

# rearrange columns and only keep tags, species, quadrats, coordinates, 
# height of measure, status, date and exact date
cnsdata = list() 

for (i in 1:n.census){
  cnsdata[[i]] <- cns[[i]][,c('treeID', 
                              "sp","dbh", "status", 'ExactDate')]
}

all <- cnsdata[[1]]
all$dbh <- all$dbh/10  # make it cm

# get rid of dead stems
all <- all[-which(all$status == 'D'), ]
all <- all[!is.na(all$dbh), ]

# divide into bins 
max(all$dbh,na.rm = TRUE)
c.x <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 
         80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 175, 200,250, 300)
c.dbh.class <- findInterval(all$dbh, vec = c.x, all.inside = TRUE)
tab.c.dbh.class <- table(c.dbh.class)
c.y <- rep(0, length(c.x)-1)
c.y[as.numeric(names(tab.c.dbh.class))] <- tab.c.dbh.class 

# Divide by 50 and by bin width
c.y <- c.y/50  # to get per ha
bins <- diff(c.x)
length(c.y)
length(bins)
c.y <- c.y/bins

write.csv(c.y, file = '../../Damage_processing/size_dist.csv', quote = FALSE)

