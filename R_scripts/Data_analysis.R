#########################################################
### Analyses using mortality and damage surveys ###
#########################################################
##############################################################################
rm(list = ls())

## LIBRARIES ##
library(stringr)

## FUNCTIONS ##
source(file = 'Source/Allometry_functions.R')

get.mort.rate.weights <- function(cens){
  N1 <- sum(cens[ ,'weights'])
  N2 <- sum(cens[which(cens[ ,'status'] == 'A'), 'weights'])
  time <- mean(cens[ ,'time'], na.rm = TRUE)
  mort <- (log(N1)-log(N2))/1
  return(mort)
}

## SET UP ##
n.census <- 4
site <- 'bci'
n.damage <- 5

load(file = 'Output/df_single_stem.RData')  # only single stem individuals
load(file = 'Output/df.RData')  # all stems
load(file = 'Output/df_main.RData')  # only the main stem of all individuals 

#####################################################################
# get just relevant info from each census
cl.st <- vector('list', n.census-1)
for(i in 1:(n.census-1)){
  cl.st[[i]] <- df.ss[ ,c(sprintf('weights.ind.%d', i), sprintf('damage.class.%d',i), sprintf('status.%d', (i+1)), 
                          sprintf('time.%d', i))]
}
cl.st <- lapply(cl.st, function(x){colnames(x) <- c('weights', 'damage.class', 'status', 'time'); return(x)})

cl.st <- do.call(rbind, cl.st)  # bind all censuses together
cl.st <- cl.st[-which(is.na(cl.st[ ,'damage.class'])), ]  # remove those that don't have a damage class in t
cl.st <- cl.st[-which(cl.st[ ,'status'] == '?'), ] # remove those with unknown status at t + 1
spl <- split(cl.st, cl.st[ ,'damage.class'])

mort.by.damage.wt <- unlist(lapply(spl, get.mort.rate.weights))[1:n.damage]

write.csv(mort.by.damage.wt, file = '../../Damage_processing/r_mort_by_damage.csv', quote = FALSE)

# Separate by census
mort <- matrix(NA, n.census-1, n.damage)
for(i in 1:(n.census-1)){
  cl.st <- df.ss[ ,c(sprintf('weights.ind.%d', i), sprintf('damage.class.%d',i), sprintf('status.%d', (i+1)), 
                     sprintf('time.%d', i))]
  colnames(cl.st) <- c('weights', 'damage.class', 'status', 'time')
  cl.st <- cl.st[-which(is.na(cl.st[ ,'damage.class'])), ] 
  cl.st <- cl.st[which(cl.st[ ,'status'] != '?'), ] # remove those with unknown status at t + 1
  spl <- split(cl.st, cl.st[ ,'damage.class'])
  mort[i, ] <- unlist(lapply(spl, get.mort.rate.weights))[1:n.damage]
}
mort
#colnames(mort) <- sapply(seq(5), function(x) paste0('Damage.', x))
#rownames(mort) <- sapply(seq(3), function(x) paste0('Census.', x))

write.csv(mort, file = '../../Damage_processing/r_mort_by_damage_by_cens.csv', quote = FALSE)

######################################################################
### And again but by canopy layer
# get just relevant info from each census
cl.st <- vector('list', n.census-1)
for(i in 1:(n.census-1)){
  cl.st[[i]] <- df.ss[ ,c(sprintf('weights.ind.%d', i), sprintf('damage.class.%d',i), sprintf('status.%d', (i+1)), 
                          sprintf('time.%d', i), sprintf('illumination.%d', i))]
}
cl.st <- lapply(cl.st, function(x){colnames(x) <- c('weights', 'damage.class', 'status', 'time', 'illumination');return(x)})

cl.st <- do.call(rbind, cl.st)  # bind all censuses together
cl.st <- cl.st[-which(is.na(cl.st[ ,'damage.class'])), ]  # remove those that don't have a damage class in t
cl.st <- cl.st[-which(cl.st[ ,'status'] == '?'), ] # remove those with unknown status at t + 1
spl <- split(cl.st, cl.st[ ,'damage.class'])

# Separate by census
u_mort <- c_mort <- matrix(NA, n.census-1, n.damage)
for(i in 1:(n.census-1)){
  cl.st <- df.ss[ ,c(sprintf('weights.ind.%d', i), sprintf('damage.class.%d',i), sprintf('status.%d', (i+1)), 
                     sprintf('time.%d', i), sprintf('illumination.%d', i))]
  colnames(cl.st) <- c('weights', 'damage.class', 'status', 'time', 'illumination')
  cl.st <- cl.st[-which(is.na(cl.st[ ,'damage.class'])), ] 
  cl.st <- cl.st[which(cl.st[ ,'status'] != '?'), ] # remove those with unknown status at t + 1
  u.df <- cl.st[which(cl.st[ ,'illumination'] %in% c(1,2)), ]
  c.df <- cl.st[which(cl.st[ ,'illumination'] %in% c(4,5)), ]
  
  u.spl <- split(u.df, u.df[ ,'damage.class'])
  u_mort[i, ] <- unlist(lapply(u.spl, get.mort.rate.weights))[1:n.damage]
  c.spl <- split(c.df, c.df[ ,'damage.class'])
  c_mort[i, ] <- unlist(lapply(c.spl, get.mort.rate.weights))[1:n.damage]
}
u_mort
c_mort
c_mort[which(is.na(c_mort))] <- 0
c_mort[which(is.infinite(c_mort))] <- 0

rownames(c_mort) <- NULL
rownames(u_mort) <- NULL

c_mort


write.csv(u_mort, file = '../../Damage_processing/r_understory_mort_by_damage_by_cens.csv', quote = FALSE)
write.csv(c_mort, file = '../../Damage_processing/r_canopy_mort_by_damage_by_cens.csv', quote = FALSE)


