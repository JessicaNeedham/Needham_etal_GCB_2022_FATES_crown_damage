####################################
### Script to read in and analyse 
### mortality survey data ###
###################################

### Note that rows correspond to STEMS and NOT INDIVIDUALS

# For more detailed descriptions of all these variables see Arellano et al., 2020
# Journal of Vegetation Science

### Some definitions: 
## STATUS 
# A - stem alive
# D - dead individuals
# X - stem entirely dead tissue
# NF - stem not found
# OK - healthy and undamaged
# ? 

## MODE
# S - standing - main axis complete (not necessarily vertical)
#     can be composed of entirely dead tissue if still continuous
# B - broken - main axis is broken or incomplete but some of it 
#     is still standing. 
# U - uprooted - roots aboveground. 
# ? - tree has died and only tag can be found

## LIVING LENGTH 
# Estimate of the proportion of remaining living tissue in the main
# axis of the stem in a B or S stem. 
# To work out biomass lost:
# Use allometry to calculate target AGB from dbh
# Then work out how much biomass has been lost above the living length

## Characterising crown loss
# Estimate remaining crown (%) within the living length
# Only applies to new damage - i.e. based on evidence of dead 
# or broken branches. 100% means no recent loss of branches.
# No branches in the first place noted by NA

# Crown illumination
# 1 = no direct light
# 2 = lateral light, < 10% vertical light
# 3 = some overhead light (10-90%)
# 4 = full overhead light
# 5 = leaves completely exposed to vertical and lateral ligth (emergent?)

# Also note that this data was collected following a spatially stratefied sampling
# design and therefore we need to correct for this to get plot level metrics
##############################################################################
rm(list = ls())

# LIBRARIES
library(stringr)

# FUNCTIONS
source(file = 'Source/Allometry_functions.R')

n.census <- 4
site <- 'bci'
n.damage <- 5

# Parameters from Gabriel and Daniel 
a <- 2.631
b <- 9.761049
prop_branches <- 1/3

# Load the data and merge it into one table
data.ls <- vector('list', n.census)
for(i in 1:n.census){
  data.ls[[i]] <- read.table(file = sprintf('Data/%s%d.subset.fates.txt', site, i))
  colnames(data.ls[[i]])[c(11:20, 26:27)] <- paste0(colnames(data.ls[[i]][c(11:20,26:27)]), '.', i)
  data.ls[[i]] <- data.ls[[i]][ ,-grep('survey', colnames(data.ls[[i]]))]
}

df <- merge(data.ls[[1]], data.ls[[2]], by = colnames(data.ls[[1]])[c(seq(9), seq(20,24,1))])
df <- merge(df, data.ls[[3]], by = colnames(df)[1:14])
df <- merge(df, data.ls[[4]], by = colnames(df)[1:14])
rm(data.ls)

# remove comments column
df <- df[ ,-grep('comments', colnames(df))]

# exclude tree ferns and palms - these all belong 
# to the families listed below
exclude.fam <- c("Arecaceae", "Thyrsopteridaceae", "Loxsomataceae",
                 "Culcitatceae", 
                 "Plagiohyriaceae", "Cibotiaceae", "Cyatheaceae",
                 "Dicksoniaceae",
                 "Meetaxyaceae", "Heliconiaceae", 'Poaceae')

cut <- which(df$family %in% exclude.fam)
df <- df[-cut, ]

# Make a wood density column
load(sprintf('Data/%s.spptable.rdata',site))
tmp <- eval(parse(text = sprintf('%s.spptable', site)))
sp.mat <- tmp
sp.mat$Genus <- word(sp.mat$Latin, 1, 1)
# add a wsg column to the main data table
wsg <- tmp[match(df[ ,'spcode'], tmp[ ,'sp']) ,'wsg']
wsg.level <- tmp[match(df[, 'spcode'], tmp[ ,'sp']), 'wsglevel']
df <- cbind(df, wsg, wsg.level)
# stems with no wsg give them mean
df$wsg <- ifelse(is.na(df$wsg), mean(sp.mat$wsg), df$wsg)
df$wsg.level <- ifelse(is.na(df$wsg.level), 'mean', df$wsg.level)

# height in m based on dbh
df$target_height <- d2h_nomax_martcano(df$dbh.full.census/10)
# get allometric AGB (kg)
df$target_agb_chave <- d2agb(df$dbh.full.census/10, df$target_height, wd = df$wsg) 

# relative crown height - (assume crown height is 1.3m )
df$Rc <- prop_branches/df$target_height

#  make living length and branch damage numeric 
df[ ,grep('living', colnames(df))] <- apply(df[ ,grep('living', colnames(df))], c(1,2), as.numeric)
df[ ,grep('branches', colnames(df))] <- apply(df[ ,grep('branches', colnames(df))], c(1,2), as.numeric)

# for each census get actual biomass
# note that branches refers only to branches lost in the current census interval. 
# So to get biomass in a given census we need to use cumulative branch loss. 
for(i in 1:n.census){
  ll <- df[ ,grep(sprintf('living_length.%d', i), colnames(df))] # living length
  # proportion alive branches in living length - cumulative
  br <- df[ ,grep('branches', colnames(df))]
  br.cumulative <- t(apply(br, 1, function(x) cumprod(as.numeric(x/100))))
  br.census <- br.cumulative[ ,i]
  # actual height as proportion of target height
  i.p <- ll/df$target_height
  # proportion of biomass in bole at the actual height, times by biomass
  trunk <- v_bole(i.p, a, b, df$Rc) * df$target_agb_chave
  # proportion of biomass in crown at actual height, times by biomass and proportion branches left
  crown <- v_crown(i.p, a, b, df$Rc) * df$target_agb_chave * br.census
  tree <- trunk + crown # total biomass 
  
  df <- cbind(df, tree)
  colnames(df)[ncol(df)] <- paste0('agb.', i)
}

# from actual biomass get change in biomass each year
agb <- df[ ,grep('agb.', colnames(df))]
damage <- t(apply(agb, 1, diff))

# remove first column - it will overestimate damage because all departures from
# target biomass would be assumed to have occurred in the first year. 
# also NA increases in biomass - that is recovery - we only want damage here
damage <- damage[ ,2:4]
damage <- apply(damage, 2, function(x){x <- ifelse(x>=0, NA, x); return(x)})
# bind new damage to dataframe
colnames(damage) <- sapply(seq(ncol(damage)), function(x) paste0('new.damage.', x))
df <- cbind(df, damage)

# based on actual v target biomass assign each stem a damage class (1-5)
for(i in 1:(n.census)){
  damage.class <- mapply(function(x,y) get_damage_class(x,y, n.damage), 
                         x = df$target_agb_chave, 
                         y = df[ ,grep(sprintf('agb.%d', i), colnames(df))])
  df <- cbind(df, damage.class)
  colnames(df)[ncol(df)] <- paste0('damage.class.', i)
}


# sort out dates
time <- matrix(NA, nrow = nrow(df), ncol = n.census-1)
for(i in 1:(n.census-1)){
  d1 <- as.Date(df[ ,grep(paste0('date.', i), colnames(df))])
  d2 <- as.Date(df[ ,grep(paste0('date.', i+1), colnames(df))])
  time[ ,i] <- as.numeric((d2-d1)/365)
}

colnames(time) <- sapply(seq(n.census-1), function(x) paste0('time.', x))
df <- cbind(df, time)

# MORTALITY
statuses <- df[ ,grep('status.', colnames(df))][ ,2:5]

# for trees that were declared dead but then came back to life... 
# make the dead year alive
statuses <- t(apply(statuses, 1, function(x){
  x[1] <- ifelse(x[2] == 'A' & x[1] == 'D', 'A', x[1]); x}))
statuses <- t(apply(statuses, 1, function(x){
  x[2] <- ifelse(x[3] == 'A' & x[2] == 'D', 'A', x[2]); x}))
statuses <- t(apply(statuses, 1, function(x){
  x[3] <- ifelse(x[4] == 'A' & x[3] == 'D', 'A', x[3]); x}))

# same for question marks
statuses <- t(apply(statuses, 1, function(x){
  x[1] <- ifelse(x[2] == 'A' & x[1] == '?', 'A', x[1]); x}))
statuses <- t(apply(statuses, 1, function(x){
  x[2] <- ifelse(x[3] == 'A' & x[2] == '?', 'A', x[2]); x}))
statuses <- t(apply(statuses, 1, function(x){
  x[3] <- ifelse(x[4] == 'A' & x[3] == '?', 'A', x[3]); x}))

# same for question marks - if in following year it was dead call it a dead
statuses <- t(apply(statuses, 1, function(x){
  x[1] <- ifelse(x[2] == 'D' & x[1] == '?', 'D', x[1]); x}))
statuses <- t(apply(statuses, 1, function(x){
  x[2] <- ifelse(x[3] == 'D' & x[2] == '?', 'D', x[2]); x}))
statuses <- t(apply(statuses, 1, function(x){
  x[3] <- ifelse(x[4] == 'D' & x[3] == '?', 'D', x[3]); x}))

df[ ,grep('status.', colnames(df))][ ,2:5] <- statuses

# if status is dead then make cd class 6
for(i in 1:n.census){
  for(j in 1:nrow(df)){
    df[j, sprintf('damage.class.%d',i)] <- ifelse(
      df[j, sprintf('status.%d', i)] == 'D', 6, df[j, sprintf('damage.class.%d',i)])
  }
}

# For mortality remove multi-stems
multis <- df[which(df$ind.id %in% names(which(table(df$ind.id) >1))), ]
# for non multi stem trees proportion of trees that died in each damage class
df.ss <- df[-which(df$ind.id %in% names(which(table(df$ind.id) >1))), ]

# All individuals but just main stem from each one!
ind.max.stem <- tapply(df$dbh.full.census, df$ind.id, which.max)
df.main <- do.call(rbind, mapply(function(x,y){x <- x[y, ]; return(x)}, 
                                 x = split(df, df$ind.id), y = ind.max.stem, SIMPLIFY = FALSE))


save(df, file = 'Output/df.RData')
save(df.ss, file = 'Output/df_single_stem.RData')
save(df.main, file = 'Output/df_main.RData')
