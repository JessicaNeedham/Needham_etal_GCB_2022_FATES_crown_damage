# clear the working directory
rm(list = ls())
# load libraries
library(lubridate)
library(plyr)
library(RColorBrewer)

# load data
cens <- c(16, 18, 20, 22, 23, 24)
n.cens <- length(cens)
bciall <- vector('list', n.cens)
for(i in 1:n.cens){
  bciall[[i]] <- read.delim(paste0('Data/Dendrometer_data/dend.bci50.c',cens[i] ,'.txt'))
  bciall[[i]] <- bciall[[i]][ ,-match(c('origorder', 'direction', 'notes', 'entrynotes', 'dendroHt', 
                              'paintDiam', 'paintHt', 'match.dendro'), 
                            colnames(bciall[[i]]))]
}

bciall <- do.call(rbind, bciall)

# remove multistem trees
n.stems <- tapply(bciall$stemtag, bciall$tag, function(x) length(unique(x)))
multis <- names(which(n.stems > 1))
bciall <- bciall[-which(bciall$tag %in% multis), ]

# split it by census
bci.ls <- split(bciall, bciall$census)

# merge based on tag, stemtag, species and dendrometer ID. 
# This means that if a tree had a new dendrometer band at any point, it will 
# essentially be treated as two different trees. 
m1 <- merge(bci.ls[[1]], bci.ls[[2]], by = c('tag', 'stemtag', 'species', 'DendroID'))
m2 <- merge(bci.ls[[3]], bci.ls[[4]], by = c('tag', 'stemtag', 'species', 'DendroID'))
m3 <- merge(bci.ls[[5]], bci.ls[[6]], by = c('tag', 'stemtag', 'species','DendroID'))

m4 <- merge(m1, m2, by = c('tag', 'stemtag', 'species','DendroID'))
m5 <- merge(m4, m3, by = c('tag', 'stemtag', 'species','DendroID'))



# Sort out column names
m5 <- m5[ ,-grep('census', colnames(m5))]
colnames(m5)
colnames(m5)[grep('stemdead', colnames(m5))] <- sapply(seq(6), function(x) paste0('stemdead.', x))
colnames(m5)[grep('lianas', colnames(m5))] <- sapply(seq(6), function(x) paste0('lianas.', x))
colnames(m5)[grep('date', colnames(m5))] <- sapply(seq(6), function(x) paste0('date.', x))
colnames(m5)[grep('illumination', colnames(m5))] <- sapply(seq(6), function(x) paste0('illumination.', x))
colnames(m5)[grep('calc.dbh', colnames(m5))] <- sapply(seq(6), function(x) paste0('calc.dbh.', x))
colnames(m5)[grep('crown', colnames(m5))] <- sapply(seq(6), function(x) paste0('crown.', x))
colnames(m5)[grep('type', colnames(m5))] <- sapply(seq(6), function(x) paste0('type.', x))
colnames(m5)[grep('measure', colnames(m5))] <- sapply(seq(6), function(x) paste0('measure.', x))
colnames(m5)[grep('dendroDiam', colnames(m5))] <- sapply(seq(6), function(x) paste0('dendroDiam.', x))
colnames(m5)[grep('code', colnames(m5))] <- sapply(seq(6), function(x) paste0('code.', x))

head(m5)
n.census <- 6

# how many stems have more than one dendrometer? 
no <- tapply(m5$DendroID, m5$tag, function(x) length(unique(x)))
multis <- names(which(no > 1))
tmp <- m5[which(m5$tag %in% multis), ]
m5 <- m5[-which(m5$tag %in% multis), ]

# these are the stems with more than one dendrometer
tmp <- split(tmp, tmp$tag)

# this function takes a stem and returns an index for which row has the most calc.dbh measurements
longest <- function(tree){
  n.cen <- apply(tree, 1, function(x) length(which(!is.na(x[grep('calc.dbh', names(x))]))))
  keep <- which.max(n.cen)
  return(keep)
}

# apply that function to the dataframe of multi dendro trees
keeps <- lapply(tmp, longest)

# then we subset the multi dendro tree dataframe by index of longest calc.dbh
treekeep <- mapply(function(x,y) x[y, ], x = tmp, y = keeps, SIMPLIFY=FALSE)
treekeep <- do.call(rbind, treekeep)

m5 <- rbind(m5, treekeep)

# Now each stem should only have one dendrometer ID

# At this point there are still some trees that have the same tag and dendro ID repeated
tab <- table(m5$tag)
which(tab > 1)
# keep only the first row
tags <- unique(m5$tag)

mat <- match(tags, m5$tag)
m5 <- m5[mat, ]

write.table(m5, file = 'Output/bciall_wide_v3.txt', row.names = FALSE, sep = ' ')

## any stem that was dead then alive, make the dead an alive
correct.zombies <- function(stem){
  status <- stem[grep('stemdead', colnames(stem))]
  status[1:(length(status)-1)] <- unlist(sapply(seq(length(status)-1), function(x) status[x] <- ifelse(status[x+1] == 'alive', 'alive', status[x])))
  stem[grep('stemdead', colnames(stem))] <- status
  return(stem)
}

for(i in 1:nrow(m5)){
  m5[i, ] <- correct.zombies(m5[i, ])
}

# make sure calc dbh is nan if stemdead is dead
for(i in 1:n.census){
  dbh <- m5[ ,grep(paste0('calc.dbh.', i), colnames(m5))]
  status <- m5[ ,grep(paste0('stemdead.', i), colnames(m5))]
  dbh <- ifelse(status  == 'dead', NA, dbh)
}

# make annual growth rate columns
dbh <- m5[ ,grep('calc.dbh', colnames(m5))]
growth <- t(apply(dbh, 1, diff))

# divide growth by time in years between censuses
for(i in 1:(n.census-1)){
  d1 <- as.Date(m5[ ,grep(paste0('date.', i), colnames(m5))])
  d2 <- as.Date(m5[ ,grep(paste0('date.', i+1), colnames(m5))])
  time <- as.numeric((d2-d1)/365)
  growth[ ,i] <- growth[ ,i]/time
}

# correct growth (can't grow more than 20 cm a year, or shrink more than 10 mm a year)
#length(which(growth < -10))
growth[which(growth < -10)] <- NA  # 14 growth rates
growth[which(growth > 200)] <- NA  # one growth rate (this is almost certainly an error - 43364.18 mm)

# can't shrink more than 10% of original size
# this only removes 4 growth rates 
growth[ ,1] <- ifelse(m5[ ,'calc.dbh.2'] < m5[ ,'calc.dbh.1']*0.9, NA, growth[ ,1])
growth[ ,2] <- ifelse(m5[ ,'calc.dbh.3'] < m5[ ,'calc.dbh.2']*0.9, NA, growth[ ,2])
growth[ ,3] <- ifelse(m5[ ,'calc.dbh.4'] < m5[ ,'calc.dbh.3']*0.9, NA, growth[ ,3])
growth[ ,4] <- ifelse(m5[ ,'calc.dbh.5'] < m5[ ,'calc.dbh.4']*0.9, NA, growth[ ,4])
growth[ ,5] <- ifelse(m5[ ,'calc.dbh.6'] < m5[ ,'calc.dbh.5']*0.9, NA, growth[ ,5])

colnames(growth) <- sapply(seq(n.census-1), function(x) paste0('growth.', x))
m5 <- cbind(m5, growth)

# Reversing damage class order so that damage class 1 is the least damaged,
# class 4 is the most; nans stay nans
# This is just so the order matches FATES
for(i in 1:n.census){
  damage <- m5[ ,grep(paste0('crown.', i), colnames(m5))]
  damage <- mapvalues(damage, from = seq(4), to = seq(4,1,-1))
  m5[ ,grep(paste0('crown.', i), colnames(m5))] <- damage
}

# how many stems in each crown illumination class
ci <- m5[ ,grep('illumination.', colnames(m5))]
apply(ci, 2, table)

# Make a binary canopy v understory column based on crown illumination
m5[ ,'canopy.1'] <- ifelse(m5[ ,'illumination.1'] >=4, 1, 0)
m5[ ,'canopy.2'] <- ifelse(m5[ ,'illumination.2'] >=4, 1, 0)
m5[ ,'canopy.3'] <- ifelse(m5[ ,'illumination.3'] >=4, 1, 0)
m5[ ,'canopy.4'] <- ifelse(m5[ ,'illumination.4'] >=4, 1, 0)
m5[ ,'canopy.5'] <- ifelse(m5[ ,'illumination.5'] >=4, 1, 0)
m5[ ,'canopy.6'] <- ifelse(m5[ ,'illumination.6'] >=4, 1, 0)

# illumination 3 is intermediate and not comparable to FATES. Make these NA
# and NA corresponding growth rates
m5[which(m5[ ,'illumination.1'] ==3) ,'canopy.1'] <- NA
m5[which(m5[ ,'illumination.2'] ==3) ,'canopy.2'] <- NA
m5[which(m5[ ,'illumination.3'] ==3) ,'canopy.3'] <- NA
m5[which(m5[ ,'illumination.4'] ==3) ,'canopy.4'] <- NA
m5[which(m5[ ,'illumination.5'] ==3) ,'canopy.5'] <- NA
m5[which(m5[ ,'illumination.6'] ==3) ,'canopy.6'] <- NA

# make sure growth is nan if illumination is 3
for(i in 1:(n.census-1)){
  growth <- m5[ ,grep(paste0('growth.', i), colnames(m5))]
  canopy <- m5[ ,grep(paste0('canopy.', i), colnames(m5))]
  growth <- ifelse(is.na(canopy), NA, growth)
  m5[ ,grep(paste0('growth.', i), colnames(m5))] <- growth
}

# Now get growth by damage by canopy layer by census
growth_can <- matrix(NA, nrow = 4, ncol = n.census-1)
growth_ustory <- matrix(NA, nrow = 4, ncol = n.census-1)

for(i in 1:(n.census-1)){
  cens <- aggregate(m5[ ,paste0('growth.', i)], 
                    by = list(m5[ ,paste0('crown.', i)], m5[ ,paste0('canopy.', i)]), 
                    mean, na.rm = TRUE)
  growth_can[ ,i] <- cens[which(cens[ ,'Group.2'] == 1), 'x']
  growth_ustory[ ,i] <- cens[which(cens[ ,'Group.2'] == 0), 'x']
}

# now the mean across all censuses
can_gr <- matrix(NA, ncol = 2)
ustr_gr <- matrix(NA, ncol = 2)

for(i in 1:(n.census-1)){
  can_grw <- m5[which(m5[ ,grep(paste0('canopy.', i), colnames(m5))] == 1),
                grep(paste0('growth.',i), colnames(m5))]
  can_dm <- m5[which(m5[ ,grep(paste0('canopy.', i), colnames(m5))] == 1),
               grep(paste0('crown.',i), colnames(m5))]
  
  can <- cbind(can_grw, can_dm)
  can_gr <- rbind(can_gr, can)
  
  ustr_grw <- m5[which(m5[ ,grep(paste0('canopy.',i), colnames(m5))] == 0),
                 grep(paste0('growth.',i), colnames(m5))]
  ustr_dm <- m5[which(m5[ ,grep(paste0('canopy.',i), colnames(m5))] == 0),
                grep(paste0('crown.',i), colnames(m5))]
  
  ustr <- cbind(ustr_grw, ustr_dm)
  ustr_gr <- rbind(ustr_gr, ustr)
  
}

canopy_growth <- tapply(can_gr[ ,'can_grw'], can_gr[, 'can_dm'], mean, na.rm = TRUE)
ustory_growth <- tapply(ustr_gr[ ,'ustr_grw'], ustr_gr[, 'ustr_dm'], mean, na.rm = TRUE)

canopy_q <- t(do.call(rbind, tapply(can_gr[ ,'can_grw'], can_gr[, 'can_dm'], quantile, c(0.25, 0.75), 
                                    na.rm = TRUE)))
ustory_q <- t(do.call(rbind, tapply(ustr_gr[ ,'ustr_grw'], ustr_gr[, 'ustr_dm'], 
                                    quantile, c(0.25, 0.75), na.rm = TRUE)))



#######################################################################################
########################################################################################
### Figues ###
#####################################################################

cols <- brewer.pal(10, 'Paired')[c(4,2,8,7)]

make.transp <- function(tmpcol, transp = 50) {
  col.out <- paste(tmpcol, transp, sep = "")
  return(col.out)
}

cols.light <- make.transp(cols)


### Presentation figure - just canopy
jpeg('Figures/Growth_fates_v_bci_canopy.jpeg', width = 480, height = 480)
par(mfrow = c(1,1), mar = c(2, 3, 2, 1), oma = c(2,2,0,0))

# Canopy low
plot(seq(12.5, 87.5, 25), canopy_growth, ylim = c(0, 7),  lwd = 2, col = cols[1],
     xlab = '', ylab = '', las = 1, xlim = c(0,100), cex.axis=1.6, cex = 1.5, pch = 1)
arrows(x0= c(0,25,50,75), x1=c(25,50,75,100),
       y0= canopy_growth, y1=canopy_growth, 
       length = 0.1, angle = 90, col = cols[1], code = 3)
arrows(x0= seq(12.5, 87.5, 25), x1=seq(12.5, 87.5, 25),
       y0= canopy_q[1,], y1=canopy_q[2,], 
       length = 0.1, angle = 90, col = cols[1], code = 3)

can_fates_low <- read.csv('/Users/JFNeedham/Damage_processing/Canopy_growth_lowrootN.csv', header = FALSE)
can_fates_low <- can_fates_low * 10  # convert cm to mm
points(seq(0, 80, 20), can_fates_low[ ,1], col = cols[3], lwd = 3, pch=1, cex = 1.5)

can_fates_high <- read.csv('/Users/JFNeedham/Damage_processing/Canopy_growth_highrootN.csv', header = FALSE)
can_fates_high <- can_fates_high * 10  # convert cm to mm
points(seq(0, 80, 20), can_fates_high[ ,1], col = cols[2], lwd = 3, pch=1, cex = 1.5)

mtext('Crown loss (%)', side = 1, line = 0.75, cex = 1.4, outer = TRUE)
mtext(expression(paste('Growth (mm yr' ^ '-1',')')),
      side = 2, line = 0, cex = 1.4, outer = TRUE)
legend('topright', legend = c('BCI', 'FATES low root R', 'FATES high root R'), 
       col = c(cols[1], cols[3], cols[2]),  cex = 1.5, pch = 1, pt.cex =2, lwd=2)
dev.off()


### Paper figure Canopy and understory
pdf('Figures/Growth_fates_v_bci_canopy.pdf', width = 12, height = 7)
par(mfrow = c(1,2), mar = c(2, 3, 3, 1), oma = c(2,2,0,0))

# Canopy low
plot(seq(12.5, 87.5, 25), canopy_growth, ylim = c(0, 7),  lwd = 2, col = cols[1],
     xlab = '', ylab = '', las = 1, xlim = c(0,100), cex.axis=1.6, cex = 1.5, pch = 1)
arrows(x0= c(0,25,50,75), x1=c(25,50,75,100),
       y0= canopy_growth, y1=canopy_growth, 
       length = 0.1, angle = 90, col = cols[1], code = 3)
arrows(x0= seq(12.5, 87.5, 25), x1=seq(12.5, 87.5, 25),
       y0= canopy_q[1,], y1=canopy_q[2,], 
       length = 0.1, angle = 90, col = cols[1], code = 3)

can_fates_low <- read.csv('/Users/JFNeedham/Damage_processing/Canopy_growth_lowrootN.csv', header = FALSE)
can_fates_low <- can_fates_low * 10  # convert cm to mm
points(seq(0, 80, 20), can_fates_low[ ,1], col = cols[3], lwd = 3, pch=1, cex = 1.5)

can_fates_high <- read.csv('/Users/JFNeedham/Damage_processing/Canopy_growth_highrootN.csv', header = FALSE)
can_fates_high <- can_fates_high * 10  # convert cm to mm
points(seq(0, 80, 20), can_fates_high[ ,1], col = cols[2], lwd = 3, pch=1, cex = 1.5)

mtext('Crown loss (%)', side = 1, line = 0.75, cex = 1.4, outer = TRUE)
mtext(expression(paste('Growth (mm yr' ^ '-1',')')),
      side = 2, line = 0, cex = 1.4, outer = TRUE)
legend('topright', legend = c('BCI', 'FATES low root R', 'FATES high root R'), 
       col = c(cols[1], cols[3], cols[2]),  cex = 1.2, pch = 1, pt.cex =2, lwd=2)
mtext('Canopy', side = 3, line = 0.75, cex = 1.4, outer =FALSE)
legend('topleft', legend = as.expression(bquote(bold('(a)'))),  bty = 'n', cex = 1.5)

# Understory low
plot(seq(12.5, 87.5, 25), ustory_growth, ylim = c(0, 7), type = 'p', lwd = 2, col = cols[1],
     xlab = '', ylab = '', las = 1, xlim = c(0,100), cex.axis=1.6, cex =1.5)
legend('topleft', legend = as.expression(bquote(bold('b)'))),  bty = 'n', cex = 1.5)
mtext('Understory', side = 3, line = 0.75, cex = 1.4)
arrows(x0= c(0,25,50,75), x1=c(25,50,75,100),
       y0= ustory_growth, y1=ustory_growth, 
       length = 0.1, angle = 90, col = cols[1], code = 3)
arrows(x0= seq(12.5, 87.5, 25), x1=seq(12.5, 87.5, 25),
       y0= ustory_q[1,], y1=ustory_q[2,], 
       length = 0.1, angle = 90, col = cols[1], code = 3)
ustory_fates_low <- read.csv('/Users/JFNeedham/Damage_processing/Understory_growth_lowrootN.csv', header = FALSE)
ustory_fates_low <- ustory_fates_low * 10  # convert cm to mm
points(seq(0, 80, 20), ustory_fates_low[ ,1], col = cols[3], lwd = 3, cex=1.5)

ustory_fates_high <- read.csv('/Users/JFNeedham/Damage_processing/Understory_growth_highrootN.csv', header = FALSE)
ustory_fates_high <- ustory_fates_high * 10  # convert cm to mm
points(seq(0, 80, 20), ustory_fates_high[ ,1], col = cols[2], lwd = 3, cex=1.5)

dev.off()


###############################
##################################################################################
##################################################################################

# Now get growth by liana load by canopy layer by census
growth_can_liana <- matrix(NA, nrow = 5, ncol = n.census-1)
growth_ustory_liana <- matrix(NA, nrow = 5, ncol = n.census-1)

for(i in 1:(n.census-1)){
  cens <- aggregate(m5[ ,paste0('growth.', i)], 
                    by = list(m5[ ,paste0('lianas.', i)], m5[ ,paste0('canopy.', i)]), 
                    mean, na.rm = TRUE)
  growth_can_liana[ ,i] <- cens[which(cens[ ,'Group.2'] == 1), 'x']
  growth_ustory_liana[ ,i] <- cens[which(cens[ ,'Group.2'] == 0), 'x']
}

###############################################################################
### Repeat but get rid of trees with very high liana loads
m5[which(m5[ ,'lianas.1'] ==4) ,'growth.1'] <- NA
m5[which(m5[ ,'lianas.2'] ==4) ,'growth.2'] <- NA
m5[which(m5[ ,'lianas.3'] ==4) ,'growth.3'] <- NA
m5[which(m5[ ,'lianas.4'] ==4) ,'growth.4'] <- NA
m5[which(m5[ ,'lianas.5'] ==4) ,'growth.5'] <- NA
m5[which(m5[ ,'lianas.6'] ==5) ,'growth.6'] <- NA

m5[which(m5[ ,'lianas.1'] ==3) ,'growth.1'] <- NA
m5[which(m5[ ,'lianas.2'] ==3) ,'growth.2'] <- NA
m5[which(m5[ ,'lianas.3'] ==3) ,'growth.3'] <- NA
m5[which(m5[ ,'lianas.4'] ==3) ,'growth.4'] <- NA
m5[which(m5[ ,'lianas.5'] ==3) ,'growth.5'] <- NA
m5[which(m5[ ,'lianas.6'] ==3) ,'growth.6'] <- NA

# Now get growth by damage by canopy layer by census
growth_can <- matrix(NA, nrow = 4, ncol = n.census-1)
growth_ustory <- matrix(NA, nrow = 4, ncol = n.census-1)

for(i in 1:(n.census-1)){
  cens <- aggregate(m5[ ,paste0('growth.', i)], 
                    by = list(m5[ ,paste0('crown.', i)], m5[ ,paste0('canopy.', i)]), 
                    mean, na.rm = TRUE)
  growth_can[ ,i] <- cens[which(cens[ ,'Group.2'] == 1), 'x']
  growth_ustory[ ,i] <- cens[which(cens[ ,'Group.2'] == 0), 'x']
}

# now the mean across all censuses
can_gr <- matrix(NA, ncol = 2)
ustr_gr <- matrix(NA, ncol = 2)

for(i in 1:(n.census-1)){
  can_grw <- m5[which(m5[ ,grep(paste0('canopy.', i), colnames(m5))] == 1),
                grep(paste0('growth.',i), colnames(m5))]
  can_dm <- m5[which(m5[ ,grep(paste0('canopy.', i), colnames(m5))] == 1),
               grep(paste0('crown.',i), colnames(m5))]
  
  can <- cbind(can_grw, can_dm)
  can_gr <- rbind(can_gr, can)
  
  ustr_grw <- m5[which(m5[ ,grep(paste0('canopy.',i), colnames(m5))] == 0),
                 grep(paste0('growth.',i), colnames(m5))]
  ustr_dm <- m5[which(m5[ ,grep(paste0('canopy.',i), colnames(m5))] == 0),
                grep(paste0('crown.',i), colnames(m5))]
  
  ustr <- cbind(ustr_grw, ustr_dm)
  ustr_gr <- rbind(ustr_gr, ustr)
  
}

canopy_growth <- tapply(can_gr[ ,'can_grw'], can_gr[, 'can_dm'], mean, na.rm = TRUE)
ustory_growth <- tapply(ustr_gr[ ,'ustr_grw'], ustr_gr[, 'ustr_dm'], mean, na.rm = TRUE)

canopy_lowerq <- tapply(can_gr[ ,'can_grw'], can_gr[, 'can_dm'], quantile, 0.25,  na.rm = TRUE)
ustory_lowerq <- tapply(ustr_gr[ ,'ustr_grw'], ustr_gr[, 'ustr_dm'], quantile, 0.25, na.rm = TRUE)

canopy_upperq <- tapply(can_gr[ ,'can_grw'], can_gr[, 'can_dm'], quantile, 0.75,  na.rm = TRUE)
ustory_upperq <- tapply(ustr_gr[ ,'ustr_grw'], ustr_gr[, 'ustr_dm'], quantile, 0.75, na.rm = TRUE)



### Paper figure Canopy and understory
pdf('Figures/Growth_fates_v_bci_canopy_no_lianas.pdf', width = 12, height = 7)
par(mfrow = c(1,2), mar = c(2, 3, 3, 1), oma = c(2,2,0,0))

# Canopy low
plot(seq(12.5, 87.5, 25), canopy_growth, ylim = c(0, 7),  lwd = 2, col = cols[1],
     xlab = '', ylab = '', las = 1, xlim = c(0,100), cex.axis=1.6, cex = 1.5, pch = 1)
arrows(x0= c(0,25,50,75), x1=c(25,50,75,100),
       y0= canopy_growth, y1=canopy_growth, 
       length = 0.1, angle = 90, col = cols[1], code = 3)
arrows(x0= seq(12.5, 87.5, 25), x1=seq(12.5, 87.5, 25),
       y0= canopy_q[1,], y1=canopy_q[2,], 
       length = 0.1, angle = 90, col = cols[1], code = 3)

can_fates_low <- read.csv('/Users/JFNeedham/Damage_processing/Canopy_growth_lowrootN.csv', header = FALSE)
can_fates_low <- can_fates_low * 10  # convert cm to mm
points(seq(0, 80, 20), can_fates_low[ ,1], col = cols[3], lwd = 3, pch=1, cex = 1.5)

can_fates_high <- read.csv('/Users/JFNeedham/Damage_processing/Canopy_growth_highrootN.csv', header = FALSE)
can_fates_high <- can_fates_high * 10  # convert cm to mm
points(seq(0, 80, 20), can_fates_high[ ,1], col = cols[2], lwd = 3, pch=1, cex = 1.5)

mtext('Crown loss (%)', side = 1, line = 0.75, cex = 1.4, outer = TRUE)
mtext(expression(paste('Growth (mm yr' ^ '-1',')')),
      side = 2, line = 0, cex = 1.4, outer = TRUE)
legend('topright', legend = c('BCI', 'FATES low root R', 'FATES high root R'), 
       col = c(cols[1], cols[3], cols[2]),  cex = 1.2, pch = 1, pt.cex =2, lwd=2)
mtext('Canopy', side = 3, line = 0.75, cex = 1.4, outer =FALSE)
legend('topleft', legend = as.expression(bquote(bold('(a)'))),  bty = 'n', cex = 1.5)

# Understory low
plot(seq(12.5, 87.5, 25), ustory_growth, ylim = c(0, 7), type = 'p', lwd = 2, col = cols[1],
     xlab = '', ylab = '', las = 1, xlim = c(0,100), cex.axis=1.6, cex =1.5)
legend('topleft', legend = as.expression(bquote(bold('b)'))),  bty = 'n', cex = 1.5)
mtext('Understory', side = 3, line = 0.75, cex = 1.4)
arrows(x0= c(0,25,50,75), x1=c(25,50,75,100),
       y0= ustory_growth, y1=ustory_growth, 
       length = 0.1, angle = 90, col = cols[1], code = 3)
arrows(x0= seq(12.5, 87.5, 25), x1=seq(12.5, 87.5, 25),
       y0= ustory_q[1,], y1=ustory_q[2,], 
       length = 0.1, angle = 90, col = cols[1], code = 3)
ustory_fates_low <- read.csv('/Users/JFNeedham/Damage_processing/Understory_growth_lowrootN.csv', header = FALSE)
ustory_fates_low <- ustory_fates_low * 10  # convert cm to mm
points(seq(0, 80, 20), ustory_fates_low[ ,1], col = cols[3], lwd = 3, cex=1.5)

ustory_fates_high <- read.csv('/Users/JFNeedham/Damage_processing/Understory_growth_highrootN.csv', header = FALSE)
ustory_fates_high <- ustory_fates_high * 10  # convert cm to mm
points(seq(0, 80, 20), ustory_fates_high[ ,1], col = cols[2], lwd = 3, cex=1.5)

dev.off()
