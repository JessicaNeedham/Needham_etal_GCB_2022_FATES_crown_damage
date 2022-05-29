#### Figure 2 - damage schematic
rm(list=ls())

library(RColorBrewer)
cols <- brewer.pal(6, 'Dark2')[2:6]

g_per_kg <- 1000.0
kg_per_Megag <- 1000.0
cm2_per_m2 <- 10000.0

get_crown_reduction <- function(cd, ncrowndamage=5){
  class_width <- 1/ncrowndamage
  return(min(1, (cd-1)*class_width))
}


# Leaf
blmax_allom <- function(d, p1=0.1266844, p2=1.281329, rho=0.8374751, dbh_maxh=200, c2b = 2){
  blmax <- (p1*min(d, dbh_maxh)^p2)/c2b
  return(blmax)
}

bleaf <- function(d,cd, ncrowndamage = 5){
  blmax <- blmax_allom(d)
  if(cd >1){
    cr <- get_crown_reduction(cd, ncrowndamage)
    bl <- blmax * (1-cr)
  }else{
    bl <- blmax
  }
  return(bl)
}

dbhs <- seq(100)

leaf1 <- sapply(dbhs, bleaf, cd = 1)
leaf2 <- sapply(dbhs, bleaf, cd = 2)
leaf3 <- sapply(dbhs, bleaf, cd = 3)
leaf4 <- sapply(dbhs, bleaf, cd = 4)
leaf5 <- sapply(dbhs, bleaf, cd = 5)

#################################################
png(file = 'Recovery_schematic.png', width = 650)
par(mfrow = c(2,1), mar = c(1.5,1,2.5,1), oma = c(3,3,1,1))

plot(dbhs, leaf1, type = 'l', col = cols[1], lwd = 2, lty = 2, 
     xlab = '', ylab = '', las =1)
points(dbhs, leaf2, type = 'l', col = cols[2], lwd = 2, lty = 2)
points(dbhs, leaf3, type = 'l', col = cols[3], lwd = 2, lty = 2)
points(dbhs, leaf4, type = 'l', col = cols[4], lwd = 2, lty = 2)
points(dbhs, leaf5, type = 'l', col = cols[5], lwd = 2, lty = 2)

points(dbhs[1:50], leaf1[1:50], type = 'l', col = cols[1], lwd = 3)
arrows(x0=50, x1=50, y0 = leaf1[50], y1=leaf5[50], col = 'black', length=0.15, lwd = 2)
points(dbhs[50:100], leaf5[50:100], type = 'l', col = cols[5], lwd = 3)
arrows(x0=70, x1=70, y0 = leaf5[70], y1=leaf4[70], col = 'black', length=0.15, lwd = 2)
points(dbhs[70:100], leaf4[70:100], type = 'l', col = cols[4], lwd = 3)
mtext('Recovery', side = 3, line = 0.2, cex = 1.3, outer = FALSE)
mtext('a)', side = 3, line = 0.2, cex = 1.3, adj = 0.02)

text(85, 2, expression(paste('(1-', italic(fr), ') * ', italic(nmax))))
text(85, 6, expression(paste(italic(fr), ' * ', italic(nmax))))
legend('topleft', col = cols, 
       legend = c('0%', '20%', '40%', '60%', '80'), lwd = 2, 
       title = 'Crown loss')


plot(dbhs, leaf1, type = 'l', col = cols[1], lwd = 2, lty = 2, 
     xlab = '', ylab = '', las = 1)
points(dbhs, leaf2, type = 'l', col = cols[2], lwd = 2, lty = 2)
points(dbhs, leaf3, type = 'l', col = cols[3], lwd = 2, lty = 2)
points(dbhs, leaf4, type = 'l', col = cols[4], lwd = 2, lty = 2)
points(dbhs, leaf5, type = 'l', col = cols[5], lwd = 2, lty = 2)

points(dbhs[1:50], leaf1[1:50], type = 'l', col = cols[1], lwd = 3)
arrows(x0=50, x1=50, y0 = leaf1[50], y1=leaf5[50], col = 'black', length=0.15, lwd = 2)
points(dbhs[50:100], leaf5[50:100], type = 'l', col = cols[5], lwd = 3)
mtext('No Recovery', side = 3, line = 0.2, cex = 1.3, outer = FALSE)
mtext('b)', side = 3, line = 0.2, cex = 1.3, adj = 0.02)

mtext('DBH (cm)', side = 1, line = 1.5, cex = 1.3, outer = TRUE)
mtext('Leaf biomass (gC)', side = 2, line = 1.5, cex = 1.3, outer = TRUE)

dev.off()


