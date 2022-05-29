rm(list = ls())

mort <- read.csv(file = '../../Damage_processing/r_mort_by_damage.csv')
mort <- mort[ ,2]

################################################
png('Figures/Damage parameters.png')
par(mfrow = c(3,3), mar = c(3,3,3,3), oma = c(2,2,0,0))
p2 <- c(0.7, 0.8, 0.9)
p1 <- c(5, 5.5, 6)
x <- seq(0, 1, length.out =100)

for(i in 1:length(p1)){
  for(j in 1:length(p2)){
    dgmort <- 1/(1+exp(-p1[i] *(x - p2[j])))
    
    plot(x, dgmort, ylim = c(0,0.9), type = 'l')
    points(seq(0, 0.8, 0.2), mort)
    mtext(paste0('p:', p2[j], '  r:', p1[i]), side = 3,
          line = -2, adj = 0.1, cex = 1.2)
  }
}
mtext('Fraction crown loss', side = 1, line = 0, outer = TRUE, cex = 1.2)
mtext(expression(paste('Damage-mortality (yr' ^ '-1', ')')),
                 side = 2, line = 0, outer = TRUE, cex = 1.2)

dev.off()

