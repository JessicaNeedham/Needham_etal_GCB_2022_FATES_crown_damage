##########################
### Allometries ###
##########################


### Diameter to height ###

d2h_martcano <- function(dbh, p1 = 58.0, p2 = 0.73, p3 = 21.8, dbh_maxh){
  h <- (p1*min(dbh, dbh_maxh)^p2)/(p3 + min(dbh, dbh_maxh)^p2)
  return(h)
}

d2h_nomax_martcano <- function(dbh, p1 = 58.0, p2 = 0.73, p3 = 21.8){
  h <- (p1*dbh^p2)/(p3 + dbh^p2)
  return(h)
}

d2h_chave <- function(dbh, p1 = 0.893, p2 = 0.76, p3 = -0.034, E){
  p1e = p1 - E
  h = exp(p1e + p2*log(dbh) + p3*log(dbh)^2)
  return(h)
}

### Diameter to AGB ###

d2agb <- function(dbh, h, p1 = 0.0673, p2 = 0.976, wd, c2b = 2.0){
  bagw <- (p1 * (wd * dbh^2 * h)^p2)/c2b
  return(bagw)
}


# Functions for the distribution of stem and crown along the tree 
#(Ver Planck & McFarlane, Forestry 2014; 87, 459–469, doi:10.1093/forestry/cpu007)
# a, b are exponents. Higher values mean more abrupt tapering
# rC is the relative position of the base of the crown, in [0, 1]
# overall.p.bole describes the overall proportion of bole vs crown

v_bole0  <- function(r, a, overall.p.bole = 2/3){
  v_bole0 <- overall.p.bole*(1 - (1-r)^a)
  return(v_bole0)
}

v_whole <- function(r, a, b, rC = 0.5){
  beta <- (1 - v_bole0(rC, a))/(1 - rC^b)
  v_whole <- ifelse(r < rC, v_bole0(r, a), 
                    v_bole0(rC, a) + beta*(1 - (1- (r - rC))^b))
  return(v_whole)
}


v_bole <- function(r, a, b, rC = 0.5){
  v_bole <- pmin(v_whole(r, a, b, rC), v_bole0(r, a))
  return(v_bole)
}

v_crown <- function(r, a, b, rC = 0.5){
  v_crown <- v_whole(r, a, b, rC) - v_bole(r, a, b, rC)
  return(v_crown)
}

get_damage_class <- function(agb.target, agb.actual, n.damage.class){
  frac <- min(agb.actual/agb.target, 1)
  class_width <- 1/n.damage.class
  cd <- max(1, ceiling((1-frac)/class_width))
  return(cd)
}


# DBH to crown area
dbh.to.CA <- function(dbh){
  CA_exponent <- 1.3
  CA <- dbh ^ CA_exponent
  return(CA)
}


