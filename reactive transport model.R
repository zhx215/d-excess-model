# developed by ZX and MJW, 2021
# run model using temperature gradient and advection transport as an example

rm(list=ls())

# Smith 1966 function
w_S66 <- function(lambda,tlcl){
  tlcl <- tlcl*9/5-459.67 # convert K to F
  exp(0.1133-log(lambda+1)+0.0393*tlcl)*10 # convert cm to mm
}

source('lcl.R')
source('raindrop re-evaporation model.R')

########################################################
# transport type
transport <- 1 # 1 = advection, 2 = eddy
raindrop <- 1 # 1 = w/ reevaporation effect, 2 = w/o

# temp and rh
ts <- seq(30+273.15,15+273.15,length.out=101) # surface temperature, list
rh <- seq(0.8,0.8,length.out=101) # relative humidity, list
tlcl <- vector() # condensation temperature, list
dm <- c(seq(1.0,1.0,length.out=101)) # drop diameter (mm), list

# fractionation functions
a18_P_lv <- function(temp) {  # function for alpha of precipitation from vapor to liquid phase (d18O) (>1)
  if (temp >= 273.15) {exp(1137/temp^2-0.4156/temp-0.00207)}
}
aD_P_lv <- function(temp){  # function for alpha of precipitation from vapor to liquid phase (dD) (>1)
  if (temp >= 273.15) {exp(24844/temp^2-76.248/temp+0.05261)}
}

# set up fractionation factors
alpha <- aD <- alpha_reevap_18 <- alpha_reevap_D <- vector()
for (i in 1:length(rh)){
  tlcl[i] <- lcl(101325,ts[i],rhl=rh[i])[2]
  alpha[i] <- a18_P_lv(tlcl[i]) # equilibrium fractionation 18O, >1, condenstation
  aD[i] <- aD_P_lv(tlcl[i]) # equilibrium fractionation D, >1, condenstation
  if (raindrop==1) {
    alpha_reevap_18[i] <- reevap(101325,ts[i],rh[i],dm[i],-12,-88)[2] # d18Ov = -12 permil; d2Hv = -88 permil; the "alpha" is not affected by vapor isotopic composition
    alpha_reevap_D[i] <- reevap(101325,ts[i],rh[i],dm[i],-12,-88)[3]
  } else {alpha_reevap_18[i] <- alpha_reevap_D[i] <- 1}
}
k18 <- 0.9723^0.8 # <1, evaporation
kD <- 0.9755^0.8 # <1, evaporation
smow18 <- 0.0020052
smowD <- 155.76*10^(-6)

# arrays
nd <- c(10e-20,0.5,1,2,5,10,20)
w <- array(dim = 101)
w <- w_S66(3,tlcl)
w <- w/w[1]

x <- 1-w # build array for distance (x')

# ET
e = seq(1,0,by=-0.1)
t = 1-e

# set up model parameters
delta_a <- de <- delta_inf <- delta_p <- delta_d_a <- de_d <- delta_d_inf <- delta_d_p <- 
  array(dim=c(length(w),length(nd), length(e)))
# d18O
delta_p[1,,] <- -2 # set initial d18O precip value
delta_a[1,,] <- (delta_p[1,,]+1000)/(alpha[1]*alpha_reevap_18[1])-1000 # set initial atmospheric vapor d18O
# dD
delta_d_p[1,,] <- 8*delta_p[1,,]+10 # set initial dD precip value
delta_d_a[1,,] <- (delta_d_p[1,,]+1000)/(aD[1]*alpha_reevap_D[1])-1000 # set initial atmospheric vapor dD
# expo term
expo <- expo_d <- array(dim=c(length(w),length(nd)))
for(j in 1:length(nd)){
  for(i in 1:length(w)){
    if (transport == 1){
      expo_d[i,j] <- (aD[i]*alpha_reevap_D[i]+nd[j]*aD[i]*alpha_reevap_D[i])-1 # advection only
      expo[i,j] <- (alpha[i]*alpha_reevap_18[i]+nd[j]*alpha[i]*alpha_reevap_18[i])-1
    } else {
      expo_d[i,j] <- (aD[i]*alpha_reevap_D[i]+nd[j]*aD[i]*alpha_reevap_D[i])^0.5-1 # eddy only
      expo[i,j] <- (alpha[i]*alpha_reevap_18[i]+nd[j]*alpha[i]*alpha_reevap_18[i])^0.5-1
    }
  }
}

# CG model
# function for alpha-evap
craig_gordon <- function(t, ts, h, delta_precip, d18 = T){
  if(d18 == T){
    r18_s <- (delta_precip/1000 + 1)*smow18
    alpha18 <- 1/a18_P_lv(ts)
    Re_18 <- (1-t)*((alpha18*k18*r18_s)/(1-h+(1-t)*k18*h)) + t*r18_s*(1/(1+(1-t)*k18*(h/(1-h))))
    Re_18/r18_s}
  else{
    rD_s <- (delta_precip/1000 + 1)*smowD
    alphaD <- 1/aD_P_lv(ts)
    Re_D <- (1-t)*((alphaD*kD*rD_s)/(1-h+(1-t)*kD*h)) + t*rD_s*(1/(1+(1-t)*kD*(h/(1-h))))
    Re_D/rD_s}
}
# integrated E in permil
isotope_e <- function(frac, delta_precip, ae){
  ((delta_precip+1000)*(1-frac)^ae-1000*(1-frac)-delta_precip)/(-frac)
}

# model iterations
for(h in 1:length(t)){
  for(j in 1:length(nd)){
    for(i in 2:length(w)){
      # calculate frac (ET/P) for isotope_e function
      frac <- 1/(1/nd[j]+1)
      # first calculate the combined equilibrium/kinetic fractionation using craig-gordon function
      # then calculate integrated ET value from location
      ae <- craig_gordon(t[h], ts[i-1], rh[i-1], delta_p[i-1,j,h])
      ae_d <- craig_gordon(t[h], ts[i-1], rh[i-1], delta_d_p[i-1,j,h], d18 = F)
      de[i-1,j,h] <- t[h]*delta_p[i-1,j,h] + e[h]*isotope_e(frac,delta_p[i-1,j,h],ae)
      de_d[i-1,j,h] <- t[h]*delta_d_p[i-1,j,h] + e[h]*isotope_e(frac,delta_d_p[i-1,j,h],ae_d)
      # Hendricks model
      # d18O
      delta_inf[i-1,j,h] <- (nd[j]*de[i-1,j,h]-(1+nd[j])*(alpha[i-1]*alpha_reevap_18[i-1]-1)*1000)/(alpha[i-1]*alpha_reevap_18[i-1]+nd[j]*alpha[i-1]*alpha_reevap_18[i-1]-1)
      delta_a[i,j,h] <- (delta_a[i-1,j,h]-delta_inf[i-1,j,h])*((w[i]/w[i-1])^expo[i-1,j])+delta_inf[i-1,j,h]
      delta_p[i,j,h] <- (delta_a[i,j,h]+1000)*alpha[i]*alpha_reevap_18[i]-1000
      # dD
      delta_d_inf[i-1,j,h] <- (nd[j]*de_d[i-1,j,h]-(1+nd[j])*(aD[i-1]*alpha_reevap_D[i-1]-1)*1000)/(aD[i-1]*alpha_reevap_D[i-1]+nd[j]*aD[i-1]*alpha_reevap_D[i-1]-1)
      delta_d_a[i,j,h] <- (delta_d_a[i-1,j,h]-delta_d_inf[i-1,j,h])*((w[i]/w[i-1])^expo_d[i-1,j])+delta_d_inf[i-1,j,h]
      delta_d_p[i,j,h] <- (delta_d_a[i,j,h]+1000)*aD[i]*alpha_reevap_D[i]-1000
    }
  }
}
dxs <- delta_d_p - 8*delta_p

# plot
colors = c("#CC79A7","#0072B2","#56B4E9","#009E73","#F0E442","#E69F00","#D55E00")

plot(-999,-999,xlim=c(-12,-1.5),ylim=c(8,20),xaxt="n",yaxt="n",xaxs="i",yaxs="i",
     xlab=expression(paste(delta^18,"O (","\u2030",")")),ylab=expression(paste("d-excess (","\u2030",")")),cex.lab=1.1)
axis(1,seq(-12,-2,2),cex.axis=1.1); axis(2,seq(8,20,2),cex.axis=1.1); title("Advection only, temperature gradient",font.main=1,adj=0,cex.main=1.1)

for(j in 1:length(nd)){
  lines(delta_p[,j,2],dxs[,j,2],col=colors[j],lty=5,lwd=1)
  lines(delta_p[,j,6],dxs[,j,6],col=colors[j],lty=1,lwd=2)
  lines(delta_p[,j,10],dxs[,j,10],col=colors[j],lty=4,lwd=1)
}

legend('topright', c('T/ET = 10%','T/ET = 50%','T/ET = 90%',
                     expression(paste(N[d]," = 0")),
                     expression(paste(N[d]," = 0.5")),
                     expression(paste(N[d]," = 1")),
                     expression(paste(N[d]," = 2")),
                     expression(paste(N[d]," = 5")),
                     expression(paste(N[d]," = 10")),
                     expression(paste(N[d]," = 20"))),
       lty = c(5,1,4,rep(1,7)),col = c('black','black','black',colors), lwd = c(1,2,1,rep(1,7)), cex=0.8)
