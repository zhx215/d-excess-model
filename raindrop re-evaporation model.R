# developed by ZX and MJW, 2021
# stand-alone function for calculating raindrop re-evaporation

setwd('C:\\Users\\Zhengyu\\Documents\\R\\github')
source('lcl.R')
reevap <- function(p0,temp0,rh0,D,d18Ov,d2Hv) {
  # function for alpha of precipitation from vapor to liquid (or ice) phase (d18O) (>1) (Majoube, 1971)
  a18_P_lv <- function(temp) {
    if (temp >= 273.15) {exp(1137/temp^2-0.4156/temp-0.00207)}
    else if (temp < 253.15) {exp(11.839/temp-0.028224)}
    else {
      lv <- exp(1137/temp^2-0.4156/temp-0.00207)
      iv <- exp(11.839/temp-0.028224)
      iv*(273.15-temp)/20+lv*(temp-253.15)/20 # aquilot mixing line
    } 
  }
  
  # function for alpha of precipitation from vapor to liquid (or ice) phase (dD) (>1) (Majoube, 1971)
  aD_P_lv <- function(temp){
    if (temp >= 273.15){exp(24844/temp^2-76.248/temp+0.05261)}
    else if (temp < 253.15){exp(16289/temp^2-0.0945)}
    else {
      lv <- exp(24844/temp^2-76.248/temp+0.05261)
      iv <- exp(16289/temp^2-0.0945)
      iv*(273.15-temp)/20+lv*(temp-253.15)/20 # aquilot mixing line
    }
  }
  
  # kinetic fractionation due to supersaturation over ice, d18O (Jouzel and Merlivat, 1984)
  a18_P_k <- function(temp){
    Si <- 1-0.004*(temp-273.15) # lambda = 0.004, Risi et al., 2010
    diff16_18 <- 1.0285 # (D/D') for 16O/18O
    a18_P_iv <- exp(11.839/temp-0.028224)
    if (temp >= 273.15) {1}
    else if (temp < 253.15) {Si/(a18_P_iv*diff16_18*(Si-1)+1)}
    else {Si/(a18_P_iv*diff16_18*(Si-1)+1)*(273.15-temp)/20+(temp-253.15)/20*1} # aquilot mixing line
  }
  
  # kinetic fractionation due to supersaturation over ice, dD (Jouzel and Merlivat, 1984)
  aD_P_k <- function(temp){
    Si <- 1-0.004*(temp-273.15) # lambda = 0.004, Risi et al., 2010
    diff1_2 <- 1.0251 # (D/D') for 1H/2H
    aD_P_iv <- exp(16289/temp^2-0.0945)
    if (temp >= 273.15) {1}
    else if (temp < 253.15) {Si/(aD_P_iv*diff1_2*(Si-1)+1)}
    else {Si/(aD_P_iv*diff1_2*(Si-1)+1)*(273.15-temp)/20+(temp-253.15)/20*1} # aquilot mixing line
  }
  
  # if RH = 1, no raindrop re-evaporation as LCL height = 0 m
  if (rh0 == 1){
    return(c(0,1,1,D))
  }
  
  # saturation vapor pressure over liquid (Romps, 2017)
  get_es <- function(temp){ # K
    Ptrip <- 611.65 # Pa
    Ttrip <- 273.16 # K
    E0v <- 2.374*10^6 # J/kg
    cvl <- 4119 # J/kg/K, specific heat capacity of liquid water
    cvv <- 1418 # J/kg/K, specific heat capacity of water vapor at constant volume
    Rv <- 461 # J/kg/K, specific gas constant of water vapor
    cpv <- cvv + Rv # J/kg/K, specific heat capacity of water vapor at constant pressure
    Ptrip*(temp/Ttrip)^((cpv-cvl)/Rv)*exp(((E0v-(cvv-cvl)*Ttrip)/Rv)*(1/Ttrip-1/temp)) # Pa
  }
  
  # pressure at a temperature of certain level, given initial P/T/RH
  get_p <- function(p0,temp0,rh0,t_z){
    Rv <- 461 # J/kg/K, specific gas constant of water vapor
    Ra <- 287.04 # J/kg/K, specific gas constant of dry air
    pv <- rh0*get_es(temp0) # Pa, partial pressure of water vapor
    qv <- Ra*pv/(Rv*p0+(Ra-Rv)*pv) # mass fraction of water vapor
    Rm <- (1-qv)*Ra+qv*Rv # J/kg/K, air parcel's specific gas constant
    cvv <- 1418 # J/kg/K, specific heat capacity of water vapor at constant volume
    cpv <- cvv + Rv # J/kg/K, specific heat capacity of water vapor at constant pressure
    cva <- 719 # J/kg/K, specific heat capacity of dry air at constant volume
    cpa <- cva + Ra # J/kg/K, heat capacity of dry air at constant pressure
    cpm <- (1-qv)*cpa + qv*cpv # J/kg/K, air parcel's specific heat capacity at constant pressure
    p0*(t_z/temp0)^(cpm/Rm)
  }
  
  # get pressure at a temperature of certain level, given initial P/T/RH
  get_rh <- function(p0,temp0,rh0,t_z){ 
    es_c <- get_es(lcl_data[2])
    p_c <- get_p(p0,temp0,rh0,lcl_data[2])
    ws_c <- 0.622*es_c/(p_c-es_c)
    es <- get_es(t_z)
    p <- get_p(p0,temp0,rh0,t_z)
    ws <-0.622*es/(p-es)
    es0 <- get_es(temp0)
    ws0 <- 0.622*es0/(p0-es0)
    w <- ws0*rh0
    rh <- w/ws
    rh_c <- w/ws_c
    (rh-rh0)/(rh_c-rh0)*(1-rh0)+rh0 # normalize final RH to 0-1 scale
  }
  
  # get effective RH for raindrops
  get_rheff <- function(temp,temp_d,rh){
    rheff <- get_es(temp)*rh/get_es(temp_d)
    if (rheff>1) {
      return(1) # set rheff as 1 if > 1
    } else {
      return(rheff)
    }
  }
  
  # air density, Salamalikis et al., 2016, Plathner and Woloszyn, 2002
  rho <- function(p,temp,rh,es){
    Rd <- 287.04 # J/kg/K, specific gas constant of dry air
    1/(Rd*temp)*(p-0.378*rh*es) # kg/m3
  }
  
  # terminal velocity, Foote and Du Toit, 1969
  tv <- function(d,rho,rho_ref){
    9.43*(1-exp(-(d/1.77)^1.147))*(rho_ref/rho)^0.4 # m/s 
  }
  
  # dynamic viscosity of air, Sutherland, 1893
  visc <- function(temp){ 
    1.458*10^(-6)*temp^1.5/(temp+110.4) # kg/m/s
  }
  
  # diffusivity of water molecules, Hall and Pruppacher, 1976
  diff <- function(p,temp){
    p0 <- 101325
    0.211*10^-4*(p0/p)*(temp/273.15)^1.94 # m2/s
  }
  
  # water density, from drop temperature, slight offset from 1000 kg/m3 or 1 g/cm3, Salamalikis et al., 2016
  rhow <- function(temp_d){ 
    t <- temp_d-273.15
    1000*(1-(t+288.9414)*(t-3.9863)^2/508929.2/(t+68.12963)) # kg/m3
  }
  
  # latent heat of evaporation, from drop temperature, Rogers and Yau, 1989
  Le <- function(temp_d){ 
    t <- temp_d-273.15
    2500800-2360*t+1.6*t^2-0.06*t^3 # J/kg
  }
  
  # thermal conductivity of air, ka is close to k, Pruppacher and Klett, 2010
  ka <- function(temp){
    418.4*(5.69+0.017*(temp-273.15))*10^-5 # J/m/s/K, the original form is in unit cal/cm/s/K (1 cal = 4.184 J)
  }
  
  # initalize the model
  lcl_data <- lcl(p0,temp0,rhl=rh0) # get lcl data
  gamma <- ((temp0-lcl_data[2])/lcl_data[1]) # K/m, lapse rate
  z_p <- lcl_data[1] # height profile
  D_p <- D # drop diameter profile
  m_p <- 4/3*pi*(D_p/2)^3*10^-3 # g, mass profile
  t_p <- td_p <- lcl_data[2] # temperature and drop temperature profile
  p_p <- get_p(p0,temp0,rh0,t_p)# pressure profile
  rh_p <- get_rh(p0,temp0,rh0,t_p) # rh profile
  rheff_p <- get_rheff(t_p,td_p,rh_p) # effective rh profile
  rho_p <- rho(p_p,t_p,rh_p,get_es(t_p)) # air density profile
  rho_ref <- rho(p0,temp0,rh0,get_es(temp0)) # reference air density
  tv_p <- tv(D_p,rho_p,rho_ref) # vertical velocity profile
  visc_p <- visc(t_p) # viscosity profile
  diff_p <- diff(p_p,t_p) # diffusivity profile
  rhow_p <- rhow(td_p) # water density profile
  Le_p <- Le(td_p) # latent heat profile
  ka_p <- ka(t_p) # conductivity profile
  Sc_p <- visc_p/(rho_p*diff_p) # Schmidt number profile
  cp <- 1004 # J/kg/K specific heat of dry air
  Pr_p <- visc_p*cp/ka_p # Prandtl number profile
  Re_p <- tv_p*(D_p/1000)*rho_p/visc_p # Reynolds number
  fv_p <- fh_p <- vector() # mass and heat ventilation coefficient
  fv_p18 <- vector() # mass ventilation coefficient for 18O
  fv_p2 <- vector() # mass ventilation coefficient for 2H
  Rw <- 461 # J/kg/K, specific gas constant of water vapor
  cw <- 4119 # J/kg/K specific heat capacity of liquid water, Romps 2017
  diff16_18 <- 1.0285 # (D/D') for 16O/18O
  diff1_2 <- 1.0251 # (D/D') for 1H/2H
  rstd18 <- 2005.2*10^-6 # VSMOW
  rstd2 <- 155.76*10^-6 # VSMOW
  rv18 <- (d18Ov/1000+1)*rstd18 # get r value for 18O/16O in vapor
  rv2 <- (d2Hv/1000+1)*rstd2 # get r value for 2H/1H in vapor
  rp18_p <- rv18*a18_P_lv(td_p)*a18_P_k(td_p) # get r value for 18O/16O in raindrop, initial value
  rp2_p <- rv2*aD_P_lv(td_p)*aD_P_k(td_p) # get r value for 2H/1H in raindrop, initial value
  m18_p <- m_p*rp18_p # get 18O mass
  m2_p <- m_p*rp2_p # get 2H mass
  i <- 1
  
  # iterate the model
  repeat{
    # mass ventilation coefficient
    if ((Sc_p[i]^(1/3)*Re_p[i]^0.5)>=1.4){
      fv_p[i] <- 0.78+0.308*Sc_p[i]^(1/3)*Re_p[i]^0.5
      fv_p18[i] <- 0.78+0.308*(Sc_p[i]*(diff16_18))^(1/3)*Re_p[i]^0.5
      fv_p2[i] <- 0.78+0.308*(Sc_p[i]*(diff1_2))^(1/3)*Re_p[i]^0.5
    } else{
      fv_p[i] <- 1+0.108*(Sc_p[i]^(1/3)*Re_p[i]^0.5)^2
      fv_p18[i] <- 1+0.108*((Sc_p[i]*(diff16_18))^(1/3)*Re_p[i]^0.5)^2
      fv_p2[i] <- 1+0.108*((Sc_p[i]*(diff1_2))^(1/3)*Re_p[i]^0.5)^2
    }
    
    # heat ventilation coefficient
    if ((Pr_p[i]^(1/3)*Re_p[i]^0.5)>=1.4){
      fh_p[i] <- 0.78+0.308*Pr_p[i]^(1/3)*Re_p[i]^0.5
    } else{
      fh_p[i] <- 1+0.108*(Pr_p[i]^(1/3)*Re_p[i]^0.5)^2
    }
    
    # mass and heat transfer equations
    if (t_p[i]<=273.15){ # if surrounding air is frozen, no evaporation
      e_rate <- e_rate_18 <- e_rate_2 <- 0
    } else {
      e_rate <- 2*pi*(D_p[i]/1000)*fv_p[i]*diff_p[i]/Rw*(rh_p[i]*get_es(t_p[i])/t_p[i]-get_es(td_p[i])/td_p[i]) # kg/s, evaporation rate
      e_rate_18 <- 2*pi*(D_p[i]/1000)*fv_p[i]*diff_p[i]/Rw*(fv_p18[i]/fv_p[i]/diff16_18)^0.58*
        (rh_p[i]*rv18*get_es(t_p[i])/t_p[i]-rp18_p[i]/a18_P_lv(td_p[i])*get_es(td_p[i])/td_p[i]) 
      e_rate_2 <- 2*pi*(D_p[i]/1000)*fv_p[i]*diff_p[i]/Rw*(fv_p2[i]/fv_p[i]/diff1_2)^0.58*
        (rh_p[i]*rv2*get_es(t_p[i])/t_p[i]-rp2_p[i]/aD_P_lv(td_p[i])*get_es(td_p[i])/td_p[i])
    }
    td_rate <- 12/((D_p[i]/1000)^2*rhow_p[i]*cw)*
      (Le_p[i]*fv_p[i]*diff_p[i]/Rw*(rh_p[i]*get_es(t_p[i])/t_p[i]-get_es(td_p[i])/td_p[i])-fh_p[i]*ka_p[i]*(td_p[i]-t_p[i])) # K/s
    
    # mass and heat changes
    t <- 0.1 # s
    de <- (e_rate*1000)*t # g
    de18 <- (e_rate_18*1000)*t
    de2 <- (e_rate_2*1000)*t
    dt <- td_rate*t # K
    
    # return results if reaching the ground
    if((z_p[i]-tv_p[i]*t)<0){
      print('hit ground')
      return(c(100-100*min(m_p)/max(m_p),tail(rp18_p,n=1)/rp18_p[1],tail(rp2_p,n=1)/rp2_p[1],tail(D_p,n=1))) # evaporation amount (%), alpha_rr for 18O, alpha_rr for 2H, final raindrop diameter (mm)
    }
    if(abs(dt)>10^(-2)){
      print('all gone (numerical model drift)')
      return(c(100,1,1,0))
    }
    if((m_p[i]+de)<0){
      print('all gone (ideally)')
      return(c(100,1,1,0))
    }
    
    # update model parameters
    z_p[i+1] <- z_p[i]-tv_p[i]*t
    m_p[i+1] <- m_p[i]+de
    D_p[i+1] <- (6*m_p[i+1]*10^3/pi)^(1/3)
    t_p[i+1] <- temp0-gamma*z_p[i+1]
    td_p[i+1] <- td_p[i]+dt
    p_p[i+1] <- get_p(p0,temp0,rh0,t_p[i+1])
    rh_p[i+1] <- get_rh(p0,temp0,rh0,t_p[i+1])
    rheff_p[i+1] <- get_rheff(t_p[i+1],td_p[i+1],rh_p[i+1])
    rho_p[i+1] <- rho(p_p[i+1],t_p[i+1],rh_p[i+1],get_es(t_p[i+1]))
    tv_p[i+1] <- tv(D_p[i+1],rho_p[i+1],rho_ref)
    visc_p[i+1] <- visc(t_p[i+1])
    diff_p[i+1] <- diff(p_p[i+1],t_p[i+1])
    rhow_p[i+1] <- rhow(td_p[i+1])
    Le_p[i+1] <- Le(td_p[i+1])
    ka_p[i+1] <- ka(t_p[i+1])
    Sc_p[i+1] <- visc_p[i+1]/(rho_p[i+1]*diff_p[i+1])
    Pr_p[i+1] <- visc_p[i+1]*cp/ka_p[i+1]
    Re_p[i+1] <- tv_p[i+1]*(D_p[i+1]/1000)*rho_p[i+1]/visc_p[i+1]
    m18_p[i+1] <- m18_p[i]+de18
    m2_p[i+1] <- m2_p[i]+de2
    rp18_p[i+1] <- m18_p[i+1]/m_p[i+1]
    rp2_p[i+1] <- m2_p[i+1]/(m_p[i+1])
    i <- i+1
  }
}