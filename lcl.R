# revised by ZX and MJW, 2021

# Version 1.0 released by David Romps on September 12, 2017.
# 
# When using this code, please cite:
# 
# @article{16lcl,
#   Title   = {Exact expression for the lifting condensation level},
#   Author  = {David M. Romps},
#   Journal = {Journal of the Atmospheric Sciences},
#   Year    = {2017},
#   Volume  = {in press},
# }
#
# This lcl function returns the height of the lifting condensation level
# (LCL) in meters.  The inputs are:
# - p in Pascals
# - T in Kelvins
# - Exactly one of rh, rhl, and rhs (dimensionless, from 0 to 1):
#    * The value of rh is interpreted to be the relative humidity with
#      respect to liquid water if T >= 273.15 K and with respect to ice if
#      T < 273.15 K. 
#    * The value of rhl is interpreted to be the relative humidity with
#      respect to liquid water
#    * The value of rhs is interpreted to be the relative humidity with
#      respect to ice
# - ldl is an optional logical flag.  If true, the lifting deposition
#   level (LDL) is returned instead of the LCL. 
# - min_lcl_ldl is an optional logical flag.  If true, the minimum of the
#   LCL and LDL is returned.

library(LambertW)

lcl <- function(p,Tx,rh=NULL,rhl=NULL,rhs=NULL,return_ldl=FALSE,return_min_lcl_ldl=FALSE) {
  
  # Parameters
  Ttrip <- 273.16     # K
  ptrip <- 611.65     # Pa
  E0v   <- 2.3740e6   # J/kg
  E0s   <- 0.3337e6   # J/kg
  ggr   <- 9.81       # m/s^2
  rgasa <- 287.04     # J/kg/K 
  rgasv <- 461        # J/kg/K 
  cva   <- 719        # J/kg/K
  cvv   <- 1418       # J/kg/K 
  cvl   <- 4119       # J/kg/K 
  cvs   <- 1861       # J/kg/K 
  cpa   <- cva + rgasa
  cpv   <- cvv + rgasv
  
  # The saturation vapor pressure over liquid water
  pvstarl <- function(Tx) {
    return( ptrip * (Tx/Ttrip)^((cpv-cvl)/rgasv) *
              exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1/Ttrip - 1/Tx) ) )
  }
  
  # The saturation vapor pressure over solid ice
  pvstars <- function(Tx) {
    return( ptrip * (Tx/Ttrip)^((cpv-cvs)/rgasv) *
              exp( (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1/Ttrip - 1/Tx) ) )
  }
  
  # Calculate pv from rh, rhl, or rhs
  rh_counter <- 0
  if (!is.null(rh )) { rh_counter <- rh_counter + 1 }
  if (!is.null(rhl)) { rh_counter <- rh_counter + 1 }
  if (!is.null(rhs)) { rh_counter <- rh_counter + 1 }
  if (rh_counter != 1) {
    stop('Error in lcl: Exactly one of rh, rhl, and rhs must be specified')
  }
  if (!is.null(rh)) {
    # The variable rh is assumed to be 
    # with respect to liquid if T > Ttrip and 
    # with respect to solid if T < Ttrip
    if (Tx > Ttrip) {
      pv <- rh * pvstarl(Tx)
    } else {
      pv <- rh * pvstars(Tx)
    }
    rhl <- pv / pvstarl(Tx)
    rhs <- pv / pvstars(Tx)
  } else if (!is.null(rhl)) {
    pv <- rhl * pvstarl(Tx)
    rhs <- pv / pvstars(Tx)
    if (Tx > Ttrip) {
      rh <- rhl
    } else {
      rh <- rhs
    }
  } else if (!is.null(rhs)) {
    pv <- rhs * pvstars(Tx)
    rhl <- pv / pvstarl(Tx)
    if (Tx > Ttrip) {
      rh <- rhl
    } else {
      rh <- rhs
    }
  }
  if (pv > p) {
    return(NA)
  }
  
  # Calculate lcl and ldl
  qv <- rgasa*pv / (rgasv*p + (rgasa-rgasv)*pv)
  rgasm <- (1-qv)*rgasa + qv*rgasv
  cpm <- (1-qv)*cpa + qv*cpv
  if (rh==0) {
    return(cpm*Tx/ggr)
  }
  al  <- -(cpv-cvl)/rgasv + cpm/rgasm
  bl  <- -(E0v-(cvv-cvl)*Ttrip)/(rgasv*Tx)
  cl  <- pv/pvstarl(Tx)*exp(-(E0v-(cvv-cvl)*Ttrip)/(rgasv*Tx))
  as  <- -(cpv-cvs)/rgasv + cpm/rgasm
  bs  <- -(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*Tx)
  cs  <- pv/pvstars(Tx)*exp(-(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*Tx))
  lcl <- cpm*Tx/ggr*( 1 - bl/(al*W(bl/al*cl^(1/al),-1)) )
  tlcl <- Tx-lcl*ggr/cpm ################ lcl temperature!!!
  ldl <- cpm*Tx/ggr*( 1 - bs/(as*W(bs/as*cs^(1/as),-1)) )
  tldl <- Tx-ldl*ggr/cpm ################ ldl temperature!!!
  
  # Return either lcl or ldl
  if (return_ldl & return_min_lcl_ldl) {
    stop('return_ldl and return_min_lcl_ldl cannot both be true')
  } else if (return_ldl) {
    return(c(ldl,tldl)) ###
  } else if (return_min_lcl_ldl) {
    return(c(min(lcl,ldl),max(tlcl,tldl))) ###
  } else {
    return(c(lcl,tlcl)) ### print lcl and lcl temp
  }
  
}