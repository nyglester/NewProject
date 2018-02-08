# ThermoPic_P1_Model
# Version: V3b (NPL: 2018-01-30)
# Combines Pgm1,2 and 3 (of Minns/Fung original code)

#------------------------------------------------------------------------------
# Builds Seasonal temperature model using lake and climate data to predicts STM parameters
# Reads Lake data from CSV files: 
#              0_User_Options  (Processing options - curretly not used for any decisisons)
#              1_Lake          (Lake data,    Primary key = Wby_Lid)
#              2_Climate       (Climate data, Primary key = Wby_Lid x Period)
# Calculates STM Parameters for each record in 2_Climate (i.e. Wby_Lid x Period combination)
# Climate records can be historical norms or future scenarios
# 'Period' variable is a text descriptor (e.g., 2001-2010, Predict_2100)
# Program ignores leap years because it assumes climate is for a peirod norm (#days in February = 28.25)
#------------------------------------------------------------------------------

# Set work directory to the ThermoPic root folder
root<-"C:/R_ThermoPic_v3"
setwd(root)
getwd()

#----------------------------------------------------------------------------
# Set file names
#----------------------------------------------------------------------------
file_in0 <- "DataIn/0_User_Options.csv"
file_in1 <- "DataIn/1_Lake.csv"
file_in2 <- "DataIn/2_Climate.csv"

file_out0<- "DataOut/0_User_Options.csv"
file_out1<- "DataOut/tmp_ClimMetrics.csv"
file_out2<- "DataOut/tmp_IceClimMetrics.csv"
file_out3<- "DataOut/3_Model_Inputs.csv"
file_out4<- "DataOut/4_STM_Parameters.csv"

#------------------------------------------------------------------------------
# Simplify code by assigning "year_fix" for calculating Day of Year and Solar elevation/angle
# Exact year not needed because because program works with Temperature norms for Specified Period
# Solar angle for a given day of year varies on 4 year cycle, but variation is very small
#------------------------------------------------------------------------------
year_fix <- 2010

#------------------------------------------------------------------------------
# Load packages required for analysis
#------------------------------------------------------------------------------
library(chron)
library(plyr)
library(XLConnect)
library(rLakeAnalyzer)

#------------------------------------------------------------------------------------
# INITIATE FUNCTION ROUTINES
#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
# Function = tdatev
#------------------
# Interpolate Julian dates for reference temperatures in spring and fall from monthly mean temperatures
# tvec - data frame with mean monthly temperatures from Dec (Y-1) through Jan-Dec (Y) to Jan(Y+1)
# tref - threshold temperature C
# tdir - value indicates whether the temperature is rising (1) or falling (-1)
# yy   - calendar year (allows day counts to be adjusted for leap years)
#-------------------------------------------------------------------------------------
tdatev <- function(tvec,tref,tdir,yy)
{
  tstore <- data.frame(postn = c(1:14),
                       dinm = c(31,31,28.25,31,30,31,30,31,31,30,31,30,31,31),
                       dcum = c(rep(0,14)), hdinm = c(rep(0,14)),hdcum = c(rep(0,14)))
  for (j in 2:14) 
  { 
    tstore$dcum[j] <- tstore$dcum[j-1] + tstore$dinm[j]
    tstore$hdcum[j] <- (tstore$dcum[j]+tstore$dcum[j-1])/2
  }
  tstore$hdcum[1] <- -tstore$dinm[1]/2
  tstore$value <- -99
  if (tdir > 0)
  {
    for (m in 2:13) 
    {
      mm <- m-1
      t1 <- tvec$tmp[mm];t2 <- tvec$tmp[m];d1 <- tstore$hdcum[mm];d2 <- tstore$hdcum[m]
      tstore$value[m] <- ifelse(t1 < tref,ifelse(t2 >= tref,d1+(d2-d1)*(tref-t1)/(t2-t1),-99),-99)
    }
  }
  else
  {
    for (m in 2:13) 
    {
      mp <- m+1
      t1 <- tvec$tmp[m];t2 <- tvec$tmp[mp];d1 <- tstore$hdcum[m];d2 <- tstore$hdcum[mp]
      tstore$value[m] <- ifelse(t1 > tref,ifelse(t2 <= tref,d1+(d2-d1)*(t1-tref)/(t1-t2),-99),-99)
    }
  }
  tdatev <- max(tstore$value)
}


#--------------------------------------------------------------------------------
# Function: DJTmps
#------------------
# Estimate Dec Temp in Year-1 from Jan in Year 
#           and Jan Temp in Y+1 from Dec in Y
#--------------------------------------------------------------------------------
DJTmps <- function(dSlice) 
{
  YRMin <- min(dSlice$period)
  YRMax <- max(dSlice$period)
  Nyr <- length(dSlice$period)
  DJd <-data.frame(TJ=dSlice$Tjan[2:Nyr],TD=dSlice$Tdec[1:(Nyr-1)])
  
  # Given J predict D with J0 and J1
  JtoD.mod <- lm(TD~TJ,DJd)
  J0 <- coef(JtoD.mod)[1]
  J1 <- coef(JtoD.mod)[2]
  
  # Given D predict J with D0 and D1
  DtoJ.mod <- lm(TJ~TD,DJd)
  D0 <- coef(DtoJ.mod)[1]
  D1 <- coef(DtoJ.mod)[2]
  return(data.frame(D0,D1,J0,J1,YRMin,YRMax))
}

#--------------------------------------------------------------------------------
# Function: SunPosition
#------------------
# Get Solar elevation and azimuth for specified date and location
# Implementation off the Astronomer's Almanac algorithm 
# (Joseph J. Michalsky. The astronomical almanacs algorithm for approximate solar position (19502050). 
# Solar Energy, 40(3):227235, 1988.) 
# It's straightforward to translate it in any other language:
# http://stackoverflow.com/questions/257717/position-of-the-sun-given-time-of-day-and-lat-long
#
#------------------------------------------------------------------
#
# http://www.esrl.noaa.gov/gmd/grad/solcalc/
# Checked against NOAA on-line calculator
#
#------------------------------------------------------------------
#
# Solar position routine to compute Solar Azimuth and Elevation
#
SunPosition <- function(year, month, day, hour=12, min=0, sec=0,lat,long) 
{
  twopi <- 2 * pi
  deg2rad <- pi / 180
  
  # Get day of the year, e.g. Feb 1 = 32, Mar 1 = 61 on leap years
  month.days <- c(0,31,28,31,30,31,30,31,31,30,31,30)
  day <- day + cumsum(month.days)[month]
  leapdays <- year %% 4 == 0 & (year %% 400 == 0 | year %% 100 != 0) & day >= 60
  day[leapdays] <- day[leapdays] + 1
  
  # Get Julian date - 2400000
  hour <- hour + min / 60 + sec / 3600 # hour plus fraction
  delta <- year - 1949
  leap <- trunc(delta / 4) # former leapyears
  jd <- 32916.5 + delta * 365 + leap + day + hour / 24
  
  # The input to the Atronomer's almanach is the difference between
  # the Julian date and JD 2451545.0 (noon, 1 January 2000)
  time <- jd - 51545.
  
  # Ecliptic coordinates
  
  # Mean longitude
  mnlong <- 280.460 + .9856474 * time
  mnlong <- mnlong %% 360
  mnlong[mnlong < 0] <- mnlong[mnlong < 0] + 360
  
  # Mean anomaly
  mnanom <- 357.528 + .9856003 * time
  mnanom <- mnanom %% 360
  mnanom[mnanom < 0] <- mnanom[mnanom < 0] + 360
  mnanom <- mnanom * deg2rad
  
  # Ecliptic longitude and obliquity of ecliptic
  eclong <- mnlong + 1.915 * sin(mnanom) + 0.020 * sin(2 * mnanom)
  eclong <- eclong %% 360
  eclong[eclong < 0] <- eclong[eclong < 0] + 360
  oblqec <- 23.429 - 0.0000004 * time
  eclong <- eclong * deg2rad
  oblqec <- oblqec * deg2rad
  
  # Celestial coordinates
  # Right ascension and declination
  num <- cos(oblqec) * sin(eclong)
  den <- cos(eclong)
  ra <- atan(num / den)
  ra[den < 0] <- ra[den < 0] + pi
  ra[den >= 0 & num < 0] <- ra[den >= 0 & num < 0] + twopi
  dec <- asin(sin(oblqec) * sin(eclong))
  
  # Local coordinates
  # Greenwich mean sidereal time
  gmst <- 6.697375 + .0657098242 * time + hour
  gmst <- gmst %% 24
  gmst[gmst < 0] <- gmst[gmst < 0] + 24.
  
  # Local mean sidereal time
  lmst <- gmst + long / 15.
  lmst <- lmst %% 24.
  lmst[lmst < 0] <- lmst[lmst < 0] + 24.
  lmst <- lmst * 15. * deg2rad
  
  # Hour angle
  ha <- lmst - ra
  ha[ha < -pi] <- ha[ha < -pi] + twopi
  ha[ha > pi] <- ha[ha > pi] - twopi
  
  # Latitude to radians
  lat <- lat * deg2rad
  
  # Azimuth and elevation
  el <- asin(sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(ha))
  az <- asin(-cos(dec) * sin(ha) / cos(el))
  elc <- asin(sin(dec) / sin(lat))
  az[el >= elc] <- pi - az[el >= elc]
  az[el <= elc & ha > 0] <- az[el <= elc & ha > 0] + twopi
  
  el <- el / deg2rad
  az <- az / deg2rad
  lat <- lat / deg2rad
  
  return(list(elevation=el, azimuth=az))
}

#--------------------------------------------------------------------------------
# Function: FindArea
#------------------

FindArea <- function(area,zmn,zmx,ztarget) {
  zrat <- zmn/zmx
  rli <- (1-zrat)/zrat
  zarea <- area*(1-ztarget/zmx)^rli
  return(zarea)
}

#----------------------------------------------------------------------------
# Function: stm
# Get temperature at given depth for given day in a lake
#----------------------------------------------------------------------------
#
stm <- function(tx,tn,js,jm,je,zm,zj,sp,jd,z)
{
  zt <- ifelse( jd > js, zm*(jd-js)/(zj+jd-js),0.0)
  pz <- ifelse( zt > 0.0,zt^sp/(zt^sp+z^sp),1.0)# Proportional scaling with depth
  jx <- (je*(jm-js)+pz*js*(je-jm))/((jm-js)+pz*(je-jm))# Julian date of peak temperature at depth z
  tp <- tn+(tx-tn)*(je-jx)/(je-jm)# Peak temperature at depth z
  spr <- (tp - tn)/(jx - js) # daily warming rate at surface in spring
  fall <- (tx - tn)/(je - jm)# daily cooling rate at all depth after peak in fall
  tmoda <- tp -spr*(jd <= jx)*(jx - jd) - fall*(jd > jx)*(jd-jx)# predicted temperature given jd and z
  tmodb <- tx - (tx - tn)*(jm-jd)/(jm-js)
  tmodc <- tx - (tx - tn)*(jd-jm)/(je-jm)
  stm <- ifelse( jd <= js,tmodb,ifelse( jd >= je,tmodc,tmoda))
}

#----------------------------------------------------------------------------
# End of Function Routines
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# I/O Files information
#----------------------------------------------------------------------------

UserOptions  <- read.csv(file_in0,header=TRUE)
head(UserOptions)

laked  <- read.csv(file_in1,header=TRUE)
head(laked)

sited <- read.csv(file_in2,header=TRUE)
head(sited)
Nsites <- length(sited$Period)

#----------------------------------------------------------------------------
# Check Lake Data for missing values (Shoreline, Secchi, DOC)
#----------------------------------------------------------------------------

# Estimate Shoreline
# OLD FORMULA: laked$elnShoreline <- 1.32844+0.67375*log(laked$Area_ha/100)+0.75298*log(laked$Depth_Max)-0.72932*log(laked$Depth_Mn)+0.3717^2/2
# NEW FORMuLA (next line) based on 721 BsM lakes (in Ontario)
laked$elnShoreline <- 1.39494+0.69596*log(laked$Area_ha/100)+0.65994*log(laked$Depth_Max)-0.68779*log(laked$Depth_Mn)+0.3772^2/2
# Alternative FORMuLA (next line) based on 9283 AHI lakes (in Ontario)
#laked$elnShoreline <- 1.6969+0.67061*log(laked$Area_ha/100)+0.38243*log(laked$Depth_Max)-0.34590*log(laked$Depth_Mn)+0.3310^2/2

laked$eShoreline <- exp(laked$elnShoreline)
laked$ShorelineFlag <- "Yes"
laked$ShorelineFlag <- ifelse(is.na(laked$Shoreline),"No",laked$ShorelineFlag)
laked$Shoreline <- ifelse(laked$ShorelineFlag == "No",laked$eShoreline,laked$Shoreline)
laked$Shoreline <- round(laked$Shoreline,digits=2)

# Estimate Secchi
laked$Secchi <- as.numeric(laked$Secchi)
# OLD FORMULAE in next 2 lines
#laked$elnSecchi <- 1.822791-0.166579*laked$DOC+0.003760*laked$DOC^2+0.186572*log(laked$Depth_Mn)+0.3324^2/2
#laked$elnSecchi <- ifelse(is.na(laked$DOC),(-77.133986-0.039255*log(laked$Area_ha/100)+3.3899*laked$Latitude-0.036865*laked$Latitude^2+0.370819*log(laked$Depth_Mn)+0.4032^2/2),(1.822791-0.166579*laked$DOC+0.003760*laked$DOC^2+0.186572*log(laked$Depth_Mn)+0.4032^2/2))
# NEW FORMULAE in next 2 lines (based on 695 bsM lakes with DOC and 717 BsM lakes without DOC)
laked$elnSecchi_DOC   <- 2.19226-0.720734*log(laked$DOC)+0.179978*log(laked$Depth_Max)+0.3617^2/2
laked$elnSecchi_noDOC <- 0.277336-0.08618*log(laked$Area_ha/100)+0.11680*log(laked$Depth_Mn)+0.36719*log(laked$Depth_Max)+0.4582^2/2
laked$elnSecchi <- ifelse(is.na(laked$DOC),laked$elnSecchi_noDOC,laked$elnSecchi_DOC)

laked$eSecchi <- exp(laked$elnSecchi)
laked$SecchiFlag <- "Yes"
laked$SecchiFlag <- ifelse(is.na(laked$Secchi),"No",laked$SecchiFlag)
laked$Secchi <- ifelse(laked$SecchiFlag == "No",laked$eSecchi,laked$Secchi)
laked$Secchi <- round(laked$Secchi,digits=2)

# Estimate DOC
# OLD FORMULA in next line - new formula proposed by NPL avoids negative DOC when Secchi is small (<0.5 m)
#laked$eDOC <- 2.0409-0.3009*log(laked$Depth_Mn)+18.7373*1/laked$Secchi-8.2407*(1/laked$Secchi)^2
# NEW FORMULA in next line (avoids negative DOC)
laked$eDOC <- exp(2.7215 - 0.1276*log(laked$Secchi) - 0.1469*laked$Secchi - 0.0734*log(laked$Depth_Max)+0.3070^2)
laked$DOCFlag <- "Yes"
laked$DOCFlag <- ifelse(is.na(laked$DOC),"No",laked$DOCFlag)
laked$DOC <- ifelse(laked$DOCFlag == "No",laked$eDOC,laked$DOC)
laked$DOC <- round(laked$DOC,digits=2)

#------------------------------------------------------------------
# Do calculations on Climate data
#------------------------------------------------------------------
sited$Tson <- sited$Tjja <- sited$Tmam <- sited$Tdjf <- NA
sited$J_Aut4 <- sited$J_Aut0 <- sited$J_Spr0 <- -99
sited$T_Aut0 <- sited$T_Spr0 <- -99
#sited$T_Spr0 <- as.numeric(sited$T_Spr0)
#sited$T_Aut0 <- as.numeric(sited$T_Aut0)

# Compute temp seasonal means
sited$Tdjf <- round((31*sited$Tdec+31*sited$Tjan+(28.25)*sited$Tfeb)/(90.25),digits=2)
sited$Tmam <- round((31*sited$Tmar+30*sited$Tapr+31*sited$Tmay)/92,digits=2)
sited$Tjja <- round((30*sited$Tjun+31*sited$Tjul+31*sited$Taug)/92,digits=2)
sited$Tson <- round((30*sited$Tsep+31*sited$Toct+30*sited$Tnov)/91,digits=2)

head(sited)

# Compute Julian day in Spring and Autumn for Zero temperature metrics (J_Spr0 and J_Aut0)
# Must Skip first 2 variables (Wby_LID, Period) when reading monthy mean temperatures (col 3 is Tjan)
# Ignores leap year issues by using using "year_fix" (e.g. 2010) - because program works with Temp norms for given Period
Ncskip <- 2
yy <- year_fix
tvec <- data.frame(tmp = c(rep(-99,14)))
for( y in 1:Nsites) 
{
  for (m in 1:12) 
  { 
    tvec$tmp[m+1] <- sited[y,m+Ncskip]
  }
  tvec$tmp[1] <- sited$Tdec[y]
  tvec$tmp[14] <- sited$Tjan[y]
  
  # Get Spring 0 C
  tref <- 0.0
  tdir <- 1
  jday <- round(tdatev(tvec,tref,tdir,yy),digits=1)
  mdy <- month.day.year(jday,c(month=12,day=31,year=year_fix-1))
  sited$J_Spr0[y] <- jday
  sited$YY_Spr0[y] <- mdy$year
  sited$MM_Spr0[y] <- mdy$month
  sited$DD_Spr0[y] <- mdy$day
  
  # Get Autumn 0 C
  tref <- 0.0
  tdir <- -1
  jday <- round(tdatev(tvec,tref,tdir,yy),digits=1)
  mdy  <- month.day.year(jday,c(month=12,day=31,year=year_fix-1))
  sited$J_Aut0[y] <- jday
  sited$YY_Aut0[y] <- mdy$year
  sited$MM_Aut0[y] <- mdy$month
  sited$DD_Aut0[y] <- mdy$day
  
  # Get utumn 4 C - do we need this
  tref <- 4.0
  tdir <- -1
  sited$J_Aut4[y] <- round(tdatev(tvec,tref,tdir,yy),digits=1)
}
sited$A40CR <- round(4.0/(sited$J_Aut0 - sited$J_Aut4),digits=3)

# Compute quarter mean temp for J_Aut0 and J_Spr0 months

nm <- data.frame(ct=c(31,28.25,31,30,31,30,31,31,30,31,30,31))

ncs <- Ncskip
for (i in 1:Nsites) 
{
  isp2 <- sited$MM_Spr0[i]
  isp1 <- ifelse(isp2 < 2,12,isp2 - 1)
  isp3 <- ifelse(isp2 > 11,1,isp2 + 1)
  nsp1 <- nm$ct[isp1];nsp2 <- nm$ct[isp2];nsp3 <- nm$ct[isp3]
  ntot <- nsp1+nsp2+nsp3;ntot
  sited$T_Spr0[i] <- round((nsp1*sited[i,ncs+isp1]+nsp2*sited[i,ncs+isp2]+nsp3*sited[i,ncs+isp3])/ntot,digits=2)
  iau2 <- sited$MM_Aut0[i]
  iau1 <- ifelse(iau2 < 2,12,iau2 - 1)
  iau3 <- ifelse(iau2 > 11,1,iau2 + 1)
  #nau1 <- nm$ct[isp1];nau2 <- nm$ct[isp2];nau3 <- nm$ct[isp3]
  nau1 <- nm$ct[iau1]
  nau2 <- nm$ct[iau2]
  nau3 <- nm$ct[iau3]
  ntot <- nau1+nau2+nau3
  sited$T_Aut0[i] <- round((nau1*sited[i,ncs+iau1]+nau2*sited[i,ncs+iau2]+nau3*sited[i,ncs+iau3])/ntot,digits=2)
}

#------------------------------------------------------------------
# CREATE A WORKSHEET WITH JUST THE ICE METRICS
#------------------------------------------------------------------
siteo <- sited[ ,c("Wby_Lid","Period","Tann","Tjja","Tson","Tmar","Tmay","Tjun","Tjul","Taug","Tsep",
                   "Paug","J_Spr0","MM_Spr0","DD_Spr0","J_Aut0","T_Aut0")]
#siteo$ssolel <- ifelse(is.na(siteo$J_Spr0),NA,siteo$ssolel)
#siteo$T_Aut0 <- ifelse(is.na(siteo$J_Aut0),NA,siteo$T_Aut0)

#------------------------------------------------------------------------
# Merge laked and climd data.frames
#------------------------------------------------------------------------

lakep <- merge(laked,siteo,by=c("Wby_Lid"))
head(lakep)

#------------------------------------------------------------------
# Compute solar elevation for J_Spr0  (Ang_Spr0)
#----------------------------------------------------------------------------

#lakep$tzone <- NULL
lakep$Ang_Spr0 <- -999.9
for (i in 1:Nsites) 
{
  lat <- lakep$Latitude[i] 
  long <- lakep$Longitude[i]
  year <- year_fix
  month <- lakep$MM_Spr0[i]
  day   <- lakep$DD_Spr0[i]
  # Add GMT offset to local noon
  tzone <- round(long/15.0,digits=0)
  locnoon <- 12 -tzone
  #lakep$tzone[i] <- tzone
  sun.values <- SunPosition (year, month, day, locnoon,0,0,lat,long)
  lakep$Ang_Spr0[i] <- sun.values$elevation
}
head(lakep)

#------------------------------------------------------------------
#  Ice Breakup and Freeze Up Predictions
#------------------------------------------------------------------
lakep$Area_km <- lakep$Area_ha/100;# Area in km2
#lakep$z1 <- as.numeric(lakep$Depth_Mn);# Mean depth in m
#
#lakep$IceBU_part1 <- round(481 + 0.0009417*lakep$Area_km + 0.73145*lakep$Longitude + 0.01477*lakep$Elevation,digits=2)
#lakep$fJ_Spr0 <- round( 0.73048*lakep$J_Spr0,digits=2)
#lakep$fJ_Aut0 <- round( - 0.73048*lakep$J_Aut0,digits=2)
#lakep$fAng_Spr0 <- round( -3.008*lakep$Ang_Spr0,digits=2)
lakep$IceBU <- 481 + 0.73048*lakep$J_Spr0 - 0.73048*lakep$J_Aut0 - 3.008*lakep$Ang_Spr0+ 0.0009417*lakep$Area_km + 0.73145*lakep$Longitude + 0.01477*lakep$Elevation
lakep$IceBU <- round(lakep$IceBU, digits = 1)

lakep$IceFU_part1 <- round(58.0924  + 7.2925*lakep$Depth_Mn^(0.5),digits=2)
#lakep$FU_Aut0 <- round(0.8303*lakep$J_Aut0,digits=2)
#lakep$FUT_Aut0 <- round(0.9435*lakep$T_Aut0,digits=2)
lakep$IceFU <- round(58.0924  + 7.2925*lakep$Depth_Mn^(0.5)+ 0.8303*lakep$J_Aut0 + 0.9435*lakep$T_Aut0,digits=1)

lakep$Icefree <- lakep$IceFU - lakep$IceBU
lakep$Icecover <- 365 - lakep$Icefree
lakep$Icethick <- round(exp(-3.968 + 0.801*log(lakep$Icecover) + 1.103*log(lakep$Latitude) + 0.224/2),digits=2)

#------------------------------------------------------------------
# Checks for odd ice in/out results
#------------------------------------------------------------------

Nsites <- nrow(lakep)
for (i in 1:Nsites) {
  if (is.na(lakep$J_Spr0[i]) || is.na(lakep$J_Aut0[i]) || lakep$J_Aut0[i] <= lakep$J_Spr0[i] || (lakep$J_Aut0[i]-365) >= lakep$J_Spr0[i]) {
    lakep$IceBU[i] <- NA
    lakep$IceFU[i] <- NA
    lakep$Icefree[i] <- 365
    lakep$Icecover[i] <- 0
    lakep$Icethick[i] <- NA
  }
}

#------------------------------------------------------------------
# OUTPUT tmp_ClimMetrics.csv, tmp_IceClimMetrics.csv and 3_Model_inputs.csv)
#------------------------------------------------------------------
write.csv(sited, file = file_out1, row.names = FALSE)
write.csv(siteo, file = file_out2, row.names = FALSE)
write.csv(lakep, file = file_out3, row.names = FALSE)

#------------------------------------------------------------------------------
# STM Parameter Calculations
#------------------------------------------------------------------------------
STM <- lakep

STM$Area_km <- as.numeric(STM$Area_km)
STM$Secchi  <- as.numeric(STM$Secchi)

# Fix for missing DOC if Secchi present
STM$DOC <- ifelse(is.na(STM$DOC),2.0409+18.7373/STM$Secchi-8.2407/STM$Secchi^2-0.3009*log(STM$Depth_Mn),STM$DOC)
#

# Note that log10 is used in some equations and log (natural log) in others 
STM$TX <- round(4.81017 -0.51763*log10(STM$Area_km)^2 +2.4337*log10(STM$Elevation) +0.25207*STM$Tann +0.55343*STM$Tjul +0.14833*STM$Taug, digits = 2)
# STM$TX <- round(4.81017 -0.09763*log(STM$Area_km)^2 +1.0569*log(STM$Elevation) +0.25207*STM$Tann +0.55343*STM$Tjul +0.14833*STM$Taug, digits = 2)
STM$TN <- round(11.9389 +1.07915*log10(STM$Area_km) +2.0226*log10(STM$Depth_Mn) -4.68475*log10(STM$Depth_Max) -0.20951*STM$Secchi +0.09426*STM$Tmar, digits=2)
#STM$TN <- round(11.9389 +0.4687*log(STM$Area_km) + 0.8784*log(STM$Depth_Mn) - 2.0357*log(STM$Depth_Max) -0.20951*STM$Secchi +0.09426*STM$Tmar, digits=2)

for (i in 1:nrow(STM)) {
  if (STM$Area_ha[i] > 800) {
# Old equation had wrong sign and included a delay of 27 days - fixed this  
#    STM$JS[i] <- round(160 + 5.14*STM$Tann[i] +5.74*log(STM$Area_km[i]/STM$Depth_Mn[i]) - 27, digits = 1)     
    STM$JS[i] <- round(160 - 5.14*STM$Tann[i] +5.74*log(STM$Area_km[i]/STM$Depth_Mn[i]), digits = 1) 
  } else {
# old equation had wrong sign and included a delay of 27 days - fixed this   
#    STM$JS[i] <- round((160 + 5.14*STM$Tann[i] +5.74*log(STM$Area_km[i]/STM$Depth_Mn[i]))/2 +(91.24 -5.87*STM$Tmay[i] -3.35*STM$DOC[i] +STM$Area_km[i] +STM$IceBU[i])/2 - 27, digits = 1) 
    STM$JS[i] <- round((160 - 5.14*STM$Tann[i] +5.74*log(STM$Area_km[i]/STM$Depth_Mn[i]))/2 +(91.24 -5.87*STM$Tmay[i] -3.35*STM$DOC[i] +STM$Area_km[i] +STM$IceBU[i])/2, digits = 1) 
  }
}

STM$JM <- round(153.592 -0.93198*STM$Longitude +3.27394*STM$Tmay -4.86477*STM$Tjul +2.83079*STM$Tsep, digits = 1) 
# check log10 in next line
STM$JE <- round(219.445 +21.221*log10(STM$Elevation) +24.5924*log10(STM$Depth_Mn) +2.43965*STM$Secchi +2.28842*STM$Taug -3.97789*STM$Tjja +5.90576*STM$Tson, digits = 1) 
#STM$JE <- round(219.445 + 9.2161*log(STM$Elevation) + 10.6803*log(STM$Depth_Mn) +2.43965*STM$Secchi +2.28842*STM$Taug -3.97789*STM$Tjja +5.90576*STM$Tson, digits = 1) 
# Old ZTH Equation
#STM$ZTH <- round(-16.6934 +0.91478*STM$Latitude +0.25176*STM$Longitude +1.82375*log10(STM$Area_km) +2.90136*log10(STM$Depth_Mn)^2 -1.36735*log10(STM$Depth_Max)^2 +0.32706*STM$Secchi, digits = 2) 
# New ZTH Equation
STM$ZTH <- round(exp(1.68062 +0.22536*log(STM$Area_km) -0.11761*log(STM$Shoreline) +0.04326*STM$Tjja +0.01575*STM$Latitude +0.02193*STM$Secchi -0.01663*STM$JM +0.00005158*STM$JM^2 +0.24518^2/2) , digits = 2)
#
#STM$ZM <- round(23.3876 +0.45338*STM$Latitude +2.29279*log10(STM$Area_km) -7.2525*log10(STM$elev) +6.66815*log10(STM$Depth_Mn) -2.32222*log10(STM$Depth_Max)^2 +0.53844*STM$Secchi -1.1782*STM$Tjul, digits = 2) 
# Use ZJ=0 form of STM model whereby ZM is a constant equal to ZTH
STM$ZM <- STM$ZTH
STM$ZJ <- 0
STM$SP <- round(-15.7148 +0.30155*STM$Latitude +0.13118*STM$Longitude -0.15883*STM$Tjun +0.50025*STM$Tjul +(0.33227/31)*STM$Paug +0.04379*STM$JM, digits = 2)

for (i in 1:nrow(STM)) 
{
  if (is.na(STM$JS[i]) || is.na(STM$JE[i])) 
  {STM$DataGood[i] <- FALSE} 
  else
  {if (STM$JM[i] > STM$JS[i] && STM$JE[i] > STM$JM[i]) 
  {STM$DataGood[i] <- TRUE} 
    else {STM$DataGood[i] <- FALSE}  
  }
  
}

STM <- STM[, c("FMZ","Group","Wby_Lid","Lake_Name","Latitude","Longitude","Elevation","Area_ha","Shoreline","Depth_Max","Depth_Mn","Secchi","DOC","Period","TX","TN","JS","JM","JE","ZTH","ZM","ZJ","SP","IceBU","IceFU","Icefree","DataGood")]
STM$LakeNumber <- NA
STM$Stratified <- NA

#------------------------------------------------------------------
# Calculating LakeNumber for lakes at JM
#------------------------------------------------------------------


for (i in 1:nrow(STM)) 
{
  Lake_ID <- STM$Wby_Lid[i]
  Year <- STM$Year[i]
  Lake_Name <- STM$Lake_Name[i]
  Zmn <- STM$Depth_Mn[i]
  Zmx <- STM$Depth_Max[i]
  A0m2 <- STM$Area_ha[i]*10^4
  TX <- STM$TX[i]
  TN <- STM$TN[i]
  JS <- STM$JS[i]
  JM <- STM$JM[i]
  JE <- STM$JE[i]
  ZM <- STM$ZM[i]
  ZJ <- STM$ZJ[i]
  SP <- STM$SP[i]
  
  jday <- round(JM,0)
  LkNum <- NA
  if (STM$DataGood[i] == TRUE) 
  {
    if (Zmx > ZM) 
    {  
      #--------------------------------------------------------------------------------
      # Generate Bathymetric curve from A0 Zmx Zmn
      #--------------------------------------------------------------------------------
      Zint <- 1;# Depth interval for hypsometric curve is 1 metre
      NZ <- round(Zmx/Zint,digits=0)
      if ( Zmx > (NZ*Zint)) {NZ <- NZ+1}
      bthD <- c(Zint*0:NZ)
      NZL <- NZ+1
      bthD[NZL] <- Zmx
      bthA <- c(0:NZ)
      for (IZ in 1:NZL) {bthA[IZ] <- round(FindArea(A0m2,Zmn,Zmx,bthD[IZ]),digits=1)}
      
      #--------------------------------------------------------------------------------
      # Generate Temperature by depth from STM parameters
      #--------------------------------------------------------------------------------
      TempD <- NULL
      for (k in 1:length(bthD)) {TempD[k] <- stm(TX,TN,JS,JM,JE,ZM,ZJ,SP,jday,bthD[k])}
      
      #--------------------------------------------------------------------------------
      # Calculate stability indicators
      #--------------------------------------------------------------------------------
      ss <- schmidt.stability(TempD,bthD,bthA,bthD, sal=TempD*0)
      #lakecond$stability[j] <- ss
      Zth <- thermo.depth(TempD,bthD,Smin=0.1,seasonal=FALSE)
      LkNum <- -99.0
      #    md <- meta.depths(TempD,bthD,slope=0.1,seasonal=FALSE,unstrat.val=-99.99)
      md <- meta.depths(TempD,bthD,slope=0.1,seasonal=FALSE)
      if (md[1] > 0.0 ) 
      {
        hd <- layer.density (md[2],Zmx,TempD,bthD,bthA,bthD,sal=TempD*0)
        ed <- layer.density (0,md[1],TempD,bthD,bthA,bthD,sal=TempD*0)
        WS <- 15;                # Wind speed km/hr, wind height 10 m std
        wms <- WS*1000/(60*60);  # wind speed m/sec
        us <- uStar(wms,10,ed);  # water friction
        # call to function lake.number (from Library = rLakeAnalyzer)
        LkNum <- lake.number(bthA,bthD,us,ss,md[1],md[2],hd)
      }
    }
    LkNum <- ifelse(is.nan(LkNum),NA,LkNum)
    STM$LakeNumber[i] <- ifelse(is.na(LkNum),NA,round(LkNum,3))
    STM$Stratified[i] <- ifelse(is.na(LkNum),FALSE,ifelse(LkNum < 1,FALSE,TRUE))
    # End if DataGood
  }
  # End of lake loop
}  

STM$Stratified <- ifelse(is.na(STM$LakeNumber),FALSE,ifelse(STM$LakeNumber < 1,FALSE,TRUE))
# Set defaults for processing by ThermoPic program
# Do_Space - calculate Thermal habitat space (volume and area)
# Do_ThermoPic - Create a ThermoPic plot
STM$Do_Space <- ifelse(STM$DataGood,TRUE,FALSE)
STM$Do_ThermoPic <- ifelse(STM$DataGood,TRUE,FALSE)


#------------------------------------------------------------------
# OUTPUT 4_STM_Parameters.csv
#------------------------------------------------------------------
write.csv(STM, file = file_out4, row.names = FALSE)


#------------------------------------------------------------------
# FINIS
#------------------------------------------------------------------

