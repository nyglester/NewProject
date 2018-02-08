# ThermoPic_P2_Report
# Version: V3b (NPL: 2018-10-19)
# ===========================================================================
# Calculates Seasonal Thermal Habitat Space (Volume and Area) in Lakes
# Produces ThermoPics if Do_ThermoPic == TRUE in Input File
# Inputs:  
#    DataIn /0_User_Options.csv
#    DataOut/4_STM_Paramaters.csv
# Ouputs:  
#    DataOut/ThermalSpace4x.csv (where x = Temperature Interval)
#    DataOut/ThermalPics/TPx_LakeNames... (*.JPG)
#
# Default is 4Degree interval, change TP_interval in 0_User_Options for other 1, 2, or 3 degree interval
# ThermoPic plots work well for 3D and 4D - not for 1D and 2D (best to set DO_Thermopics='NO')
# Production of ThermoPic for each lake is controlled by Do_ThermoPic in 4_STM_Parameters.csv
# Adapted from Pgm5v3: ThermoPic_Pgm4Rev_NormLakeThermalCalc&Plot_4D20160906 (Minns)
#----------------------------------------------------------------------------

# Set work directory to the ThermoPic root folder
root<-"C:/R_ThermoPic_v3"
setwd(root)
getwd()

#==============================================================================
# Set file names
# For ThermoPic plots, set location, prefix for name and format
#----------------------------------------------------------------------------
file_in0<- "DataIn/0_User_Options.csv"
file_in1<- "DataOut/4_STM_Parameters.csv"

file_out5<- "DataOut/5_ThermalSpace4D.csv"

TP_folder  <- "DataOut/ThermoPics"
TP_root  <- paste(TP_folder,"/TP4_",sep="")
#TP_format <- "TIFF"
TP_format <-"JPEG"

#============================================================================
# Set switches for testing and other options
#----------------------------------------------------------------------------
#To block plotting of any ThermoPics, set ThermoPic_on = FALSE
ThermoPic_on <- TRUE
#ThermoPic_on <- FALSE

# To make 2D (2 degree interval) report, set Report_2D = TRUE
# Note:  2D ThermoPics only show first 5 Temperature Groups, but Report has all intervals
Report_2D <- FALSE
#Report_2D <- TRUE
if(Report_2D == TRUE)
{
  # 2 Degree intervals from 8 to 28
  file_outsheet4<- "5_ThermalSpace2D"
  TP_root  <- paste(TP_folder,"/TP2_",sep="")
}

# Test code using small number of lakes (If Nlakes_test = 0, do all lakes)
Nlakes_test <- 2

#============================================================================
# Required Libraries
#----------------------------------------------------------------------------
library(rootSolve)
library(XLConnect)
library(plyr)
library(chron)
library(RColorBrewer)

#============================================================================
# Function Sub-Routines
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
# Function: FindDay
# Calculate a julian day for a given temperature based on a known
# day with known temperature and known slope of temperature change
#----------------------------------------------------------------------------
#
FindDay <- function (endtemp,startday,starttemp,tslope)
{
  dtemp <- endtemp-starttemp
  dday <- dtemp/tslope
  endday <- startday+dday
  return(endday)
}

#----------------------------------------------------------------------------
# Function: FindTemp
# Calculate a temperature for a given day based on a known day
# with known temperature and known slope of temperature change
#----------------------------------------------------------------------------
#
FindTemp <- function (endday,startday,starttemp,tslope)
{
  dday <- endday-startday
  dtemp <- dday*tslope
  endtemp <- starttemp+dtemp
  return(endtemp)
}

#----------------------------------------------------------------------------
# Function: stm
# Find the water temperature at a given depth on a given day in a lake
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
# Function: FindVol
# Find the volume between two depths in a lake
#----------------------------------------------------------------------------
#
FindVol <- function(area,zmn,zmx,zstart,zend) {
  zrat <- zmn/zmx
  rli <- (1-zrat)/zrat
  integrand <- function(z) {area*(1-z/zmx)^rli}
  vol <- integrate(integrand,zstart,zend)$val
  return(vol)
}

#----------------------------------------------------------------------------
# Function: FindArea
# Find the area between two depths in a lake
#----------------------------------------------------------------------------
#
FindArea <- function(area,zmn,zmx,ztarget) {
  zrat <- zmn/zmx
  rli <- (1-zrat)/zrat
  zarea <- area*(1-ztarget/zmx)^rli
  return(zarea)
}

#----------------------------------------------------------------------------
# Function: FindDepth
# Find the depth at a particular temperature in a lake
#----------------------------------------------------------------------------
#
FindDepth <- function(zmx,tx,tn,js,jm,je,zm,zj,sp,jd,targettemp) {
  topdepth <- 0
  bottomdepth <- zmx
  middepth <- ifelse( jd > js, zm*(jd-js)/(zj+jd-js),(topdepth+bottomdepth)/2)
  while (bottomdepth - topdepth > 0.005) {
    toptemp <- stm(tx,tn,js,jm,je,zm,zj,sp,jd,topdepth)
    midtemp <- stm(tx,tn,js,jm,je,zm,zj,sp,jd,middepth)
    bottomtemp <- stm(tx,tn,js,jm,je,zm,zj,sp,jd,bottomdepth)
    if (targettemp >= midtemp && targettemp <= toptemp) {
      topdepth <- topdepth
      bottomdepth <- middepth
      middepth <- (topdepth+bottomdepth)/2}
    if (targettemp < midtemp && targettemp >= bottomtemp) {
      topdepth <- middepth
      bottomdepth <- bottomdepth
      middepth <- (topdepth+bottomdepth)/2}
    if (targettemp > toptemp || targettemp < bottomtemp) {
      middepth <- NA
      break}}
  finaldepth <- ifelse (middepth < 0.5,0,middepth)
  return(finaldepth)}

#----------------------------------------------------------------------------
# Function: FindHabitat
# Find the volume and area between two temperatures in a lake
#----------------------------------------------------------------------------
#
FindHabitat <- function(tx,tn,js,jm,je,zm,zj,sp,jd,zmx,area,zmn,lakevol,upperlimit,lowerlimit){
  maxtemp <- stm(tx,tn,js,jm,je,zm,zj,sp,jd,0)
  mintemp <- stm(tx,tn,js,jm,je,zm,zj,sp,jd,zmx)
  if (lowerlimit <= mintemp) {upperdepth <- NA
  } else {
    if (upperlimit >= maxtemp) {upperdepth <- 0
    } else {upperdepth <- FindDepth(zmx,tx,tn,js,jm,je,zm,zj,sp,jd,upperlimit)}}
  if (lowerlimit >= maxtemp) {lowerdepth <- NA
  } else {
    if (lowerlimit <= mintemp) {lowerdepth <- zmx
    } else {lowerdepth <- FindDepth(zmx,tx,tn,js,jm,je,zm,zj,sp,jd,lowerlimit)}}
  HabitatVol <- ifelse (is.na(upperdepth) || is.na(lowerdepth),0,FindVol(area,zmn,zmx,upperdepth,lowerdepth))
  HabitatSpace <- round(HabitatVol/lakevol, digits = 3)
  return(c(HabitatSpace,upperdepth,lowerdepth))}


#----------------------------------------------------------------------------
# Function: FindHabZ
# Find the depth for a given temperature in a lake
#----------------------------------------------------------------------------
#
FindHabZ <- function(tx,tn,js,jm,je,zm,zj,sp,jd,zmx,area,zmn,T1,T2){
  maxtemp <- stm(tx,tn,js,jm,je,zm,zj,sp,jd,0)
  mintemp <- stm(tx,tn,js,jm,je,zm,zj,sp,jd,zmx)
  if ( mintemp < T1 && maxtemp < T1) {Tdepth <- NA}
  if ( mintemp >= T2 && maxtemp >= T2) {Tdepth <- NA}
  if ( mintemp < T1 && maxtemp >= T1) {Tdepth <- FindDepth(zmx,tx,tn,js,jm,je,zm,zj,sp,jd,T1)}
  if ( mintemp < T2 && maxtemp >= T2) {Tdepth <- zmx}
  if ( mintemp < T1 && maxtemp >= T2) {Tdepth <- FindDepth(zmx,tx,tn,js,jm,je,zm,zj,sp,jd,T1)}
  if ( mintemp >= T1 && maxtemp < T2) {Tdepth <- zmx}
  return(Tdepth)
}

#----------------------------------------------------------------------------
# Function: SpStats
# Calculate Annual Summary Habitat Space Statistics
# based on Daily Predictions from the STM model for the period
# Spring isothermal at 4C to Fall/Autumn isothermal at 4C
# Detects Presence of One Continuous  or Two Seasons of Availability
# For Stratified and Unstratified Lakes Respectively
#----------------------------------------------------------------------------
SpStats <-  function(hssub,JM){
  # Seasonal Component
  Ct <- length(hssub$Doy)
  WinV  <-  ifelse(length(hssub$Doy) <= 1,0,round(length(hssub$Doy)/365, digits = 3))
  WinSt <-  hssub$Doy[1]
  WinEn <-  hssub$Doy[Ct]
  hssub$Seas <- c(0:(Ct-1))
  hssub$Seas <- hssub$Seas +WinSt - hssub$Doy
  hssub$Seas <- ifelse(hssub$Seas < 0,2,1)
  ts <- as.data.frame(table(hssub$Seas))
  Seas <- nrow(ts)
  PSpr <- round(ts$Freq[1]/Ct,digits=3)
  SprEn <- NA
  AutSt <- NA
  if(Seas > 1) {
    SprEn <- hssub$Doy[ts$Freq[1]]
    AutSt <- hssub$Doy[ts$Freq[1]+1]
  }
  # Annual Volume Stats
  Vmean <-  round(mean(hssub$V), digits = 2)
  Vsd   <-  round(sd(hssub$V), digits = 2)
  Vmax  <-  round(max(hssub$V),digits = 2)
  Vmin  <-  round(min(hssub$V),digits = 2)
  # Annual Area Stats
  Amean <-  round(mean(hssub$A), digits = 2)
  Asd   <-  round(sd(hssub$A), digits = 2)
  Amax  <-  round(max(hssub$A),digits = 2)
  Amin  <-  round(min(hssub$A),digits = 2)
  # Midsummer Metrics
  if( Seas == 1 && round(JM) > min(hssub$Doy) && round(JM) < max(hssub$Doy)) {
    VJM   <-  round(hssub$V[hssub$Doy == round(JM)],digits=2)
    AJM   <-  round(hssub$A[hssub$Doy == round(JM)],digits=2)
  } else {
    VJM <- NA
    AJM <- NA
  }
  list(WinV,WinSt,WinEn,Vmean,Vsd,Vmax,Vmin,VJM,Amean,Asd,Amax,Amin,AJM,Seas,PSpr,SprEn,AutSt)
}

#----------------------------------------------------------------------------
# End of Function Sub-Routines
#============================================================================


#===========================================================================
# Load Data File with Lake Characteristics, STM parameter estimates, and
# flags for Do_ThermoPic and Stratified (TRUE or FALSE)
#
#---------------------------------------------------------------------------

Options <- read.csv(file_in0, header = TRUE)
head(Options)
spaced <- read.csv(file_in1, header = TRUE)
head(spaced)

#==================================================================
#Process User Options
#------------------------------------------------------------------
Nlakes_test<-as.numeric(as.character(Options$User[5]))


# Switch for turning off all plotting
ThermoPic_on <- ifelse(Options$User[1]=="Yes",TRUE,FALSE)
if (!ThermoPic_on) {spaced$Do_ThermoPic <- FALSE}

TP_interval  <- as.numeric(as.character(Options$User[2]))
TP_text      <- as.character(TP_interval)
TP_interval  <- as.numeric(TP_interval)
TP_interval  <- ifelse(TP_interval>4,4,TP_interval)
TP_interval  <- ifelse(TP_interval<1,4,TP_interval)

# Choose format = TIFF or JPEG
TP_format    <- Options$User[3]

# Change default folder for ThermoPics
TP_folder    <- Options$User[4]

#----------------------------------------------------------------------------
# Set Paramaters and Labels for Space Calculations
#----------------------------------------------------------------------------
#
# Temperature Boundaries (Only Inner Temperature Ranges are Analyzed)
# Default TP_interval is 4 degrees

TmpBnds <- c(0,8,12,16,20,24,28,32,50)

if (TP_interval<4)
{
  if (TP_interval == 1) {TmpBnds <- c(0,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,50)}
  if (TP_interval == 2) {TmpBnds <- c(0,8,10,12,14,16,18,20,22,24,26,28,30,32,50)}
  if (TP_interval == 3) {TmpBnds <- c(0,8,11,15,18,21,24,27,30,33,50)}
  file_outsheet4<- paste("5_ThermalSpace",TP_text,"D",sep="")
  TP_root  <- paste(TP_folder,"/TP",TP_interval,"_",sep="")
}
NTmps <- length(TmpBnds)
NTmpI <- NTmps-1
NTInt <- NTmpI-2
TmpCHARS <- as.character(TmpBnds)
TmpCHARS <- ifelse(nchar(TmpCHARS) < 2,paste("0",TmpCHARS,sep=""),TmpCHARS)

# The generation of these col names could be automated
Ztype <- paste("Z",TmpCHARS,sep="")
Ttype <- rep("Txxxx",NTmpI)
Vtype <- rep("V:xx-xxC",NTmpI)
Atype <- rep("V:xx-xxC",NTmpI)
for (i in 1:NTmpI)
{
  Ttype[i] <- paste("T",TmpCHARS[i],TmpCHARS[i+1],sep="")
  Vtype[i] <- paste("V:",TmpCHARS[i],"-",TmpCHARS[i+1],"C",sep="")
  Atype[i] <- paste("A:",TmpCHARS[i],"-",TmpCHARS[i+1],"C",sep="")
}
Vtype[1] <- paste("V:LT ",TmpCHARS[2],"C",sep="")
Vtype[NTmpI] <- paste("V:GT ",TmpCHARS[NTmpI],"C",sep="")
Atype[1] <- paste("A:LT ",TmpCHARS[2],"C",sep="")
Atype[NTmpI] <- paste("A:GT ",TmpCHARS[NTmpI],"C",sep="")

Speccols <- rev(brewer.pal(n = NTmpI-2, name = 'Spectral'))
Tmpcol <-c("grey75",Speccols,"grey25")

#Tmpcol <-c("grey75","blue2","green2","yellow2","orange2","red2","grey25")
#Tmpcol <-c("grey75","blue2","xx","green2","xx","yellow2","xx","orange2","xx","red2","xx","grey25")

#
#--------------------------------------------------------
# Computing Daily Thermal Habitat (Volume and Area)
# Calculate Annual Summary Statistics
#--------------------------------------------------------
#
ii <- 0;# Internal Counter For Tracking Accumulation of Summary Statistics
HdSpSum.Accum <- NULL  #Initialize summary data collecting dataframe

Nlakes <- length(spaced$Lake_Name)
if (Nlakes_test > 0)
{
  Nlakes <- Nlakes_test
}

for (i in 1:Nlakes)
  # Begin: Lake loop
{
  if (spaced$Do_Space[i] == TRUE )
    # Begin: If Do_Space
  {
    ii <- ii+1
    # Prepare DF for descriptive data
    Lake.info <- data.frame(FMZ=spaced$FMZ[i], Wby_Lid=spaced$Wby_Lid[i],Lake_Name=spaced$Lake_Name[i],Area_ha=spaced$Area_ha[i],Depth_Max=spaced$Depth_Max[i],Depth_Mn=spaced$Depth_Mn[i],Stratified=spaced$Stratified[i],Period=spaced$Period[i],TRange=NA)
    # Set column number of WinSt
    cstart <- 11

    # Extract Lake and STM Characteristics
    Area_ha <- spaced$Area_ha[i]
    Depth_Mn <- spaced$Depth_Mn[i]
    Depth_Max <- spaced$Depth_Max[i]
    TX <- spaced$TX[i]
    TN <- spaced$TN[i]
    JS <- spaced$JS[i]
    JM <- spaced$JM[i]
    JE <- spaced$JE[i]
    ZM <- spaced$ZM[i]
    ZJ <- spaced$ZJ[i]
    SP <- spaced$SP[i]
    # Compute 4C dates
    JS4 <- round(JS - (JM - JS)*(TN - 4)/(TX - TN),0)
    JE4 <- round(JE + (JE-JM)*(TN - 4)/(TX - TN),0)
    # Extract Predicted Ice Dates
    IBU <- spaced$IceBU[i]
    IFU <- spaced$IceFU[i]
    # Prepare variables and DF for computation
    duration <- JE4-JS4+1
    hs <- matrix(data=NA,nrow = duration,ncol = (1+NTmps+2*NTmpI), byrow = FALSE, dimnames = NULL)
    habsp <- as.data.frame(hs)
    colnames(habsp) <- c("Doy",Ztype,Vtype,Atype)
    habsp$Doy <- c(JS4:JE4)

    # Determine total lake volume
    LakeV <- FindVol(Area_ha,Depth_Mn,Depth_Max,0,Depth_Max)

    if (spaced$Stratified[i] == TRUE)
    {
      # Begin: Stratified lake
      for (j in 1:duration)
        # Begin: Day loop (Strat=Yes)
      {
        # Loop to find habitat suitability from JE4 to JS4 for each Doy
        # Begin: Temperature loop
        for (k in 1:NTmpI)
        {
          # Finds depths for temperatures
          DoyTz <- FindHabZ(TX,TN,JS,JM,JE,ZM,ZJ,SP,habsp$Doy[j],Depth_Max,Area_ha,Depth_Mn,TmpBnds[k],TmpBnds[k+1])
          habsp[j,(1+k)] <- DoyTz
        }
        # Compute volumes and areas
        for (l in 1:NTmpI)
        {
          DoyHVol <- NA
          DoyHArea <- NA
          if (is.na(habsp[j,(1+l)])|| habsp[j,(1+l)] <= 0)
          {
            if (is.na(habsp[j,(2+l)]) || habsp[j,(2+l)] <= 0)
            {
              DoyHVol <- NA
              DoyHArea <- NA
            }
            else
            {
              DoyHVol <- FindVol(Area_ha,Depth_Mn,Depth_Max,habsp[j,(l+2)],Depth_Max)
              DoyHArea <- FindArea(Area_ha,Depth_Mn,Depth_Max,habsp[j,(l+2)])
            }
          }
          else
          {
            if (is.na(habsp[j,(2+l)]))
            {
              DoyHVol <- FindVol(Area_ha,Depth_Mn,Depth_Max,0,habsp[j,(l+1)])
              DoyHArea <- Area_ha-FindArea(Area_ha,Depth_Mn,Depth_Max,habsp[j,(l+1)])
            }
            else
            {
              DoyHVol <- FindVol(Area_ha,Depth_Mn,Depth_Max,habsp[j,(l+2)],habsp[j,(l+1)])
              DoyHArea <- FindArea(Area_ha,Depth_Mn,Depth_Max,habsp[j,(l+2)])-FindArea(Area_ha,Depth_Mn,Depth_Max,habsp[j,(l+1)])
            }
          }
          habsp[j,NTmps+l+1] <- round(100.0*DoyHVol/LakeV,digits=2)
          habsp[j,NTmps+NTmpI+l+1] <- round(100.0*DoyHArea/Area_ha,digits=2)
        }
        # End: Temperature loop
      }
      # End: Stratified lake
    }
    else
    {
      # Begin: Non-Stratified lake
      for (j in 1:duration)
      {
        # Begin: Day loop (Non-Stratified)
        # Loop to find habitat suitability from JE4 to JS4 for each Doy
        DJ <- habsp$Doy[j]
        TJ <- ifelse( DJ <= JM, 4+(TX-4)*(DJ-JS4)/(JM-JS4),4+(TX-4)*(JE4-DJ)/(JE4-JM))

        for (l in 1:NTmpI)
        {
          if ( TJ >= TmpBnds[l] && TJ < TmpBnds[l+1])
          {
            # Assign volumes and areas if present
            habsp[j,NTmps+l+1] <- round(100.0,digits=2)
            habsp[j,NTmps+NTmpI+l+1] <- round(100.0,digits=2)
          }
        }
        # End: Non-Stratified lake
      }
      # End: If/Else for Stratified (back 73 lines)
    }

    #----------------------------------------------------------------------------
    # Compute Lake Summary Space Statistics
    #----------------------------------------------------------------------------

    # Temperature range loop begins
    for (k in 2:(NTmps-2))
    {
      # Select Temp Range Results
      Lake.info$TRange <- as.character(Ttype[k])
      hssub <- na.omit(habsp[ ,c(1,(1+NTmps+k),(1+NTmps+NTmpI+k))])
      if (length(hssub$Doy)> 0)
      {
        colnames(hssub) <- c("Doy","V","A")
        # Compute Statistics
        SpSum <- SpStats(hssub,JM)
        SpSum <- as.data.frame(SpSum)
        colnames(SpSum) <- c("WinV","WinSt","WinEn","Vmean","Vsd","Vmax","Vmin","VJM","Amean","Asd","Amax","Amin","AJM","Seas","PSpr","SprEn","AutSt")
      }
      else
      {
        SpSum <- data.frame(WinV=c(0),WinSt=NA,WinEn=NA,Vmean=c(0),Vsd=NA,Vmax=NA,Vmin=NA,VJM=c(0),Amean=c(0),Asd=NA,Amax=NA,Amin=NA,AJM=c(0),Seas=c(0),
                            PSpr=NA,SprEn=NA,AutSt=NA)
      }
      # Add Header
      HdSpSum <- data.frame(Lake.info,SpSum)
      #--------------------------------------------------------
      # Accumulate Results

      if (ii == 1 && k == 2)
      {HdSpSum.Accum <- HdSpSum}
      else {HdSpSum.Accum <- rbind(HdSpSum.Accum,HdSpSum)}
    }

    #----------------------------------------------------------------------------
    # Plotting Occupancy Polygons
    #----------------------------------------------------------------------------
    # Only do ThermoPics for cases with Do_ThermoPic = TRUE

    if (spaced$Do_ThermoPic[i] == TRUE )
      # Begin: if Do_Thermopic
    {
      # Create Label for a lake
      FMZ<- spaced$FMZ[i]
      Wby_Lid <- spaced$Wby_Lid[i]
      Lake_Name <- spaced$Lake_Name[i]
      Period <- spaced$Period[i]
      figure_label <- paste(FMZ," ",Lake_Name," (",Wby_Lid,"): ",Period,sep="")
      #TP_root <- "Data/TP4_"
      # Write TIFF or JPEG File for a Lake (Or plot on screen)
      #--------------------------------------------------------
      if(TP_format == "TIFF")
      {TIFFName <- paste(TP_root,FMZ,"_",Lake_Name,"_",Wby_Lid,"_P",Period,".tiff",sep="")
      #tiff(filename = TIFFName, width = 5, height = 8, units = "in", compression= "lzw",,bg = "white",res=400)
      tiff(filename = TIFFName, width = 5, height = 7, units = "in", pointsize=18, compression= "lzw",bg = "white",res=400)
      }
      if(TP_format == "JPEG")
      {JPEGName <- paste(TP_root,FMZ,"_",Lake_Name,"_",Wby_Lid,"_P",Period,".jpeg",sep="")
      jpeg(filename = JPEGName, width = 400, height = 600, units = "px", pointsize=20, bg="white", res=NA, family="")
      }

      # Plot Instructions Begin
      opar <- par
      par(mfrow=c(1,1), mar=c(2.5,3.0,0.5,0.5),mgp=c(2,0.5,0),cex=0.8)

      for (ifish in 2:(NTmpI-1))
      {
        fhcol <- Tmpcol[ifish]
        #Doy.St <- HdSpSum.Accum[(ii-1)*NTInt+ifish-1,11]
        #Doy.En <- HdSpSum.Accum[(ii-1)*NTInt+ifish-1,12]
        # cstart is 11 to access column with WinSt
        Doy.St <- HdSpSum.Accum[(ii-1)*NTInt+ifish-1,cstart]
        Doy.En <- HdSpSum.Accum[(ii-1)*NTInt+ifish-1,cstart+1]
        if (is.na(Doy.St) == "FALSE" )
        {
          x <- c(Doy.St:Doy.En)
          Doy.min <- min(habsp$Doy)
          x.rows <- c((Doy.St-Doy.min+1):(Doy.En-Doy.min+1))
          y.low <- rep((0+(ifish-2)*100),(Doy.En-Doy.St+1))
          # base is y.low or NA
          y.high <- ifelse(is.na(habsp[x.rows ,NTmps+1+ifish]),y.low,y.low+habsp[x.rows ,NTmps+1+ifish])
          pldt <- na.omit(data.frame(Doy=x,Vlo=y.low,Vhi=y.high))

          if (ifish == 2)
          {
            plot(Vhi~Doy,pldt,type = 'n', ylim = c(-25,650),xlim=c(1,365),las=1,
                 ylab = "Percentage of lake volume", xlab = NA,axes=FALSE)
            xticks = c(0, 32, 60, 91, 121, 152, 182,213,244,274,305,335,366)
            xtlabs = c("  Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan")
            axis(side = 1, at = xticks,labels=xtlabs,cex.axis=0.8)
            yticks = c(0,50,100,150,200,250,300,350,400,450,500,550,600)
            ytlabs = c("  0","50","100","50","100","50","100","50","100","50","100", "50", "100")
            axis(side = 2,at = yticks,labels=ytlabs,cex.axis=0.8,las=1)
          }

          if(HdSpSum.Accum[(ii-1)*NTInt+ifish-1,cstart+12]== 1)
          {
            #One Season
            lines(pldt$Doy,pldt$Vlo, col = fhcol)
            lines(pldt$Doy,pldt$Vhi, col = fhcol)
            polygon(c(pldt$Doy, rev(pldt$Doy)), c(pldt$Vhi, rev(pldt$Vlo)),
                    col = fhcol, border = NA)
          }
          else
          {
            #Two Seasons
            #Spring
            #Doy.SprEn <- HdSpSum.Accum[(ii-1)*NTInt+ifish-1,25]
            Doy.SprEn <- HdSpSum.Accum[(ii-1)*NTInt+ifish-1,cstart+14]
            x.rows <- c(1:(Doy.SprEn-Doy.St+1))
            lines(pldt$Doy[x.rows],pldt$Vlo[x.rows], col = fhcol)
            lines(pldt$Doy[x.rows],pldt$Vhi[x.rows], col = fhcol)
            polygon(c(pldt$Doy[x.rows], rev(pldt$Doy[x.rows])), c(pldt$Vhi[x.rows], rev(pldt$Vlo[x.rows])),
                    col = fhcol, border = NA)
            #Autumn
            Doy.AutSt <- HdSpSum.Accum[(ii-1)*NTInt+ifish-1,cstart+15]
            x.rows <- c((Doy.AutSt-Doy.St+1):(Doy.En-Doy.St+1))
            lines(pldt$Doy[x.rows],pldt$Vlo[x.rows], col = fhcol)
            lines(pldt$Doy[x.rows],pldt$Vhi[x.rows], col = fhcol)
            polygon(c(pldt$Doy[x.rows], rev(pldt$Doy[x.rows])), c(pldt$Vhi[x.rows], rev(pldt$Vlo[x.rows])),
                    col = fhcol, border = NA)
          }
        }
        if (ifish > 1 && ifish < NTmpI)
        {
          #text(50,(ifish-1)*100-7.5,substr(Vtype[ifish],3,8),adj=0,col=fhcol,cex=0.7)
          text(20,(ifish-2)*100+15,substr(Vtype[ifish],3,8),adj=0,col="black",cex=0.7)
        }
        xx <- c(0,365)
        yy <- rep((ifish-2)*100,2)
        lines(xx,yy,col="grey50",lwd=1,lty="dotted")
      }

      lines(xx,c(600,600),col="grey50",lwd=1,lty="dotted")
      lines(c(JM,JM),c(0,600),col="grey50",lwd=1,lty="dotted")
      text(JM,615,label="JM",col="grey50",adj=0.5,cex=0.7)
      lines(c(JS4,JS4),c(0,600),col="darkblue",lwd=2)
      lines(c(JE4,JE4),c(0,600),col="darkblue",lwd=2)
      text(JS4-7,615,label="4C",col="darkblue",adj=0,cex=0.7)
      text(JE4-7,615,label="4C",col="darkblue",adj=0,cex=0.7)
      lines(c(IBU,IBU),c(0,600),col="darkgreen",lwd=2)
      lines(c(IFU,IFU),c(0,600),col="darkgreen",lwd=2)
      text(IBU-7,-15,label="BU",col="darkgreen",adj=0,cex=0.7)
      text(IFU-7,-15,label="FU",col="darkgreen",adj=0,cex=0.7)
      text(186.5,640,label=figure_label,col="black",adj=0.5,cex=0.9)

      box()
      par <- opar
      dev.off()
    }
    # End: If Do_ThermoPic
  }
  # End: If Do_Space
}
# End: Lake loop

#============================================================================
# Prepare for creating CSV outputs
#----------------------------------------------------------------------------
rownames(HdSpSum.Accum) <- c(1:(nrow(HdSpSum.Accum)))

# ADD Annual Yield Index Values as %Year*%Spacemean
HdSpSum.Accum$YVol <- round(HdSpSum.Accum$WinV*HdSpSum.Accum$Vmean,digits=2)
HdSpSum.Accum$YArea <- round(HdSpSum.Accum$WinV*HdSpSum.Accum$Amean,digits=2)

#Rename Thermal Habitat Variables
habnew <- data.frame(PD_year=c(0),TSeasons=c(0),PD_season1=NA,Jstart_Spr=NA,Jend_Spr=NA,Jstart_Aut=NA,Jend_Aut=NA,PV_JM=c(0),PV_mean=c(0),PV_sd=NA,PV_min=NA,PV_max=NA,PV_year=c(0),PA_JM=c(0),PA_mean=c(0),PA_sd=NA,PA_min=NA,PA_max=NA,PA_year=c(0))
habitat <- data.frame(HdSpSum.Accum[,c(1:9)],habnew)
habitat$PD_year <- HdSpSum.Accum$WinV
habitat$TSeasons <- HdSpSum.Accum$Seas
# Rename PD_icefree to PD_season1 (PD_icefree was wrong interpretation)
habitat$PD_season1 <- HdSpSum.Accum$PSpr
habitat$Jstart_Spr <- HdSpSum.Accum$WinSt
habitat$Jend_Spr <- HdSpSum.Accum$SprEn
habitat$Jstart_Aut <- HdSpSum.Accum$AutSt
habitat$Jend_Aut <- HdSpSum.Accum$WinEn
habitat$PV_JM <- HdSpSum.Accum$VJM
habitat$PV_mean <- HdSpSum.Accum$Vmean
habitat$PV_sd <- HdSpSum.Accum$Vsd
habitat$PV_min <- HdSpSum.Accum$Vmin
habitat$PV_max <- HdSpSum.Accum$Vmax
habitat$PV_year<- HdSpSum.Accum$YVol
habitat$PA_JM <- HdSpSum.Accum$AJM
habitat$PA_mean <- HdSpSum.Accum$Amean
habitat$PA_sd <- HdSpSum.Accum$Asd
habitat$PA_min <- HdSpSum.Accum$Amin
habitat$PA_max <- HdSpSum.Accum$Amax
habitat$PA_year<- HdSpSum.Accum$YArea
head(habitat)

# Clean up output
#spaced$DOC <- round(spaced$DOC,digits=1)

write.csv(habitat, file = file_out5, row.names = FALSE)

#----------------------------------------------------------------------------
# FINIS
#----------------------------------------------------------------------------

