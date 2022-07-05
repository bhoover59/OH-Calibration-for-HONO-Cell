###############################################################
### OH Calibration and LOD Calculations                     ###
###############################################################
# Edit extractor constants before running
# Update calibrator characterizations if necessary
# Change directory of user for output
# Txt file must find and replace space with , before running this code
###############################################################

# Clear memory
rm(list = ls()) 

library(readr)
library(dplyr)
library(zoo)
library(data.table)
library(ggplot2)
library(ggpmisc)
library(minpack.lm) # allows nlsLM 
library(grid)
library(gtable)
library(gridExtra)

# Default Extractor Constants, user needs to input if different
P308_start <- 3.05 # mW
P355_start <- 18.0 # W
UVout1_start <- 0.0058 # V
UVout3_start <- 0.0045 # V

# Functions used within Extractor
char_to_time <- function(TimeX){
  formats = c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y %I:%M:%S %p", "%Y-%m-%d %H:%M:%S", "%m/%d/%Y %H:%M", "%m/%d/%Y", "%Y-%m-%d", "%H:%M:%S")
  #Add extra formats to character vector above
  #formats tested in order, so go from most specific to least specific, 
  #as the first one that works will be used
  for (i in 1:length(formats)){
    test = as.POSIXct(TimeX, format = formats[i])
    if (length(na.omit(test)) == length(test)) {
      TimeX = as.POSIXct(TimeX, format = formats[i])
      break
    } else if (i == length(formats)){
      print("Could not find adequate format, char_to_time function needs to be edited")
    }
  }
  return(TimeX)
}

# Input data file
  # Txt file must find and replace space with , before running this code
alldata = read.table(file.choose(), header = TRUE, sep = ",")
# alldata = read.table(file.choose(), header = TRUE, sep = " ")
#alldata = alldata[,!(names(alldata) %in% c("X"))] # removes column X
alldata = alldata[,!(names(alldata) %in% c("TimeX", "dttm"))]
setDT(alldata) # converts list to data table
colnames(alldata)[28] <- "Notes" # Renames column X to Notes 
# Not sure why R sees X instead of Notes from output. File output by LabView shows Notes as column title

# Calibrator Constants @ 10 LPM, needs periodically updated
Cal1_ozone <- 131.49
Cal1_O2_cross_section <- 1.12e-20 
Cal1_ROH <- 0.483 # Wall loss
Cal2_ozone <- 389.23
Cal2_O2_cross_section <- 7.805e-21
Cal2_ROH <- 0.601

# Getting user input for calibrator
calibrator = readline(prompt = "Calibrator (enter 1 or 2): ")
if(calibrator == 1){
  O3 <- Cal1_ozone
  O2_cross_section <- Cal1_O2_cross_section
  WallLoss <- Cal1_ROH
  HONO_PE <- 0.0041
  } else{
  O3 <- Cal2_ozone
  O2_cross_section <- Cal2_O2_cross_section
  WallLoss <- Cal2_ROH
  HONO_PE <- 0.0041
  }

# Getting input from user on water source
water_source = readline(prompt = "Water source, enter 1 for box monitor or 2 for probe: ")
  # If using box monitor (1), RH and Temp must be calculated from box values. 
  # Must remove rows without notes due to incorrect RH during transitions
if(water_source == 1){
  # alldata <- subset(alldata, alldata$Notes != NA) # removes row if there are no notes
  pos_del <- which(alldata$Notes == delete | alldata$Notes == Delete)
  len_pos_del <- length(pos_del) # might need to replace length with nrow
  alldata <- alldata[-(1:len_pos_del), ] # delete rows equal to and before delete or Delete in notes
  alldata <- alldata[!(alldata$Notes == "" | is.na(alldata$Notes)), ] # deletes rows with no notes or NA
  alldata$Temp1 <- (alldata$Notes + 40) * 5 - 40
  alldata$Temp2 <- alldata$Temp1
  alldata$RH1 <- alldata$Notes * 5
  alldata$RH2 <- alldata$RH1
  } # if select 2, then nothing needs to change. Use direct probe output without corrections
exp <- 10**((7.591386 * alldata$Temp1) / (240.7263 + alldata$Temp1)) # separated to make equation easier to read
alldata$Water <- 6.117823 * (alldata$RH1/100) * exp / (1013.25 - (6.116441 * exp * (alldata$RH1/100)))

# Edit classes of variables
if (!("DateTime" %in% names(alldata))){
  alldata$DateTime = paste(alldata$Date, alldata$Time)
}
if ("Time" %in% names(alldata)){
  alldata$Time = char_to_time(alldata$Time)
  if (!class(alldata$Time)[1] == "POSIXct"){
    print("Error in time conversion")
    break
  }
}
if ("Date" %in% names(alldata)){
  alldata$Date = char_to_time(alldata$Date)
  if (!class(alldata$Date)[1] == "POSIXct"){
    print("Error in time conversion")
    break
  }
}

alldata$DateTime = char_to_time(alldata$DateTime)
if (!class(alldata$DateTime)[1] == "POSIXct"){
  print("Error in time conversion")
  break
}

# Separate out notes, as the rest of the data is numeric
# Still need to find a way to incorporate notes back into final data files
notes = alldata$Notes # store notes for later
alldata = alldata[,-c("Notes")] # delete notes. They have been used by now to calculate RH and Temp
# Gets replaced by DateTime somehow

#Change all classes from generic "factor" to numeric
for (i in colnames(alldata)){
  if (class(alldata[[i]])[1]=="factor"){
    alldata[[i]] = as.numeric(as.character(alldata[[i]]))
  } else if (class(alldata[[i]])[1]=="logical"){
    alldata[[i]] = as.numeric(alldata[[i]])
  }
}

# Cycle column indicates one online+offline cycle, could speed extractor if put in program output

if (!"Cycle" %in% colnames(alldata)){
  #Create new column of zeros
  alldata$Cycle = 0
  #creates new column that's the difference between the Online value for each row and the previous
    # mutate creates new variable newcol
    # lag finds previous value
    # %>% translates to "and then" meaning it passes value onto next thing
  alldata = alldata %>% mutate(newcol = Online - lag(Online)) 
  # Finds rows that go from online to offline
    # which returns position of element meeting criteria
  test = which(alldata$newcol == 1)
  # Change data in Cycle column (takes a long time)
  for (i in 2:length(test)){
    alldata$Cycle[test[i-1]:test[i]-1] = i-1
  }
  # Remove newcol used to get cycles
  alldata = alldata[,-c("newcol")]
}

# Add locks (at least 2 before lock indicators and 1 after) to remove all points where laser is moving
# Find rows that currently have lock indicator
test = which(alldata$Lock == 1) # returns positions when lock is on
delete_before <- readline(prompt = "points to delete before lock (at least 2): ")
delete_after <- readline(prompt = "points to delete after lock (at least 1): ")
# Define initial variable
test2 = test[1] # starting with first lock
for (i in 1:length(test)){
  test2 = append(test2, (test[i]-delete_before):(test[i]+delete_after))
} # adds points previous to, equal to, and after lock to test2
test2 = unique(test2) # removes duplicates
alldata$Lock[test2] = 1 # when lock is online, changes to 1
rawdata = alldata[alldata$Lock == 0, -c("Lock")] # data frame without lock column when lock off
rm(test)
rm(test2)

# Laser Scatter (LS) Test
if (sum(rawdata$LS) > 0){
  LS = rawdata[rawdata$LS != 0, ]
  rawdata = rawdata[rawdata$LS == 0, ] # removes points where LS is not 0
} 
# Photo-diode (PD) Test
if (sum(rawdata$PD) > 0){
  PD = rawdata[rawdata$PD != 0, ]
  rawdata = rawdata[rawdata$PD == 0, ] # removes points where PD is not 0
}
rawdata = rawdata[,-c("LS", "PD")] # remove LS and PD columns


# Remove first point of cycle where shutter hasn't opened yet
# Get new column that's the difference between shutter status for a row and the row before it
rawdata = rawdata %>% mutate(newcol = Shutter - lag(Shutter))
# Find first point of opened shutter (start of HONO cycles)
test = which(rawdata$newcol == -1)
# Row before is last point of closed shutter
test = test-1
# Remove points from rawdata and newcol from dataframe
rawdata = rawdata[-test, -c("newcol")]
rm(test)

# Get averages for each cycle for the variables that don't depend on cycling
  # Aggregate returns a dataframe, applies function to dataframe by grouping via list
rawavg = aggregate(rawdata, list(rawdata$Cycle), mean)
rawsd = aggregate(rawdata, list(rawdata$Cycle), sd)
rawsd$Cycle = rawavg$Cycle
rawsd = rawsd[ ,!(names(rawsd) %in% c("Group.1"))]

# Separate and average onlines and offlines
onl = rawdata[rawdata$Online == 1, -c("Online")]
off = rawdata[rawdata$Online == 0, -c("Online")]

# Averages based on Cycle number
onlavg = aggregate(onl, list(onl$Cycle), mean)
offavg = aggregate(off, list(off$Cycle), mean)
# Counts within each cycle
onlct = aggregate(onl[,1], list(onl$Cycle), length)
offct = aggregate(off[,1], list(off$Cycle), length)
# sd within each cycle
onlsd = aggregate(onl, list(onl$Cycle), sd)
offsd = aggregate(off, list(off$Cycle), sd)
# Append counts and signal sd to rest of data
onlavg$count = onlct[[2]]
offavg$count = offct[[2]]
onlavg$sigsd = onlsd$Signal
offavg$sigsd = offsd$Signal

onlavg = onlavg[ , !(names(onlavg) %in% c("Group.1"))]
offavg = offavg[ , !(names(offavg) %in% c("Group.1"))]

# If data set starts on an offline, removes the first offline to match online
# if either on or offline has an extra cycle at the beginning, then the rawdata
# will also have an extra

while(!offavg$Cycle[1] == onlavg$Cycle[1]){
  if(offavg$Cycle[1] < onlavg$Cycle[1]){
    offavg = offavg[-1, ]
    rawavg = rawavg[-1, ]
    rawsd = rawsd[-1, ]
  } else if (offavg$Cycle[1]>onlavg$Cycle[1]){
    onlavg = onlavg[-1, ]
    rawavg = rawavg[-1, ]
    rawsd = rawsd[-1, ]
  }
  print("In while loop, stop code if you see this too many times")
}

if (!(length(offavg$UVout1) == length(onlavg$UVout1))){
  if(length(offavg$UVout1)>length(onlavg$UVout1)){
    notinboth = setdiff(offavg$Cycle, onlavg$Cycle)
    offavg = offavg[!(offavg$Cycle %in% notinboth), ]
    rawavg = rawavg[!(rawavg$Cycle %in% notinboth), ]
    rawsd = rawsd[!(rawsd$Cycle %in% notinboth), ]
  } else {
    notinboth = setdiff(onlavg$Cycle, offavg$Cycle)
    onlavg = onlavg[!(onlavg$Cycle %in% notinboth), ]
    rawavg = rawavg[!(rawavg$Cycle %in% notinboth), ]
    rawsd = rawsd[!(rawsd$Cycle %in% notinboth), ]
  }
}

#Get cycle averages
col_names = colnames(onlavg)
cycles = data.frame(matrix(ncol = length(col_names), nrow = length(onlavg$UVout1)))
names(cycles) = col_names
for (i in col_names){
  if (i == "UVout1"){
    cycles[[i]] = (onlavg[[i]]*onlavg$count+offavg[[i]]*offavg$count)/(onlavg$count+offavg$count)
  } else if(i == "Signal") {
    cycles[[i]] = ((onlavg[[i]]/onlavg$UVout1)-(offavg[[i]]/offavg$UVout1))*((onlavg$count*onlavg$UVout1)+(offavg$count*offavg$UVout1))/(onlavg$count+offavg$count)
  } else if(i == "sigsd"){
    cycles[[i]] = onlavg[[i]]
  } else {
    cycles[[i]] = rawavg[[i]]
  } 
} # Is this net weighted signal? I think so. Will preserve as is and reassign as WeightedNetSignal
cycles$WeightedNetSignal <- cycles$Signal

cycles$bkgsd = offavg$sigsd

# Correct UVout1 for 355 interference
cycles$Shutter[cycles$Shutter == 0 & cycles$UVout3 < 0.01] = 1 # return columns where shutter off and UVout3 less than 0.01
cycles_shutter = cycles[cycles$Shutter == 1, ] # returns rows where shutter is on
cycles_shutter$newcol = 0
cycles_shutter = cycles_shutter %>%
  mutate(newcol = zoo::rollmean(UVout1, k = 10, fill = NA)) #average every 10 points? 10 min resolution?
  # change 10 to variable so it can be 15 or 30 min?
test2 = which(is.na(cycles_shutter$newcol))
for (i in 1:length(test2)){
  if (i < (length(cycles_shutter$newcol)/2)){
    cycles_shutter$newcol[test2[i]]=mean(cycles_shutter$UVout1[1:10])
  } else {
    cycles_shutter$newcol[test2[i]] = mean(cycles_shutter$UVout1[length((cycles_shutter$UVout1)-10):length(cycles_shutter$UVout1)])
  }
}
cycles_shutter = cycles_shutter[,c("Cycle", "newcol")]
setDT(cycles_shutter)
setDT(cycles)
setkey(cycles_shutter, "Cycle")
setkey(cycles, "Cycle")
cycles = cycles_shutter[cycles, roll = TRUE]
cycles$UVout1[cycles$Shutter == 0] = cycles$newcol[cycles$Shutter == 0]
test2 = which(is.na(cycles$newcol))
for (i in 1:length(test2)){
  if (i < length(cycles$newcol)/2){
    cycles$newcol[test2[i]]=mean(cycles$newcol[1:10], na.rm = TRUE)
  } else {
    cycles$newcol[test2[i]] = mean(cycles$newcol[length((cycles$newcol)-10):length(cycles$newcol)], na.rm = TRUE)
  }
}

# Default Dry Sensitivity and Water Dependence
  # need to vary these values to minimize sum of square of differences
dry_sensitivity <- 1.77e-8
water_dependence <- 103.3 

# Calculations for Sensitivities and LODs
if(calibrator == 1){
  cycles$Water <- cycles$Water
  } else{
  cycles$Water <- cycles$Water * 100
}
water_background <- 0 # always 0 in Excel extractor
UVout2_background <- 0 # always 0 in Excel extractor
setkey(cycles, DateTime)
cycles$P308 <- cycles$UVout1 * (P308_start / UVout1_start)
cycles$P355 <- cycles$UVout3 * (P355_start / UVout3_start)
cycles <- cycles[ , -c("newcol")] # not found. When is it created and deleted?
cycles$EstROH <- dry_sensitivity /(1 + (water_dependence * (cycles$Water)))
cycles$PE <- HONO_PE
cycles$NetWater <- cycles$Water - water_background  
cycles$NetUVout2 <- cycles$UVout2 - UVout2_background
cycles$O3 <- cycles$NetUVout2 * O3
cycles$NetOH <- cycles$O3* 1e-9 * cycles$NetWater * 2.4e19 * 7.1e-20 / 2 / 0.2 / O2_cross_section
cycles$ROH <- cycles$WeightedNetSignal / cycles$NetOH / cycles$P308 / WallLoss 
cycles$ROH_scaled <- cycles$ROH * 1e8
cycles$ROH2_scaled <- cycles$ROH * 1e8
cycles$EstROH2 <- dry_sensitivity / (1 + (cycles$Water * water_dependence))
cycles$EstROH_scaled <- cycles$EstROH * 1e8
cycles$EstROH2_scaled <- cycles$EstROH2 * 1e8
cycles$DiffSquare <- (cycles$EstROH_scaled - cycles$ROH_scaled)^2
cycles$DiffSquare2 <- (cycles$EstROH2_scaled - cycles$ROH2_scaled)^2

sumSquare <- sum(cycles$DiffSquare) # not needed
sumSquare2 <- sum(cycles$DiffSquare2) # not needed

# Delete row if ROH below 5e-10 which indicates error
cycles <- subset(cycles, cycles$ROH > 3e-12)

# Linear Fit
fit <- lm(formula = cycles$ROH ~ cycles$NetWater)
cf <- coef(fit)
OH_sensitivity_linear <- cf["(Intercept)"]
slope_linear <- cf["cycles$NetWater"]
summary(fit) # Use to check that the coefficients (cf) are correct
plot(cycles$NetWater, cycles$ROH)
points(cycles$NetWater, cycles$EstROH, col="red")
abline(fit)
title("OH Calibration")
legend("topright", bty="n", legend=paste("R^2 = ", format(summary(fit)$adj.r.squared, digits=4))) # n suppresses box

# Non Linear Fit 
EstROH_Fit <- function(x,a,b){
  a / (1 + b * x) 
} # x is net water, a is dry sensitivity, b is water dependence
y <- cycles$ROH
x <- cycles$NetWater
# CANNOT use nls function. Will get "singular gradient" error
EstROH_nlsLM_fit <- nlsLM(y ~ EstROH_Fit(x,a,b), start = list(a=1e-10, b=100))
Fit_coeff <- coef(EstROH_nlsLM_fit)
dry_sensitivity <- coef(EstROH_nlsLM_fit)[["a"]] 
water_dependence <- coef(EstROH_nlsLM_fit)[["b"]] 
# Need to apply this correction to change EstROH values
cycles$EstROH <- EstROH_Fit(cycles$NetWater, dry_sensitivity, water_dependence)
plot(cycles$NetWater, cycles$ROH)
points(cycles$NetWater, cycles$EstROH, col="red")
fit2 <- lm(formula = cycles$EstROH ~ cycles$NetWater)
cf2 <- coef(fit2)
OH_sensitivity_SV <- cf2["(Intercept)"]
slope_SV <- cf2["cycles$NetWater"]
summary(fit2)

ggplot(cycles, aes(x = NetWater, y = ROH)) +
  geom_smooth(method='lm', formula= y ~ x) +
  geom_point(size=2, shape=16) +
  theme_bw() +
  xlab("Net Water") +
  ylab("ROH") +
  ggtitle("OH Calibration") +
  geom_point(data = cycles, mapping = aes(x = NetWater, y = EstROH), color = "red", shape = 23)


# Background
background <- subset(off, off$RH1 < 10) # reassigning off cycles to background data frame. Only want offline when RH = 0 and laser off. RH = 0 impossible so filter with less than 10
N_background <- readline(prompt = "Time for background (15 min default): ")
N_background <- N_background * 60 # point per second
sd_background <- sd(background$Signal)
avg_background <- mean(background$Signal)

# Calculate OH LOD (molec/cm^3)
P308_avg <- mean(cycles$P308)
RH <- 0.01 # test at relative humidity of 1%
SNR <- 3 # signal to noise ratio (1-3)
ROH_1perc_linear <- slope_linear * RH + OH_sensitivity_linear
ROH_1perc_SV <- dry_sensitivity / (1 + water_dependence * RH)
cycles$P308_avg <- mean(cycles$P308)
cycles$OH_LOD_linear <- (SNR * sd_background) / (sqrt(N_background) * ROH_1perc_linear * P308_avg)
cycles$OH_LOD_SV <- (SNR * sd_background) / (sqrt(N_background) * ROH_1perc_SV * P308_avg)

# Calculate HONO LOD (ppb)
cycles$HONO_LOD_linear <- cycles$OH_LOD_linear/ HONO_PE / 2.4606e10
cycles$HONO_LOD_SV <- cycles$OH_LOD_SV / HONO_PE / 2.4606e10

# For visual ouput
OH_LOD_linear <- (SNR * sd_background) / (sqrt(N_background) * ROH_1perc_linear * P308_avg)
OH_LOD_SV <- (SNR * sd_background) / (sqrt(N_background) * ROH_1perc_SV * P308_avg)
HONO_LOD_linear <- OH_LOD_linear/ HONO_PE / 2.4606e10
HONO_LOD_SV <- OH_LOD_SV / HONO_PE / 2.4606e10

# Output (create table for output at some point)
sprintf("OH Sensitivity = %f", OH_sensitivity_SV)
sprintf("OH LOD = %f", OH_LOD_SV)
sprintf("HONO LOD (ppb) = %f", HONO_LOD_SV)

# Export data
file.directory <- "C:/Users/bodehoov/Desktop/HONO Extracted Files/"
write.csv(cycles, paste(file.directory, "Calibrations/", format(Sys.time(), "%Y%m%d %H%M%S"),
                        "OH Calibration for ", csv.names[1], sep = ""))

