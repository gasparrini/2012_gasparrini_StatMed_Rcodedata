################################################################################
# Updated version of the code for the analysis in:
#
#   "Multivariate meta-analysis for non-linear and other 
#     multi-parameter associations"
#   Gasparrini, Armstrong and Kenward
#   Statistics in Medicine 2012
#   http://www.ag-myresearch.com/statmed2012.html
#
# Update: 26 May 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
################################################################################

################################################################################
# NB: THE EXAMPLE IS DIFFERENT IF COMPARED TO THE PUBLICATION, AS THE ORIGINAL
#     DATA ARE NOT AVAILABLE ANY MORE THROUGH THE PACKAGE NMMAPSlite, NOW
#     ARCHIVED
################################################################################

# LOAD THE PACKAGE
library(dlnm) ; library(splines) ; library(xtable)

# CHECK VERSION OF THE PACKAGE
  if(packageVersion("dlnm")<"2.2.0")
    stop("update dlnm package to version >= 2.2.0")

# LOAD THE DATA
data <- read.csv("regEngWales.csv",row.names=1)

####################################################################
# REGIONS

regions <- as.character(unique(data$regnames))

####################################################################
# LIST OF DATAFRAMES FOR 10 REGIONS

datalist <- lapply(regions, function(region) data[data$regnames==region,])
names(datalist) <- regions

####################################################################
# CITY-LEVEL META-PREDICTORS

lat <- c(54.84815,53.58832,53.72352,52.85539,52.53304,52.03734,51.50583,
  51.24213,51.05361,52.02615)
perclat <- round(quantile(lat,c(1,3)*0.25),1)

####################################################################
# ADDITIONAL INFO

m <- length(datalist)

# MOVING AVERAGE OF TMEAN OVER LAG 0-6
for(i in seq(datalist)) datalist[[i]]$tmean05 <- 
  filter(datalist[[i]]$tmean,rep(1,6)/6,side=1)

# TEMPERATURE RANGES (FOR LAG 0-5)
ranges <- t(sapply(datalist,function(x) range(x$tmean05,na.rm=T)))

# COMPUTE 25TH-75TH PERCENTILES OF META-VARIABLES

# DEFINE THE AVERAGE RANGE, CENTERING POINT, DEGREE AND TYPE OF THE SPLINE
# (THESE PARAMETERS CAN BE CHANGED BY THE USER FOR ADDITIONAL ANALYSES)
cen <- 17
bound <- colMeans(ranges)
degree <- 2
type <- "bs"
df <- 6

# DEFINE THE KNOTS AT TEMPERATURE CORRESPONDING TO AVERAGE PERCENTILES
knotperc <- c(5,35,65,95)
knots <- rowMeans(sapply(datalist,function(x) 
  quantile(x$tmean05,knotperc/100,na.rm=T)))

# SAVE
save.image("data.RData")
	
#
