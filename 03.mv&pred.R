################################################################################
# Updated version of the code for the analysis in:
#
#   "Multivariate meta-analysis for non-linear and other 
#     multi-parameter associations"
#   Gasparrini, Armstrong and Kenward
#   Statistics in Medicine 2012
#   http://www.ag-myresearch.com/2012_gasparrini_statmed.html
#
# Update: 15 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2012_gasparrini_StatMed_Rcodedata
################################################################################

################################################################################
# NB: THE EXAMPLE IS DIFFERENT IF COMPARED TO THE PUBLICATION, AS THE ORIGINAL
#     DATA ARE NOT AVAILABLE ANY MORE THROUGH THE PACKAGE NMMAPSlite, NOW
#     ARCHIVED
################################################################################

####################################################################
# 1) RUN THE MODELS WITH mvmeta
#
# 2) CREATE BASIS VARIABLES USING onebasis, USED FOR PREDICTION
#
# 3) PREDICT THE OUTCOME PARAMETERS OVER SPECIFIC VALUES OF STUDY-LEVEL
#   COVARIATES USING predict (mvmeta),THEN RE-BUILD THE PREDICTED CURVE AT
#   THOSE VALUES USING crosspred AGAIN
#
# NOTE: THE USE OF dlnm FUNCTIONS FACILITATES PREDICTION AND PLOTTING
#
####################################################################

####################################################################
# PERFORM MULTIVARIATE META-ANALYSIS
####################################################################

# LOAD THE PACKAGES (mvmeta PACKAGE IS ASSUMED TO BE INSTALLED)
library(mvmeta)

# MULTIVARIATE META-ANALYSIS
mv <- mvmeta(ymat,Slist,method="ml")
summary(mv)

# MULTIVARIATE META-REGRESSION
(mvlat <- mvmeta(ymat~lat,Slist,method="ml"))

# NB: IN VERSION 0.4.1, CONVERGENCE MAY BE INSPECTED USING THE ARGUMENT:
#   control=list(showiter=T)
# NB: LESS STRICT CONVERGENCE CRITERIA, USEFUL FOR HIGH DIMENSIONAL
#   MODELS, MAY BE SELECTED BY ADDING A reltol ARGUMENT, FOR EXAMPLE:
#   control=list(showiter=T,reltol=10^-3)

####################################################################
# CREATE BASIS FOR PREDICTION
####################################################################

# BASIS USED TO PREDICT TEMPERATURE, EQUAL TO THAT USED FOR ESTIMATION
#   NOTE: INTERNAL AND BOUNDARY KNOTS PLACED AT SAME VALUES AS IN ESTIMATION
tmean <- seq(bound[1],bound[2],length=30)
btmean <- onebasis(tmean,fun=type,degree=degree,knots=knots,Bound=bound)
 
####################################################################
# PREDICTION FROM MODELS
####################################################################

# USE OF crosspred TO PREDICT THE EFFECTS FOR THE CHOSEN VALUES

# PREDICTION FROM SIMPLE META-ANALYSES WITH NO PREDICTORS
# CENTERED TO SPECIFIC VALUE
cp <- crosspred(btmean,coef=coef(mv),vcov=vcov(mv),model.link="log",
  by=0.1,cen=cen)

# COMPUTE PREDICTION FOR MULTIVARIATE META-REGRESSION MODELS
#   1ST STEP: PREDICT THE OUTCOME PARAMETERS FOR SPECIFIC VALUES OF META-PREDICTOR
#   2ND STEP: PREDICT THE RELATIONSHIP AT CHOSEN VALUES GIVEN THE PARAMETERS

predlat <- predict(mvlat,data.frame(lat=perclat),vcov=T)
cplat25 <- crosspred(btmean,coef=predlat[[1]]$fit,vcov=predlat[[1]]$vcov,
  model.link="log",by=0.1,cen=cen)
cplat75 <- crosspred(btmean,coef=predlat[[2]]$fit,vcov=predlat[[2]]$vcov,
  model.link="log",by=0.1,cen=cen)

#
