################################################################################
# Updated version of the code for the analysis in:
#
#   "Multivariate meta-analysis for non-linear and other 
#     multi-parameter associations"
#   Gasparrini, Armstrong and Kenward
#   Statistics in Medicine 2012
#   http://www.ag-myresearch.com/2012_gasparrini_statmed.html
#
# Update: 05 December 2017
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
# RUN THE FIRST STAGE MODEL
#   COMPUTING TIME IS ~12SEC (IN A 2.66GHz-4GBRAM PC UNDER WINDOWS)
####################################################################

# BUILT OBJECTS WHERE RESULTS WILL BE STORED:
#   ymat IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   Slist IS THE LIST WITH (CO)VARIANCE MATRICES
ymat <- matrix(NA,m,df,dimnames=list(regions,paste("spl",seq(df),sep="")))
Slist <- vector("list",m)
names(Slist) <- regions

####################################################################
# RUN THE FIRST-STAGE ANALYSIS

# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
options(warn=-1)

system.time(
for(i in seq(m)) {
  
  # PRINT ITERATION
  cat(i,"")

  # LOAD
  data <- datalist[[i]]

  # CREATE THE SPLINE
  # NB: KNOTS AND BOUNDARIES FIXED AT SAME VALUES
  btmean05 <- onebasis(data$tmean05,fun=type,degree=degree,knots=knots,
    Bound=bound)

  # RUN THE MODEL
  model <- glm(death ~ btmean05 + dow + ns(time,7*14),family=quasipoisson(),data)
	
  # EXTRACT AND SAVE THE RELATED COEF AND VCOV
  predtmean05 <- crosspred(btmean05,model,cen=cen)
  ymat[i,] <- predtmean05$coef
  Slist[[i]] <- predtmean05$vcov
})

# RESET WARNING
options(warn=0)

#
