####################################################################
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
####################################################################

################################################################################
# NB: THE EXAMPLE IS DIFFERENT IF COMPARED TO THE PUBLICATION, AS THE ORIGINAL
#     DATA ARE NOT AVAILABLE ANY MORE THROUGH THE PACKAGE NMMAPSlite, NOW
#     ARCHIVED
################################################################################


####################################################################
# KNOTS POSITION
####################################################################
knots
bound

####################################################################
# MINIMUM MORTALITY POINT
####################################################################

# ON ABSOLUTE SCALE
cp$predvar[which.min(cp$allRRfit)]

# ON RELATIVE SCALE
names(perctmean)[which.min(cprel$allRRfit)]

####################################################################
# PREDICTED AVERAGE RR FOR 0 AND 20C
####################################################################

with(cp,cbind(allRRfit,allRRlow,allRRhigh)["0",])
with(cp,cbind(allRRfit,allRRlow,allRRhigh)["20",])


####################################################################
# PREDICTED AVERAGE RR FOR 29C, AT 25TH-75TH PERCENTILES OF LATITUDE
####################################################################

with(cplat25,cbind(allRRfit,allRRlow,allRRhigh)["20",])
with(cplat75,cbind(allRRfit,allRRlow,allRRhigh)["20",])


####################################################################
# COMPARISON ML-REML
####################################################################

# RUN THE MODEL WITH REML
(mvreml <- mvmeta(ymat,Slist,method="reml"))
cpreml <- crosspred(btmean,coef=coef(mvreml),vcov=vcov(mvreml),
  model.link="log",by=0.5,cen=cen)

# COMPARISON OF STANDARD ERRORS
(sqrt(diag(vcov(mvreml))) - sqrt(diag(vcov(mv)))) / sqrt(diag(vcov(mv))) * 100

# COMPARISON OF CURVES AND CI
par(mar=c(5,4,2,1)+0.1,cex.axis=0.7)
layout(1)
plot(cp,"overall",ci="lines",col=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
  xlab="Temperature (C)")
lines(cpreml,"overall",ci="lines",col=4)
legend("top",c("ML","REML"),col=c(2,4),lty=1)

####################################################################
# CORRELATIONS
####################################################################

# IN FIRST-STAGE ESTIMATES
cor(ymat)

# BETWEEN-STUDY CORRELATIONS FROM MODEL
cov2cor(mv$Psi)

#
