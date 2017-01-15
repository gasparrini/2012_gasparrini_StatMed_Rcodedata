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
# RUN THE TWO-STAGE ANALYSIS FOR RELATIVE SCALE
#   TOTAL COMPUTING TIME IS ~15SEC (IN A 2.66GHz-4GBRAM PC UNDER WINDOWS)
####################################################################

# BUILT OBJECTS WHERE RESULTS WILL BE STORED:
#   ymat IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   Slist IS THE LIST WITH (CO)VARIANCE MATRICES
ymatrel <- matrix(NA,m,df,dimnames=list(regions,paste("spl",seq(df),sep="")))
Slistrel <- vector("list",m)
names(Slistrel) <- regions

####################################################################
# RUN THE FIRST-STAGE ANALYSIS 
####################################################################

system.time(
for(i in seq(m)) {
  
  # PRINT ITERATION
  cat(i,"")

  # LOAD
  data <- datalist[[i]]
  
  # DEFINE THE KNOTS AT STUDY-SPECIFIC PERCENTILES
  tempknots <- quantile(data$tmean05,knotperc/100,na.rm=TRUE)
  cenrel <- quantile(data$tmean05,0.75,na.rm=TRUE)

  # CREATE THE QUADRATIC SPLINE: STUDY-SPECIFIC KNOTS
  bperc <- onebasis(data$tmean05,fun=type,degree=degree,knots=tempknots)

  # RUN THE MODEL
  model <- glm(death ~ bperc + dow + ns(time,7*14),family=quasipoisson(),data)
	
  # EXTRACT AND SAVE THE RELATED COEF AND VCOV
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  predperc <- crosspred(bperc,model)
  ymatrel[i,] <- predperc$coef
  Slistrel[[i]] <- predperc$vcov
})

####################################################################
# MULTIVARIATE META-ANALYSIS AND PREDICTION
####################################################################

# MULTIVARIATE META-ANALYSIS
(mvrel <- mvmeta(ymatrel,Slistrel,method="ml"))
# MULTIVARIATE META-REGRESSION MODELS
(mvrellat <- mvmeta(ymatrel~lat,Slistrel,method="ml"))

# BASIS USED TO PREDICT TEMPERATURE, EQUAL TO THAT USED FOR ESTIMATION
#   NOTE: INTERNAL AND BOUNDARY KNOTS PLACED AT SAME VALUES AS IN ESTIMATION
tperc <- seq(bound[1],bound[2],length=30)
perctmean <- rowMeans(sapply(datalist,function(x) {
  quantile(x$tmean05,c(0:10/10,1:99,990:1000/10)/100,na.rm=T)
}))
bperc <- onebasis(tperc,fun=type,degree=degree,knots=knots,Bound=bound)

# PREDICTION FROM SIMPLE META-ANALYSES WITH NO PREDICTORS
cprel <- crosspred(bperc,coef=coef(mvrel),vcov=vcov(mvrel),model.link="log",
  at=c(perctmean,knots),cen=perctmean["75.0%"])

# COMPUTE PREDICTION FOR MULTIVARIATE META-REGRESSION MODELS
#   1ST STEP: PREDICT THE OUTCOME PARAMETERS FOR SPECIFIC VALUES OF PREDICTOR
#   2ND STEP: PREDICT THE RELATIONSHIP AT CHOSEN VALUES GIVEN THE PARAMETERS
predrellat <- predict(mvrellat,data.frame(lat=perclat),vcov=T)
cprellat25 <- crosspred(bperc,coef=predrellat[[1]]$fit,
  vcov=predrellat[[1]]$vcov,model.link="log",by=0.5,cen=perctmean["75.0%"])
cprellat75 <- crosspred(bperc,coef=predrellat[[2]]$fit,
  vcov=predrellat[[2]]$vcov,model.link="log",by=0.5,cen=perctmean["75.0%"])

####################################################################
# GENERATE SECOND PART OF TABLE 2 IN THE MANUSCRIPT
####################################################################

tab2[3,] <- ftab(mvrel)
tab2[4,] <- ftab(mvrellat,mvrel)

tab2
xtable(tab2,digits=c(0,1,0,3,1,1,1,1,0,3,1,0,3),align="lccccccrccrcc")

####################################################################
# PLOT
####################################################################

# VALUES LABEL FOR TEMPERATURE
labperc <- paste(c(0,1,5,25,50,75,95,100),sep="")

pdf("relscale.pdf",width=11.2,height=3.8)
# SET par OPTIONS AND MULTIPANEL
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(matrix(1:2,1,2))

plot(cprel,"overall",type="n",ylab="RR",ylim=c(.9,1.3),
  xlab="Temperature percentiles",xaxt="n")
axis(1,at=perctmean[paste(labperc,".0%",sep="")],
  labels=paste(labperc,"%",sep=""))

for(i in seq(m)) {
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  lines(crosspred(bperc,coef=ymatrel[i,],vcov=Slistrel[[i]],model.link="log",
    cen=perctmean["75.0%"]),col=grey(0.8),lty=5)
}
lines(cprel,"overall",col=1,lwd=2)
points(round(knots,1),cprel$allRRfit[as.character(knots)],pch=19,cex=0.6)
abline(h=1)

plot(cprellat25,"overall",type="n",ylab="RR",ylim=c(.9,1.3),
  xlab="Temperature percentiles",ci.arg=list(density=20,col=grey(0.7)),xaxt="n")
lines(cprellat25,"overall",col=1,lty=4,lwd=2)
lines(cprellat75,"overall",col=1,lty=5,lwd=2,ci="area",
  ci.arg=list(density=20,angle=-45,col=grey(0.7)))
abline(h=1)
axis(1,at=perctmean[paste(labperc,".0%",sep="")],
  labels=paste(labperc,"%",sep=""))
legend("top",paste(perclat),lty=c(4,5),lwd=2,cex=1,inset=0.1,bty="n",
  title="Latitude (degree North)")

dev.off()

#
