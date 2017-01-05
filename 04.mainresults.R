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
# TEMPERATURE RANGES
####################################################################

par(mar=c(5,6,4,1)+0.1)
layout(1)
plot(seq(-10,30,length=m),seq(m),type="n",yaxt="n",
  ylab="",xlab="Temperature range",main="Temperature ranges (C) - Lag 0-5")
axis(2,at=seq(m),labels=regions,las=1,cex.axis=0.7)
arrows(ranges[,1],seq(m),ranges[,2],seq(m),angle=90,length=0.05,code=3)
title(ylab="Cities",mgp=c(5,1,0))
abline(v=knots)
abline(v=bound,lty=2)
mtext("Vertical lines identify internal (continuous) and boundary (dashed) knots",
  cex=0.8)

####################################################################
# GENERATE FIRST PART OF TABLE SIMILAR TO TABLE 2 IN THE MANUSCRIPT
####################################################################

tab2 <- matrix(NA,4,12)
colnames(tab2) <- c("Q","df","p","I-square","AIC","BIC","stat","df","p",
  "stat","df","p")
rownames(tab2) <- c("intercept-only","with latitude","rel intercept-only",
  "rel with latitude")

# FUNCTION TO COMPUTE THE STATISTICS FOR EACH MODEL
# THESE COMPUTATIONS ARE LIKELY TO BE REPLACED BY PROPER anova METHODS IN
#   FUTURE RELEASES OF THE PACKAGE
ftab <- function(m,mref=NULL) {
  # HETEROGENEITY AND IC STATS
  q <- qtest(m)
  het <- c(q$Q[1],q$df[1],q$pvalue[1],(q$Q[1]-q$df[1])/q$Q[1]*100,AIC(m),BIC(m))
  # LR TEST (ONLY FOR META-REGRESSION)
  if(!is.null(mref)) {
    lrstat <- -2*(logLik(mref)-logLik(m))[1]
    df <- attr(logLik(m),"df")-attr(logLik(mref),"df")
    pvalue <- 1-pchisq(lrstat,df)
    lr <- c(lrstat,df,pvalue)
  }
  # WALD TEST (ONLY FOR META-REGRESSION)
  if(!is.null(mref)) {
    coef <- coef(m)[-grep("Int",names(coef(m)))]
    vcov <- vcov(m)[-grep("Int",names(coef(m))),-grep("Int",names(coef(m)))]
    waldstat <- coef%*%solve(vcov)%*%coef
    df <- length(coef)
    pvalue <- 1-pchisq(waldstat,df)
    wald <- c(waldstat,df,pvalue)
  }
  # RESULTS
  if(!is.null(mref)) {
    return(c(het,lr,wald))
  } else return(c(het,rep(NA,6))) 
}

tab2[1,] <- ftab(mv)
tab2[2,] <- ftab(mvlat,mv)

# THE TABLE WILL BE COMPLETED LATER

####################################################################
# POOLED RELATIONSHIP
####################################################################

pdf("pooled.pdf",width=5.6,height=3.8)
# SET par OPTIONS
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(1)

plot(cp,"overall",col=1,lwd=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
  xlab="Temperature (C)")
points(round(knots,1),cp$allRRfit[as.character(round(knots,1))],pch=19,cex=0.6)

dev.off()

####################################################################
# BLUP
####################################################################

pdf("bluptot.pdf",width=11.2,height=3.8)
# SET par OPTIONS AND MULTIPANEL
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(matrix(1:2,1,2))

# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
options(warn=-1)

# PLOT OF AVERAGE AND CITY-SPECIFIC ESTIMATES
plot(cp,"overall",type="n",ci="n",ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
  xlim=range(ranges),xlab="Temperature (C)")
for(i in seq(m)) {
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  lines(crosspred(btmean,coef=ymat[i,],vcov=Slist[[i]],model.link="log",
    from=ranges[i,1],to=ranges[i,2],cen=cen),col=grey(0.8),lty=5)
}
lines(cp,"overall",col=1,lwd=2)
abline(h=1)
mtext("Study-specific",cex=1.3)

# PLOT OF AVERAGE AND BLUP ESTIMATES
blup <- blup(mv,vcov=TRUE)
plot(cp,"overall",type="n",ci="n",ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
  xlim=range(ranges),xlab="Temperature (C)")
for(i in seq(m)) {
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  lines(crosspred(btmean,coef=blup[[i]]$blup,vcov=blup[[i]]$vcov,
    model.link="log",from=ranges[i,1],to=ranges[i,2],cen=cen),col=grey(0.8),
    lty=5)
}
lines(cp,"overall",col=1,lwd=2)
abline(h=1)
mtext("BLUP",cex=1.3)

# RESET WARNING
options(warn=0)

dev.off()


pdf("blupregion.pdf",width=11.2,height=3.8)
# SET par OPTIONS AND MULTIPANEL
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(matrix(1:2,1,2))

plot(cp,"overall",col=1,lwd=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),ci="n",
  xlab="Temperature (C)")
lines(crosspred(btmean,coef=ymat["London",],vcov=Slist[["London"]],
  model.link="log",by=0.1,cen=cen),lty=4)
lines(crosspred(btmean,coef=blup[["London"]]$blup,vcov=blup[["London"]]$vcov,
  model.link="log",by=0.1,cen=cen),lty=5)
mtext("London",cex=1.3)
legend("top",c("Population-average","First-stage","BLUP"),lty=c(1,4,5),
  lwd=c(2,1,1),cex=0.8,inset=0.1,bty="n")

plot(cp,"overall",col=1,lwd=2,ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),ci="n",
  xlab="Temperature (C)")
lines(crosspred(btmean,coef=ymat["N-East",],vcov=Slist[["N-East"]],
  model.link="log",by=0.1,cen=cen),lty=4)
lines(crosspred(btmean,coef=blup[["N-East"]]$blup,vcov=blup[["N-East"]]$vcov,
  model.link="log",by=0.1,cen=cen),lty=5)
mtext("N-East",cex=1.3)
legend("top",c("Population-average","First-stage","BLUP"),lty=c(1,4,5),
  lwd=c(2,1,1),cex=0.8,inset=0.1,bty="n")

dev.off()

####################################################################
# META-REGRESSION
####################################################################

pdf("metareg.pdf",width=6,height=4)
# SET par OPTIONS AND MULTIPANEL
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9)
layout(1)

plot(cplat25,"overall",type="n",ylab="RR",ylim=c(.9,1.4),xlim=c(-5,25),
  xlab="Temperature (C)",ci.arg=list(density=20,col=grey(0.7)))
lines(cplat25,"overall",col=1,lty=4,lwd=2)
lines(cplat75,"overall",col=1,lty=5,lwd=2,ci="area",
  ci.arg=list(density=20,angle=-45,col=grey(0.7)))
abline(h=1)
legend("top",paste(perclat),lty=c(4,5),lwd=2,cex=1,inset=0.1,bty="n",
  title="Latitude (degree North)")

dev.off()

#
