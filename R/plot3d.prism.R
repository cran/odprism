## M.Maechler:
## requiring & loading fields (and spam!) just for image.plot()
## is definitely an overkill

plot3d.prism <-
function(x, variable, prec=2, absol=TRUE, zmax1=0, zmax2=0) {
  if (!inherits(x, "odprism")) {stop("Only works with \"odprism\" objects") }
  if(prec!=1 & prec!=2) { stop("Please choose either prec=1  or prec=2") }
  if(!(variable %in% c("b0", "bX", "bX2", "VI", "VS", "VR", "C")))
      stop("Please choose variable to be  b0, bX, bX2, VI, VS, VR or C")
  if(variable=="b0") {
    a1<-x$EstInt
    a2<-x$Int
    b1<-x$Q1Int
    b2<-x$Q2Int
    b3<-x$Q3Int
    b4<-x$Q4Int
    d<-x$PowerFixInt
  }
  if(variable=="bX") {
    a1<-x$EstSlope
    a2<-x$Slope
    b1<-x$Slope
    b2<-x$Q2Slope
    b3<-x$Q3Slope
    b4<-x$Q4Slope
    d<-x$PowerFixSlope
  }
   if(variable=="bX2") {
    a1<-x$EstSlope2
    a2<-x$Slope2
    b1<-x$Slope2
    b2<-x$Q2Slope2
    b3<-x$Q3Slope2
    b4<-x$Q4Slope2
    d<-x$PowerFixSlope2
  }
  if(variable=="VI") {
    a1<-x$EstVI
    a2<-x$VI
    b1<-x$Q1VI
    b2<-x$Q2VI
    b3<-x$Q3VI
    b4<-x$Q4VI
    d<-x$PowerRandInt
  }
  if(variable=="VS") {
    a1<-x$EstVS
    a2<-x$VS
    b1<-x$Q1VS
    b2<-x$Q2VS
    b3<-x$Q3VS
    b4<-x$Q4VS
    d<-x$PowerRandSlope
  }
  if(variable=="C" &  length(which(names(x)=="CorIS"))==1) {
    a1<-x$EstCorIS
    a2<-x$CorIS
    b1<-x$Q1CorIS
    b2<-x$Q2CorIS
    b3<-x$Q3CorIS
    b4<-x$Q4CorIS
    d<-x$PowerRandSlope
  }
    if(variable=="C" & length(which(names(x)=="CovIS"))==1) {
    a1<-x$EstCovIS
    a2<-x$CovIS
    b1<-x$Q1CovIS
    b2<-x$Q2CovIS
    b3<-x$Q3CovIS
    b4<-x$Q4CovIS
    d<-x$PowerRandSlope
  }
  if(variable=="VR") {
    a1<-x$EstVR
    a2<-x$VR
    b1<-x$Q1VR
    b2<-x$Q2VR
    b3<-x$Q3VR
    b4<-x$Q4VR
    d<-x$PowerRandInt
  }

  dev.new() # rather than x11()
  if(is.na(d[1]))
      par(mfrow =c(1,2), mar=c(10,5,5,5))
  else par(mfrow=c(1,3), mar=c(10,5,5,5))
  # plot 1 accuracy
  if(absol) a<-abs(a1-a2) else a<-a1-a2
  if(zmax1!=0) { a[a>zmax1]<-zmax1 }
  dim(a)<-c(length(unique(x$R)),length(unique(x$G)))
  a<-t(a)
  revgradient<-palette(rev(c("#FF0000FF","#FF2400FF","#FF4900FF","#FF6D00FF","#FF9200FF","#FFB600FF","#FFDB00FF","#FFFF00FF","#FFFF40FF","#FFFFBFFF","#00FFB2","#0080FFFF")))
  gradient<-palette(c("#FF0000FF","#FF2400FF","#FF4900FF","#FF6D00FF","#FF9200FF","#FFB600FF","#FFDB00FF","#FFFF00FF","#FFFF40FF","#FFFFBFFF","#00FFB2","#0080FFFF"))
  if(x$FLAG[1]) {
      graphics::image(x=sort(unique(x$G)), y=sort(unique(x$R)), z=a,
                      col = gradient, xlab="Population size monitored",
                      ylab="Study period (# sampling occasions)",
                      xlim=c(min(x$G),max(x$G)), ylim=c(min(x$R),max(x$R)),
                      axes=FALSE)
  } else { # ..FLAG is FALSE
      graphics::image(x=sort(unique(x$G)), y=sort(unique(x$R)), z=a,
                      col = gradient, xlab="Number of grouping units",
                      ylab="Replicates per grouping unit",
                      xlim=c(min(x$G),max(x$G)), ylim=c(min(x$R),max(x$R)),
                      axes=FALSE)
  }
  axis(1,las=2, at=sort(unique(x$G)))
  axis(2,las=1, at=sort(unique(x$R)))
  image.plot(a, col = gradient, legend.only=TRUE)
  if(absol==TRUE) {title(main = c("Bias: |Deviation estimates from 'real' value|",variable,"=",a2[1]), font.main = 4)}
  if(absol==FALSE) {title(main = c("Bias: Deviation estimates from 'real' value",variable,a2[1]), font.main = 4)}
  # plot 2 precision
  if(prec==1) b<-b4-b1 else b<-b3-b2
  if(zmax2!=0) { b[b>zmax2]<-zmax2 }
  dim(b)<-c(length(unique(x$R)),length(unique(x$G)))
  b<-t(b)
  if(x$FLAG[1]==TRUE) {   graphics::image(x=sort(unique(x$G)), y=sort(unique(x$R)), z=b, col = gradient, xlab="Population size monitored",ylab="Study period (# sampling occasions)", xlim=c(min(x$G),max(x$G)), ylim=c(min(x$R),max(x$R)), axes=FALSE)}
  if(x$FLAG[1]==FALSE) {   graphics::image(x=sort(unique(x$G)), y=sort(unique(x$R)), z=b, col = gradient, xlab="Number of grouping units",ylab="Replicates per grouping unit", xlim=c(min(x$G),max(x$G)), ylim=c(min(x$R),max(x$R)), axes=FALSE)}
  axis(1,las=2, at=sort(unique(x$G)))
  axis(2,las=1, at=sort(unique(x$R)))
  image.plot(b, col = gradient, legend.only=TRUE)
  if(prec==1) { title(main = c("Imprecision: Quantile 4 - Quantile 1 of", variable), font.main = 4) }
  if(prec==2) { title(main = c("Imprecision: Quantile 3 - Quantile 2 of", variable), font.main = 4) }
  #plot 3 statistical power
  if(is.na(d[1])==FALSE) {
    dim(d)<-c(length(unique(x$R)),length(unique(x$G)))
    d<-t(d)
    if(x$FLAG[1]==TRUE) {  graphics::image(x=sort(unique(x$G)), y=sort(unique(x$R)), z=d, col = revgradient, xlab="Population size monitored",ylab="Study period (# sampling occasions)", xlim=c(min(x$G),max(x$G)), ylim=c(min(x$R),max(x$R)), axes=FALSE)}
    if(x$FLAG[1]==FALSE) {  graphics::image(x=sort(unique(x$G)), y=sort(unique(x$R)), z=d, col = revgradient, xlab="Number of grouping units",ylab="Replicates per grouping unit", xlim=c(min(x$G),max(x$G)), ylim=c(min(x$R),max(x$R)), axes=FALSE)}
    axis(1,las=2, at=sort(unique(x$G)))
    axis(2,las=1, at=sort(unique(x$R)))
    image.plot(d, col = revgradient, legend.only=TRUE)
    title(main = "Power to detect random intercepts", font.main = 4)
  }
}

