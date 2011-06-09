plot2d.prism <-
function(x, variable, xvar="G", cons, ylim1=c(0,0)) {
  if (!inherits(x, "odprism")) {stop("Only works with \"odprism\" objects") }
  if(variable!="b0" & variable!="bX" & variable!="bX2" & variable!="VI" & variable!="VS" & variable!="VR" & variable!="C" )   { stop("Please choose varaible to be b0, bX, VI, VS, VR or C") }
  if(xvar!="G" & xvar!="R")   { stop("Please choose xvar to be either G or R") }
  x11()
 ifelse(is.na(x$PowerRandInt[1])==TRUE, par(mfrow=c(1,1), mar=c(10,5,5,5)), par(mfrow=c(1,2), mar=c(10,5,5,5)))
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
    b1<-x$Q1Slope
    b2<-x$Q2Slope
    b3<-x$Q3Slope
    b4<-x$Q4Slope
    d<-x$PowerFixInt
  }   
  if(variable=="bX2") {
    a1<-x$EstSlope2
    a2<-x$Slope2
    b1<-x$Q1Slope2
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
  if(variable=="C" & length(which(names(x)=="CorIS"))==1) {
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
  if(xvar=="G") {
    # put G on x-axis
    selection<-which(x$R==cons)
    if(length(selection)==0) {
      cons<-max(x$R)
      selection<-which(x$R==cons)
      if(x$FLAG[1]==TRUE) { print(c("You have not selected the study period that needs to be kept constant, therefore the maximum is chosen",cons)) }
      if(x$FLAG[1]==FALSE) { print(c("You have not selected the number of replicates that needs to be kept constant, therefore the maximum is chosen",cons)) }
    }
    if(ylim1[1]==0 & ylim1[2]==0) { ylim1=c(min(b1[selection]),max(b4[selection])) }
    if(x$FLAG[1]==TRUE) {  plot(x$G[selection],b4[selection], type="l", col="red", ylim=ylim1, xlab="Population size monitored", ylab="parameter value") }
    if(x$FLAG[1]==FALSE) {  plot(x$G[selection],b4[selection], type="l", col="red", ylim=ylim1, xlab="Number of grouping units", ylab="parameter value") }
    points(x$G[selection],a1[selection], col="blue")
    lines(x$G[selection],b1[selection], col="red")
    lines(x$G[selection],b2[selection], col="green")
    lines(x$G[selection],b3[selection], col="green")
    lines(x$G[selection],a1[selection], col="blue")
    lines(x$G[selection],a2[selection])
    if(x$FLAG[1]==TRUE) { title(main = c("Distribution of estimated values from simulations for",variable, "and study period kept constant at",cons), font.main = 4)}
    if(x$FLAG[1]==FALSE) { title(main = c("Distribution of estimated values from simulations for",variable, "and replicates kept constant at",cons), font.main = 4)}
    if(is.na(d[1])==FALSE) {    
      if(x$FLAG[1]==TRUE) {plot(x$G[selection],d[selection], type="l", ylim=c(0,1), xlab="Population size monitored", ylab="Statistical power") }
      if(x$FLAG[1]==FALSE) {plot(x$G[selection],d[selection], type="l", ylim=c(0,1), xlab="Number of grouping units", ylab="Statistical power") }
      points(x$G[selection],d[selection])
      if(variable=="b0" | variable=="bX" | variable=="bX2") { title(main = c("test for fixed parameter", variable, "included in model"), font.main = 4) }
      if(variable=="VI" | variable=="VR" ) { title(main = c("test for random intercept included in model"), font.main = 4) }
      if(variable=="VS"| variable=="C") { title(main = c("test for random slopes included in model"), font.main = 4) }
    }
  }
  if(xvar=="R") {
    # put R on x-axis
    selection<-which(x$G==cons)
    if(length(selection)==0) {
      cons<-max(x$G)
      selection<-which(x$G==cons)
      if(x$FLAG[1]==TRUE) { print(c("You have not selected the monitored population size that needs to be kept constant, therefore the maximum is chosen",cons)) }
      if(x$FLAG[1]==FALSE) { print(c("You have not selected the number of individuals sampled that needs to be kept constant, therefore the maximum is chosen",cons)) }
    }
    if(ylim1[1]==0 & ylim1[2]==0) { ylim1=c(min(b1[selection]),max(b4[selection])) }
    if(x$FLAG[1]==TRUE) {  plot(x$R[selection],b4[selection], type="l", col="red", ylim=ylim1, xlab="Study period (# sampling occasions)", ylab="parameter value") }
    if(x$FLAG[1]==FALSE) {  plot(x$R[selection],b4[selection], type="l", col="red", ylim=ylim1, xlab="Number of replicates", ylab="parameter value") }
    points(x$R[selection],a1[selection], col="blue")
    lines(x$R[selection],b1[selection], col="red")
    lines(x$R[selection],b2[selection], col="green")
    lines(x$R[selection],b3[selection], col="green")
    lines(x$R[selection],a1[selection], col="blue")
    lines(x$R[selection],a2[selection])
    if(x$FLAG[1]==TRUE) { title(main = c("Distribution of estimated values from simulations for",variable, "and population size kept constant at",cons), font.main = 4)}
    if(x$FLAG[1]==FALSE) { title(main = c("Distribution of estimated values from simulations for",variable, "and number of individuals kept constant at",cons), font.main = 4)}
    if(is.na(d[1])==FALSE) {    
      if(x$FLAG[1]==TRUE) { plot(x$R[selection],d[selection], type="l", ylim=c(0,1), xlab="Study period (# sampling occasions)", ylab="Statistical power")}
      if(x$FLAG[1]==FALSE) { plot(x$R[selection],d[selection], type="l", ylim=c(0,1), xlab="Number of replicates", ylab="Statistical power")} 
      points(x$R[selection],d[selection])
      if(variable=="b0" | variable=="bX" | variable=="bX2") { title(main = c("test for fixed parameter", variable, "included in model"), font.main = 4) }
      if(variable=="VI" | variable=="VR" ) { title(main = c("test for random intercept included in model"), font.main = 4) }
      if(variable=="VS"| variable=="C") { title(main = c("test for random slopes included in model"), font.main = 4) }
    }
  }
}

