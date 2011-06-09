odprismB <-
function(indiv, repl, sims=1000, fixed=c(0,0), correl=TRUE, random, Xvar=c(1,0,FALSE), alpha=0, Quant=c(0.025,0.25,0.75,0.975), Verbal=TRUE) {
  VI<-random[1]
  VS<-random[2]
  if(correl==TRUE) {
    CorIS<-random[3]
    CovIS<-CorIS*sqrt(VI)*sqrt(VS)  
  }
  if(correl==FALSE) {
    CovIS<-random[3]
    CorIS<-CovIS/(sqrt(VI)*sqrt(VS))
    if(abs(CorIS)>1) { stop("Please select covariance value that is realistic i.e. abs(CovIS/(sqrt(VI)*sqrt(VS)))<=1") }
  }
  ifelse((length(random)<4), VR<-1-VI, VR<-random[4]) 
  beta0<-fixed[1]
  betaX<-fixed[2]
  VX<-Xvar[1]
  Ar1X<-Xvar[2]
  SameEnv<-Xvar[3]
  if(SameEnv==FALSE && Ar1X>0) { stop("Sorry, autocorrelation not possible if SameEnv==FALSE") }
  Nrepl<-rep(repl,length(indiv))  
  Nindiv<-sort(rep(indiv,length(repl))) 
  modeloutput<-matrix(data = NA, nrow = sims, ncol = 10)    
  g<-numeric(length(Nindiv))
  r<-numeric(length(Nindiv))
  estInt<-numeric(length(Nindiv))
  estSlope<-numeric(length(Nindiv))
  estVR<-numeric(length(Nindiv))
  estVI<-numeric(length(Nindiv))
  estVS<-numeric(length(Nindiv))
  estCorIS<-numeric(length(Nindiv))
  estCovIS<-numeric(length(Nindiv))
  q1Int<-numeric(length(Nindiv))
  q2Int<-numeric(length(Nindiv))
  q3Int<-numeric(length(Nindiv))
  q4Int<-numeric(length(Nindiv))
  q1Slope<-numeric(length(Nindiv))
  q2Slope<-numeric(length(Nindiv))
  q3Slope<-numeric(length(Nindiv))
  q4Slope<-numeric(length(Nindiv))
  q1VR<-numeric(length(Nindiv))
  q2VR<-numeric(length(Nindiv))
  q3VR<-numeric(length(Nindiv))
  q4VR<-numeric(length(Nindiv))
  q1VI<-numeric(length(Nindiv))
  q2VI<-numeric(length(Nindiv))
  q3VI<-numeric(length(Nindiv))
  q4VI<-numeric(length(Nindiv))
  q1VS<-numeric(length(Nindiv))
  q2VS<-numeric(length(Nindiv))
  q3VS<-numeric(length(Nindiv))
  q4VS<-numeric(length(Nindiv))
  if(correl==TRUE) {
    q1CorIS<-numeric(length(Nindiv))
    q2CorIS<-numeric(length(Nindiv))
    q3CorIS<-numeric(length(Nindiv))
    q4CorIS<-numeric(length(Nindiv))
  }
  if(correl==FALSE) {
    q1CovIS<-numeric(length(Nindiv))
    q2CovIS<-numeric(length(Nindiv))
    q3CovIS<-numeric(length(Nindiv))
    q4CovIS<-numeric(length(Nindiv))
  }  
  powerFixInt<-numeric(length(Nindiv))
  powerFixSlope<-numeric(length(Nindiv))
  powerRandInt<-numeric(length(Nindiv))
  powerRandSlope<-numeric(length(Nindiv))   
  Scounter<-0
  for (m in 1:(length(Nindiv))) {
    Scounter<-Scounter+1
    Nindi<-Nindiv[m]
    Nrep<-Nrepl[m]
    if (Verbal==TRUE) { print(c("Now running individuals",Nindi,"and replicates",Nrep))}
    for (i in 1:sims) {
      dataset<-matrix(data = NA, nrow =Nindi*Nrep, ncol = 4)
      alive<-matrix(data = NA, nrow = Nindi, ncol = 3) 
      alive[,1]<-seq(1,Nindi)
      IndHet <- rmvnorm(Nindi, c(0, 0), sigma=matrix(c(VI,CovIS, CovIS, VS), ncol=2), method = "svd")  
      alive[,2]<-IndHet[,1]   # elevation of individuals
      alive[,3]<-IndHet[,2]   # slope of individuals
      IDcounter<-Nindi
      for (j in 1:Nrep) {
        # select values for X in each sampling occasion      
        if (SameEnv==TRUE) { ifelse (j==1, sX<-rnorm(1,0,sqrt(VX)), sX<-sX*Ar1X+sqrt(1-Ar1X*Ar1X)*rnorm(1,0,sqrt(VX)) ) }
        for (k in 1:Nindi) {
          if (SameEnv==FALSE) { sX<-rnorm(1,0,sqrt(VX)) } # no autocorrelation here      
          dataset[(Nindi*(j-1)+k),1]<-j               # store time     
          dataset[(Nindi*(j-1)+k),2]<-alive[k,1]      # save individual's ID
          dataset[(Nindi*(j-1)+k),3]<-sX              # save individual's X                                                                     
          dataset[(Nindi*(j-1)+k),4]<-beta0+alive[k,2]+(betaX+alive[k,3])*sX+rnorm(1,0,sqrt(VR))      # generate Y values   
          }
        }
      dataset<-as.data.frame(dataset)  
      model <- lmer(V4 ~ V3 + (V3 | V2), data=dataset) 
      if (alpha!=0) {
        model1 <- lm(V4 ~ V3, data=dataset)
        model2 <- lmer(V4 ~ V3 + (1 | V2), data=dataset) 
        modeloutput[i,1]<-(pchisq(-2 * (logLik(model1, REML = TRUE) - logLik(model2, REML = TRUE))[[1]], 1, lower.tail = FALSE))<alpha # P value random intercept
        modeloutput[i,2] <- anova(model, model2)[2, "Pr(>Chisq)"]<alpha    # P value random slope
        z <- fixef(model)/sqrt(diag(vcov(model, useScale = FALSE)))
        Pfixed <- 2 * (1 - pnorm(abs(z)))
        modeloutput[i,9] <- Pfixed[1]<alpha   # P value fixed intercept
        modeloutput[i,10] <- Pfixed[2]<alpha  # P value fixed slope
      }
      VC<-VarCorr(model)
      modeloutput[i,3]<-fixef(model)[1]   # estimate population intercept
      modeloutput[i,4]<-fixef(model)[2]   # estimate population slope
      modeloutput[i,5]<-attr(VC, "sc")^2    # estimate VR
      modeloutput[i,6]<-VC[[1]][1]          # estimate VI
      modeloutput[i,7]<-VC[[1]][4]          # estimate VS
      ifelse(correl==TRUE,modeloutput[i,8]<-attr(VC$V2,"correlation")[1,2],modeloutput[i,8]<-VC[[1]][2]) # estimate covariance/correlation intercept and slopes 
    }
    g[Scounter]<-Nindi
    r[Scounter]<-Nrep
    estInt[Scounter]<-(quantile(modeloutput[,3],0.5, na.rm=TRUE))
    estSlope[Scounter]<-(quantile(modeloutput[,4],0.5, na.rm=TRUE))
    estVR[Scounter]<-(quantile(modeloutput[,5],0.5, na.rm=TRUE))
    estVI[Scounter]<-(quantile(modeloutput[,6],0.5, na.rm=TRUE))
    estVS[Scounter]<-(quantile(modeloutput[,7],0.5, na.rm=TRUE))
    q1Int[Scounter]<-quantile(modeloutput[,3],Quant[1], na.rm=TRUE)
    q2Int[Scounter]<-quantile(modeloutput[,3],Quant[2], na.rm=TRUE)
    q3Int[Scounter]<-quantile(modeloutput[,3],Quant[3], na.rm=TRUE)
    q4Int[Scounter]<-quantile(modeloutput[,3],Quant[4], na.rm=TRUE)
    q1Slope[Scounter]<-quantile(modeloutput[,4],Quant[1], na.rm=TRUE)
    q2Slope[Scounter]<-quantile(modeloutput[,4],Quant[2], na.rm=TRUE)
    q3Slope[Scounter]<-quantile(modeloutput[,4],Quant[3], na.rm=TRUE)
    q4Slope[Scounter]<-quantile(modeloutput[,4],Quant[4], na.rm=TRUE)
    q1VR[Scounter]<-quantile(modeloutput[,5],Quant[1], na.rm=TRUE)
    q2VR[Scounter]<-quantile(modeloutput[,5],Quant[2], na.rm=TRUE)
    q3VR[Scounter]<-quantile(modeloutput[,5],Quant[3], na.rm=TRUE)
    q4VR[Scounter]<-quantile(modeloutput[,5],Quant[4], na.rm=TRUE)
    q1VI[Scounter]<-quantile(modeloutput[,6],Quant[1], na.rm=TRUE)
    q2VI[Scounter]<-quantile(modeloutput[,6],Quant[2], na.rm=TRUE)
    q3VI[Scounter]<-quantile(modeloutput[,6],Quant[3], na.rm=TRUE)
    q4VI[Scounter]<-quantile(modeloutput[,6],Quant[4], na.rm=TRUE)
    q1VS[Scounter]<-quantile(modeloutput[,7],Quant[1], na.rm=TRUE)
    q2VS[Scounter]<-quantile(modeloutput[,7],Quant[2], na.rm=TRUE)
    q3VS[Scounter]<-quantile(modeloutput[,7],Quant[3], na.rm=TRUE)
    q4VS[Scounter]<-quantile(modeloutput[,7],Quant[4], na.rm=TRUE)
    if(correl==TRUE) {
      estCorIS[Scounter]<-(quantile(modeloutput[,8],0.5, na.rm=TRUE))
      q1CorIS[Scounter]<-quantile(modeloutput[,8],Quant[1], na.rm=TRUE)
      q2CorIS[Scounter]<-quantile(modeloutput[,8],Quant[2], na.rm=TRUE)
      q3CorIS[Scounter]<-quantile(modeloutput[,8],Quant[3], na.rm=TRUE)
      q4CorIS[Scounter]<-quantile(modeloutput[,8],Quant[4], na.rm=TRUE)
      }
    if(correl==FALSE) {
      q1CovIS[Scounter]<-quantile(modeloutput[,8],Quant[1], na.rm=TRUE)
      q2CovIS[Scounter]<-quantile(modeloutput[,8],Quant[2], na.rm=TRUE)
      q3CovIS[Scounter]<-quantile(modeloutput[,8],Quant[3], na.rm=TRUE)
      q4CovIS[Scounter]<-quantile(modeloutput[,8],Quant[4], na.rm=TRUE)
      estCovIS[Scounter]<-(quantile(modeloutput[,8],0.5, na.rm=TRUE))
      }  
    ifelse(alpha!=0, powerRandInt[Scounter]<-mean(modeloutput[,1],na.rm=TRUE), powerRandInt[Scounter]<-NA)
    ifelse(alpha!=0, powerRandSlope[Scounter]<-mean(modeloutput[,2],na.rm=TRUE), powerRandSlope[Scounter]<-NA)  
    ifelse(alpha!=0, powerFixInt[Scounter]<-mean(modeloutput[,9],na.rm=TRUE), powerFixInt[Scounter]<-NA) 
    ifelse(alpha!=0, powerFixSlope[Scounter]<-mean(modeloutput[,10],na.rm=TRUE), powerFixSlope[Scounter]<-NA) 
  }
  if (correl==TRUE) {
    sims.sum <- data.frame(G=g,R=r,
    Q1Int=q1Int,Q2Int=q2Int,EstInt=estInt,Q3Int=q3Int,Q4Int=q4Int,    
    Q1Slope=q1Slope,Q2Slope=q2Slope,EstSlope=estSlope,Q3Slope=q3Slope,Q4Slope=q4Slope,  
    Q1VI=q1VI,Q2VI=q2VI,EstVI=estVI,Q3VI=q3VI,Q4VI=q4VI,  
    Q1VS=q1VS,Q2VS=q2VS,EstVS=estVS,Q3VS=q3VS,Q4VS=q4VS,  
    Q1CorIS=q1CorIS,Q2CorIS=q2CorIS,EstCorIS=estCorIS,Q3CorIS=q3CorIS,Q4CorIS=q4CorIS,  
    Q1VR=q1VR,Q2VR=q2VR,EstVR=estVR,Q3VR=q3VR,Q4VR=q4VR,
    PowerFixInt=powerFixInt, PowerFixSlope=powerFixSlope, PowerRandInt=powerRandInt, PowerRandSlope=powerRandSlope, 
    Int=beta0,Slope=betaX,VI=VI,VS=VS,CorIS=CorIS,VR=VR,FLAG=FALSE)
  }
  if (correl==FALSE) {
    sims.sum <- data.frame(G=g,R=r,
    Q1Int=q1Int,Q2Int=q2Int,EstInt=estInt,Q3Int=q3Int,Q4Int=q4Int,    
    Q1Slope=q1Slope,Q2Slope=q2Slope,EstSlope=estSlope,Q3Slope=q3Slope,Q4Slope=q4Slope,  
    Q1VI=q1VI,Q2VI=q2VI,EstVI=estVI,Q3VI=q3VI,Q4VI=q4VI,  
    Q1VS=q1VS,Q2VS=q2VS,EstVS=estVS,Q3VS=q3VS,Q4VS=q4VS,  
    Q1CovIS=q1CovIS,Q2CovIS=q2CovIS,EstCovIS=estCovIS,Q3CovIS=q3CovIS,Q4CovIS=q4CovIS,  
    Q1VR=q1VR,Q2VR=q2VR,EstVR=estVR,Q3VR=q3VR,Q4VR=q4VR,
    PowerFixInt=powerFixInt, PowerFixSlope=powerFixSlope, PowerRandInt=powerRandInt, PowerRandSlope=powerRandSlope,
    Int=beta0,Slope=betaX,VI=VI,VS=VS,CovIS=CovIS,VR=VR,FLAG=FALSE)
  }
  class(sims.sum) = c("odprism", "data.frame")
  sims.sum
}

