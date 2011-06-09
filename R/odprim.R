odprim <-
function(indiv, repl, sims=1000, fixed, random, Xvar=c(1,0,FALSE), Xvar2=1, alpha=0, Quant=c(0.025,0.25,0.75,0.975), Verbal=TRUE) {
  VI<-random[1]
  ifelse((length(random)<2), VR<-1-VI, VR<-random[2]) 
  beta0<-fixed[1]
  betaX<-fixed[2]
  if(length(fixed)==3) { betaX2<-fixed[3]  }
  VX<-Xvar[1]
  Ar1X<-Xvar[2]
  SameEnv<-Xvar[3]
  VX2<-Xvar2
  if(SameEnv==FALSE && Ar1X>0) { stop("Sorry, autocorrelation not possible if SameEnv==FALSE") }
  Nrepl<-rep(repl,length(indiv))  
  Nindiv<-sort(rep(indiv,length(repl))) 
  modeloutput<-matrix(data = NA, nrow = sims, ncol = 9)    
  g<-numeric(length(Nindiv))
  r<-numeric(length(Nindiv))
  estInt<-numeric(length(Nindiv))
  estSlope<-numeric(length(Nindiv))
  estSlope2<-numeric(length(Nindiv))  
  estVR<-numeric(length(Nindiv))
  estVI<-numeric(length(Nindiv))
  q1Int<-numeric(length(Nindiv))
  q2Int<-numeric(length(Nindiv))
  q3Int<-numeric(length(Nindiv))
  q4Int<-numeric(length(Nindiv))
  q1Slope<-numeric(length(Nindiv))
  q2Slope<-numeric(length(Nindiv))
  q3Slope<-numeric(length(Nindiv))
  q4Slope<-numeric(length(Nindiv))
  q1Slope2<-numeric(length(Nindiv))
  q2Slope2<-numeric(length(Nindiv))
  q3Slope2<-numeric(length(Nindiv))
  q4Slope2<-numeric(length(Nindiv))
  q1VR<-numeric(length(Nindiv))
  q2VR<-numeric(length(Nindiv))
  q3VR<-numeric(length(Nindiv))
  q4VR<-numeric(length(Nindiv))
  q1VI<-numeric(length(Nindiv))
  q2VI<-numeric(length(Nindiv))
  q3VI<-numeric(length(Nindiv))
  q4VI<-numeric(length(Nindiv))
  powerFixInt<-numeric(length(Nindiv))
  powerFixSlope<-numeric(length(Nindiv))
  powerFixSlope2<-numeric(length(Nindiv))
  powerRandInt<-numeric(length(Nindiv))
  Scounter<-0
  for (m in 1:(length(Nindiv))) {
    Scounter<-Scounter+1
    Nindi<-Nindiv[m]
    Nrep<-Nrepl[m]
    if (Verbal==TRUE) { print(c("Now running individuals",Nindi,"and replicates",Nrep))}
    for (i in 1:sims) {
      dataset<-matrix(data = NA, nrow =Nindi*Nrep, ncol = 5)
      alive<-matrix(data = NA, nrow = Nindi, ncol = 3) 
      alive[,1]<-seq(1,Nindi)
      alive[,2]<-rnorm(Nindi,0,sqrt(VI))   # elevation of individuals
      if(length(fixed)==3) { alive[,3]<-rnorm(Nindi,0,sqrt(VX2))}   # values of X for each grouping unit (e.g. birthweight)
      IDcounter<-Nindi
      for (j in 1:Nrep) {
        # select values for X in each sampling occasion      
        if (SameEnv==TRUE) { ifelse (j==1, sX<-rnorm(1,0,sqrt(VX)), sX<-sX*Ar1X+sqrt(1-Ar1X*Ar1X)*rnorm(1,0,sqrt(VX)) ) }
        for (k in 1:Nindi) {
          if (SameEnv==FALSE) { sX<-rnorm(1,0,sqrt(VX)) } # no autocorrelation here      
          dataset[(Nindi*(j-1)+k),1]<-j               # store time     
          dataset[(Nindi*(j-1)+k),2]<-alive[k,1]      # save individual's ID
          dataset[(Nindi*(j-1)+k),3]<-sX              # save individual's X                                                                     
          dataset[(Nindi*(j-1)+k),4]<-alive[k,3]       # save individual's X2                                                                     
          if(length(fixed)<3)  { dataset[(Nindi*(j-1)+k),5]<-beta0+alive[k,2]+betaX*sX+rnorm(1,0,sqrt(VR))}   # generate Y values   
          if(length(fixed)==3) { dataset[(Nindi*(j-1)+k),5]<-beta0+alive[k,2]+betaX*sX+betaX2*alive[k,3]+rnorm(1,0,sqrt(VR)) } # generate Y values   
        }
      }
      dataset<-as.data.frame(dataset)  
      if(length(fixed)<3) {
        model <- lmer(V5 ~ V3 + (1 | V2), data=dataset) 
        if (alpha!=0) {
          model1 <- lm(V5 ~ V3, data=dataset)
          modeloutput[i,1]<-(pchisq(-2 * (logLik(model1, REML = TRUE) - logLik(model, REML = TRUE))[[1]], 1, lower.tail = FALSE))<alpha # P value random intercept
          z <- fixef(model)/sqrt(diag(vcov(model, useScale = FALSE)))
          Pfixed <- 2 * (1 - pnorm(abs(z)))
          modeloutput[i,2]  <- Pfixed[1]<alpha   # P value fixed intercept beta0
          modeloutput[i,3] <- Pfixed[2]<alpha  # P value fixed slope betaX
        }
        VC<-VarCorr(model)
        modeloutput[i,4]<-fixef(model)[1]   # estimate population intercept beta0
        modeloutput[i,5]<-fixef(model)[2]   # estimate population slope betaX
        modeloutput[i,6]<-attr(VC, "sc")^2  # estimate VR
        modeloutput[i,7]<-VC[[1]][1]        # estimate VI
      }
      if(length(fixed)==3) {
        model <- lmer(V5 ~ V3 + V4 + (1 | V2), data=dataset) 
        if (alpha!=0) {
          model1 <- lm(V5 ~ V3 + V4, data=dataset)
          modeloutput[i,1]<-(pchisq(-2 * (logLik(model1, REML = TRUE) - logLik(model, REML = TRUE))[[1]], 1, lower.tail = FALSE))<alpha # P value random intercept
          z <- fixef(model)/sqrt(diag(vcov(model, useScale = FALSE)))
          Pfixed <- 2 * (1 - pnorm(abs(z)))
          modeloutput[i,2]  <- Pfixed[1]<alpha   # P value fixed intercept beta0
          modeloutput[i,3] <- Pfixed[2]<alpha  # P value fixed slope betaX
          modeloutput[i,8] <- Pfixed[3]<alpha  # P value fixed slope betaX2
        }
        VC<-VarCorr(model)
        modeloutput[i,4]<-fixef(model)[1]   # estimate population intercept beta0
        modeloutput[i,5]<-fixef(model)[2]   # estimate population slope betaX
        modeloutput[i,9]<-fixef(model)[3]   # estimate population slope betaX2
        modeloutput[i,6]<-attr(VC, "sc")^2  # estimate VR
        modeloutput[i,7]<-VC[[1]][1]        # estimate VI
      }
    }
    g[Scounter]<-Nindi
    r[Scounter]<-Nrep
    estInt[Scounter]<-(quantile(modeloutput[,4],0.5, na.rm=TRUE))
    estSlope[Scounter]<-(quantile(modeloutput[,5],0.5, na.rm=TRUE))
    estSlope2[Scounter]<-(quantile(modeloutput[,9],0.5, na.rm=TRUE))
    estVR[Scounter]<-(quantile(modeloutput[,6],0.5, na.rm=TRUE))
    estVI[Scounter]<-(quantile(modeloutput[,7],0.5, na.rm=TRUE))
    q1Int[Scounter]<-quantile(modeloutput[,4],Quant[1], na.rm=TRUE)
    q2Int[Scounter]<-quantile(modeloutput[,4],Quant[2], na.rm=TRUE)
    q3Int[Scounter]<-quantile(modeloutput[,4],Quant[3], na.rm=TRUE)
    q4Int[Scounter]<-quantile(modeloutput[,4],Quant[4], na.rm=TRUE)
    q1Slope[Scounter]<-quantile(modeloutput[,5],Quant[1], na.rm=TRUE)
    q2Slope[Scounter]<-quantile(modeloutput[,5],Quant[2], na.rm=TRUE)
    q3Slope[Scounter]<-quantile(modeloutput[,5],Quant[3], na.rm=TRUE)
    q4Slope[Scounter]<-quantile(modeloutput[,5],Quant[4], na.rm=TRUE)
    q1Slope2[Scounter]<-quantile(modeloutput[,9],Quant[1], na.rm=TRUE)
    q2Slope2[Scounter]<-quantile(modeloutput[,9],Quant[2], na.rm=TRUE)
    q3Slope2[Scounter]<-quantile(modeloutput[,9],Quant[3], na.rm=TRUE)
    q4Slope2[Scounter]<-quantile(modeloutput[,9],Quant[4], na.rm=TRUE)
    q1VR[Scounter]<-quantile(modeloutput[,6],Quant[1], na.rm=TRUE)
    q2VR[Scounter]<-quantile(modeloutput[,6],Quant[2], na.rm=TRUE)
    q3VR[Scounter]<-quantile(modeloutput[,6],Quant[3], na.rm=TRUE)
    q4VR[Scounter]<-quantile(modeloutput[,6],Quant[4], na.rm=TRUE)
    q1VI[Scounter]<-quantile(modeloutput[,7],Quant[1], na.rm=TRUE)
    q2VI[Scounter]<-quantile(modeloutput[,7],Quant[2], na.rm=TRUE)
    q3VI[Scounter]<-quantile(modeloutput[,7],Quant[3], na.rm=TRUE)
    q4VI[Scounter]<-quantile(modeloutput[,7],Quant[4], na.rm=TRUE)
    ifelse(alpha!=0, powerRandInt[Scounter]<-mean(modeloutput[,1],na.rm=TRUE), powerRandInt[Scounter]<-NA) 
    ifelse(alpha!=0, powerFixInt[Scounter]<-mean(modeloutput[,2],na.rm=TRUE), powerFixInt[Scounter]<-NA) 
    ifelse(alpha!=0, powerFixSlope[Scounter]<-mean(modeloutput[,3],na.rm=TRUE), powerFixSlope[Scounter]<-NA) 
    ifelse(alpha!=0, powerFixSlope2[Scounter]<-mean(modeloutput[,8],na.rm=TRUE), powerFixSlope2[Scounter]<-NA) 
  }
  if (length(fixed)==3) {
    sims.sum <- data.frame(G=g,R=r,
    Q1Int=q1Int,Q2Int=q2Int,EstInt=estInt,Q3Int=q3Int,Q4Int=q4Int,    
    Q1Slope=q1Slope,Q2Slope=q2Slope,EstSlope=estSlope,Q3Slope=q3Slope,Q4Slope=q4Slope,  
    Q1Slope2=q1Slope2,Q2Slope2=q2Slope2,EstSlope2=estSlope2,Q3Slope2=q3Slope2,Q4Slope2=q4Slope2, 
    Q1VI=q1VI,Q2VI=q2VI,EstVI=estVI,Q3VI=q3VI,Q4VI=q4VI,  
    Q1VR=q1VR,Q2VR=q2VR,EstVR=estVR,Q3VR=q3VR,Q4VR=q4VR,
    PowerFixInt=powerFixInt, PowerFixSlope=powerFixSlope, PowerFixSlope2=powerFixSlope2, PowerRandInt=powerRandInt, 
    Int=beta0,Slope=betaX,Slope2=betaX2,VI=VI,VR=VR,FLAG=FALSE)
  }
  if (length(fixed)<3) {
    sims.sum <- data.frame(G=g,R=r,
    Q1Int=q1Int,Q2Int=q2Int,EstInt=estInt,Q3Int=q3Int,Q4Int=q4Int,    
    Q1Slope=q1Slope,Q2Slope=q2Slope,EstSlope=estSlope,Q3Slope=q3Slope,Q4Slope=q4Slope,  
    Q1VI=q1VI,Q2VI=q2VI,EstVI=estVI,Q3VI=q3VI,Q4VI=q4VI,  
    Q1VR=q1VR,Q2VR=q2VR,EstVR=estVR,Q3VR=q3VR,Q4VR=q4VR,
    PowerFixInt=powerFixInt, PowerFixSlope=powerFixSlope, PowerRandInt=powerRandInt, 
    betaX=betaX,Int=beta0,Slope=betaX,VI=VI,VR=VR,FLAG=FALSE)
  }
  class(sims.sum) = c("odprism", "data.frame")
  sims.sum
}

