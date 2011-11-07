odprismU <-
function(pops, years, sims=1000, fixed=c(0,0), correl=TRUE, random, Survival, Xvar=c(1,0,FALSE), alpha=0, Quant=c(0.025,0.25,0.75,0.975), ViabilitySelection=0, Verbal=TRUE) {
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
  Nyears<-rep(years,length(pops))  
  Npops<-sort(rep(pops,length(years))) 
  Sbeta0<-(-log(1/Survival-1))
  modeloutput<-matrix(data = NA, nrow = sims, ncol = 10)    
  g<-numeric(length(Npops))
  r<-numeric(length(Npops))
  estInt<-numeric(length(Npops))
  estSlope<-numeric(length(Npops))
  estVR<-numeric(length(Npops))
  estVI<-numeric(length(Npops))
  estVS<-numeric(length(Npops))
  estCorIS<-numeric(length(Npops))
  estCovIS<-numeric(length(Npops))
  q1Int<-numeric(length(Npops))
  q2Int<-numeric(length(Npops))
  q3Int<-numeric(length(Npops))
  q4Int<-numeric(length(Npops))
  q1Slope<-numeric(length(Npops))
  q2Slope<-numeric(length(Npops))
  q3Slope<-numeric(length(Npops))
  q4Slope<-numeric(length(Npops))
  q1VR<-numeric(length(Npops))
  q2VR<-numeric(length(Npops))
  q3VR<-numeric(length(Npops))
  q4VR<-numeric(length(Npops))
  q1VI<-numeric(length(Npops))
  q2VI<-numeric(length(Npops))
  q3VI<-numeric(length(Npops))
  q4VI<-numeric(length(Npops))
  q1VS<-numeric(length(Npops))
  q2VS<-numeric(length(Npops))
  q3VS<-numeric(length(Npops))
  q4VS<-numeric(length(Npops))
  if(correl==TRUE) {
    q1CorIS<-numeric(length(Npops))
    q2CorIS<-numeric(length(Npops))
    q3CorIS<-numeric(length(Npops))
    q4CorIS<-numeric(length(Npops))
  }
  if(correl==FALSE) {
    q1CovIS<-numeric(length(Npops))
    q2CovIS<-numeric(length(Npops))
    q3CovIS<-numeric(length(Npops))
    q4CovIS<-numeric(length(Npops))
  }  
  powerFixInt<-numeric(length(Npops))
  powerFixSlope<-numeric(length(Npops))
  powerRandInt<-numeric(length(Npops))
  powerRandSlope<-numeric(length(Npops))   
  dev.new() # rather than x11()
  plot(seq(-4*(VR+VI),4*(VR+VI),0.1),(1/(1+exp(-1*(Sbeta0+ViabilitySelection*seq(-4*(VR+VI),4*(VR+VI),0.1))))), type="l", xlim=c(-3*(VR+VI),3*(VR+VI)), ylim=c(0,1), xlab="Y value", ylab="annual survival probability", main="viability selection on trait Y" ) 
  Scounter<-0
  for (m in 1:(length(Npops))) {
    Scounter<-Scounter+1
    Npop<-Npops[m]
    Nyear<-Nyears[m]
    if (Verbal==TRUE) { print(c("Now running population size ",Npop,"and study period",Nyear))}
    for (i in 1:sims) {
      dataset<-matrix(data = NA, nrow =Npop*Nyear, ncol = 4)
      alive<-matrix(data = NA, nrow = Npop, ncol = 3) 
      alive[,1]<-seq(1,Npop)
      IndHet <- rmvnorm(Npop, c(0, 0), sigma=matrix(c(VI,CovIS, CovIS, VS), ncol=2), method = "svd")  
      alive[,2]<-IndHet[,1]   # elevation of individuals
      alive[,3]<-IndHet[,2]   # slope of individuals
      IDcounter<-Npop
      for (j in 1:Nyear) {
        # select values for X in each year       
        if (SameEnv==TRUE) { ifelse (j==1, sX<-rnorm(1,0,sqrt(VX)), sX<-sX*Ar1X+sqrt(1-Ar1X*Ar1X)*rnorm(1,0,sqrt(VX)) ) }
        for (k in 1:Npop) {
          if (SameEnv==FALSE) { sX<-rnorm(1,0,sqrt(VX)) } # no autocorrelation here      
          dataset[(Npop*(j-1)+k),1]<-j               # store time     
          dataset[(Npop*(j-1)+k),2]<-alive[k,1]      # save individual's ID
          dataset[(Npop*(j-1)+k),3]<-sX              # save individual's X                                                                     
          dataset[(Npop*(j-1)+k),4]<-beta0+alive[k,2]+(betaX+alive[k,3])*sX+rnorm(1,0,sqrt(VR))      # generate Y values   
          if(Survival<1) {
           SurvProb<-1/(1+exp(-1*(Sbeta0+ViabilitySelection*dataset[k,4]))) # viability selection function
            if (rbinom(1,1,SurvProb)==0) { 
              IDcounter<-IDcounter+1
              IndHet2 <- rmvnorm(1, c(0, 0), sigma=matrix(c(VI,CovIS, CovIS, VS), ncol=2), method = "svd")       
              alive[k,1]<-IDcounter   
              alive[k,2]<-IndHet2[,1]   # elevation of individual
              alive[k,3]<-IndHet2[,2]   # slope of individual    
            }
          }
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
    g[Scounter]<-Npop
    r[Scounter]<-Nyear
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
    Int=beta0,Slope=betaX,VI=VI,VS=VS,CorIS=CorIS,VR=VR,FLAG=TRUE)
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
    Int=beta0,Slope=betaX,VI=VI,VS=VS,CovIS=CovIS,VR=VR,FLAG=TRUE)
  }
  class(sims.sum) = c("odprism", "data.frame")
  sims.sum
}

