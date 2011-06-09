\name{DataExample}
\alias{DataExample}
\docType{data}
\title{
Example of results that can be obtained with package odprism
}
\description{
Example of results from a simulation obtained with the code:
   \code{DataExample<-odprismB( indiv=c(10,20,30,40,50,75,100,150,200), 
          repl=c(2,3,4,5,6,7,8,9,10), random=c(0.2,0.1,0.5), 
          alpha=0.05, Quant=c(0.025,0.25,0.75,0.975) )}.
}
\usage{data(DataExample)}
\format{
  A data frame with 81 observations on the following 43 variables.
  \describe{
    \item{\code{G}}{a numeric vector describing the number of grouping units sampled}
    \item{\code{R}}{a numeric vector describing the number of replicates per grouping unit sampled}
    \item{\code{Q1Int}}{a numeric vector describing the lowest quantile of the distribution of the parameter estimate of intercept beta0 (e.g. the 2.5-percentile)}  
    \item{\code{Q2Int}}{a numeric vector describing the lower quantile of the distribution of the parameter estimate of intercept beta0 (e.g. the 25-percentile)}
    \item{\code{EstInt}}{a numeric vector describing the median of the parameter estimate of intercept beta0 accross all simulated datasets (e.g. the 50-percentile)}
    \item{\code{Q3Int}}{a numeric vector describing the upper quantile of the distribution of the parameter estimate of intercept beta0 (e.g. the 75-percentile)}
    \item{\code{Q4Int}}{a numeric vector describing the highest quantile of the distribution of the parameter estimate of intercept beta0 (e.g. the 97.5-percentile)}
    \item{\code{Q1Slope}}{a numeric vector describing the lowest quantile of the distribution of the parameter estimate of slope betaX (e.g. the 2.5-percentile)}  
    \item{\code{Q2Slope}}{a numeric vector describing the lower quantile of the distribution of the parameter estimate of slope betaX (e.g. the 25-percentile)}
    \item{\code{EstSlope}}{a numeric vector describing the median of the parameter estimate of slope betaX accross all simulated datasets (e.g. the 50-percentile)}
    \item{\code{Q3Slope}}{a numeric vector describing the upper quantile of the distribution of the parameter estimate of slope betaX (e.g. the 75-percentile)}
    \item{\code{Q4Slope}}{a numeric vector describing the highest quantile of the distribution of the parameter estimate of slope betaX (e.g. the 97.5-percentile)}
     \item{\code{Q1VI}}{a numeric vector describing the lowest quantile of the distribution of the parameter estimate of VI (e.g. the 2.5-percentile)}   
    \item{\code{Q2VI}}{a numeric vector describing the lower quantile of the distribution of the parameter estimate of VI (e.g. the 25-percentile)}
    \item{\code{EstVI}}{a numeric vector describing the median of the parameter estimate of VI accross all simulated datasets (e.g. the 50-percentile)}
    \item{\code{Q3VI}}{a numeric vector describing the upper quantile of the distribution of the parameter estimate of VI (e.g. the 75-percentile)}
    \item{\code{Q4VI}}{a numeric vector describing the highest quantile of the distribution of the parameter estimate of VI (e.g. the 97.5-percentile)}
    \item{\code{Q1VS}}{a numeric vector describing the lowest quantile of the distribution of the parameter estimate of VS (e.g. the 2.5-percentile)}  
    \item{\code{Q2VS}}{a numeric vector describing the lower quantile of the distribution of the parameter estimate of VS (e.g. the 25-percentile)}
    \item{\code{EstVS}}{a numeric vector describing the median of the parameter estimate of VS accross all simulated datasets (e.g. the 50-percentile)}
    \item{\code{Q3VS}}{a numeric vector describing the upper quantile of the distribution of the parameter estimate of VS (e.g. the 75-percentile)}
    \item{\code{Q4VS}}{a numeric vector describing the highest quantile of the distribution of the parameter estimate of VS (e.g. the 97.5-percentile)}
    \item{\code{Q1CorIS}}{a numeric vector describing the lowest quantile of the distribution of the parameter estimate of CORIS (e.g. the 2.5-percentile)}   
    \item{\code{Q2CorIS}}{a numeric vector describing the lower quantile of the distribution of the parameter estimate of CORIS (e.g. the 25-percentile)}
    \item{\code{EstCorIS}}{a numeric vector describing the median of the parameter estimate of CORIS accross all simulated datasets (e.g. the 50-percentile)}
    \item{\code{Q3CorIS}}{a numeric vector describing the upper quantile of the distribution of the parameter estimate of CORIS (e.g. the 75-percentile)}
    \item{\code{Q4CorIS}}{a numeric vector describing the highest quantile of the distribution of the parameter estimate of CORIS (e.g. the 97.5-percentile)}
    \item{\code{Q1VR}}{a numeric vector describing the lowest quantile of the distribution of the parameter estimate of VR (e.g. the 2.5-percentile)}   
    \item{\code{Q2VR}}{a numeric vector describing the lower quantile of the distribution of the parameter estimate of VR (e.g. the 25-percentile)}   
    \item{\code{EstVR}}{a numeric vector describing the median of the parameter estimate of VR accross all simulated datasets (e.g. the 50-percentile)}
    \item{\code{Q3VR}}{a numeric vector describing the upper quantile of the distribution of the parameter estimate of VR (e.g. the 75-percentile)}
    \item{\code{Q4VR}}{a numeric vector describing the highest quantile of the distribution of the parameter estimate of VR (e.g. the 97.5-percentile)}
    \item{\code{PowerFixInt}}{a numeric vector describing the proportion of simulation that rejected the null hypothesis of no effect for beta0}
    \item{\code{PowerFixSlope}}{a numeric vector describing the proportion of simulation that rejected the null hypothesis of no effect for betaX}
    \item{\code{PowerRandInt}}{a numeric vector describing the proportion of simulation that rejected the null hypothesis of no effect for VI}
    \item{\code{PowerRandSlope}}{a numeric vector describing the proportion of simulation that rejected the null hypothesis of no effect for VS=CorIS=0}
    \item{\code{Int}}{a numeric vector describing the value of intercept beta0 used to generate the datasets}
    \item{\code{Slope}}{a numeric vector describing the value of slope betaX used to generate the datasets}
    \item{\code{VI}}{a numeric vector describing the value of the among-grouping units variance in intercepts VI used to generate the datasets}
    \item{\code{VS}}{a numeric vector describing the value of the among-grouping units variance in slopes VS used to generate the datasets}
    \item{\code{CorIS}}{a numeric vector describing the value of the correlation between intercepts and slopes CorIS used to generate the datasets}
    \item{\code{VR}}{a numeric vector describing the value of the within-grouping units (residual) variance VR used to generate the datasets}
    \item{\code{FLAG}}{a logical vector coding for balanced datasets \code{FLAG=FALSE} and unbalanced datasets \code{FLAG=TRUE}}
  }
}

\source{
van de Pol, M. (submitted) Quantifying individual reaction norms: How study design affects the accuracy, precision and power of random regression models.
}
\references{
van de Pol, M. (submitted) Quantifying individual reaction norms: How study design affects the accuracy, precision and power of random regression models.
}
\examples{
## plot DataExample in various ways
DataExample<-odprismB(indiv=c(10,25,50,75,100), repl=c(2,4,6,8,10), 
  random=c(0.2,0.1,0.5), sims=10, alpha=0.05)
plot2d.prism(x=DataExample, variable="VS", xvar="G", cons=10)
plot3d.prism(x=DataExample, variable="VS")
}
\keyword{datasets}
