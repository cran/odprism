\name{odprismU}
\alias{odprismU}
 \concept{GLMM}
\title{Simulation function to asses the optimal design (unbalanced datasets) for studies planning to use random intercept and slopes models to analyze their data consiting of annually measured traits. }
\description{
In cases where traits are only measured or expressed once a year, the number of replicates per individual one can obtain is limited by their lifespan. odprismU is a simulation function that allows users to vary the sampling design in terms of population size being monitored and the number of years (occasions) this population is being sampled. The mean annual survival of individuals can be varied too (dead individuals are replaced), which together with the population size and sampling duration directly affects how many individuals are being sampled and how many replicates per individuals (i.e. unbalanced data).
Finally, one can choose to expose the trait of interest Y to viability selection, to determine whether this will bias the outcome of parameter estimates obtianed using standard random regression models. 
}
\usage{
odprismU(pops, years, sims = 1000, fixed = c(0, 0), correl = TRUE, 
  random, Survival, Xvar = c(1, 0, FALSE), alpha = 0, 
  Quant = c(0.025, 0.25, 0.75, 0.975), ViabilitySelection = 0, 
  Verbal = TRUE)
}
\arguments{
  \item{pops}{ Population size  being sampled at each given occasion. }
  \item{years}{ The number of sampling occasions ('years') the population is being sampled }
  \item{sims}{  Number of different simulated dataset that are generated. Default: \code{sims=1000} }
  \item{fixed}{ Vector of parameter values b0 and bX, which are used to generate trait values (i.e. mean intercept and mean slope). Default \code{fixed=c(0,0)}.}
  \item{correl}{If \code{correl=TRUE} then the third entry in the vector of argument random is taken to be CorIS (correlation intercept and slopes), if \code{FALSE} then it is taken to be CovIS (covariance intercept and slopes). Default: \code{correl=TRUE}.}
  \item{random}{ Vector of parameter values VI, VS, CorIS/CovIS and VR, which are used to generate trait values (residual varaince and (co)variance in intercepts and slopes). If no value is given for VR (e.g. when \code{random=c(0.2,0.1,0.5)}) it is assumed VR=1-VI.}
  \item{Survival}{ The mean annual survival of individuals, determining how many of the individuals measured at a given occasion can be sampled again at the next occasion. All dead individuals are replaced with new indviduals (without inheritance of trait values). }
  \item{Xvar}{Vector describing properties used to generate X values. The first entry gives the value of VX, the second value the autocorrelation in X between subsequent samping occasions, if the third entry is \code{TRUE} all individuals sampled at a specific occasion are given the same X value, while if \code{FALSE} all individuals are given different values at a specific occasion.  Default: \code{Xvar=c(1,0,FALSE)}.}
  \item{alpha}{ If \code{alpha!=0} then the statistical power is calculated, defined as the proportion of analysed simulated dataset that resulted in P-value < \code{alpha}. P-values of fixed effects (b0 & bX) were generated using Wald-tests, while P-value of random effects were generated using likelihood rate tests. Specifically, the P-value for the random slopes term was obtained by comparing the deviance of full model \code{Y~X+(X|Individual)} with the reduced model \code{Y~X+(1||Individual)}, while the P-value for the random intercept term was obtained by comparing the deviance of full model \code{Y~X+(1|Individual)} with the reduced model \code{Y~X}.  Default: \code{alpha=0}.}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  \item{Quant}{ Each sampling design is simulated sims times. The argument \code{Quant} is used to describe which percentiles of the distribution of parameter estimates should be stored. Default: \code{Quant=c(0.025,0.25,0.75,0.975)}.}
  \item{ViabilitySelection}{ Determines the strength of directional viability selection on trait Y, using the function P(survival)=1/(1+exp(-1(g0+ViabilitySelection*Y)). 
    When \code{ViabilitySelection=0} there is no selection on Y. The parameter g0 is determined by the argument \code{Survival}, such that P(survival) at Y=0 equals the chosen values of \code{Survival}, see also the automatically generated plotted survival function.}
  \item{Verbal}{ If \code{TRUE}, the function prints to the screen which combination of indiv and repl it is currently evaluating. Default: \code{Verbal=TRUE}}
}
\details{Please be aware that statistical inference in terms of P-values for fixed and especially for random parameters is a contentious topic. The methods used to calculate P-values here (Wald for fixed effects and LRT test for random effects) are widely used, but not necessarily the best. More generally, instead of focussing on statisical power I would advice to focus more on the acuracy and precision of parameter estimates. Therefore the default options (\code{alpha=0}) do not generate descriptions of statistical power.}
\value{
Dataframe of class 'odprism' that contains the chosen values used to generate the data as the distribution of estimates derived from the random regression models applied to these datasets.
For example the column 'VR' gives the value used to generate residual noise, while the columns 'Q1VR', 'Q2VR', 'EstVR', 'Q3VR' and 'Q4VR' give respectively the 2.5, 25, 50 (i.e. median), 75 and 97.5 percentiles of the parameters estimates of VR in all simulated dataset.  
See the help file of DataExample for more details.

}
\references{
van de Pol, M. 2012 Quantifying individual reaction norms: How study design affects the accuracy, precision and power of random regression models. Methods in Ecology & Evolution, In press.
}
\author{
Martijn van de Pol 
}

\seealso{
\code{\link{odprim}}, \code{\link{odprismU}}, \code{\link{plot2d.prism}}, \code{\link{plot3d.prism}}, \code{\link{DataExample}}
}
\examples{
## Note example uses only few sims to speed things up, normally sims>1000
results<-odprismU(pops=c(10,25,50,75,100), years=c(2,4,6,8,10), 
  random=c(0.2,0.1,0.5), Survival=0.7, sims=10, alpha=0.05)
results
plot2d.prism(x=results, variable="VS", xvar="G", cons=10)

## Alternatively look at an example datafile DataExample, 
## which  is available with the package as an example.
## ensure that class(DataExample) = c("odprism", "data.frame")
## then run: plot3d.prism(x=DataExample, variable="VS") 
}
\keyword{ models}
\keyword{ methods }
