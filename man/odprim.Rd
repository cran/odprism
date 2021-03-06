\name{odprim}
\alias{odprim}
 \concept{GLMM}
\title{
Simulation function to asses the optimal design (balanced datasets) for studies planning to use random intercept models to analyze their data.    
}
\description{
\code{odprim} is a simulation functions that allows users to vary the sampling design in terms of number of grouping units sampled (e.g. individuals, schools, populations, etc.) and replicates per grouping unit (limited to balanced datasets; all individuals have the same number of replicates) and also allow users to vary the parameter conditions used to generate the data. Subsequently, the performance of random intecept models of the form \code{Y~X+X2+(1|Individual)} fitted on these datasets is assessed in terms of the accuracy and the precision of estimates of fixed and random parameter, as well as the statistical power. 
Y values are generated using the formula Y=b0+bX*X+bX2*X2+u+e, where residual noise e is generated by a Gaussian random variable with mean zero and residual variance VR and the variation between grouping units in their intercept u is generated by a normal distribution with mean zero and varaince VI. X values are generated by a random normal variable with zero mean and variance VX, where X is a variable that varies both within and among grouping units, while X2 only varies among grouping units.  
}
\usage{
odprim(indiv, repl, sims = 1000, fixed, random, 
  Xvar = c(1, 0, FALSE), Xvar2 = 1, alpha = 0, 
  Quant = c(0.025, 0.25, 0.75, 0.975), Verbal = TRUE)
}
\arguments{
  \item{indiv}{ Number of individuals (or any other type of grouping units) sampled }
  \item{repl}{ Number of replicates per individual }
  \item{sims}{ Number of different simulated dataset that are generated. Default: \code{sims=1000} }
  \item{fixed}{ Vector of parameter values for b0, bX, and optionally bX2 which are used to generate trait values (i.e. mean intercept and mean slopes of X and X2). If fixed[3] is not provided (e.g. when \code{fixed=c(0,1)}) then the varaible X2 is removed from the model. }
  \item{random}{Vector of parameter values VI and VR, which are used to generate trait values (residual varaince and variance in intercepts). If no value is given for VR (e.g. when \code{random=0.2}) it is assumed VR=1-VI.}
  \item{Xvar}{ Vector describing properties used to generate X values. The first entry gives the value of VX, the second value the autocorrelation in X between subsequent samping occasions, if the third entry is \code{TRUE} all individuals sampled at a specific occasion are given the same X value, while if \code{FALSE} all individuals are given different values at a specific occasion.  Default: \code{Xvar=c(1,0,FALSE)}. }
  \item{Xvar2}{ Vector describing properties used to generate X2 values. The first entry gives the value of VX2, such that X2~Normal(0,VX2).}
  \item{alpha}{ If \code{alpha!=0} then the statistical power is calculated, defined as the proportion of analysed simulated dataset that resulted in P-value < \code{alpha}. P-values of fixed effects (b0 & bX) were generated using Wald-tests, while P-value of random effects were generated using likelihood rate tests. Specifically, the P-value for the random slopes term was obtained by comparing the deviance of full model \code{Y~X+(X|Individual)} with the reduced model \code{Y~X+(1||Individual)}, while the P-value for the random intercept term was obtained by comparing the deviance of full model \code{Y~X+(1|Individual)} with the reduced model \code{Y~X}.  Default: \code{alpha=0}.}
  \item{Quant}{ Each sampling design is simulated sims times. The argument \code{Quant} is used to describe which percentiles of the distribution of parameter estimates should be stored. Default: \code{Quant=c(0.025,0.25,0.75,0.975)}.}
  \item{Verbal}{ If \code{TRUE}, the function prints to the screen which combination of indiv and repl it is currently evaluating. Default: \code{Verbal=TRUE}}
}
\details{ Please be aware that statistical inference in terms of P-values for fixed and especially for random parameters is a contentious topic. The methods used to calculate P-values here (Wald for fixed effects and LRT test for random effects) are widely used, but not necessarily the best. More generally, instead of focussing on statisical power I would advice to focus more on the acuracy and precision of parameter estimates. Therefore the default options (\code{alpha=0}) do not generate descriptions of statistical power.  }
\value{ Dataframe of class 'odprism' that contains the chosen values used to generate the data as the distribution of estimates derived from the random regression models applied to these datasets.
For example the column 'VR' gives the value used to generate residual noise, while the columns 'Q1VR', 'Q2VR', 'EstVR', 'Q3VR' and 'Q4VR' give respectively the 2.5, 25, 50 (i.e. median), 75 and 97.5 percentiles of the parameters estimates of VR in all simulated dataset.  
See the help file of DataExample for more details.}
\references{ van de Pol, M. 2012 Quantifying individual reaction norms: How study design affects the accuracy, precision and power of random regression models. Methods in Ecology & Evolution, In press.}
\author{ Martijn van de Pol }

\seealso{
\code{\link{odprismB}}, \code{\link{odprismU}}, \code{\link{plot2d.prism}}, \code{\link{plot3d.prism}}, \code{\link{DataExample}}
}
\examples{
## Note that this example (model is Y~X+(1|Individual)) has only few sims   
results<-odprim(indiv=c(10,25,50,75,100), repl=c(2,4,6,8,10), 
  fixed=c(0,1), random=0.2, sims=10, alpha=0.05)
results
plot2d.prism(x=results, variable="bX", xvar="G", cons=2)

## The difference wht the above model is that here X2 is included, 
## which is a covariate that varies only among grouping units.   
results<-odprim(indiv=c(10,25,50,75,100), repl=c(2,4,6,8,10), 
  fixed=c(0,1,1), random=0.2, sims=10, alpha=0.05)
results
plot3d.prism(x=results, variable="VI")
}
\keyword{ models }
\keyword{ methods }