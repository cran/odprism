\name{odprism-package}
\alias{odprism-package}
\alias{odprism}
\docType{package}
\title{ The package 'odprism' can be used to determine the Optimal Design and Performane of Random Intercept and Slope Models }
\description{Simulation functions to assess the optimal design for studies planning to use random intercept and/or slopes models to analyze their data, such as reaction norm studies. The functions 'odprim', 'odprismB' and odprismU' allow users to simulate many dataset while varying the sampling design in terms of number of grouping units sampled (e.g. individuals, schools, populations, etc.) and replicates per grouping unit (balanced as well as unbalanced designs) and also allow users to vary the parameter conditions used to generate the data. Subsequently, the performance of mixed models (based on package "lme4" )fitted on these datasets is assessed in terms of the accuracy and the precision of estimates of fixed and random parameter, as well as the statistical power. Results can be visualised using the functions 'plot2d.prism' and 'plot3d.prism'.}
\details{   
  \tabular{ll}{
  Package: \tab odprism\cr
  Type: \tab Package\cr
  Version: \tab 1.0\cr
  Date: \tab 2011-06-11\cr
  License: \tab GPL 2.00\cr
  LazyLoad: \tab yes\cr
  }
}
\author{
Martijn van de Pol
}
\references{
van de Pol, M. 2012 Quantifying individual reaction norms: How study design affects the accuracy, precision and power of random regression models. Methods in Ecology & Evolution, In press.
}
\keyword{ package }
\seealso{ \code{\link{odprim}}, \code{\link{odprismB}}, \code{\link{odprismU}}, \code{\link{plot2d.prism}}, \code{\link{plot3d.prism}},  \code{\link{DataExample}}}

