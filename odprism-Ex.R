pkgname <- "odprism"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('odprism')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("DataExample")
### * DataExample

flush(stderr()); flush(stdout())

### Name: DataExample
### Title: Example of results that can be obtained with package odprism
### Aliases: DataExample
### Keywords: datasets

### ** Examples

## plot DataExample in various ways
DataExample<-odprismB(indiv=c(10,25,50,75,100), repl=c(2,4,6,8,10), 
  random=c(0.2,0.1,0.5), sims=10, alpha=0.05)
plot2d.prism(x=DataExample, variable="VS", xvar="G", cons=10)
plot3d.prism(x=DataExample, variable="VS")



cleanEx()
nameEx("odprim")
### * odprim

flush(stderr()); flush(stdout())

### Name: odprim
### Title: Simulation function to asses the optimal design (balanced
###   datasets) for studies planning to use random intercept models to
###   analyze their data.
### Aliases: odprim
### Keywords: models methods

### ** Examples

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



cleanEx()
nameEx("odprismB")
### * odprismB

flush(stderr()); flush(stdout())

### Name: odprismB
### Title: Simulation function to asses the optimal design (balanced
###   datasets) for studies planning to use random intercept and slopes
###   models to analyze their data.
### Aliases: odprismB
### Keywords: models methods

### ** Examples

## Note that this example has only few sims 
results<-odprismB(indiv=c(10,25,50,75,100), repl=c(2,4,6,8,10), 
  random=c(0.2,0.1,0.5), sims=10, alpha=0.05)
results
plot2d.prism(x=results, variable="VS", xvar="R", cons=10)

## Alternatively look at the example datafile DataExample, 
## which  is available with the package as an example)
## ensure that class(DataExample) = c("odprism", "data.frame")
## then run: plot3d.prism(x=DataExample, variable="VS") 



cleanEx()
nameEx("odprismU")
### * odprismU

flush(stderr()); flush(stdout())

### Name: odprismU
### Title: Simulation function to asses the optimal design (unbalanced
###   datasets) for studies planning to use random intercept and slopes
###   models to analyze their data consiting of annually measured traits.
### Aliases: odprismU
### Keywords: models methods

### ** Examples

## Note example uses only few sims to speed things up, normally sims>1000
results<-odprismU(pops=c(10,25,50,75,100), years=c(2,4,6,8,10), 
  random=c(0.2,0.1,0.5), Survival=0.7, sims=10, alpha=0.05)
results
plot2d.prism(x=results, variable="VS", xvar="G", cons=10)

## Alternatively look at an example datafile DataExample, 
## which  is available with the package as an example.
## ensure that class(DataExample) = c("odprism", "data.frame")
## then run: plot3d.prism(x=DataExample, variable="VS") 



cleanEx()
nameEx("plot2d.prism")
### * plot2d.prism

flush(stderr()); flush(stdout())

### Name: plot2d.prism
### Title: Function to plot the performance (accuracy, precision and
###   statistical power) as a function of the number of individuals
###   (x-axis) and replicates (y-axis) sampled.
### Aliases: plot2d.prism
### Keywords: models methods

### ** Examples

## Note example uses only few sims to speed things up, normally sims>1000
results<-odprismB(indiv=c(10,25,50,75,100), repl=c(2,4,6,8,10), 
  random=c(0.2,0.1,0.5), sims=10, alpha=0.05)
results
plot2d.prism(x=results, variable="VI", xvar="G", cons=10)

## Alternatively look at an example datafile DataExample, 
## which is available with the package as an example.
## ensure that class(DataExample) = c("odprism", "data.frame")
## then run: plot2d.prism(x=DataExample, variable="VI", xvar="R", cons=100)



cleanEx()
nameEx("plot3d.prism")
### * plot3d.prism

flush(stderr()); flush(stdout())

### Name: plot3d.prism
### Title: Function to plot the performance (accuracy, precision and
###   statistical power) as a function of the number of individuals
###   (x-axis) and replicates (y-axis) sampled.
### Aliases: plot3d.prism
### Keywords: models methods

### ** Examples

## Note example uses only few sims to speed things up, normally sims>1000
results<-odprismB(indiv=c(10,25,50,75,100), repl=c(2,4,6,8,10), 
  random=c(0.2,0.1,0.5), sims=10, alpha=0.05)
results
plot3d.prism(x=results, variable="C")

## Alternatively look at an example datafile DataExample, 
## which is available with the package as an example.
## ensure that class(DataExample) = c("odprism", "data.frame")
## then run plot3d.prism(x=DataExample, variable="C") 



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
