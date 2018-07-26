\name{plot.swgee}
\alias{plot.swgee}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.swgee
}
\description{
Produce the plot of the quadratic extrapolation curve for any covariables with measurement error in the swgee model}
\usage{
## S3 method for class 'swgee'
\method{plot}{swgee}(x, covariate, ...)
}

\arguments{
  \item{x}{
object of class 'swgee'
}
  \item{covariate}{
covariates specified in the formula
}
  \item{\dots}{
further arguments passed to or from other functions.
}
}


\value{
Plot the simulation and extrapolation step
}

\references{
Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. \emph{Journal of the American Statistical Association}, 89, 1314-1328.

Carrol, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu, C. (2006) Measurement error in nonlinear models: A modern perspective., Second Edition. London: Chapman and Hall.

Yi, G. Y. (2008) A simulation-based marginal method for longitudinal data with dropout and mismeasured covariates. \emph{Biostatistics}, 9, 501-512.
}

\author{
Juan Xiong<jxiong@szu.edu.cn>, Grace Y. Yi<yyi@uwaterloo.ca>
}


\examples{
require(gee)
require(mvtnorm)
data(BMI)
bmidata <- BMI

rho <- 0
sigma1 <- 0.5
sigma2 <- 0.5

sigma <- matrix(0,2,2)
sigma[1,1] <- sigma1*sigma1
sigma[1,2] <- rho*sigma1*sigma2
sigma[2,1] <- sigma[1,2]
sigma[2,2] <- sigma2*sigma2

set.seed(1000)
##swgee method ##########
output2 <- swgee(bbmi~sbp+chol+age, data = bmidata, id = id, 
            family = binomial(link="logit"),corstr = "independence", 
            missingmodel = O~bbmi+sbp+chol+age, SIMEXvariable = c("sbp","chol"), 
            SIMEX.err = sigma, repeated = FALSE, B = 20, lambda = seq(0, 2, 0.5))

summary(output2)

plot(output2,"sbp")


}





