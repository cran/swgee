\name{swgee}
\alias{swgee}
\title{
Simulation Extrapolation Inverse Probability Weighted Generalized Estimating Equations 
}
\description{
Implementation of the SIMEX inverse probability weighted GEE method for longitudinal data with missing observations and measurement error in covariates
}
\usage{
swgee(formula, data = parent.frame(), id, family = family, 
    corstr = "independence", missingmodel, SIMEXvariable, SIMEX.err, 
    repeated = FALSE, repind = NULL, B = 50, lambda = seq(0, 2, 0.5))
}
\arguments{
  \item{formula}{
specifies the model to be fitted, with the variables coming with data. This argument has the same format as the formula argument in the function \code{geeglm} from
package \pkg{geepack}, of the form \code{response ~ predictors}. See documentation of \code{geeglm} and \code{formula} for details.}
  \item{data}{an optional data frame in which to interpret the variables occurring in the formula, along with the id variable.}
  \item{id}{
a vector which identifies the clusters. The length of id should be the same as the number of observations. Data are assumed to be sorted so that observations 
on a cluster are contiguous rows for all entities in the formula.}
  \item{family}{
a family object as the family argument in the function \code{gee} from package \pkg{gee}. Families supported in \code{swgee} are \code{gaussian}, 
\code{binomial}, \code{poisson}, \code{Gamma}, and \code{quasi}. See documentation of \code{gee} and \code{family} for details.}
  \item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{"independence"}, \code{"fixed"}, \code{"stat_M_dep"}, 
\code{"non_stat_M_dep"}, \code{"exchangeable"}, \code{"AR-M"} and \code{"unstructured"}.}
  \item{missingmodel}{specifies the misisng model to be fitted, of the form \code{O~ predictors}, where O is the missing data indicator.}
  \item{SIMEXvariable}{a vector of characters containing the name of the covariates subject to measurement error.}
  \item{SIMEX.err}{specifies the covariance matrix of measurement errors in error model.}  
  \item{repeated}{This is the indicator if there are repeated measurements for the covariates with measurement error. The default value is FALSE.}
  \item{repind}{This is the index of the repeated measurement variables for each covariate with measurement error. It has an R list form. If repeated = TRUE, repind must be specified.}
  \item{B}{the number of simulated samples for the simulation step. The default is set to be 50.}
  \item{lambda}{a vector of lambdas for which the simulation step should be done.}
}
\details{
The quadratic extrapolation method is implemented as described in Cook and Stefanski
}

\value{
    \item{call}{the function call}
    \item{family}{family}
    \item{corstr}{correlation structure}
    \item{SIMEXvariable}{a vector of characters containing the name of the covariates subject to measurement error}
    \item{B}{the number of iterations}
    \item{beta}{the coefficients associated with the response process}
    \item{alpha}{the coefficients associated with the missing process}
    \item{simex.plot}{the estimates for every B and lambda}
}

\references{
Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. \emph{Journal of the American Statistical Association}, 89, 1314-1328.

Carrol, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu, C. (2006) Measurement error in nonlinear models: A modern perspective., Second Edition. London: Chapman and Hall.

Yi, G. Y. (2008) A simulation-based marginal method for longitudinal data with dropout and mismeasured covariates. \emph{Biostatistics}, 9, 501-512.
}

\author{
Juan Xiong<jxiong@szu.edu.cn>, Grace Y. Yi<yyi@uwaterloo.ca>
}

\seealso{\code{\link[geepack]{geeglm}}}

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
##naive method, ignore missingness and measurement error 
output1 <- gee(bbmi~sbp+chol+age, id = id, data = bmidata,
            family = binomial(link="logit"), corstr = "independence")

##swgee method ##########
output2 <- swgee(bbmi~sbp+chol+age, data = bmidata, id = id, 
            family = binomial(link="logit"),corstr = "independence", 
            missingmodel = O~bbmi+sbp+chol+age, SIMEXvariable = c("sbp","chol"), 
            SIMEX.err = sigma, repeated = FALSE, B = 20, lambda = seq(0, 2, 0.5))

summary(output2)
}





