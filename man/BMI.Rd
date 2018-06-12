\name{BMI}
\alias{BMI}
\docType{data}
\title{
BMI dataset
}
\description{
a subset of the the Framingham Heart Study Data. The data set consists of  measurements of 100 patients from a series of exams with 5 assessments for each individual. 
Individual's obesity status, age, systolic blood pressure (SBP) and cholesterol level (CHOL) are collected at each assessment. }
\usage{data("BMI")}
\format{
  A data frame with 500 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{a numeric vector for subject id}
    \item{\code{visit}}{a numeric vector for assessment time}
    \item{\code{age}}{a numeric vector of age}
    \item{\code{sbp}}{a numeric vector of systolic blood pressure}
    \item{\code{chol}}{a numeric vector of cholesterol level}
    \item{\code{bbmi}}{an indicator of obesity status(1=yes, 0=no)}
    \item{\code{O}}{an indicator of observed measurement(1=yes, 0=no)}
  }
}
\details{
The author thanks Boston University and the National Heart, Lung, and Blood Institute (NHLBI) for providing the data set from the Framingham Heart
Study (No. N01-HC-25195) in the illustration. The Framingham Heart Study is conducted and supported by the NHLBI in collaboration with 
Boston University. This package was not prepared in collaboration with investigators of the Framingham Heart Study and 
does not necessarily reflect the opinions or views of the Framingham Heart Study, Boston University, or NHLBI. }

\source{
Strug, L., Sun, L. and Corey, M. (2003). The genetics of cross-sectional and longitudinal body mass index. \emph{BMC Genetics} 4 (Suppl 1), S14

Yoo, Y. J., Huo, Y., Ning, Y., Gordon, D., Finch, S. and Mendell, N. R. (2003). Power of maximum HLOD tests to detect linkage to obesity genes. \emph{BMC Genetics} 4 (Suppl 1), S16.
}

\references{
Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. \emph{Journal of the American Statistical Association}, 89, 1314-1328.

Carrol, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu, C. (2006) Measurement error in nonlinear models: A modern perspective., Second Edition. London: Chapman and Hall.

Yi, G. Y. (2008) A simulation-based marginal method for longitudinal data with dropout and mismeasured covariates. \emph{Biostatistics}, 9, 501-512.
}

\examples{
data(BMI)
}

