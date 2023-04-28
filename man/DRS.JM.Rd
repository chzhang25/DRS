\name{DRS.JM}
\alias{DRS.JM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Dynamic Risk Score Modeling for Multiple Longitudinal Risk Factors and Survival
}
\description{
  Dynamic Risk Score Modeling for Multiple Longitudinal Risk Factors and Survival Subject to Both Independent and Dependent Censoring
}
\usage{
DRS.JM(coxForm, jmFixedForm, jmRandomForm, jmCoxForm1, jmCoxForm2, timeVar,
       data.long, data.id, kType = "epan",
       coxControl = list(), jmControl = list())
}
\arguments{
  \item{coxForm}{a formula to integrate the multiple longitudinal risk factors}
  \item{jmFixedForm}{a formula to specify the fixed effect of the longitudinal submodel in the joint model}
  \item{jmRandomForm}{a formula to specify the random effect of the longitudinal submodel in the joint model}
  \item{jmCoxForm1}{a formula to specify the main event submodel in the joint model}
  \item{jmCoxForm2}{a formula to specify the dependent censoring event submodel in the joint model}
  \item{timeVar}{a character to specify the time variable for the longitudinal risk factors}
  \item{data.long}{a data.frame for the longitudinal data}
  \item{data.id}{a data.frame for the survival data, i.e., one row per subject}
  \item{kType}{the smoothing kernel used in the estimation function. Must be one of {"epan", "gauss"}, where "epan" is the Epanechnikov kernel and "gauss" is the Gaussian kernel.}
  \item{coxControl}{control parameters for convergence criteria}
  \item{jmControl}{control parameters for the joint model}
}
\value{
A list is returned with the following components:
  \item{coef}{the estimated coefficients}
  \item{se.coef}{the standard errors of the estimated coefficient}
  \item{convergence}{two logic values to show the absolute and relative convergence of the algorithm}
  \item{coef.iter}{the estimated coefficients at each iteration}
}
\references{
Rizopoulos, D. (2010). JM: An R Package for the Joint Modelling of Longitudinal and Time-to-Event Data. Journal of Statistical Software, 35(9), 1–33. https://doi.org/10.18637/jss.v035.i09

Hongyuan Cao, Mathew M. Churpek, Donglin Zeng & Jason P. Fine (2015) Analysis of the Proportional Hazards Model With Sparse Longitudinal Covariates, Journal of the American Statistical Association, 110:511, 1187-1196, DOI: 10.1080/01621459.2014.957289
}
\author{
Cuihong Zhang, Ruosha Li <Ruosha.Li@uth.tmc.edu>
}
\examples{
## load required data
data(simda)
head(da.long)
da.id <- subset(da.long, v.time==0)

## fit the proposed model
out <- DRS.JM(coxForm = ~ Lt1 + Lt2 + Lt3 + Lt4 + Lt5,
              jmFixedForm = ~ v.time + W,
              jmRandomForm = ~ v.time | id,
              jmCoxForm1 = Surv(obsT, delta==1) ~ W,
              jmCoxForm2 = Surv(obsT, delta==2) ~ W,
              timeVar = "v.time", data.long = da.long, data.id = da.id,
              coxControl = list(n.iter = 30))

est <- out$coefficients
se <- out$se.coef
print(round(cbind(est, se)[-grep("gammas.bs", names(est)), ], 3))
}
