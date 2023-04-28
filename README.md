README
================

## Introduction

Modelling disease risk and survival using longitudinal risk factor
trajectories is of interest in various clinical scenarios. We propose a
dynamic risk score modeling framework for multiple longitudinal risk
factors and survival in the presence of dependent censoring, where both
events depend on participants’ post-baseline clinical progression. Our
model requires relatively few random effects regardless of the number of
longitudinal risk factors and can therefore accommodate multiple
longitudinal risk factors in a parsimonious manner.

## Installation

You can install the development version of DRS using:

``` r
library(devtools) 
install_github("chzhang25/DRS", quiet = TRUE) 
```

## Example

``` r
library(DRS) 
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

est <- out$coef
se <- out$se.coef
print(round(cbind(est, se)[-grep("gammas.bs", names(est)), ], 3)) 
```

## Reference

Reference will be added later if published
