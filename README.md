README
================

## Introduction

XXX

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

est <- out$coefficients
se <- out$se.coef
print(round(cbind(est, se)[-grep("gammas.bs", names(est)), ], 3)) 
```

## Reference
