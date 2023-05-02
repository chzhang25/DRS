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

The development version of package `DRS` can be installed using:

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
da.id <- subset(da.long, v.time == 0)

## fit the proposed model
out <- DRS.JM(coxForm = ~ Lt1 + Lt2 + Lt3 + Lt4 + Lt5,
              jmFixedForm = ~ v.time + W,
              jmRandomForm = ~ v.time | id,
              jmCoxForm1 = Surv(obsT, delta == 1) ~ W,
              jmCoxForm2 = Surv(obsT, delta == 2) ~ W,
              timeVar = "v.time", data.long = da.long, data.id = da.id,
              coxControl = list(n.iter = 30))
EST <- out$coef
SE  <- out$se.coef
print(round(cbind(EST, SE)[-grep("gammas.bs|D", names(EST)), ], 3)) 

#                        EST    SE
# betas.(Intercept)_1  0.225 0.089
# betas.v.time_1       0.195 0.024
# betas.W_1           -0.089 0.050
# betas.(Intercept)_2 -0.233 0.090
# betas.v.time_2       0.414 0.030
# betas.W_2            0.059 0.049
# log.sigma1          -1.238 0.090
# log.sigma2          -1.142 0.115
# gammas1.W           -0.239 0.076
# alpha1               1.402 0.106
# gammas2.W            0.064 0.077
# alpha2               1.151 0.118
# Lt1_1                0.316 0.055
# Lt2_1                0.517 0.029
# Lt3_1               -0.409 0.053
# Lt4_1                0.477 0.021
# Lt5_1                0.487 0.030
# Lt1_2                0.310 0.075
# Lt2_2                0.495 0.036
# Lt3_2               -0.423 0.061
# Lt4_2                0.254 0.033
# Lt5_2                0.645 0.034
```

`betas.(Intercept)_1`, `betas.v.time_1` and `betas.W_1` are the
estimates for fixed effects $\boldsymbol{\beta_1}$ of the longitudinal
risk score of Event 1, whereas `betas.(Intercept)_2`, `betas.v.time_2`
and `betas.W_2` are the estimates for fixed effects
$\boldsymbol{\beta_2}$ of longitudinal risk score of Event 2.

`log.sigma1` and `log.sigma2` are the estimates of $log(\sigma_1)$ and
$log(\sigma_2)$ separately.

`gammas1.W` and `gammas2.W` are the estimates for
$\boldsymbol{\gamma_1}$ and $\boldsymbol{\gamma_2}$ separately. `alpha1`
and `alpha2` are the estimates for $\boldsymbol{\alpha_1}$ and
$\boldsymbol{\alpha_2}$ separately.

`Lt1_1`, `Lt2_1`, `Lt3_1`, `Lt4_1` and `Lt5_1` are the estimates of the
components in $\boldsymbol{\xi_1}$, whereas `Lt1_2`, `Lt2_2`, `Lt3_2`,
`Lt4_2` and `Lt5_2` are the estimates of the components in
$\boldsymbol{\xi_2}$.

`gammas.bs.11 - gammas.bs.19` (not shown) are estimates for
approximating log baseline hazards $h_{01}(t)$, whereas
`gammas.bs.21 - gammas.bs.29` (not shown) are estimates for
approximating log baseline hazards $h_{02}(t)$.

`D1- D10` (not shown) are estimates for Choleski factors of the random
effects covariance matrix ($\textbf{D}$).

## Reference

Cuihong Zhang, Jing Ning, Jianwen Cai, James E. Squires, Steven H. Belle
& Ruosha Li (2023). Dynamic Risk Score Modeling for Multiple
Longitudinal Risk Factors and Survival \[Manuscript submitted for
publication\]
