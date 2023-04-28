\name{DRS-package}
\alias{DRS-package}
\docType{package}
\title{
  Dynamic Risk Score Modeling for Multiple Longitudinal Risk Factors and Survival
}
\description{
  Dynamic Risk Score Modeling for Multiple Longitudinal Risk Factors and Survival Subject to Both Independent and Dependent Censoring
}

\details{
  \tabular{ll}{
    Package: \tab DRS\cr
    Type: \tab Package\cr
    Version: \tab 0.1.0\cr
    Date: \tab 2023-04-25\cr
    License: \tab GPL (>= 2)\cr
    LazyLoad: \tab yes\cr
  }
  The joint modeling framework associates the longitudinal risk factors to the survival outcome through joining two or more submodels: a longitudinal submodel for each longitudinally measured risk factor and a survival submodel for the event outcome. Subject-specific random effects are introduced to bridge the longitudinal model(s) and the survival model, enabling dynamic risk prediction based on the trajectories of longitudinal risk factors. However, the number of parameters and model complexity grow quickly with the number of longitudinal risk factors. Further, it is desirable to define a longitudinal, summary risk score that can be easily calculated by clinicians, for monitoring changes in disease status and incorporation into patient care.

Motivated by these considerations, we propose a new modeling framework by combining multiple longitudinal risk factors into two longitudinal summary risk scores,
to capture participants’ risk of the main event and the event which censors the main event. Instead of formulating a submodel for each longitudinal risk factor, we pose two submodels for the two risk scores. Thus, the main thrust of the proposed model is the reduced dimension of random effects regardless of the number of longitudinal risk factors. Furthermore, the proposed model is more parsimonious and scalable when the number of longitudinal risk factors increases.
}

\author{
Cuihong Zhang, Ruosha Li <Ruosha.Li@uth.tmc.edu>
}

\keyword{package}
