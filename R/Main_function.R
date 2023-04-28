## Main functions for dynamic risk score modeling (DRS.JM) with two underlying risk scores.
## Iteratively update estimators from two parts.
## coxBetasEst (Part 1), modified estimation function from SurvLong package. The original estimation function is betaEst with score function scoreHalf. Name the normalized estimators as coxbetas.
## JM_modified (Part 2), modified estimation function from JM package to accommodate the setting with two longitudinal risk scores and two events. Name the estimators as jmpars.
## Score2Y, score functions.

DRS.JM <- function(coxForm, jmFixedForm, jmRandomForm, jmCoxForm1, jmCoxForm2, timeVar, data.long, data.id, kType = "epan", coxControl = list(), jmControl = list()){

  ## Control parameters
  coxCtr <- list(n.iter = 4, tol0 = 0.005, tol1 = 0.0015, tol2 = 0.005)
  coxCtr[names(coxControl)] <- coxControl

  ## individual marker matrix
  Ltmat <- model.matrix(coxForm, data.long)[, -1]     ## remove intercept
  WVar1 <- attr(terms(jmCoxForm1), "term.labels")     ## W variable name
  WVar2 <- attr(terms(jmCoxForm2), "term.labels")     ## W variable name
  if(!identical(WVar1, WVar2)){
    stop("The RHS must be the same for the two jmCoxForm")
  }
  WVar <- WVar1
  if(identical(WVar, character(0))) WVar <- NULL    ## set as null if no W
  idVar <- all.vars(jmRandomForm)[length(all.vars(jmRandomForm))]
  TVar <- all.vars(jmCoxForm1)[1]                    ## obsT
  deltaVar <- all.vars(jmCoxForm1)[2]                ## delta
  jmcoxfit1 <- coxph(jmCoxForm1, data = data.id, x = TRUE)
  jmcoxfit2 <- coxph(jmCoxForm2, data = data.id, x = TRUE)

  ## calculate hn
  survDat1 <- jmcoxfit1$y
  effn1 <- sum(survDat1[, 2])
  hn1 <- 2 * IQR(data.long[[timeVar]]) * effn1^(-0.4)
  survDat2 <- jmcoxfit2$y
  effn2 <- sum(survDat2[, 2])
  hn2 <- 2 * IQR(data.long[[timeVar]]) * effn2^(-0.4)

  ## updating coxbetas and jmpars
  ThetasIter <- result.cox <- jm.iter <- NULL # store all parameters for each iteration
  for(iter in 1:coxCtr$n.iter){
    if(iter == 1){
      # init.coxbetas1
      coxfit <- tryCatch(halfKernel(X = cbind(data.id[, idVar], survDat1), Z = cbind(data.long[, c(idVar, timeVar, WVar)], Ltmat),
                                    tau = max(survDat1[,1], survDat2[,1]), kType = kType, bw = hn1, verbose = F), warning = function(w) NA, error = function(e) NA)
      coxbetas1 <- coxfit$betaHat
      coxbetas1 <- norm1fun(coxbetas1[, !colnames(coxbetas1) %in% WVar]) ## remove estimator for baseline covariate, and normalize the estimators for individual risk factors

      # init.coxbetas2
      coxfit <- tryCatch(halfKernel(X = cbind(data.id[, idVar], survDat2), Z = cbind(data.long[, c(idVar, timeVar, WVar)], Ltmat),
                                    tau = max(survDat1[,1], survDat2[,1]), kType = kType, bw = hn2, verbose = F), warning = function(w) NA, error = function(e) NA)
      coxbetas2 <- coxfit$betaHat
      coxbetas2 <- norm1fun(coxbetas2[, !colnames(coxbetas2) %in% WVar]) ## remove estimator for baseline covariate, and normalize the estimators for individual risk factors

      conv1 <- conv2 <- 0
    } else {
      for(iter.cox in 1:3){
        # update coxbetas1
        coxfit1 <- tryCatch(coxBetasEst(delta = 1, X = cbind(data.id[, idVar,drop=F], survDat1[, 1:2], jmcoxfit1$x), Lt = cbind(data.long[, c(idVar, timeVar)], Ltmat),
                                        knowbetas = coxbetas2, tau = max(survDat1[,1], survDat2[,1]), h = hn1, kType = kType, thetas = list.jmpars, init.betas = coxbetas1, WVar = WVar), warning = function(w) w, error = function(e) NA)
        coxbetas1 <- coxfit1$betaHat
        # update coxbetas2
        coxfit2 <- tryCatch(coxBetasEst(delta = 2, X = cbind(data.id[, idVar,drop=F], survDat2[, 1:2], jmcoxfit2$x), Lt = cbind(data.long[, c(idVar, timeVar)], Ltmat),
                                        knowbetas = coxbetas1, tau = max(survDat1[,1], survDat2[,1]), h = hn2, kType = kType, thetas = list.jmpars, init.betas = coxbetas2, WVar = WVar), warning = function(w) w, error = function(e) NA)
        coxbetas2 <- coxfit2$betaHat
        result.cox <- rbind(result.cox, c(iter.cox, coxbetas1, coxbetas2))
      }
    }
    names(coxbetas1) <- paste0(colnames(Ltmat), "_1")
    names(coxbetas2) <- paste0(colnames(Ltmat), "_2")

    ### change EM iterations to be faster for the first few iterations
    if(iter <= 2){
      jmCtr <- list(tol2 = 0.010, only.EM = TRUE, iter.EM = 15)        ## for the first two iters, relax the JM convergence criteria to save time
    } else if (iter > 2 & iter <= 10) {
      jmCtr <- list(tol2 = 0.005)
    } else {
      jmCtr <- list(tol2 = 0.005)
      coxCtr$tol2 <- 0.010                                             ## after 10 times of iters, relax the final convergence criteria
    }

    data.long$Rt <- Ltmat %*% coxbetas1
    data.long$Ut <- Ltmat %*% coxbetas2
    jmlmefit1 <- eval(substitute(lme(form, random = jmRandomForm, data = data.long, control = lmeControl(opt='optim')),
                                 list(form = formula(paste0("Rt~", as.character(jmFixedForm)[2])))))
    jmlmefit2 <- eval(substitute(lme(form, random = jmRandomForm, data = data.long, control = lmeControl(opt='optim')),
                                 list(form = formula(paste0("Ut~", as.character(jmFixedForm)[2])))))
    # lme enviroment wired environment bug https://stat.ethz.ch/pipermail/r-help/2003-January/029199.html

    jmfit <- JM_modified(jmlmefit1, jmlmefit2, jmcoxfit1, jmcoxfit2, timeVar = timeVar, method = "spline-PH-aGH", control = jmCtr)
    list.jmpars <- jmfit$coefficients
    if(is.null(WVar)){ # above line will extract gammas.bs for gammas if W is null
      list.jmpars$gammas1 <- NULL
      list.jmpars$gammas2 <- NULL
    }

    ## h0 & cumulative baseline hazards H0, check the baseline hazard estimators together
    h0bs1 <- function(tt) {
      kn <- jmfit$control$knots[[1]]
      W2s <- splineDesign(kn, tt, ord = jmfit$control$ord, outer.ok = TRUE)
      exp(unlist(W2s) %*% list.jmpars$gammas.bs1)
    }
    h0bs2 <- function(tt) {
      kn <- jmfit$control$knots[[1]]
      W2s <- splineDesign(kn, tt, ord = jmfit$control$ord, outer.ok = TRUE)
      exp(unlist(W2s) %*% list.jmpars$gammas.bs2)
    }
    ch0bs1 <- integrate(h0bs1, 0, max(survDat1[, 1]))$value       ## integral of h0 from 0 to C.max
    ch0bs2 <- integrate(h0bs2, 0, max(survDat2[, 1]))$value       ## integral of h0 from 0 to C.max

    ThetasIter <- rbind(ThetasIter, c(unlist(list.jmpars), coxbetas1, coxbetas2, ch0bs1, ch0bs2))
    gammsbs1idx <- grep("gammas.bs1", colnames(ThetasIter))
    gammsbs2idx <- grep("gammas.bs2", colnames(ThetasIter))
    jm.iter <- rbind(jm.iter, jmfit$iters)

    if(iter > 1){
      conv1 <- max(abs(ThetasIter[iter, -c(gammsbs1idx, gammsbs2idx)] - ThetasIter[iter - 1, -c(gammsbs1idx, gammsbs2idx)])) < coxCtr$tol2
      conv2 <- max(abs(ThetasIter[iter, -c(gammsbs1idx, gammsbs2idx)] - ThetasIter[iter - 1, -c(gammsbs1idx, gammsbs2idx)]) / abs(ThetasIter[iter - 1, -c(gammsbs1idx, gammsbs2idx)] + coxCtr$tol1)) < coxCtr$tol0
    }
    if(conv1 | conv2){
      break
    }
  }
  ThetasIter <- ThetasIter[, -c(ncol(ThetasIter)-1, ncol(ThetasIter))]  ## delete the last column about ch0bs1 and ch0bs2
  if(sum(conv1, conv2) == 0) {
    warning(paste0("Parameter estimates did not converge within ",  coxCtr$n.iter, " iterations."))
  }

  list.Full <- list(coxfit1 = coxfit1, coxfit2 = coxfit2, jmfit = jmfit)
  list.Thetas <- c(list.jmpars, list(coxbetas1 = coxbetas1[-length(coxbetas1)]), list(coxbetas2 = coxbetas2[-length(coxbetas2)]))
  Thetas <- unlist(as.relistable(list.Thetas)) ## must use relistable

  environment(Score2Y) <- environment()
  Score <- tryCatch(Score2Y(Thetas, list.Full), warning = function(w) NA, error = function(e) NA)
  dUdpar <- tryCatch(fd.vec1(Thetas, Score2Y, list.Full = list.Full, eps = 1e-05), warning = function(w) NA, error = function(e) NA)
  if(!any(is.na(dUdpar))){
    varCov <-  tryCatch(solve(dUdpar) %*% Score$VarScore %*% t(solve(dUdpar)), warning = function(w) NA, error = function(e) NA)
  } else {
    varCov <- NA
  }

  # delta method to find the variance of the last individual risk factor
  coxbetas1.n <- coxbetas1[-length(coxbetas1)]  ## extract the first p-1 normalized coxbetas1
  coxbetas2.n <- coxbetas2[-length(coxbetas2)]
  dBetadNorm1 <- rbind(diag(length(coxbetas1.n)), c(-coxbetas1.n/c(sqrt(1-sum(coxbetas1.n^2)))))
  dBetadNorm2 <- rbind(diag(length(coxbetas2.n)), c(-coxbetas2.n/c(sqrt(1-sum(coxbetas2.n^2)))))
  if(!any(is.na(varCov))){
    varCov.coxbetas1.n <- varCov[(ncol(varCov)-2*length(coxbetas1.n)+1):(ncol(varCov)-length(coxbetas1.n)), (ncol(varCov)-2*length(coxbetas1.n)+1):(ncol(varCov)-length(coxbetas1.n))]
    varCov.coxbetas1 <- dBetadNorm1 %*% varCov.coxbetas1.n %*% t(dBetadNorm1)

    varCov.coxbetas2.n <- varCov[(ncol(varCov)-length(coxbetas2.n)+1):ncol(varCov), (ncol(varCov)-length(coxbetas2.n)+1):ncol(varCov)]
    varCov.coxbetas2 <- dBetadNorm2 %*% varCov.coxbetas2.n %*% t(dBetadNorm2)
    se.coef <- tryCatch(sqrt(c(diag(varCov)[1:(ncol(varCov)-2*length(coxbetas1.n))], diag(varCov.coxbetas1), diag(varCov.coxbetas2))), warning = function(w) NA, error = function(e) NA)
  } else {
    se.coef <- NA
  }

  return(list(coefficients = ThetasIter[nrow(ThetasIter), ], varCov = varCov, se.coef = se.coef, convergence = c(conv1, conv2), coefs.iter = cbind(ThetasIter, jm.iter), coxfit1 = coxfit1, coxfit2 = coxfit2, jmfit = jmfit))
}

## coxBetasEst (Part 1), modified estimation function (betaEST) from SurvLong R package to estimate normalized coxbetas
coxBetasEst <- function(Lt, delta, X, tau, h, thetas, init.betas = NULL, WVar, knowbetas,
                        kType = "epan", tol = 0.001, maxiter = 100){
  pre <- preprocessInputs(data.X = X, data.Lt = Lt)
  X <- pre$data.X
  Lt <- pre$data.Lt

  if( is.null(init.betas) ) {
    coxbetas <- numeric(ncol(Lt) - 2L)
    coxbetas[] <- 0.01
  } else {
    coxbetas <- init.betas
  }
  coxbetas <- coxbetas/sqrt(sum(coxbetas^2)) # to satisfy the ||coxbetas|| = 1 constrain

  iter <- 0L
  while( TRUE ) {
    ################## Main modification of this part to calculate new parameter values, start ###############
    Lvec <- halfKernalScore(coxbetas, Lt, delta, X, tau, h, kType, thetas, knowbetas, WVar)
    Ldet <- det(Lvec$dUdBeta)
    if( (Ldet < 1.5e-8 && Ldet > -1.5e-8) ) {
      stop("singular matrix encountered in Newton-Raphson")
    }
    change <- solve(cbind(rbind(Lvec$dUdBeta, 2*coxbetas), c(2*coxbetas, 0))) %*% c(-Lvec$U, 1-sum(coxbetas^2))
    change <- change[-length(change)]
    (coxbetas.hat <- coxbetas + change)
    ################# Main modification of this part to calculate new parameter values, end ##################


    # Determine if parameter estimates have converged.
    test <- TRUE
    for( i in 1L:length(coxbetas) ) {
      if( abs(coxbetas[i]) > 0.001 ) {
        if( abs(change[i])/abs(coxbetas[i]) > tol ) test <- FALSE
      } else {
        if( abs(change[i]) > tol ) test <- FALSE
      }
    }
    if( test ) break
    coxbetas <- coxbetas.hat

    iter <- iter + 1L
    if(iter >= maxiter) {
      warning(paste("Parameter estimates did not converge within ",
                    maxiter, "iterations.", sep=""))
      break
    }
  }
  coxbetas <- as.vector(coxbetas)
  names(coxbetas) <- colnames(Lt)[c(-1, -2)]
  return( list("betaHat" = coxbetas, "X" = X, "delta" = delta, "Lt" = Lt, "tau" = tau, "kType" = kType, "h" = h, "knowbetas" = knowbetas, "thetas" = thetas, "U" = Lvec$U, "dUdBeta" = Lvec$dUdBeta) )
}


## halfKernalScore, modified score function (scoreHalf) from SurvLong R package
halfKernalScore <- function(coxbetas, Lt, delta, X, tau, h, kType, thetas, knowbetas, WVar){

  if(!is.list(thetas)){
    thetas <- relist(thetas)
  }
  betas <- thetas$betas
  sigma <- exp(thetas$log.sigma)
  D <- JM:::chol.transf(thetas$D)
  gammas1 <- thetas$gammas1
  alpha1 <- thetas$alpha1
  gammas2 <- thetas$gammas2
  alpha2 <- thetas$alpha2

  D11idx <- 1:(nrow(D) / 2)
  D11 <- D[ D11idx,  D11idx]
  D12 <- D[ D11idx, -D11idx]
  D22 <- D[-D11idx, -D11idx]

  p <- ncol(Lt) - 2L
  Ltp <- data.matrix(Lt[,3L:{p+2L}, drop=FALSE])
  Ltid <- Lt[,1L, drop=FALSE]
  Lttime <- Lt[,2L]

  rho <- t(apply(cbind(1, Lttime, 1, Lttime), 1, function(x) sig12inv22(x[1:2], x[3:4], D11, D12, D22, sigma[1], sigma[2])))
  if (delta == 1){
    Ltp <- alpha1 * rho[, 1] * Ltp
  } else if (delta == 2){
    Ltp <- alpha2 * rho[, 4] * Ltp
  }

  dUdBeta <- matrix(data = 0.0, nrow = p, ncol = p)
  Mmatrix <- matrix(data = 0.0, nrow = p, ncol = p)
  Uvec <- numeric(length = p)

  outerOnce <- matrix(data = apply(Ltp,1L,function(x){x %*% t(x)}),
                      nrow = nrow(Ltp),
                      ncol = p*p,
                      byrow = TRUE)

  n <- nrow(X)
  U.sub <- matrix(0, n, p)   ## subject-specific Score
  for( i in 1L:n ) {
    if( X[i,3L] < 0.5 ) next
    time <- X[i,2L]
    if( time > tau ) next
    use <- (Lttime <= time) & Ltid == X[i,1L]
    if( !any(use) ) next
    kern <- SurvLong:::local_kernel(t = {time - Lttime}, h = h, kType = kType )
    ptIDs <- X[time <= X[,2L],1L]
    LtptIDs <- (unlist(Ltid) %in% ptIDs)
    LtptIDs <- LtptIDs & (Lttime <= time)
    cova <- Ltp[LtptIDs,,drop = FALSE]

    ################## Main modification of this part to accommodate current model need, start ##############################
    if(!is.null(WVar)){
      if (delta == 1){
        gamma0 <- (alpha1 * rho[, 2] * (Ltp %*% knowbetas))[LtptIDs]  ## given the coxbetas of Rt/Ut
        Wp <- merge(Ltid[LtptIDs, , drop = FALSE], X[, c(-2, -3)], by = names(X)[1])[, -1, drop = FALSE]  ## long format Wi
        gammap <- sweep(alpha1 * ((1 - rho[, 1][LtptIDs]) %*% t(betas[paste0(WVar, "_1")]) - (rho[, 2][LtptIDs]) %*% t(betas[paste0(WVar, "_2")])), 2L, gammas1, "+")
        gammapWp <- rowSums(gammap * Wp)
      } else if (delta == 2){
        gamma0 <- (alpha2 * rho[, 3] * (Ltp %*% knowbetas))[LtptIDs]  ## given the coxbetas of Rt/Ut
        Wp <- merge(Ltid[LtptIDs, , drop = FALSE], X[, c(-2, -3)], by = names(X)[1])[, -1, drop = FALSE]  ## long format Wi
        gammap <- sweep(alpha2 * ((1 - rho[, 4][LtptIDs]) %*% t(betas[paste0(WVar, "_2")]) - (rho[, 3][LtptIDs]) %*% t(betas[paste0(WVar, "_1")])), 2L, gammas2, "+")
        gammapWp <- rowSums(gammap * Wp)
      }
    } else {
      if (delta == 1){
        gamma0 <- (alpha1 * rho[, 2] * (Ltp %*% knowbetas))[LtptIDs]  ## given the coxbetas of Rt/Ut
        gammapWp <- 0
      } else if (delta == 2){
        gamma0 <- (alpha2 * rho[, 3] * (Ltp %*% knowbetas))[LtptIDs]  ## given the coxbetas of Rt/Ut
        gammapWp <- 0
      }
    }
    ################## Main modification of this part to accommodate current model need, end ##############################

    prod <- kern[LtptIDs]*exp(cova %*% coxbetas)*exp(gamma0)*exp(gammapWp)
    s0 <- sum(prod)

    if( (s0 > -1.5e-8) && (s0 < 1.5e-8) ) next
    s1 <- colSums(prod[,1L]*cova)
    s2 <- matrix(data = colSums(outerOnce[LtptIDs,,drop=FALSE]*prod[,1L]),
                 nrow = p,
                 ncol = p)

    LtmLt <- sweep(x = Ltp[use,,drop = FALSE],
                   MARGIN = 2L,
                   STATS = s1/s0,
                   FUN = "-")

    tmp <- colSums( kern[use] * LtmLt )
    Mmatrix <- Mmatrix + tmp %*% t(tmp)
    Uvec <- Uvec + tmp
    U.sub[i, ] <- tmp
    dUdBeta <- dUdBeta + sum(kern[use])*(s1 %*% t(s1) - s2*s0)/(s0*s0)
  }
  return(list("U" = Uvec, "dUdBeta" = dUdBeta, "U.sub" = U.sub, "mMatrix" = Mmatrix))
}



## Score2Y, score function for all jmpars and the first p-1 normalized coxbetas1 and coxbetas2.
Score2Y <- function (thetas, list.Full) {

  if(!is.list(thetas)){
    thetas <- relist(thetas)
  }
  list.thetas <- thetas
  betas <- thetas$betas
  sigma <- exp(thetas$log.sigma)
  gammas1 <- thetas$gammas1
  alpha1 <- thetas$alpha1
  gammas.bs1 <- thetas$gammas.bs1
  gammas2 <- thetas$gammas2
  alpha2 <- thetas$alpha2
  gammas.bs2 <- thetas$gammas.bs2
  D <- JM:::chol.transf(thetas$D)
  coxbetas1 <- thetas$coxbetas1
  coxbetas2 <- thetas$coxbetas2

  ################################################################################################
  ######## Score function for jmpars, modified from JM R package #################################
  ######## To accommodate the current setting with two underlying risk scores and two events #####
  ######## Also save the scores of each jmpars for each subject ##################################
  ################################################################################################

  ############ Extract model information (JM::splinePHGH.fit) ############
  xDat1 <- list.Full$jmfit$xDat1
  xDat2 <- list.Full$jmfit$xDat2
  xDat1$y <- c(Ltmat %*% c(coxbetas1, sqrt(1-sum(coxbetas1^2))))  ## replace y to Ltmat %*% standard(coxbetas1)
  xDat2$y <- c(Ltmat %*% c(coxbetas2, sqrt(1-sum(coxbetas2^2))))  ## replace here since later function will use xDat as input, such as gr.longSplinePH
  control <- list.Full$jmfit$control
  parameterization <- list.Full$jmfit$parameterization

  # response vectors and design matrix
  logT <- xDat1$logT    # all T
  idT <- xDat1$idT
  y <- c(xDat1$y, xDat2$y)
  id <- c(xDat1$id, xDat2$id)
  X <- as.matrix(bdiag(xDat1$X, xDat2$X))
  Xtime <- as.matrix(bdiag(xDat1$Xtime, xDat2$Xtime))
  Xs <- as.matrix(bdiag(xDat1$Xs, xDat2$Xs))
  Z <- as.matrix(bdiag(xDat1$Z, xDat2$Z))
  Ztime <- as.matrix(bdiag(xDat1$Ztime, xDat2$Ztime))
  Zs <- as.matrix(bdiag(xDat1$Zs, xDat2$Zs))
  id.GK <- xDat1$id.GK

  # sample size settings
  ncx <- ncol(X)
  ncz <- ncol(Z)
  N <- length(y) / 2
  n <- length(logT)
  nGK <- n * control$GKk

  # crossproducts and others
  XtX <- crossprod(X) # 2*2 t(X)%*%X
  ZtZ <- t(sapply(split(Z, id), function (x) crossprod(matrix(x, ncol = ncz)))) # 100*4 = 100 list with 2*2

  # Gauss-Hermite quadrature rule components
  GH <- JM:::gauher(control$GHk)
  b <- as.matrix(expand.grid(rep(list(GH$x), ncz))) # 25*2 b_t sum_(t1,..,tq)
  k <- nrow(b)
  wGH <- as.matrix(expand.grid(rep(list(GH$w), ncz)))
  wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b)) # 25
  if (control$typeGH == "simple") {
    b <- sqrt(2) * t(control$inv.chol.VC %*% t(b))  # r_t, why not + bi?
    wGH <- wGH * control$det.inv.chol.VC
  } else {
    b <- sqrt(2) * b
    VCdets <- control$det.inv.chol.VCs
  }
  dimnames(b) <- NULL
  b2 <- t(apply(b, 1, function (x) x %o% x)) # 25*4
  Zb <- Z %*% t(b)    # 779*25 each line of long data, every b point
  if (parameterization %in% c("value", "both")) {
    Ztimeb <- Ztime %*% t(b) # 100*25
    Zsb <- Zs %*% t(b) # 1500*25
  }

  # pseudo-adaptive Gauss-Hermite
  if (control$typeGH != "simple") {
    lis.b <- vector("list", n)
    for (i in 1:n) {
      lis.b[[i]] <- t(control$inv.chol.VCs[[i]] %*% t(b)) +
        rep(control$ranef[i, ], each = k) # lis.b is r_t
      Zb[id == i, ] <- Z[id == i, , drop = FALSE] %*% t(lis.b[[i]])
    }
    lis.b2 <- lapply(lis.b, function (b) if (ncz == 1) b * b else
      t(apply(b, 1, function (x) x %o% x)))
    for (i in seq_along(logT)) {
      if (parameterization %in% c("value", "both")) {
        bb <- t(lis.b[[idT[i]]])
        Ztimeb[c(idT, idT) == i, ] <- Ztime[c(idT, idT) == i, , drop = FALSE] %*% bb
        Zsb[id.GK == i, ] <- Zs[id.GK == i, ] %*% bb
      }
    }
  }

  # calculate separate survival log-likelihood: used in EM and LogLik function
  logSurv <- function(gammas, alpha, gammas.bs, Y, Ys, xDat) {
    # Req Data: logT, d, idT, W, W2, W2s, WintF.vl, Ws.intF.vl, P, wk, (Y, Ys, id.GK)
    # Req par: gammas, alpha, gammas.bs
    list2env(xDat, environment())
    eta.tw1 <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, length(logT)) # 100
    eta.tw2 <- as.vector(W2 %*% gammas.bs) # 100
    eta.ws <- as.vector(W2s %*% gammas.bs) # 1500
    if (parameterization %in% c("value", "both")) {
      eta.t <- eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y # 100*25
      eta.s <- c(Ws.intF.vl %*% alpha) * Ys #1500*25
    }
    log.hazard <- eta.t # 100*25
    log.survival <- - exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + eta.s), id.GK, reorder = FALSE) # 100*25
    dimnames(log.survival) <- NULL
    log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE) # 100*25
    return(listNamed(eta.tw1, eta.ws, log.p.tb))
  }

  environment(Score.splineGH) <- environment()
  Score.JM.sub <- Score.splineGH(thetas, indScore = TRUE)

  ####################################################################
  ######## Score function for the first p-1 normalized coxbetas ######
  ####################################################################
  WVar <- colnames(xDat1$W)

  ULt1.sub <- NormScore( coxbetas = thetas$coxbetas1,
                         Lt = list.Full$coxfit1$Lt,
                         delta = list.Full$coxfit1$delta,
                         X = list.Full$coxfit1$X,
                         tau = list.Full$coxfit1$tau,
                         h = list.Full$coxfit1$h,
                         kType = list.Full$coxfit1$kType,
                         knowbetas = c(thetas$coxbetas2, sqrt(1-sum(thetas$coxbetas2^2))),
                         thetas = thetas, WVar = WVar)$Score.sub

  ULt2.sub <- NormScore(coxbetas = thetas$coxbetas2,
                        Lt = list.Full$coxfit2$Lt,
                        delta = list.Full$coxfit2$delta,
                        X = list.Full$coxfit2$X,
                        tau = list.Full$coxfit2$tau,
                        h = list.Full$coxfit2$h,
                        kType = list.Full$coxfit2$kType,
                        knowbetas = c(thetas$coxbetas1, sqrt(1-sum(thetas$coxbetas1^2))),
                        thetas = thetas, WVar = WVar)$Score.sub

  Score.sub <- cbind(Score.JM.sub, -ULt1.sub, -ULt2.sub)  ## Combine the score vectors from JM and Cox parts
  colnames(Score.sub) <- names(unlist(thetas))
  Score <- colSums(Score.sub, na.rm = TRUE)
  VarScore <- t(Score.sub) %*% Score.sub                 ## variance matrix of combined score function
  colnames(VarScore) <- names(unlist(thetas))
  rownames(VarScore) <- names(unlist(thetas))

  list.score <- list(Score.sub = Score.sub, Score = Score, VarScore = VarScore)
}


## Score function for the first p-1 normalized coxbetas (coxbetas.n)
NormScore <- function(coxbetas.n, Lt, delta, X, tau, h, kType, thetas, knowbetas, WVar){

  if(!is.list(thetas)){
    thetas <- relist(thetas)
  }
  coxbetas <- c(coxbetas.n, sqrt(1-sum(coxbetas.n^2)))  ## express normalized coxbetas as the function of the first p-1 normalized coxbetas
  ################# similar as halfKernalScore, start here ######################
  betas <- thetas$betas
  sigma <- exp(thetas$log.sigma)
  D <- JM:::chol.transf(thetas$D)
  gammas1 <- thetas$gammas1
  alpha1 <- thetas$alpha1
  gammas2 <- thetas$gammas2
  alpha2 <- thetas$alpha2

  D11idx <- 1:(nrow(D) / 2)
  D11 <- D[ D11idx,  D11idx]
  D12 <- D[ D11idx, -D11idx]
  D22 <- D[-D11idx, -D11idx]

  p <- ncol(Lt) - 2L
  Ltp <- data.matrix(Lt[,3L:{p+2L}, drop=FALSE])
  Ltid <- Lt[,1L, drop=FALSE]
  Lttime <- Lt[,2L]

  rho <- t(apply(cbind(1, Lttime, 1, Lttime), 1, function(x) sig12inv22(x[1:2], x[3:4], D11, D12, D22, sigma[1], sigma[2])))
  if (delta == 1){
    Ltp <- alpha1 * rho[, 1] * Ltp
  } else if (delta == 2){
    Ltp <- alpha2 * rho[, 4] * Ltp
  }

  Uvec <- numeric(length = p-1)

  n <- nrow(X)
  U.sub <- matrix(0, n, p-1)   ## subject-specific Score

  for( i in 1L:n ) {
    if( X[i,3L] < 0.5 ) next
    time <- X[i,2L]
    if( time > tau ) next
    use <- (Lttime <= time) & Ltid == X[i,1L]
    if( !any(use) ) next
    kern <- SurvLong:::local_kernel(t = {time - Lttime}, h = h, kType = kType )
    ptIDs <- X[time <= X[,2L],1L]
    LtptIDs <- (unlist(Ltid) %in% ptIDs)
    LtptIDs <- LtptIDs & (Lttime <= time)
    cova <- Ltp[LtptIDs,,drop = FALSE]
    if(!is.null(WVar)){
      if (delta == 1){
        gamma0 <- (alpha1 * rho[, 2] * (Ltp %*% knowbetas))[LtptIDs]
        Wp <- merge(Ltid[LtptIDs, , drop = FALSE], X[, c(-2, -3)], by = names(X)[1])[, -1, drop = FALSE]  ## long format Wi
        gammap <- sweep(alpha1 * ((1 - rho[, 1][LtptIDs]) %*% t(betas[paste0(WVar, "_1")]) - (rho[, 2][LtptIDs]) %*% t(betas[paste0(WVar, "_2")])), 2L, gammas1, "+")
        gammapWp <- rowSums(gammap * Wp)
      } else if (delta == 2){
        gamma0 <- (alpha2 * rho[, 3] * (Ltp %*% knowbetas))[LtptIDs]
        Wp <- merge(Ltid[LtptIDs, , drop = FALSE], X[, c(-2, -3)], by = names(X)[1])[, -1, drop = FALSE]  ## long format Wi
        gammap <- sweep(alpha2 * ((1 - rho[, 4][LtptIDs]) %*% t(betas[paste0(WVar, "_2")]) - (rho[, 3][LtptIDs]) %*% t(betas[paste0(WVar, "_1")])), 2L, gammas2, "+")
        gammapWp <- rowSums(gammap * Wp)
      }
    } else {
      if (delta == 1){
        gamma0 <- (alpha1 * rho[, 2] * (Ltp %*% knowbetas))[LtptIDs]
        gammapWp <- 0
      } else if (delta == 2){
        gamma0 <- (alpha2 * rho[, 3] * (Ltp %*% knowbetas))[LtptIDs]
        gammapWp <- 0
      }
    }

    prod <- kern[LtptIDs]*exp(cova %*% coxbetas)*exp(gamma0)*exp(gammapWp)
    s0 <- sum(prod)

    if( (s0 > -1.5e-8) && (s0 < 1.5e-8) ) next
    s1 <- colSums(prod[,1L]*cova)

    LtmLt <- sweep(x = Ltp[use,,drop = FALSE],
                   MARGIN = 2L,
                   STATS = s1/s0,
                   FUN = "-")

    tmp <- colSums( kern[use] * LtmLt )
    ################# similar to halfKernalScore, end here ######################

    dBetadNorm <- rbind(diag(p-1), c(-coxbetas.n/c(sqrt(1-sum(coxbetas.n^2)))))     ## chain rule to find score for the first p-1 normalized coxbetas (coxbetas.n)
    tmp.n <- t(tmp) %*% dBetadNorm
    Uvec <- Uvec + tmp.n
    U.sub[i, ] <- tmp.n
  }
  return(list(Score = Uvec, Score.sub = U.sub))
}


## Other supportive functions
## preprocessInputs, slight modification from preprocessInputs function in SurvLong R package
preprocessInputs <- function(data.X, data.Lt) {
  # nc <- ncol(data.x)
  # if( nc != 3L ) stop("data.x must include {ID, time, delta}.") ## modification due to different data structure

  ncz <- ncol(data.Lt)
  if( ncz < 3L ) stop("data.Lt must include {ID, time, measurement}.")

  if( !is.integer(data.Lt[,1L]) ) {
    data.Lt[,1L] <- as.integer(round(data.Lt[,1L],0))
    # message("Patient IDs in data.Lt were coerced to integer.\n")
  }

  if( !is.integer(data.X[,1L]) ) {
    data.X[,1L] <- as.integer(round(data.X[,1L],0))
    # message("Patient IDs in data.X were coerced to integer.\n")
  }

  rmRow <- apply(data.Lt, 1, function(x){all(is.na(x))})
  data.Lt <- data.Lt[!rmRow,]

  tst <- is.na(data.Lt)
  data.Lt[tst] <- 0.0

  tst <- is.na(data.X[,2L])
  data.X <- data.X[!tst,]

  if( any(data.Lt[,2L] < {-1.5e-8}) ) {
    stop("Time is negative in data.Lt.")
  }

  if( any(data.X[,2L] < {-1.5e-8}) ) {
    stop("Time is negative in data.X.")
  }

  return(list(data.X = data.X, data.Lt = data.Lt))

}


# normalize the coefficients
norm1fun <- function(x) {x / sqrt(sum(x ^ 2))}


# calculate Sigma12 %*% inv(Sigma22)
sig12inv22 <- function(z1t, z2t, D11, D12, D22, sigma1, sigma2){
  sigmaR12 <- c(t(z1t) %*% D11 %*% z1t, t(z1t) %*% D12 %*% z2t)
  sigmaU12 <- c(t(z1t) %*% D12 %*% z2t, t(z2t) %*% D22 %*% z2t)

  varR <- t(z1t) %*% D11 %*% z1t + sigma1^2
  varU <- t(z2t) %*% D22 %*% z2t + sigma2^2
  RUcov <- t(z1t) %*% D12 %*% z2t
  invSig22 <- solve(matrix(c(varR, RUcov, RUcov, varU), 2, 2))

  return <- c(sigmaR12 %*% invSig22, sigmaU12 %*% invSig22)
}

