## JM_modified, modified functions from JM package to accommodate the setting with two longitudinal risk scores and two events.

JM_modified <- function (lmeObject1, lmeObject2, survObject1, survObject2, timeVar, parameterization = c("value"),
                         method = c("spline-PH-aGH", "spline-PH-GH"), control = list(), ...) {

  ############## simplify scenario and prepare required dataset to accommodate current need, start ##################################
  cl <- match.call()
  method. <- match.arg(method)
  method <- switch(method.,
                   "spline-PH-GH" =,
                   "spline-PH-aGH" = "spline-PH-GH")
  # control values
  ind.noadapt <- method. %in% c("spline-PH-GH")
  con <- list(only.EM = FALSE, iter.EM = if (method == "spline-PH-GH") 120 else 50,
              iter.qN = 350, optimizer = "optim", tol1 = 1e-03, tol2 = 1e-04,
              tol3 = sqrt(.Machine$double.eps), numeriDeriv = "fd", eps.Hes = 1e-06,
              parscale = NULL, step.max = 0.1,
              knots = NULL, ObsTimes.knots = TRUE,
              lng.in.kn = 5, ord = 4,
              equal.strata.knots = TRUE, typeGH = if (ind.noadapt) "simple" else "adaptive",
              GHk = if (nrow(lmeObject1$data) < 2000) 15 else 9,
              GKk = 15, verbose = FALSE)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (con$typeGH != "simple" && !"GHk" %in% namc) {
    con$GHk <- if (nrow(lmeObject1$data) < 2000) 5 else 3
  }
  if (length(noNms <- namc[!namc %in% namC]) > 0)
    warning("unknown names in 'control': ", paste(noNms, collapse = ", "))

  # prepare data
  for(obj in 1:2){
    lmeObject <- switch(obj, lmeObject1, lmeObject2)
    survObject <- switch(obj, survObject1, survObject2)

    if (!inherits(lmeObject, "lme"))
      stop("\n'lmeObject' must inherit from class lme.")
    if (length(lmeObject$group) > 1)
      stop("\nnested random-effects are not allowed in lme().")
    if (!is.null(lmeObject$modelStruct$corStruct))
      warning("correlation structure in 'lmeObject' is ignored.\n")
    if (!is.null(lmeObject$modelStruct$varStruct))
      warning("variance structure in 'lmeObject' is ignored.\n")
    if (!inherits(survObject, "coxph"))
      stop("\n'survObject' must inherit from class coxph.")
    if (!is.matrix(survObject$x))
      stop("\nuse argument 'x = TRUE' in coxph()")
    if (length(timeVar) != 1 || !is.character(timeVar))
      stop("\n'timeVar' must be a character string.")

    # survival process
    formT <- formula(survObject)
    W <- survObject$x
    keepW <- suppressWarnings(!is.na(survObject$coefficients))
    W <- W[, keepW, drop = FALSE]
    surv <- survObject$y
    if (attr(surv, "type") == "right") {
      LongFormat <- FALSE
      Time <- survObject$y[, 1]
      d <- survObject$y[, 2]
    }
    idT <- seq_along(Time)
    idT <- match(idT, unique(idT))
    nT <- length(unique(idT))
    if (LongFormat && is.null(survObject$model$cluster))
      stop("\nuse argument 'model = TRUE' and cluster() in coxph().")
    if (!length(W))
      W <- NULL
    if (sum(d) < 5)
      warning("\nmore than 5 events are required.")
    WintF.vl <- as.matrix(rep(1, length(Time)))

    # longitudinal process
    id <- lmeObject$groups[[1]]
    id <- match(id, unique(id))
    b <- data.matrix(ranef(lmeObject)) # !!!
    dimnames(b) <- NULL
    nY <- nrow(b)
    if (nY != nT)
      stop("sample sizes in the longitudinal and event processes differ; ",
           "maybe you forgot the cluster() argument.\n")
    TermsX <- lmeObject$terms
    data <- lmeObject$data[all.vars(TermsX)]
    data <- data[complete.cases(data), ]
    formYx <- formula(lmeObject)
    mfX <- model.frame(TermsX, data = data)
    X <- model.matrix(formYx, mfX)
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    Z <- model.matrix(formYz, mfZ)
    y.long <- model.response(mfX, "numeric")
    data.id <- data[!duplicated(id), ]
    data.id <- data.id[idT, ]
    if (!timeVar %in% names(data))
      stop("\n'timeVar' does not correspond to one of the columns in the model.frame of 'lmeObject'.")

    # check if there are any longitudinal measurements after the event times
    max.timeY <- tapply(data[[timeVar]], factor(id, unique(id)), max)
    max.timeT <- tapply(Time, factor(idT, unique(idT)), max)
    if (!isTRUE(all(max.timeT >= max.timeY))) {
      idnams <- factor(lmeObject$groups[[1]])
      stop("\nit seems that there are longitudinal measurements taken after the event times for some subjects ",
           "(i.e., check subject(s): ", paste(levels(idnams)[(max.timeT < max.timeY)], collapse = ", "), ").")
    }

    # extra design matrices for the longitudinal part
    data.id[[timeVar]] <- Time
    mfX.id <- model.frame(TermsX, data = data.id)
    mfZ.id <- model.frame(TermsZ, data = data.id)
    Xtime <- model.matrix(formYx, mfX.id)
    Ztime <- model.matrix(formYz, mfZ.id)

    # Gauss-Kronrod for 'method = "spline-PH-GH"'
    wk <- JM:::gaussKronrod(con$GKk)$wk
    sk <- JM:::gaussKronrod(con$GKk)$sk
    if (LongFormat) {
      Time0 <- survObject$y[, 1]
      P <- (Time - Time0) / 2
      P1 <- (Time + Time0) / 2
      st <- outer(P, sk) + P1
    } else {
      P <- as.vector(Time)/2
      st <- outer(P, sk + 1)
    }
    dimnames(st) <- names(P) <- NULL
    id.GK <- rep(seq_along(Time), each = con$GKk)
    data.id2 <- data.id[id.GK, , drop = FALSE]
    data.id2[[timeVar]] <- c(t(st))
    mfX <- model.frame(TermsX, data = data.id2)
    mfZ <- model.frame(TermsZ, data = data.id2)
    Xs <- model.matrix(formYx, mfX)
    Zs <- model.matrix(formYz, mfZ)
    Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]

    # prepare W2 and W2s: baseline spline
    strt <- if (is.null(survObject$strata)) gl(1, length(Time)) else survObject$strata
    nstrt <- length(levels(strt))
    split.Time <- split(Time, strt)
    ind.t <- if (LongFormat) {
      unlist(tapply(idT, idT,
                    function (x) c(rep(FALSE, length(x) - 1), TRUE)))
    } else {
      rep(TRUE, length(Time))
    }
    kn <- if (con$equal.strata.knots) { # spline knots
      pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
      pp <- tail(head(pp, -1), -1)
      kk <- quantile(Time[ind.t], pp, names = FALSE)
      kk <- kk[kk < max(Time)]
      rr <- rep(list(sort(c(rep(range(Time, st), con$ord), kk))), nstrt)
      names(rr) <- names(split.Time)
      rr
    }
    con$knots <- kn
    W2 <- mapply(function (k, t) splineDesign(k, t, ord = con$ord), kn,
                 split.Time, SIMPLIFY = FALSE)
    if (any(sapply(W2, colSums) == 0))
      stop("\nsome of the knots of the B-splines basis are set outside the range",
           "\n   of the observed event times for one of the strata; refit the model",
           "\n   setting the control argument 'equal.strata.knots' to FALSE.")
    W2 <- mapply(function (w2, ind) {
      out <- matrix(0, length(Time), ncol(w2))
      out[strt == ind, ] <- w2
      out
    }, W2, levels(strt), SIMPLIFY = FALSE)
    W2 <- do.call(cbind, W2)
    strt.s <- rep(strt, each = con$GKk)
    split.Time <- split(c(t(st)), strt.s)
    W2s <- mapply(function (k, t) splineDesign(k, t, ord = con$ord),
                  kn, split.Time, SIMPLIFY = FALSE)
    W2s <- mapply(function (w2s, ind) {
      out <- matrix(0, length(Time) * con$GKk, ncol(w2s))
      out[strt.s == ind, ] <- w2s
      out
    }, W2s, levels(strt), SIMPLIFY = FALSE)
    W2s <- do.call(cbind, W2s)
    # response vectors and design matrices: assign as xDat1, xDat2
    assign(paste0("xDat", obj), listNamed(y = y.long, logT = log(Time), d, strata = strt, formYx = formYx, X, formYz = formYz, Z, Xtime, Ztime,
                                          Xs, Zs, W, WintF.vl, Ws.intF.vl, W2, W2s, id, idT, id.GK, P, wk,
                                          ncx = ncol(X), ncz = ncol(Z), XtX = crossprod(X), Xnm = colnames(X), Wnm = colnames(W)))
  }
  ############## simplify scenario and prepare required dataset to accommodate current need, end ##################################


  ############## main modification for initial values, start ###########################################
  # initial values: Longitudinal data, one Big model for two markers
  iid <- c(xDat1$id, xDat2$id)
  yy <- c(xDat1$y, xDat2$y)
  XX <- as.matrix(bdiag(xDat1$X, xDat2$X))
  ZZ <- as.matrix(bdiag(xDat1$Z, xDat2$Z))
  yyg <- rep(1:2, each = length(xDat1$id)) # yy group: marker1 and marker2
  lmeObject <- nlme::lme(yy ~ XX - 1, random = ~ ZZ - 1 | iid, weights = varIdent(form = ~ 1 | yyg),
                   control = list(maxIter = 200, opt = "optim"))
  sigma <- lmeObject$sigma * c(1, coef(lmeObject$modelStruct$varStruct, uncons=FALSE)) # two sigmas
  VC <- getVarCov(lmeObject) # initial for D matrix
  VC <- JM:::dropAttr(VC)
  if (all(VC[upper.tri(VC)] == 0)) VC <- diag(VC)
  if (con$typeGH != "simple") {
    Vs <- vector("list", nY)  # D_i: for each subject
    inv.VC <- solve(VC)
    for (i in 1:nY) {
      Z.i1 <- xDat1$Z[xDat1$id == i, , drop = FALSE]
      Z.i2 <- xDat2$Z[xDat2$id == i, , drop = FALSE]
      ZtZ.i <- as.matrix(bdiag(crossprod(Z.i1) / sigma[1]^2, crossprod(Z.i2) / sigma[2]^2))
      Vs[[i]] <- solve(ZtZ.i + inv.VC)
    }
    con$inv.chol.VCs <- lapply(Vs, function (x) solve(chol(solve(x))))
    con$det.inv.chol.VCs <- sapply(con$inv.chol.VCs, det)
  }
  con$inv.chol.VC <- solve(chol(solve(VC)))  # inv.B: H = solve(D)
  con$det.inv.chol.VC <- det(con$inv.chol.VC)
  b <- data.matrix(ranef(lmeObject))
  dimnames(b) <- NULL
  con$ranef <- b

  # initial values: survival data, separate models
  betas <- fixef(lmeObject)
  long <- c(xDat1$X %*% betas[1:xDat1$ncx]) + rowSums(xDat1$Z * b[xDat1$id, 1:xDat1$ncz]) # used for initial alpha
  init.surv1 <- JM:::initial.surv(exp(xDat1$logT), xDat1$d, xDat1$W, xDat1$WintF.vl, xDat1$WintF.sl, xDat1$id,
                                  times = lmeObject1$data[[timeVar]], method, parameterization, long = long,
                                  extra = list(W2 = xDat1$W2, control = con, ii = xDat1$idT, strata = survObject1$strata),
                                  LongFormat = LongFormat)
  names(init.surv1) <- paste0(names(init.surv1), "1")
  long <- c(xDat2$X %*% betas[-(1:xDat1$ncx)]) + rowSums(xDat2$Z * b[xDat2$id, -(1:xDat1$ncz)]) # used for initial alpha
  init.surv2 <- JM:::initial.surv(exp(xDat2$logT), xDat2$d, xDat2$W, xDat2$WintF.vl, xDat2$WintF.sl, xDat2$id,
                                  times = lmeObject2$data[[timeVar]], method, parameterization, long = long,
                                  extra = list(W2 = xDat2$W2, control = con, ii = xDat2$idT, strata = survObject2$strata),
                                  LongFormat = LongFormat)
  names(init.surv2) <- paste0(names(init.surv2), "2")

  initial.values <- c(list(betas = betas, sigma = sigma, D = VC), init.surv1, init.surv2)
  initial.values

  # remove objects
  rm(list=setdiff(ls(), c("xDat1", "xDat2", "initial.values", "parameterization", "timeVar", "con", lsf.str())))
  # con$iter.EM <- 10
  ############## main modification for initial values, end #############################################




  # joint model fit, modified splinePHGH.fit
  out <- splinePHGH.fit(xDat1, xDat2, initial.values, parameterization, con)

  # check for problems with the Hessian at convergence
  H <- out$Hessian
  if (any(is.na(H) | !is.finite(H))) {
    warning("infinite or missing values in Hessian at convergence.\n")
  } else {
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (!all(ev >= -1e-06 * abs(ev[1])))
      warning("Hessian matrix at convergence is not positive definite.\n")
  }
  # manage output
  VarCov <- JM:::vcov.jointModel(out)
  out$Score <- NULL
  out$Hessian <- NULL
  out$coefficients <- out$coefficients[!sapply(out$coefficients, is.null)]
  out$VarCov <- VarCov
  out$xDat1 <- xDat1
  out$xDat2 <- xDat2
  out$control <- con
  out$parameterization <- parameterization
  out$timeVar <- timeVar
  # out$call <- cl
  return(out)
}




###### modified splinePHGH.fit, modification including data structure, initial values, EM iterations, etc ##############
splinePHGH.fit <- function (xDat1, xDat2, initial.values, parameterization, control) {
  # control = con

  # response vectors and design matrix
  logT <- xDat1$logT    # all T
  idT <- xDat1$idT
  y <- c(xDat1$y, xDat2$y) # all 2Y
  id <- c(xDat1$id, xDat2$id)
  X <- as.matrix(bdiag(xDat1$X, xDat2$X))
  Xtime <- as.matrix(bdiag(xDat1$Xtime, xDat2$Xtime))
  Xs <- as.matrix(bdiag(xDat1$Xs, xDat2$Xs))
  Z <- as.matrix(bdiag(xDat1$Z, xDat2$Z))
  Ztime <- as.matrix(bdiag(xDat1$Ztime, xDat2$Ztime))
  Zs <- as.matrix(bdiag(xDat1$Zs, xDat2$Zs))
  id.GK <- xDat1$id.GK

  # X <- JM:::dropAttr(X); Z <- JM:::dropAttr(Z); WW <- JM:::dropAttr(WW)
  # WintF.vl <- JM:::dropAttr(WintF.vl); WintF.sl <- JM:::dropAttr(WintF.sl)
  # Ws.intF.vl <- JM:::dropAttr(Ws.intF.vl); Ws.intF.sl <- JM:::dropAttr(Ws.intF.sl)
  # if (parameterization == "value") {
  #     Xtime <- JM:::dropAttr(Xtime); Ztime <- JM:::dropAttr(Ztime); Xs <- JM:::dropAttr(Xs); Zs <- JM:::dropAttr(Zs)
  # }

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

  # initial values
  list2env(initial.values, environment()) # betas, sigma, D, (gammas, gammas.bs, alpha) * 2
  gammas1 <- if (!is.null(xDat1$W)) gammas1 else NULL
  gammas2 <- if (!is.null(xDat2$W)) gammas2 else NULL

  # fix environments to use the variables calculated in this splinePHGH.fit function
  environment(opt.survSplinePH) <- environment(gr.survSplinePH) <- environment()
  environment(gr.longSplinePH) <- environment(H.longSplinePH) <- environment()
  environment(LogLik.splineGH) <- environment(Score.splineGH) <- environment()
  # old <- options(warn = (-1))
  # on.exit(options(old))

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

  # EM iterations
  iter <- control$iter.EM
  Y.mat <- matrix(0, iter + 1, ncx + 2) # betas, sigma
  B.mat <- matrix(0, iter + 1, ncz * ncz) # D
  T.mat1 <- matrix(0, iter + 1, length(c(gammas1, alpha1, gammas.bs1))) # gammas.bs, gammas, alpha
  T.mat2 <- matrix(0, iter + 1, length(c(gammas2, alpha2, gammas.bs2)))
  lgLik <- rep(NA, iter + 1)
  conv <- TRUE
  for (it in 1:iter) {
    # save parameter values in matrix
    Y.mat[it, ] <- c(betas, sigma)
    B.mat[it, ] <- D
    T.mat1[it, ] <- c(gammas1, alpha1, gammas.bs1)
    T.mat2[it, ] <- c(gammas2, alpha2, gammas.bs2)

    # E-step: longitudinal log-likelihood
    eta.yx <- as.vector(X %*% betas) # 779
    if (parameterization %in% c("value", "both")) {
      Y <- as.vector(Xtime %*% betas) + Ztimeb # 100*25
      Ys <- as.vector(Xs %*% betas) + Zsb # 1500*25
    }
    mu.y <- eta.yx + Zb # 779*25
    logNorm <- dnorm(y, mu.y, rep(sigma, each = N), TRUE) # 779*25
    log.p.yb <- rowsum(logNorm, id, reorder = FALSE) # 100*25
    dimnames(log.p.yb) <- NULL

    log.p.b <- if (control$typeGH == "simple") {
      rep(JM:::dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
    } else {
      matrix(JM:::dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE), n, k, byrow = TRUE)
    } # 100*25

    # E-step: survival log-likelihood
    logLikSurv1 <- logSurv(gammas1, alpha1, gammas.bs1, Y[1:n, ], Ys[1:nGK, ], xDat1)
    logLikSurv2 <- logSurv(gammas2, alpha2, gammas.bs2, Y[-(1:n), ], Ys[-(1:nGK), ], xDat2)

    # compute total log-likelihood
    p.ytb <- exp(log.p.yb + logLikSurv1$log.p.tb + logLikSurv2$log.p.tb + log.p.b) # 100*25
    if (control$typeGH != "simple")
      p.ytb <- p.ytb * VCdets
    p.yt <- c(p.ytb %*% wGH) # 100 P(y, T), not related to b, so one vector
    log.p.yt <- log(p.yt) # 100
    lgLik[it] <- sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)

    # posterior distribution of b
    p.byt <- p.ytb / p.yt    # 100*25 P(b | y, T)
    # E(b|y, T)
    post.b <- if (control$typeGH == "simple") {
      p.byt %*% (b * wGH)
    } else {
      sapply(seq_len(ncz), function (i)
        (p.byt * t(sapply(lis.b, "[", seq_len(k), i))) %*% wGH)
    } # 100*2 Eb
    # Var(b|y, T)
    post.vb <- if (control$typeGH == "simple") {
      if (ncz == 1) {
        c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
      } else {
        (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
      }
    } else {
      dd <- sapply(seq_len(ncz^2), function (i)
        (p.byt * t(sapply(lis.b2, "[", seq_len(k), i))) %*% wGH)
      bb <- apply(post.b, 1, function (x) x %o% x)
      dd - if (ncz == 1) c(bb) else t(bb)
    } # 100*4

    # check convergence
    if (it > 5 && lgLik[it] > lgLik[it - 1]) {
      thets1 <- c(Y.mat[it - 1, ], T.mat1[it - 1, ], T.mat2[it - 1, ], B.mat[it - 1, ])
      thets2 <- c(Y.mat[it, ], T.mat1[it, ], T.mat2[it, ], B.mat[it, ])
      check1 <- max(abs(thets2 - thets1) / (abs(thets1) + control$tol1)) < control$tol2
      check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * (abs(lgLik[it - 1]) + control$tol3)
      if (check1 || check2) {
        conv <- FALSE
        if (control$verbose)
          cat("\nEM converged!\n")
        break
      }
    }
    if (iter == 0) break

    # M-step
    ZEb <- rowSums(Z * post.b[id, ], na.rm = TRUE)
    mu <- y - eta.yx
    # update D
    Dn <- if (control$typeGH == "simple") {
      matrix(colMeans(p.byt %*% (b2 * wGH), na.rm = TRUE), ncz, ncz)
    } else {
      matrix(colMeans(dd, na.rm = TRUE), ncz, ncz)
    }
    Dn <- 0.5 * (Dn + t(Dn))
    # update sigma
    tr.idx <- matrix(1:ncz^2, ncz, ncz) # to extract separate block of Z'ZvarB
    tr.idx1 <- c(tr.idx[1:xDat1$ncz, 1:xDat1$ncz])
    tr.tZZvarb <- sum((ZtZ * post.vb)[, tr.idx1], na.rm = TRUE)
    sigman1 <- sqrt(c(crossprod(mu[1:N], mu[1:N] - 2 * ZEb[1:N]) + crossprod(ZEb[1:N]) + tr.tZZvarb) / N)
    tr.idx2 <- c(tr.idx[-(1:xDat1$ncz), -(1:xDat1$ncz)])
    tr.tZZvarb <- sum((ZtZ * post.vb)[, tr.idx2], na.rm = TRUE)
    sigman2 <- sqrt(c(crossprod(mu[-(1:N)], mu[-(1:N)] - 2 * ZEb[-(1:N)]) + crossprod(ZEb[-(1:N)]) + tr.tZZvarb) / N)
    sigman <- c(sigman1, sigman2)
    # update betas
    longExtra1 <- c(listNamed(ZEb = ZEb[1:N], Zsb = Zsb[1:nGK, ], alpha = alpha1, sigma = sigma[1]), logLikSurv1[c("eta.tw1", "eta.ws")])
    Hbetas1 <- JM:::nearPD(H.longSplinePH(betas[1:xDat1$ncx], xDat1, longExtra1))
    scbetas1 <- gr.longSplinePH(betas[1:xDat1$ncx], xDat1, longExtra1)
    longExtra2 <- c(listNamed(ZEb = ZEb[-(1:N)], Zsb = Zsb[-(1:nGK), ], alpha = alpha2, sigma = sigma[2]), logLikSurv2[c("eta.tw1", "eta.ws")])
    Hbetas2 <- JM:::nearPD(H.longSplinePH(betas[-(1:xDat1$ncx)], xDat2, longExtra2))
    scbetas2 <- gr.longSplinePH(betas[-(1:xDat1$ncx)], xDat2, longExtra2)
    betasn <- betas - c(solve(Hbetas1, scbetas1), solve(Hbetas2, scbetas2))
    # update surv parameters: gammas, gammas.bs, alpha
    optSurv <- function(gammas, alpha, gammas.bs, xDat, survExtra){
      list.survpars <- listNamed(gammas, alpha, gammas.bs)
      list.survpars <- list.survpars[!sapply(list.survpars, is.null)]
      survpars <- unlist(as.relistable(list.survpars))
      optim(par = survpars, fn = opt.survSplinePH, gr = gr.survSplinePH, xDat = xDat,
            list.survpars = list.survpars, survExtra = survExtra, method = "BFGS",
            control = list(maxit = ifelse(it < 5, 20, 4), parscale = rep(ifelse(it < 5, 0.01, 0.1), length(survpars)))
      )
    }
    survExtra1 <- list(Y = Y[1:n, ], Ys = Ys[1:nGK, ])
    optz.surv1 <- optSurv(gammas1, alpha1, gammas.bs1, xDat1, survExtra1)
    survExtra2 <- list(Y = Y[-(1:n), ], Ys = Ys[-(1:nGK), ])
    optz.surv2 <- optSurv(gammas2, alpha2, gammas.bs2, xDat2, survExtra2)

    # update parameter values
    betas <- betasn
    sigma <- sigman
    D <- Dn
    thetasn1 <- relist(optz.surv1$par)
    gammas1 <- if(!is.null(xDat1$W)) thetasn1$gammas else NULL
    alpha1 <- thetasn1$alpha
    gammas.bs1 <- thetasn1$gammas.bs
    thetasn2 <- relist(optz.surv2$par)
    gammas2 <- if(!is.null(xDat2$W)) thetasn2$gammas else NULL
    alpha2 <- thetasn2$alpha
    gammas.bs2 <- thetasn2$gammas.bs
  }

  list.thetas <- listNamed(betas, log.sigma = log(sigma), D = JM:::chol.transf(D),
                           gammas1, alpha1, gammas.bs1, gammas2, alpha2, gammas.bs2)
  list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
  thetas <- unlist(as.relistable(list.thetas))

  # if not converged, start quasi-Newton iterations
  if (conv && !control$only.EM) {
    if (is.null(control$parscale))
      control$parscale <- rep(0.01, length(thetas))
    if (control$verbose)
      cat("\n\nquasi-Newton iterations start.\n\n")
    out <- if (control$optimizer == "optim") {
      optim(thetas, LogLik.splineGH, Score.splineGH, method = "BFGS",
            control = list(maxit = control$iter.qN, parscale = control$parscale,
                           trace = 10 * control$verbose))
    } else {
      nlminb(thetas, LogLik.splineGH, Score.splineGH, scale = control$parscale,
             control = list(iter.max = control$iter.qN, trace = 1 * control$verbose))
    }
    if ((conv <- out$convergence) == 0 || - out[[2]] > lgLik[iter]) {
      lgLik[iter + 1] <- - out[[2]]
      list.thetas <- relist(out$par, skeleton = list.thetas)
      list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
      thetas <- unlist(as.relistable(list.thetas))
      it <- it + if(control$optimizer == "optim") out$counts[1] else out$iterations
    }
  }
  # calculate Score vector
  Score <- Score.splineGH(thetas, indScore = FALSE)

  # calculate Hessian matrix
  Hessian <- if (control$numeriDeriv == "fd") {
    JM:::fd.vec(unlist(thetas), Score.splineGH, eps = control$eps.Hes)
  } else {
    JM:::cd.vec(unlist(thetas), Score.splineGH, eps = control$eps.Hes)
  }
  # coeficients names
  names(list.thetas$betas) <- c(paste0(xDat1$Xnm, "_1"), paste0(xDat2$Xnm, "_2"))
  names(list.thetas$gammas1) <- xDat1$Wnm
  names(list.thetas$gammas2) <- xDat2$Wnm
  names(list.thetas$alpha1) <- NULL
  names(list.thetas$alpha2) <- NULL
  return(list(coefficients = list.thetas, logLik = lgLik[!is.na(lgLik)], iters = it,
              convergence = conv, Score = Score, Hessian = Hessian, post.b = post.b))
}

## Modified gr.longSplinePH function from JM R package
gr.longSplinePH <- function (betas, xDat, longExtra, indScore = FALSE) {
  list2env(xDat, environment())
  list2env(longExtra, environment())
  eta.yx <- as.vector(X %*% betas)
  Ys <- as.vector(Xs %*% betas) + Zsb
  WintF.vl.alph <- c(WintF.vl %*% alpha)
  Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
  eta.s <- Ws.intF.vl.alph * Ys
  exp.eta.tw.P <- exp(eta.tw1) * P
  Int <- wk * exp(eta.ws + eta.s)
  if (!indScore){
    sc1 <- - crossprod(X, y - eta.yx - ZEb) / sigma^2
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
      ki <- exp.eta.tw.P * rowsum(Int * Ws.intF.vl.alph * Xs[, i], id.GK, reorder = FALSE)
      ki <- c(rowsum(ki, idT, reorder = FALSE))
      kii <- c((p.byt * ki) %*% wGH)
      ddd <- tapply(d * WintF.vl.alph * Xtime[, i], idT, sum)
      sc2[i] <- - sum(ddd - kii, na.rm = TRUE)
    }
    return(c(sc1 + sc2))
  } else {
    sc1.data <- split.data.frame(cbind(y - eta.yx - ZEb, X), id)
    sc1.sub <- lapply(sc1.data, function(x) t(x[, -1, drop = FALSE]) %*% x[, 1] / (-sigma^2))   ## sc1 for each subject
    sc1.sub <- matrix(unlist(sc1.sub), ncol = ncx, byrow = TRUE)

    sc2.sub <- matrix(0, n, ncx)  ## sc2 for each subject
    for (i in 1:ncx) {
      ki <- exp.eta.tw.P * rowsum(Int * Ws.intF.vl.alph * Xs[, i], id.GK, reorder = FALSE)
      ki <- c(rowsum(ki, idT, reorder = FALSE))
      kii <- c((p.byt * ki) %*% wGH)
      ddd <- tapply(d * WintF.vl.alph * Xtime[, i], idT, sum)
      sc2.sub[, i] <- - (ddd - kii)
    }
    return(sc1.sub + sc2.sub)
  }
}


## Modified gr.survSplinePH from JM R package
gr.survSplinePH <- function (survpars, list.survpars, xDat, survExtra, indScore = FALSE) {
  list2env(xDat, environment())
  list2env(survExtra, environment())
  survpars <- relist(survpars, skeleton = list.survpars)
  gammas <- survpars$gammas
  alpha <- survpars$alpha
  gammas.bs <- survpars$gammas.bs
  nk <- length(gammas.bs)
  eta.tw1 <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
  eta.tw2 <- as.vector(W2 %*% gammas.bs)
  eta.t <- eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y
  eta.s <- c(Ws.intF.vl %*% alpha) * Ys
  eta.ws <- as.vector(W2s %*% gammas.bs)
  exp.eta.tw.P <- exp(eta.tw1) * P
  Int <- wk * exp(eta.ws + eta.s)

  if (!indScore){
    scgammas1 <- if (!is.null(W)) {
      ki <- exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
      scg1 <- numeric(ncol(W))
      for (jj in seq_along(scg1)) {
        tt <- rowsum(W[, jj] * ki, idT, reorder = FALSE)
        scg1[jj] <- sum(c((p.byt * tt) %*% wGH), na.rm = TRUE)
      }
      - colSums(W * d, na.rm = TRUE) + scg1
    } else {
      NULL
    }

    scgammas2 <- numeric(nk)
    for (i in 1:nk) {
      kk <- exp.eta.tw.P * rowsum(Int * W2s[, i], id.GK, reorder = FALSE)
      kk <- rowsum(kk, idT, reorder = FALSE)
      scgammas2[i] <- - sum(W2[, i] * d) + sum(c((p.byt * kk) %*% wGH))
    }

    scalpha <- if (parameterization %in% c("value", "both")) {
      rr <- numeric(ncol(WintF.vl))
      for (k in seq_along(rr)) {
        rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.vl[, k] * Ys,
                                     id.GK, reorder = FALSE)
        rrr <- rowsum(rrr, idT, reorder = FALSE)
        rr[k] <- - sum((p.byt * (rowsum(d * WintF.vl[, k] * Y, idT,
                                        reorder = FALSE) - rrr)) %*% wGH, na.rm = TRUE)
      }
      rr
    } else {
      NULL
    }
    return(c(scgammas1, scalpha, scgammas2))
  } else {
    scgammas1.sub <- if (!is.null(W)) {
      ki <- exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
      scg1 <- matrix(0, n, ncol(W))
      for (jj in 1:ncol(W)) {
        tt <- rowsum(W[, jj] * ki, idT, reorder = FALSE)
        scg1[, jj] <- c((p.byt * tt) %*% wGH)
      }
      - (W * d) + scg1
    } else {
      NULL
    }

    scgammas2.sub <- matrix(0, n, nk)
    for (i in 1:nk) {
      kk <- exp.eta.tw.P * rowsum(Int * W2s[, i], id.GK, reorder = FALSE)
      kk <- rowsum(kk, idT, reorder = FALSE)
      scgammas2.sub[,i] <- - W2[, i] * d + c((p.byt * kk) %*% wGH)
    }

    scalpha.sub <- if (parameterization %in% c("value", "both")) {
      rr <- matrix(0, n, ncol(WintF.vl))
      for (l in 1:ncol(rr)) {
        rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.vl[, l] *
                                       Ys, id.GK, reorder = FALSE)
        rrr <- rowsum(rrr, idT, reorder = FALSE)
        rr[,l] <- - (p.byt * (rowsum(d * WintF.vl[, l] * Y, idT,
                                     reorder = FALSE) - rrr)) %*% wGH
      }
      rr
    } else {
      NULL
    }
    return(cbind(scgammas1.sub, scalpha.sub, scgammas2.sub))
  }
}


## Modified H.longSplinePH function from JM R package
H.longSplinePH <- function (betas, xDat, longExtra) {
  list2env(xDat, environment())
  list2env(longExtra, environment())
  eta.yx <- as.vector(X %*% betas)
  Ys <- as.vector(Xs %*% betas) + Zsb
  Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
  eta.s <- Ws.intF.vl.alph * Ys

  exp.eta.tw.P <- exp(eta.tw1) * P
  H1 <- XtX / sigma^2
  Int <- wk * exp(eta.ws + eta.s) #* alpha^2
  H2 <- matrix(0, ncx, ncx)
  for (i in 1:ncx) {
    for (j in i:ncx) {
      XX <- Ws.intF.vl.alph^2 * Xs[, i] * Xs[, j]
      ki <- exp.eta.tw.P * rowsum(Int * XX, id.GK, reorder = FALSE)
      ki <- rowsum(ki, idT, reorder = FALSE)
      kii <- c((p.byt * ki) %*% wGH)
      H2[i, j] <- sum(kii, na.rm = TRUE)
    }
  }
  H2[lower.tri(H2)] <- t(H2)[lower.tri(H2)]
  H1 + H2
}

## Modified LogLik.splineGH function from JM R package
LogLik.splineGH <- function (thetas) {
  thetas <- relist(thetas, skeleton = list.thetas)

  list2env(thetas, environment())
  sigma <- exp(log.sigma)
  D <- if (length(D) == 1) exp(D) else JM:::chol.transf(D)

  # E-step: longitudinal log-likelihood
  eta.yx <- as.vector(X %*% betas) # 779
  if (parameterization %in% c("value", "both")) {
    Y <- as.vector(Xtime %*% betas) + Ztimeb # 100*25
    Ys <- as.vector(Xs %*% betas) + Zsb # 1500*25
  }
  mu.y <- eta.yx + Zb # 779*25
  logNorm <- dnorm(y, mu.y, rep(sigma, each = N), TRUE) # 779*25
  log.p.yb <- rowsum(logNorm, id, reorder = FALSE) # 100*25

  log.p.b <- if (control$typeGH == "simple") {
    rep(JM:::dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
  } else {
    matrix(JM:::dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE), n, k, byrow = TRUE)
  } # 100*25

  # E-step: survival log-likelihood
  logLikSurv1 <- logSurv(gammas1, alpha1, gammas.bs1, Y[1:n, ], Ys[1:nGK, ], xDat1)
  logLikSurv2 <- logSurv(gammas2, alpha2, gammas.bs2, Y[-(1:n), ], Ys[-(1:nGK), ], xDat2)
  # compute total log-likelihood
  p.ytb <- exp(log.p.yb + logLikSurv1$log.p.tb + logLikSurv2$log.p.tb + log.p.b) # 100*25
  if (control$typeGH != "simple")
    p.ytb <- p.ytb * VCdets
  p.yt <- c(p.ytb %*% wGH) # 100 P(y, T), not related to b, so one vector
  log.p.yt <- log(p.yt) # 100
  negLogLik <- - sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
  return(negLogLik)
}


## Modified opt.survSplinePH function from JM R package
opt.survSplinePH <- function (survpars, list.survpars, xDat, survExtra) {
  list2env(xDat, environment())
  list2env(survExtra, environment())
  survpars <- relist(survpars, skeleton = list.survpars)
  gammas <- survpars$gammas
  alpha <- survpars$alpha
  gammas.bs <- survpars$gammas.bs
  eta.tw1 <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
  eta.tw2 <- as.vector(W2 %*% gammas.bs)
  eta.t <- eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y
  eta.s <- c(Ws.intF.vl %*% alpha) * Ys
  eta.ws <- as.vector(W2s %*% gammas.bs)
  log.hazard <- eta.t
  log.survival <- - exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + eta.s),
                                              id.GK, reorder = FALSE)
  dimnames(log.survival) <- NULL
  log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE)
  p.bytn <- p.byt * log.p.tb
  -sum(p.bytn %*% wGH, na.rm = TRUE)
}

## Modified Score.splineGH from JM R package
Score.splineGH <- function (thetas, indScore = FALSE) {
  environment(gr.survSplinePH) <- environment(gr.longSplinePH) <- environment()
  if(!is.list(thetas)){
    thetas <- relist(thetas, skeleton = list.thetas)
  }
  list2env(thetas, environment())
  sigma <- exp(log.sigma)
  D <- if (length(D) == 1) exp(D) else JM:::chol.transf(D)

  # log-likelihood: same as EM
  eta.yx <- as.vector(X %*% betas)
  if (parameterization %in% c("value", "both")) {
    Y <- as.vector(Xtime %*% betas) + Ztimeb
    Ys <- as.vector(Xs %*% betas) + Zsb
  }
  mu.y <- eta.yx + Zb
  logNorm <- dnorm(y, mu.y, rep(sigma, each = N), TRUE)
  log.p.yb <- rowsum(logNorm, id, reorder = FALSE)

  log.p.b <- if (control$typeGH == "simple") {
    rep(JM:::dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
  } else {
    matrix(JM:::dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE),
           n, k, byrow = TRUE)
  }

  logLikSurv1 <- logSurv(gammas1, alpha1, gammas.bs1, Y[1:n, ], Ys[1:nGK, ], xDat1)
  logLikSurv2 <- logSurv(gammas2, alpha2, gammas.bs2, Y[-(1:n), ], Ys[-(1:nGK), ], xDat2)

  p.ytb <- exp(log.p.yb + logLikSurv1$log.p.tb + logLikSurv2$log.p.tb + log.p.b)
  if (control$typeGH != "simple")
    p.ytb <- p.ytb * VCdets
  dimnames(p.ytb) <- NULL
  p.yt <- c(p.ytb %*% wGH)
  p.byt <- p.ytb / p.yt
  post.b <- if (control$typeGH == "simple") {
    p.byt %*% (b * wGH)
  } else {
    sapply(seq_len(ncz), function (i)
      (p.byt * t(sapply(lis.b, "[", seq_len(k), i))) %*% wGH)
  }
  post.vb <- if (control$typeGH == "simple") {
    if (ncz == 1) {
      c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
    } else {
      (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
    }
  } else {
    dd <- sapply(seq_len(ncz^2), function (i)
      (p.byt * t(sapply(lis.b2, "[", seq_len(k), i))) %*% wGH)
    bb <- apply(post.b, 1, function (x) x %o% x)
    dd - if (ncz == 1) c(bb) else t(bb)
  }

  ZEb <- rowSums(Z * post.b[id, ], na.rm = TRUE)
  mu <- y - eta.yx

  if (!indScore){
    # Scores of log(sigma)
    tr.idx <- matrix(1:ncz^2, ncz, ncz) # to extract separate block of Z'ZvarB
    tr.idx1 <- c(tr.idx[1:xDat1$ncz, 1:xDat1$ncz])
    tr.tZZvarb <- sum((ZtZ * post.vb)[, tr.idx1], na.rm = TRUE)
    sclogsigma1 <- - sigma[1] * (- N / sigma[1] + c(crossprod(mu[1:N], mu[1:N] - 2 * ZEb[1:N]) + crossprod(ZEb[1:N]) + tr.tZZvarb) / sigma[1]^3)
    tr.idx2 <- c(tr.idx[-(1:xDat1$ncz), -(1:xDat1$ncz)])
    tr.tZZvarb <- sum((ZtZ * post.vb)[, tr.idx2], na.rm = TRUE)
    sclogsigma2 <- - sigma[2] * (- N / sigma[2] + c(crossprod(mu[-(1:N)], mu[-(1:N)] - 2 * ZEb[-(1:N)]) + crossprod(ZEb[-(1:N)]) + tr.tZZvarb) / sigma[2]^3)

    # Scores of betas
    longExtra1 <- c(listNamed(ZEb = ZEb[1:N], Zsb = Zsb[1:nGK, ], alpha = alpha1, sigma = sigma[1]), logLikSurv1[c("eta.tw1", "eta.ws")])
    scbetas1 <- gr.longSplinePH(betas[1:xDat1$ncx], xDat1, longExtra1)
    longExtra2 <- c(listNamed(ZEb = ZEb[-(1:N)], Zsb = Zsb[-(1:nGK), ], alpha = alpha2, sigma = sigma[2]), logLikSurv2[c("eta.tw1", "eta.ws")])
    scbetas2 <- gr.longSplinePH(betas[-(1:xDat1$ncx)], xDat2, longExtra2)

    score.y <- c(scbetas1, scbetas2, sclogsigma1, sclogsigma2)

    # Scores of survival parameters
    survExtra1 <- list(Y = Y[1:n, ], Ys = Ys[1:nGK, ])
    list.survpars <- list(gammas = gammas1, alpha = alpha1, gammas.bs = gammas.bs1)
    list.survpars <- list.survpars[!sapply(list.survpars, is.null)]
    survpars <- unlist(as.relistable(list.survpars))
    score.t1 <- gr.survSplinePH(survpars, list.survpars, xDat1, survExtra1)

    survExtra2 <- list(Y = Y[-(1:n), ], Ys = Ys[-(1:nGK), ])
    list.survpars <- list(gammas = gammas2, alpha = alpha2, gammas.bs = gammas.bs2)
    list.survpars <- list.survpars[!sapply(list.survpars, is.null)]
    survpars <- unlist(as.relistable(list.survpars))
    score.t2 <- gr.survSplinePH(survpars, list.survpars, xDat2, survExtra2)

    # Scores of D
    score.b <- if (!is.matrix(D)) {
      svD <- 1 / D
      svD2 <- svD^2
      cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
      dim(cS.postVB) <- c(ncz, ncz)
      D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - colSums(as.matrix(post.b^2), na.rm = TRUE) * svD2)
    } else {
      svD <- solve(D)
      dD <- JM:::deriv.D(D)
      ndD <- length(dD)
      D1 <- sapply(dD, function (x) sum(svD * x))
      D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
      cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
      out <- numeric(ndD)
      for (i in seq_along(dD)) {
        D.mat <- D2[i, ]
        dim(D.mat) <- c(ncz, ncz)
        out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) + sum((post.b %*% D.mat) * post.b, na.rm = TRUE)
      }
      J <- JM:::jacobian2(attr(D, "L"), ncz)
      drop(0.5 * (n * D1 - out) %*% J)
    }
    return(c(score.y, score.b, score.t1, score.t2))
  } else {
    # Scores of log(sigma)
    tr.idx <- matrix(1:ncz^2, ncz, ncz) # to extract separate block of Z'ZvarB
    tr.idx1 <- c(tr.idx[1:xDat1$ncz, 1:xDat1$ncz])
    tr.tZZvarb <- rowSums((ZtZ * post.vb)[, tr.idx1], na.rm = TRUE)
    sc.sigma.data <- split.data.frame(cbind(mu[1:N], ZEb[1:N]), xDat1$id)
    ni <- as.vector(tapply(xDat1$id, xDat1$id, length))
    sc.sigma.sub <- lapply(sc.sigma.data, function(x) crossprod(x[,1], (x[,1]-2*x[,2])) + crossprod(x[,2]))
    sclogsigma1 <- - sigma[1] * (- ni / sigma[1] + (unlist(sc.sigma.sub) + tr.tZZvarb) / sigma[1]^3)

    tr.idx2 <- c(tr.idx[-(1:xDat1$ncz), -(1:xDat1$ncz)])
    tr.tZZvarb <- rowSums((ZtZ * post.vb)[, tr.idx2], na.rm = TRUE)
    sc.sigma.data <- split.data.frame(cbind(mu[-(1:N)], ZEb[-c(1:N)]), xDat2$id)
    ni <- as.vector(tapply(xDat2$id, xDat2$id, length))
    sc.sigma.sub <- lapply(sc.sigma.data, function(x) crossprod(x[,1], (x[,1]-2*x[,2])) + crossprod(x[,2]))
    sclogsigma2 <- - sigma[2] * (- ni / sigma[2] + (unlist(sc.sigma.sub) + tr.tZZvarb) / sigma[2]^3)

    # Scores of betas
    longExtra1 <- c(listNamed(ZEb = ZEb[1:N], Zsb = Zsb[1:nGK, ], alpha = alpha1, sigma = sigma[1]), logLikSurv1[c("eta.tw1", "eta.ws")])
    scbetas1 <- gr.longSplinePH(betas[1:xDat1$ncx], xDat1, longExtra1, indScore = TRUE)
    longExtra2 <- c(listNamed(ZEb = ZEb[-(1:N)], Zsb = Zsb[-(1:nGK), ], alpha = alpha2, sigma = sigma[2]), logLikSurv2[c("eta.tw1", "eta.ws")])
    scbetas2 <- gr.longSplinePH(betas[-(1:xDat1$ncx)], xDat2, longExtra2, indScore = TRUE)

    score.y.sub <- cbind(scbetas1, scbetas2, sclogsigma1, sclogsigma2)

    # Scores of survival parameters
    survExtra1 <- list(Y = Y[1:n, ], Ys = Ys[1:nGK, ])
    list.survpars <- list(gammas = gammas1, alpha = alpha1, gammas.bs = gammas.bs1)
    list.survpars <- list.survpars[!sapply(list.survpars, is.null)]
    survpars <- unlist(as.relistable(list.survpars))
    score.t1.sub <- gr.survSplinePH(survpars, list.survpars, xDat1, survExtra1, indScore = TRUE)

    survExtra2 <- list(Y = Y[-(1:n), ], Ys = Ys[-(1:nGK), ])
    list.survpars <- list(gammas = gammas2, alpha = alpha2, gammas.bs = gammas.bs2)
    list.survpars <- list.survpars[!sapply(list.survpars, is.null)]
    survpars <- unlist(as.relistable(list.survpars))
    score.t2.sub <- gr.survSplinePH(survpars, list.survpars, xDat2, survExtra2, indScore = TRUE)

    # Scores of D
    diag.D <- !is.matrix(D)
    score.b.sub <- if (diag.D) {
      svD <- 1 / D
      svD2 <- svD^2
      dim(cS.postVB) <- c(ncz, ncz)
      D * 0.5 * (svD - diag(post.vb) * svD2 - as.matrix(post.b^2) * svD2)
    } else {
      svD <- solve(D)
      dD <- JM:::deriv.D(D)
      ndD <- length(dD)
      D1 <- sapply(dD, function (x) sum(svD * x))
      D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
      out <- matrix(0, n, ndD)
      for (i in seq_along(dD)) {
        D.mat <- D2[i, ]
        dim(D.mat) <- c(ncz, ncz)
        out[,i] <- rowSums(matrix(rep(D2[i, ], n), n, byrow=T) * post.vb, na.rm = TRUE) + rowSums((post.b %*% D.mat) * post.b, na.rm = TRUE)
      }
      J <- JM:::jacobian2(attr(D, "L"), ncz)
      drop(0.5 * (matrix(rep(D1, n), n, byrow = TRUE) - out) %*% J)
    }
    return(cbind(score.y.sub, score.b.sub, score.t1.sub, score.t2.sub))
  }
}


listNamed <- function(...){
  dots <- list(...)
  inferred <- sapply(substitute(list(...)), function(x) deparse(x)[1])[-1]
  if(is.null(names(inferred))){
    names(dots) <- inferred
  } else {
    names(dots)[names(inferred) == ""] <- inferred[names(inferred) == ""]
  }
  dots
}

## fd.vec1, slight modification from fd.vec function in JM R package
fd.vec1 <- function (x, f, ..., eps = 1e-05) { # forward numerical differentiation (extract score)
  n <- length(x)
  res <- matrix(0, n, n)
  ex <- pmax(abs(x), 1)
  f0 <- f(x, ...)
  for (i in 1:n) {
    x1 <- x
    x1[i] <- x[i] + eps * ex[i]
    diff.f <- c(f(x1, ...)$Score - f0$Score)  ## modified due to different data structure
    diff.x <- x1[i] - x[i]
    res[, i] <- diff.f / diff.x
  }
  res
}




#############################################################################################
## survfitJM: Predict conditional probability of main event for one id at one given time
## **************** In fact, this is 1-SurvProb, NOT survival prob as the source code
## Input:
##    object: fitted joint model with competing risk
##    newdata: longitudinal information for the predicted subject (only one subject is allowed)
##    last.time: the last time when the subject is known to be alive; if NULL then the recent measurement time
##    time.to.pred: time to predict
##    simulate: TRUE for proposed estimator; FALSE for empirical Bayes estimator [See Rizopoulos(2011)]
##    M: number of Monte Carlo samples
##    scale: a numeric scalar that controls the acceptance rate of the Metropolis-Hastings algorithm
## Output: conditional probability of main event prob under dependent and independent censoring
##############################################################################################
survfitJM_modified <- function (object, newdata, idVar = "id", deltaVar = "delta",
                                last.time = NULL, time.to.pred,
                                simulate = FALSE, M = 200, scale = 1.6) {

  ########################### Extract model info #############################
  ############################################################################
  timeVar <- object$timeVar
  parameterization <- object$parameterization ## "value"/"both"

  ########################### Extract longitudinal info ######################
  ############################################################################
  formYx1 <- object$xDat1$formYx
  formYx2 <- object$xDat2$formYx
  formYz <- object$xDat1$formYz               ## Z: ~time
  y1 <- newdata[, all.vars(formYx1)[1]]       ## y data
  y2 <- newdata[, all.vars(formYx2)[1]]       ## y data
  X1 <- model.matrix(formYx1, newdata)        ## X design matrix
  X2 <- model.matrix(formYx2, newdata)        ## X design matrix
  Z <- model.matrix(formYz, newdata)          ## Z design matrix

  id <- as.numeric(unclass(newdata[[idVar]]))  ## long id
  id <- match(id, unique(id))                  ## change id to 1, 2, 3...

  ########################### Extract surv info ##############################
  ############################################################################
  data.id <- newdata[!duplicated(id), ]          ## first row
  delta <- data.id[[deltaVar]]

  idT <- data.id[[idVar]]
  idT <- match(idT, unique(idT))

  TermsT <- object$xDat1$Wnm
  if(!length(TermsT)){TermsT <- "1"}
  formT <- reformulate(TermsT)
  W <- model.matrix(formT, data.id)              ## W design matrix

  ########################## Extract time.to.predict #########################
  ############################################################################
  obs.times <- split(newdata[[timeVar]], id)             ## Y measurements times
  last.time <- if(is.null(last.time)) sapply(obs.times, max) else last.time
  time.to.pred <- if(length(time.to.pred) == 1) rep(time.to.pred, length(last.time))
  if(any(last.time > time.to.pred) | length(last.time) != length(idT) | length(time.to.pred) != length(idT))
    stop("time.to.pred should be larger than last.time")

  ########################### Extract parameter info #########################
  ############################################################################
  ncx1 <- object$xDat1$ncx
  ncx2 <- object$xDat2$ncx
  ncz <- object$xDat1$ncz + object$xDat2$ncz
  ncww <- ncol(W)
  if (ncww == 1) {
    W <- NULL
    ncww <- 0
  } else {
    W <- W[,-1, drop = FALSE]
    ncww <- ncww - 1
  }

  betas1 <- object$coefficients[['betas']][1:(length(object$coefficients[['betas']])/2)]
  betas2 <- object$coefficients[['betas']][(length(object$coefficients[['betas']])/2+1):length(object$coefficients[['betas']])]
  sigma1 <- exp(object$coefficients[['log.sigma']])[1]
  sigma2 <- exp(object$coefficients[['log.sigma']])[2]
  D <- JM:::chol.transf(object$coefficients[['D']])
  diag.D <- ncol(D) == 1 & nrow(D) > 1
  D <- if (diag.D) diag(c(D)) else D
  gammas1 <- object$coefficients[['gammas1']]
  gammas2 <- object$coefficients[['gammas2']]
  alpha1 <- object$coefficients[['alpha1']]
  alpha2 <- object$coefficients[['alpha2']]
  gammas.bs1 <- object$coefficients[['gammas.bs1']]
  gammas.bs2 <- object$coefficients[['gammas.bs2']]
  list.thetas <- list(betas1 = betas1,  betas2 = betas2, log.sigma1 = log(sigma1), log.sigma2 = log(sigma2), gammas1 = gammas1, gammas2 = gammas2,
                      alpha1 = alpha1,  alpha2 = alpha2, gammas.bs1 = gammas.bs1, gammas.bs2 = gammas.bs2,
                      D = if (diag.D) log(diag(D)) else JM:::chol.transf(D))

  list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
  thetas <- unlist(as.relistable(list.thetas))
  Var.thetas <- object$VarCov

  environment(log.posterior.b_modified) <- environment(S.b_modified) <- environment(ModelMats_modified) <- environment()
  ############################################################################
  ### construct model matrices to calculate the survival functions
  n.tp <- length(time.to.pred)                  ## n.to.predict
  survMats <- survMats.last <- vector("list", n.tp)
  for(i in 1:n.tp){
    survMats[[i]] <- ModelMats_modified(time.to.pred[i], timeVar, data.id[i, ], formYx1, formYx2, formYz)
    survMats.last[[i]] <- ModelMats_modified(last.time[i], timeVar, data.id[i, ], formYx1, formYx2, formYz)
  }

  ### calculate the Empirical Bayes estimates and their (scaled) variance;
  ### ".new" are used in log.posterior.b() & S.b()
  modes.b <- matrix(0, n.tp, ncz)       ## modes.b for each subject to be predicted
  Vars.b <- vector("list", n.tp)
  betas1.new <- betas1; betas2.new <- betas2
  sigma1.new <- sigma1; sigma2.new <- sigma2
  D.new <- D
  gammas1.new <- gammas1; gammas2.new <- gammas2
  alpha1.new <- alpha1; alpha2.new <- alpha2
  gammas.bs1.new <- gammas.bs1; gammas.bs2.new <- gammas.bs2
  for(i in 1:n.tp){
    ff <- function (b, y1, y2, Mats, CR, ii) -log.posterior.b_modified(b, y1, y2, Mats, CR, ii)  ## use info upto last.time
    CR <- if(delta[i] == 2) TRUE else if(delta[i] == 0) FALSE
    opt <- optim(rep(0, ncz), ff, y1 = y1[id==i],  y2 = y2[id==i], Mats = survMats.last[[i]], CR = CR, ii = i,
                 method = "BFGS", hessian = TRUE)
    modes.b[i, ] <- opt$par
    Vars.b[[i]] <- scale * solve(opt$hessian)
  }

  if(!simulate) {
    v12 <- rep(NA, n.tp)
    for(i in 1:n.tp){
      S.last <- S.b_modified(last.time[i], modes.b[i, ], survMats.last[[i]], i)  ## use longitudinal upto last.time
      S.pred <- S.b_modified(time.to.pred[i], modes.b[i, ], survMats[[i]], i)
      v12[i] <- 1 - S.pred/S.last
    }
    return(list("out" = v12))
  } else {
    out <- matrix(NA, nrow = M, ncol = n.tp)     ## output to store v1, v2 (M*n.tp)
    success.rate <- matrix(FALSE, nrow = M, ncol = n.tp)
    b.old <- b.new <- modes.b
    if (n.tp == 1) dim(b.old) <- dim(b.new) <- c(1, ncz)
    thetas.all <- matrix(NA, nrow = M, ncol = dim(object$VarCov)[1])
    for (m in 1:M) {
      # Step 1: simulate new parameter values
      if(!is.positive.definite(Var.thetas)){
        Var.thetas <- JM:::nearPD(Var.thetas)
      }                                                        ### updated 6/21/2021
      thetas.new <- mvrnorm(1, thetas, Var.thetas)
      thetas.all[m, ] <- thetas.new
      thetas.new <- relist(thetas.new, skeleton = list.thetas)
      betas1.new <- thetas.new$betas1; betas2.new <- thetas.new$betas2
      sigma1.new <- exp(thetas.new$log.sigma1); sigma2.new <- exp(thetas.new$log.sigma2)
      gammas1.new <- thetas.new$gammas1; gammas2.new <- thetas.new$gammas2
      alpha1.new <- thetas.new$alpha1; alpha2.new <- thetas.new$alpha2
      D.new <- thetas.new$D
      D.new <- if (diag.D) exp(D.new) else JM:::chol.transf(D.new)
      gammas.bs1.new <- thetas.new$gammas.bs1; gammas.bs2.new <- thetas.new$gammas.bs2

      v12 <- rep(NA, n.tp)
      for(i in 1:n.tp){
        # Step 2: simulate new random effects values
        proposed.b <- JM:::rmvt(1, modes.b[i, ], Vars.b[[i]], 4)
        dmvt.old <- JM:::dmvt(b.old[i, ], modes.b[i, ], Vars.b[[i]], 4, TRUE)
        dmvt.proposed <- JM:::dmvt(proposed.b, modes.b[i, ], Vars.b[[i]], 4, TRUE)
        CR <- if(delta[i] == 2) TRUE else if(delta[i] == 0) FALSE
        a <- min(exp(log.posterior.b_modified(proposed.b, y1[id==i], y2[id==i], survMats.last[[i]], CR, i) + dmvt.old -
                       log.posterior.b_modified(b.old[i, ], y1[id==i], y2[id==i], survMats.last[[i]], CR, i) - dmvt.proposed), 1)
        ind <- (runif(1) <= a)
        success.rate[m, 1] <- ind
        if (!is.na(ind) && ind) b.new[i, ] <- proposed.b

        # Step 3: compute Pr(T1* <= t| theta.new, b.new)
        S.last <- S.b_modified(last.time[i], b.new[i, ], survMats.last[[i]], i)
        S.pred <- S.b_modified(time.to.pred[i], b.new[i, ], survMats[[i]], i)
        v12[i] <- 1 - S.pred/S.last
      }
      b.old <- b.new
      out[m, ] <- v12
    }
    return(list("out" = t(out), "success.rate" = colMeans(success.rate)))
  }
}


#############################################################################################
## log.posterior.b: Calculate log of joint p(T, delta, Y, b)
##    require extra theta.hat: not an arg in this function
## Input:
##    b: random effects
##    y: observed longitudinal biomarkers
##    Mats: output from ModelMats() to approximate integral in survival function
##    CR: TRUE if competing risk observed; FALSE if independent censor observed
## Output: joint probability
##############################################################################################
log.posterior.b_modified <- function (b, y1, y2, Mats, CR = T, ii) {
  ###############################################################

  id.i <- id %in% ii
  idT.i <- idT %in% ii
  X1.i <- X1[id.i, , drop = FALSE]
  X2.i <- X2[id.i, , drop = FALSE]
  Z.i <- Z[id.i, , drop = FALSE]
  W.i <- W[idT.i, , drop = FALSE]

  mu.y1 <- as.vector(X1.i %*% betas1.new) + rowSums(Z.i * rep(b[1:2], each = nrow(Z.i)))
  mu.y2 <- as.vector(X2.i %*% betas2.new) + rowSums(Z.i * rep(b[3:4], each = nrow(Z.i)))
  logNorm1 <- dnorm(y1, mu.y1, sigma1.new, TRUE)
  logNorm2 <- dnorm(y2, mu.y2, sigma2.new, TRUE)
  log.p.yb <- sum(logNorm1) + sum(logNorm2)                                ## log{Pr(Y|b)}

  log.p.b <- JM:::dmvnorm(b, rep(0, ncz), D.new, TRUE)     ## log{Pr(b)}

  st <- Mats$st         ## use the whole history of Y; st = t(1+sq)/2
  wk <- Mats$wk         ## Gauss-Kronrod weight
  P <- Mats$P           ## P = time/2
  X1s <- Mats$X1s       ## X history
  X2s <- Mats$X2s       ## X history
  Zs <- Mats$Zs
  if (parameterization %in% c("value", "both")){   ## Ys = m_i(s), trajectory
    Y1s <- as.vector(X1s %*% betas1.new + rowSums(Zs * rep(b[1:2], each = nrow(Zs))))
    Y2s <- as.vector(X2s %*% betas2.new + rowSums(Zs * rep(b[3:4], each = nrow(Zs))))
  }

  gamma.t.w1 <- ifelse(!is.null(W.i), as.vector(W.i %*% gammas1.new), 0)  ## eta.tw = gamma %*% w
  gamma.t.w2 <- ifelse(!is.null(W.i), as.vector(W.i %*% gammas2.new), 0)  ## eta.tw = gamma %*% w

  ## baseline hazard: method == "spline-PH-GH": JM_R_book P53
  Bspline1 <- splineDesign(object$control$knots[[1]], st, ord = object$control$ord, outer.ok = TRUE)       ## ord = 4: cubic spline
  Bspline2 <- splineDesign(object$control$knots[[1]], st, ord = object$control$ord, outer.ok = TRUE)       ## ord = 4: cubic spline

  logS.st1 <- exp(c(Bspline1 %*% gammas.bs1.new) + alpha1.new * Y1s)    ## time term: h01(st) + alpha1*m_i(st)
  logS.st2 <- exp(c(Bspline2 %*% gammas.bs2.new) + alpha2.new * Y2s)    ## h02(st) + alpha2*m_i(st)
  log.survival1 <- -exp(gamma.t.w1) * P * sum(wk * logS.st1)
  log.survival2 <- -exp(gamma.t.w2) * P * sum(wk * logS.st2)
  if (!CR) {
    return(log.survival1 + log.survival2 + log.p.yb + log.p.b)
  } else{
    log.haz2 <- log(logS.st2[length(logS.st2)]) + gamma.t.w2 ## log of competing risk hazard
    return(log.haz2 + log.survival1 + log.survival2 + log.p.yb + log.p.b)
  }
}


#############################################################################################
## S.b: Calculate survival prob S(t|biomarker history, thetas) for the main event (e.g. death)
##    require extra theta.hat: not an arg in this function
## Input:
##    t: time to predict
##    b: random effects
##    Mats: output from ModelMats() to approximate integral in survival function
## Output: survival probability
##############################################################################################
S.b_modified <- function (t, b, Mats, ii) {
  ###############################################################
  idT.i <- idT %in% ii
  st <- Mats$st
  wk <- Mats$wk
  P <- Mats$P
  X1s <- Mats$X1s
  X2s <- Mats$X2s
  Zs <- Mats$Zs
  if (parameterization %in% c("value", "both")){   ## Ys = m_i(s), trajectory
    Y1s <- as.vector(X1s %*% betas1.new + rowSums(Zs * rep(b[1:2], each = nrow(Zs))))
  }

  W.i <- W[idT.i, , drop = FALSE]
  gamma.t.w1 <- ifelse(!is.null(W.i), as.vector(W.i %*% gammas1.new), 0)  ## eta.tw = gamma %*% w

  ## baseline hazard: method == "spline-PH-GH": JM_R_book P53
  Bspline <- splineDesign(object$control$knots[[1]], st, ord = object$control$ord,
                          outer.ok = TRUE)                             ## ord = 4: cubic spline
  logS.st1 <- exp(c(Bspline %*% gammas.bs1.new) + alpha1.new * Y1s)    ## time term: h01(st) + alpha1*m_i(st)
  log.survival1 <- -exp(gamma.t.w1) * P * sum(wk * logS.st1)
  return(exp(log.survival1))
}


#############################################################################################
## ModelMats: Prepare for Gauss-Kronrod to approximate integral in survival function,
##    but not engage survival function
## Input:
##    time: the integral upper bound, i.e., use longitudinal history upto "time" in survival
##    timeVar: time variable for the longitudinal measurements
##    data.id: survival data to build time-dependent X & Z history
##    formYx: formula of X in longitudinal sub-model, to build X history
##    formYz: formula of Z in longitudinal sub-model, to build Z history
## Output: wk, sk, Xs & Zs matrix
##############################################################################################
ModelMats_modified <- function (time, timeVar, data.id, formYx1, formYx2, formYz) {
  id.GK <- rep(1, each = object$control$GKk + 1)        ## add 1 for log(h_02)

  wk <- c(JM:::gaussKronrod(object$control$GKk)$wk, 0)  ## add one 0 weight, a trick to calculate log(h_02)
  sk <- JM:::gaussKronrod(object$control$GKk)$sk
  P <- time / 2
  st <- c(P * (sk + 1), time)    ## st = t(1+sk)/2, add "time" to calculate log(h_02)

  data.id2 <- data.id[id.GK, ]   ## data.id: not an argument in this function
  data.id2[[timeVar]] <- st

  out <- list(st = st, wk = wk, P = P)
  out$X1s <- model.matrix(formYx1, data.id2)
  out$X2s <- model.matrix(formYx2, data.id2)
  out$Zs <- model.matrix(formYz, data.id2)
  out
}
