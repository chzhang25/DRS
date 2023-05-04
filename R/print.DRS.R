
print.DRS <- function(x, ...){
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  EST <- x$coef
  SE  <- x$se.coef
  result <- cbind(EST, SE, EST/SE, 2*(1-pnorm(abs(EST/SE))))
  colnames(result) <- c("Estimate", "StdErr", "z value", "Pr(>|z|)")
  print(round(result[-grep("gammas.bs|D", names(EST)), ], 4))
  invisible(x)
}
