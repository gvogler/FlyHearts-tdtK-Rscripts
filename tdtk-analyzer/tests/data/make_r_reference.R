# Regenerates r_reference.json. Needs baseline.rollingBall.R from the CRAN "baseline" sources.
source("baseline.rollingBall.R")
set.seed(1)
out <- list()
for (n in c(60, 137, 500)) {
  x <- round(cumsum(rnorm(n)) + 5*sin(seq_len(n)/7), 2)
  out[[length(out)+1]] <- list(n=n, x=x, b=c(baseline.rollingBall(rbind(x), 20, 20)$baseline))
}
for (n in c(12, 30, 49, 50, 120, 400, 1000)) {
  t <- seq(0, by=0.0077, length.out=n)
  y <- sin(t*9) + rnorm(n, sd=0.1)
  s <- smooth.spline(t, y)
  g <- seq(0, max(s$x), by=0.0002)
  p <- predict(s, g)$y
  out[[length(out)+1]] <- list(n=n, t=t, y=y, spar=s$spar, pred=p[seq(1, length(p), length.out=25)], idx=seq(1, length(p), length.out=25))
}
xs <- c(3,5,7,9,NA,11); ys <- 1:6
z <- summary(lm(ys ~ xs)); out[[length(out)+1]] <- list(slope=z$coefficients[2], adj=z$adj.r.squared)
xs <- c(4,4,4); ys <- 1:3
z <- suppressWarnings(summary(lm(ys ~ xs))); out[[length(out)+1]] <- list(slope=z$coefficients[2], adj=z$adj.r.squared)
tojson <- function(v) {
  if (is.list(v)) {
    if (!is.null(names(v))) return(paste0("{", paste0('"', names(v), '":', sapply(v, tojson), collapse=","), "}"))
    return(paste0("[", paste(sapply(v, tojson), collapse=","), "]"))
  }
  f <- ifelse(is.na(v), "null", sprintf("%.17g", v))
  if (length(v) == 1) f else paste0("[", paste(f, collapse=","), "]")
}
writeLines(tojson(out), "r_reference.json")
