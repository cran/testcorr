########################################################################################################################
### UNIVARIATE STATISTICS
########################################################################################################################

#' Testing zero autocorrelation
#'
#' The function ac.test computes the test statistics for examining the null hypothesis
#' of zero autocorrelation for univariate time series given in Dalla, Giraitis and Phillips (2019).
#'
#' @usage ac.test(x, max.lag, alpha = 0.05, lambda = 2.576, plot = TRUE,
#'         table = TRUE, var.name = NULL, scale.font = 1)
#' @param x A numeric vector or a univariate numeric time series object or a data frame.
#' @param max.lag Maximum lag at which to calculate the test statistics.
#' @param alpha Significance level for hypothesis testing used in the plots. Default is 0.05.
#' @param lambda Threshold in \eqn{\widetilde{Q}}{Q-tilde} test statistics. Default is 2.576.
#' @param plot Logical. If TRUE the sample autocorrelations with their confidence bands and the cumulative statistics with their critical values are plotted. Default is TRUE.
#' @param table Logical. If TRUE the sample autocorrelations, the confidence bands, the test statistics and their p-values are printed out. Default is TRUE.
#' @param var.name NULL or a character string specifying the variable name. If NULL and x has name, the name of x is used. If NULL and x has no name, the string "x" is used. Default is NULL.
#' @param scale.font A positive number indicating the scaling of the font size in the plots. Default is 1.
#' @details
#' The standard \eqn{t} and robust \eqn{\widetilde{t}}{t-tilde} statistics are for testing the null hypothesis \eqn{H_0:\rho_k=0}{H[0]:\rho[k]=0} at lags \eqn{k=1,...,max.lag},
#' and the standard \eqn{LB} and robust \eqn{\widetilde{Q}}{Q-tilde} statistics are for testing the null hypothesis \eqn{H_0:\rho_1=...=\rho_m=0}{H[0]:\rho[1]=...=\rho[m]=0} at lags \eqn{m=1,...,max.lag},
#' where \eqn{\rho_k}{\rho[k]} denotes the autocorrelation of \eqn{x_t}{x[t]} at lag \eqn{k}.
#' @return An object of class "ac.test", which is a list with the following components:
#' \item{lag}{The lags used.}
#' \item{ac}{The sample autocorrelations.}
#' \item{scb}{The lower and upper limit of the confidence bands based on the standard test statistics.}
#' \item{rcb}{The lower and upper limit of the confidence bands based on the robust test statistics.}
#' \item{t}{The \eqn{t}{t} test statistics.}
#' \item{pvt}{The p-values for the \eqn{t}{t} test statistics.}
#' \item{ttilde}{The \eqn{\widetilde{t}}{t-tilde} test statistics.}
#' \item{pvttilde}{The p-values for the \eqn{\widetilde{t}}{t-tilde} test statistics.}
#' \item{lb}{The \eqn{LB} test statistics.}
#' \item{pvlb}{The p-values for the \eqn{LB} test statistics.}
#' \item{qtilde}{The \eqn{\widetilde{Q}}{Q-tilde} test statistics.}
#' \item{pvqtilde}{The p-values for the \eqn{\widetilde{Q}}{Q-tilde} test statistics.}
#' @note
#' Missing values are not allowed.
#' @author
#' Violetta Dalla, Liudas Giraitis and Peter C. B. Phillips
#' @references
#' Dalla, V., Giraitis, L. and Phillips, P. C. B. (2019). "Robust Tests for White Noise and Cross-Correlation". Cowles Foundation, Discussion Paper No. 2194, \url{https://cowles.yale.edu/sites/default/files/files/pub/d21/d2194.pdf}.
#' @examples
#' x <- rnorm(100)
#' ac.test(x, max.lag = 10)
#' @importFrom stats acf pchisq pnorm qchisq qnorm is.ts
#' @importFrom assertthat is.string
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom knitr kable
#' @importFrom methods show
#' @export

ac.test <- function(x, max.lag, alpha = 0.05, lambda = 2.576, plot = TRUE, table = TRUE, var.name = NULL, scale.font = 1) {

  if (is.null(x)) stop()
  if (is.ts(x) | is.numeric(x) | is.data.frame(x)) {}
  else{
    stop('argument "x" must be numeric time series object or numeric vector or data frame')
  }
  if (!is.numeric(max.lag)) stop('argument "max.lag" must be numeric')
  if (!is.numeric(alpha)) stop('argument "alpha" must be numeric')
  if (!is.numeric(lambda)) stop('argument "lambda" must be numeric')
  if (!is.logical(plot)) stop('argument "plot" must be logical')
  if (!is.logical(table)) stop('argument "table" must be logical')
  if (!is.null(var.name) & NROW(var.name) == 1 & !is.string(var.name)) stop('argument "var.names" must be NULL or string')
  if (!is.numeric(scale.font)) stop('argument "scale.font" must be numeric')
  if (NCOL(x) != 1) stop('argument "x" must be univariate')
  if (any(is.na(x))) stop('argument "x" must not contain missing values')
  if (max.lag != ceiling(max.lag) | max.lag < 1 | max.lag > (NROW(x) - 1)) stop('argument "max.lag" must be an integer value greater or equal to one and less than the sample size')
  if (alpha < 0 | alpha > 1) stop('argument "alpha" must be between 0 and 1')
  if (lambda < 0) stop('argument "lambda" must be non-negative')
  if (is.null(var.name) & !is.null(colnames(x))) {
    if (NROW(colnames(x)) != NCOL(x)) stop('argument "x" must have one name')
  }
  if (!is.null(var.name)) {
    if (NROW(var.name) != NCOL(x)) stop('argument "var.names" must contain one name')
  }
  if (scale.font <= 0) stop('argument "scale.font" must be positive')

  if (is.null(var.name)) {
    if (!is.null(colnames(x))) {
      my.name <- colnames(x)
    }
    if (is.null(colnames(x))) {
      my.name <- "x"
    }
  }
  if (!is.null(var.name)) {
    my.name <- var.name
  }

  x <- as.matrix(x)
  x <- as.vector(x)

  n <- NROW(x)

  x.tilde <- x - mean(x)

  seq1maxlag <- seq(1, max.lag, 1)

  matrix.e <- x.tilde[2:n] * x.tilde[1:(n - 1)]
  for (kk in 2:max.lag) {
    matrix.e <- cbind(matrix.e, c(matrix(0, (kk - 1), 1), x.tilde[(kk + 1):n] * x.tilde[1:(n - kk)]))
  }

  matrix.esq <- matrix.e ^ 2

  t.tilde <- colSums(matrix.e) / sqrt(colSums(matrix.esq))

  big.rstar <- diag(max.lag) / 2
  for (jj in 1:(max.lag - 1)) {
    for (kk in (jj + 1):max.lag) {
      max.jk <- max(jj, kk)
      num.r <- t(matrix.e[, jj]) %*% matrix.e[, kk]
      rjkstar.ind <- num.r / sqrt(t(matrix.esq[max.jk:(n - 1), jj]) %*% matrix.esq[max.jk:(n - 1), kk])
      if (abs(rjkstar.ind) > lambda) {
        big.rstar[jj, kk] <- num.r / sqrt(sum(matrix.esq[max.jk:(n - 1), jj]) * sum(matrix.esq[max.jk:(n - 1), kk]))
      }
    }
  }

  big.rstar <- big.rstar + t(big.rstar)

  q.tilde <- matrix(nrow = max.lag, ncol = 1)
  for (mm in 1:max.lag) {
    q.tilde[mm] <- t(t.tilde[1:mm]) %*% solve(big.rstar[1:mm, 1:mm]) %*% t.tilde[1:mm]
  }

  is.na(q.tilde) <- q.tilde <= 0

  ac <- stats::acf(x, max.lag, plot = FALSE)$acf[2:(max.lag + 1)]

  t <- sqrt(n) * ac

  sq.ac <- ac ^ 2
  seq.overnmink <- (n - seq1maxlag) ^ (-1)
  lb <- n * (n + 2) * cumsum(sq.ac * seq.overnmink)

  z.cv <- qnorm(1 - alpha / 2)
  vec1 <- matrix(1, nrow = max.lag, ncol = 1)
  s.cb <- z.cv / sqrt(n) * vec1
  r.cb <- z.cv * ac / t.tilde

  lu.s.cb <- cbind(-s.cb, s.cb)
  lu.r.cb <- cbind(-r.cb, r.cb)
  rownames(lu.r.cb) <- rownames(lu.s.cb)
  colnames(lu.r.cb) <- colnames(lu.s.cb)

  pv.lb <- 1 - pchisq(lb, seq1maxlag)
  pv.qtilde <- 1 - pchisq(q.tilde, seq1maxlag)

  if (plot == TRUE) {

    chisq.cv <- qchisq(alpha, seq1maxlag, lower.tail = FALSE)

    plotcorr(max.lag, seq1maxlag, ac, s.cb, r.cb, alpha, n, 1, my.name, scale.font)
    plotstat(max.lag, seq1maxlag, lb, q.tilde, alpha, chisq.cv, n, 1, my.name, "LB", expression(tilde(Q)), 1, scale.font)
  }

  if (table == TRUE) {
    options(max.print = 1000000)

    results.cb <- cbind(ac, s.cb, r.cb)
    results.cb <- format(round(results.cb, 3), nsmall = 3)
    results.t <- cbind(t, 2 * (1 - pnorm(abs(t))), t.tilde, 2 * (1 - pnorm(abs(t.tilde))))
    results.t <- format(round(results.t, 3), nsmall = 3)
    results.q <- cbind(lb, pv.lb, q.tilde, pv.qtilde)
    results.q <- format(round(results.q, 3), nsmall = 3)

    results.tq <- cbind(seq1maxlag, results.cb, seq1maxlag, results.t, seq1maxlag, results.q)
    results.tq <- data.frame(results.tq, row.names = NULL)
    names(results.tq) <- c("Lag", "AC", paste("Stand. CB(", 100 * (1 - alpha), "%)", sep = ""), paste("Robust CB(", 100 * (1 - alpha), "%)", sep = ""), "Lag", "t", "p-value", "t-tilde", "p-value", "Lag", "LB", "p-value", "Q-tilde", "p-value")

    results.cb.lu <- cbind(s.cb, -s.cb, r.cb, -r.cb)
    results.cb.lu <- format(round(results.cb.lu, 3), nsmall = 3)
    results.tq[, 3] <- paste(results.cb.lu[, 2], paste(",", results.cb.lu[, 1], sep = ""), sep = "")
    results.tq[, 4] <- paste(results.cb.lu[, 4], paste(",", results.cb.lu[, 3], sep = ""), sep = "")
    results.tq[, 3] <- paste("(", paste(results.tq[, 3], ")", sep = ""), sep = "")
    results.tq[, 4] <- paste("(", paste(results.tq[, 4], ")", sep = ""), sep = "")

    table.tq <- kable(results.tq, align = "r")
    cat("\n")
    cat("Tests for zero autocorrelation of ", my.name, sep = "")
    print(table.tq)
    cat("\n")
  }

  colnames(seq1maxlag) <- NULL
  rownames(seq1maxlag) <- NULL
  colnames(ac) <- NULL
  rownames(ac) <- NULL
  colnames(lu.s.cb) <- NULL
  rownames(lu.s.cb) <- NULL
  colnames(lu.r.cb) <- NULL
  rownames(lu.r.cb) <- NULL
  colnames(t) <- NULL
  rownames(t) <- NULL
  colnames(t.tilde) <- NULL
  rownames(t.tilde) <- NULL
  colnames(lb) <- NULL
  rownames(lb) <- NULL
  colnames(q.tilde) <- NULL
  rownames(q.tilde) <- NULL
  colnames(pv.lb) <- NULL
  rownames(pv.lb) <- NULL
  colnames(pv.qtilde) <- NULL
  rownames(pv.qtilde) <- NULL

  invisible(structure(list(lag = as.numeric(seq1maxlag), ac = as.numeric(ac), scb = as.matrix(lu.s.cb), rcb = as.matrix(lu.r.cb), t = as.numeric(t), pvt = as.numeric(2 * (1 - pnorm(abs(t)))), ttilde = as.numeric(t.tilde), pvttilde = as.numeric(2 * (1 - pnorm(abs(t.tilde)))), lb = as.numeric(lb), pvlb = as.numeric(pv.lb), qtilde = as.numeric(q.tilde), pvqtilde = as.numeric(pv.qtilde)), class = "ac.test"))

}

########################################################################################################################
### UNIVARIATE STATISTICS: iid
########################################################################################################################

#' Testing iid property
#'
#' The function iid.test computes the test statistics for examining the null hypothesis
#' of i.i.d. property for univariate series given in Dalla, Giraitis and Phillips (2019).
#'
#' @usage iid.test(x, max.lag, alpha = 0.05, plot = TRUE, table = TRUE,
#'          var.name = NULL, scale.font = 1)
#' @param x A numeric vector or a univariate numeric time series object or a data frame.
#' @param max.lag Maximum lag at which to calculate the test statistics.
#' @param alpha Significance level for hypothesis testing used in the plots. Default is 0.05.
#' @param plot Logical. If TRUE the test statistics and their critical values are plotted. Default is TRUE.
#' @param table Logical. If TRUE the test statistics and their p-values are printed out. Default is TRUE.
#' @param var.name NULL or a character string specifying the variable name. If NULL and x has name, the name of x is used. If NULL and x has no name, the string "x" is used. Default is NULL.
#' @param scale.font A positive number indicating the scaling of the font size in the plots. Default is 1.
#' @details
#' The \eqn{J_{x,|x|}}{J[x,|x|]} and \eqn{J_{x,x^2}}{J[x,x^2]} statistics are for testing the null hypothesis of i.i.d. at lag \eqn{k}, \eqn{k=1,...,max.lag},
#' and the \eqn{C_{x,|x|}}{C[x,|x|]} and \eqn{C_{x,x^2}}{C[x,x^2]} statistics are for testing the null hypothesis of i.i.d. at lags \eqn{1,...,m}, \eqn{m=1,...,max.lag}.
#' @return An object of class "iid.test", which is a list with the following components:
#' \item{lag}{The lags used.}
#' \item{jab}{The \eqn{J_{x,|x|}}{J[x,|x|]} test statistics.}
#' \item{pvjab}{The p-values for the \eqn{J_{x,|x|}}{J[x,|x|]} test statistics.}
#' \item{jsq}{The \eqn{J_{x,x^2}}{J[x,x^2]} test statistics.}
#' \item{pvjsq}{The p-values for the \eqn{J_{x,x^2}}{J[x,x^2]} test statistics.}
#' \item{cab}{The \eqn{C_{x,|x|}}{C[x,|x|]} test statistics.}
#' \item{pvcab}{The p-values for the \eqn{C_{x,|x|}}{C[x,|x|]} test statistics.}
#' \item{csq}{The \eqn{C_{x,x^2}}{C[x,x^2]} test statistics.}
#' \item{pvcsq}{The p-values for the \eqn{C_{x,x^2}}{C[x,x^2]} test statistics.}
#' @note
#' Missing values are not allowed.
#' @author
#' Violetta Dalla, Liudas Giraitis and Peter C. B. Phillips
#' @references
#' Dalla, V., Giraitis, L. and Phillips, P. C. B. (2019). "Robust Tests for White Noise and Cross-Correlation". Cowles Foundation, Discussion Paper No. 2194, \url{https://cowles.yale.edu/sites/default/files/files/pub/d21/d2194.pdf}.
#' @examples
#' x <- rnorm(100)
#' iid.test(x, max.lag = 10)
#' @importFrom stats acf pchisq pnorm qchisq qnorm
#' @importFrom assertthat is.string
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom knitr kable
#' @importFrom methods show
#' @export

iid.test <- function(x, max.lag, alpha = 0.05, plot = TRUE, table = TRUE, var.name = NULL, scale.font = 1) {

  if (is.null(x)) stop()
  if (is.ts(x) | is.numeric(x) | is.data.frame(x)) {}
  else{
    stop('argument "x" must be numeric time series object or numeric vector or data frame')
  }
  if (!is.numeric(max.lag)) stop('argument "max.lag" must be numeric')
  if (!is.numeric(alpha)) stop('argument "alpha" must be numeric')
  if (!is.logical(plot)) stop('argument "plot" must be logical')
  if (!is.logical(table)) stop('argument "table" must be logical')
  if (!is.null(var.name) & NROW(var.name) == 1 & !is.string(var.name)) stop('argument "var.names" must be NULL or string')
  if (!is.numeric(scale.font)) stop('argument "scale.font" must be numeric')
  if (NCOL(x) != 1) stop('argument "x" must be univariate')
  if (any(is.na(x))) stop('argument "x" must not contain missing values')
  if (max.lag != ceiling(max.lag) | max.lag < 1 | max.lag > (NROW(x) - 1)) stop('argument "max.lag" must be an integer value greater or equal to one and less than the sample size')
  if (alpha < 0 | alpha > 1) stop('argument "alpha" must be between 0 and 1')
  if (is.null(var.name) & !is.null(colnames(x))) {
    if (NROW(colnames(x)) != NCOL(x)) stop('argument "x" must have one name')
  }
  if (!is.null(var.name)) {
    if (NROW(var.name) != NCOL(x)) stop('argument "var.names" must contain one name')
  }
  if (scale.font <= 0) stop('argument "scale.font" must be positive')

  if (is.null(var.name)) {
    if (!is.null(colnames(x))) {
      my.name <- colnames(x)
    }
    if (is.null(colnames(x))) {
      my.name <- "x"
    }
  }
  if (!is.null(var.name)) {
    my.name <- var.name
  }

  x <- as.matrix(x)
  x <- as.vector(x)

  n <- NROW(x)

  x.tilde <- x - mean(x)

  seq1maxlag <- seq(1, max.lag, 1)

  ac <- acf(x, max.lag, plot = FALSE)$acf[2:(max.lag + 1)]
  ac.abs <- acf(abs(x.tilde), max.lag, plot = FALSE)$acf[2:(max.lag + 1)]
  ac.sq <- acf(x.tilde ^ 2, max.lag, plot = FALSE)$acf[2:(max.lag + 1)]

  j.abs <- (n ^ 2) / (n - seq1maxlag) * (ac ^ 2 + ac.abs ^ 2)
  j.sq <- (n ^ 2) / (n - seq1maxlag) * (ac ^ 2 + ac.sq ^ 2)

  c.abs <- cumsum(j.abs)
  c.sq <- cumsum(j.sq)

  if (plot == TRUE) {
    chisq2.cv <- matrix(qchisq(alpha, 2, lower.tail = FALSE), max.lag, 1)
    chisq.cv <- qchisq(alpha, 2 * seq1maxlag, lower.tail = FALSE)

    fvar.name <- substr(my.name, 1, 1)
    sub.abs <- bquote(.(fvar.name) * ",|" * .(fvar.name) * "|")
    sub.sq <- bquote(.(fvar.name) * "," * .(fvar.name)^2)

    name.j.abs <- bquote(J[.(sub.abs)])
    name.j.sq <- bquote(J[.(sub.sq)])
    name.c.abs <- bquote(C[.(sub.abs)])
    name.c.sq <- bquote(C[.(sub.sq)])

    plotstat(max.lag, seq1maxlag, j.abs, j.sq, alpha, chisq2.cv, n, 1, my.name, name.j.abs, name.j.sq, 2, scale.font)
    plotstat(max.lag, seq1maxlag, c.abs, c.sq, alpha, chisq.cv, n, 1, my.name, name.c.abs, name.c.sq, 3, scale.font)
  }

  if (table == TRUE) {
    options(max.print = 1000000)

    results.j <- cbind(j.abs, 1 - pchisq(j.abs, 2), j.sq, 1 - pchisq(j.sq, 2))
    results.j <- format(round(results.j, 3), nsmall = 3)
    results.c <- cbind(c.abs, 1 - pchisq(c.abs, 2 * seq1maxlag), c.sq, 1 - pchisq(c.sq, 2 * seq1maxlag))
    results.c <- format(round(results.c, 3), nsmall = 3)

    results.jc <- cbind(seq1maxlag, results.j, seq1maxlag, results.c)
    results.jc <- data.frame(results.jc, row.names = NULL)

    fvar.name <- substr(my.name, 1, 1)
    name.j.abs <- paste("J[", fvar.name, ",", "\u01C0", fvar.name, "\u01C0]", sep = "")
    name.j.sq <- paste("J[", fvar.name, ",", fvar.name, "\u00B2]", sep = "")
    name.c.abs <- paste("C[", fvar.name, ",", "\u01C0", fvar.name, "\u01C0]", sep = "")
    name.c.sq <- paste("C[", fvar.name, ",", fvar.name, "\u00B2]", sep = "")

    names(results.jc) <- c("Lag", name.j.abs, "p-value", name.j.sq, "p-value", "Lag", name.c.abs, "p-value", name.c.sq, "p-value")

    table.jc <- kable(results.jc, align = "r")
    cat("\n")
    cat("Tests for i.i.d. property of ", my.name, sep = "")
    print(table.jc, justify = "right")
    cat("\n")
  }

  colnames(seq1maxlag) <- NULL
  rownames(seq1maxlag) <- NULL
  colnames(j.abs) <- NULL
  rownames(j.abs) <- NULL
  colnames(j.sq) <- NULL
  rownames(j.sq) <- NULL
  colnames(c.abs) <- NULL
  rownames(c.abs) <- NULL
  colnames(c.sq) <- NULL
  rownames(c.sq) <- NULL

  invisible(structure(list(lag = as.numeric(seq1maxlag), jabs = as.numeric(j.abs), pvjabs = as.numeric(1 - pchisq(j.abs, 2)), jsq = as.numeric(j.sq), pvjsq = as.numeric(1 - pchisq(j.sq, 2)), cabs = as.numeric(c.abs), pvcabs = as.numeric(1 - pchisq(c.abs, 2 * seq1maxlag)), csq = as.numeric(c.sq), pvcsq = as.numeric(1 - pchisq(c.sq, 2 * seq1maxlag))), class = "iid.test"))

}

########################################################################################################################
### BIVARIATE STATISTICS
########################################################################################################################

#' Testing zero cross-correlation
#'
#' The function cc.test computes the test statistics for examining the null hypothesis
#' of zero cross-correlation for bivariate time series given in Dalla, Giraitis and Phillips (2019).
#'
#' @usage cc.test(x, y, max.lag, alpha = 0.05, lambda = 2.576, plot = TRUE,
#'         table = TRUE, var.names = NULL, scale.font = 1)
#' @param x A numeric vector or a univariate numeric time series object or a data frame.
#' @param y A numeric vector or a univariate numeric time series object or a data frame.
#' @param max.lag Maximum lag at which to calculate the test statistics.
#' @param alpha Significance level for hypothesis testing used in the plots. Default is 0.05.
#' @param lambda Threshold in \eqn{\widetilde{Q}}{Q-tilde} test statistics. Default is 2.576.
#' @param plot Logical. If TRUE the sample cross-correlations with their confidence bands and the cumulative statistics with their critical values are plotted. Default is TRUE.
#' @param table Logical. If TRUE the sample cross-correlations, the confidence bands, the test statistics and their p-values are printed out. Default is TRUE.
#' @param var.names NULL or a character string specifying the variable names. If NULL and x,y have names, the names of x,y are used. If NULL and x,y have no names, the string c("x","y") is used. Default is NULL.
#' @param scale.font A positive number indicating the scaling of the font size in the plots. Default is 1.
#' @details
#' The standard \eqn{t} and robust \eqn{\widetilde{t}}{t-tilde} statistics are for testing the null hypothesis \eqn{H_0:\rho_k=0}{H[0]:\rho[k]=0} at lags \eqn{k=-max.lag,...,-1,0,1,max.lag},
#' and the standard \eqn{HB} and robust \eqn{\widetilde{Q}}{Q-tilde} statistics are for testing the null hypothesis \eqn{H_0:\rho_0=...=\rho_m=0}{H[0]:\rho[0]=...=\rho[m]=0} at lags \eqn{m=-max.lag,...,-1,0,1,max.lag},
#' where \eqn{\rho_k}{\rho[k]} denotes the cross-correlation of \eqn{x_t}{x[t]} and \eqn{y_{t-k}}{y[t-k]} at lag \eqn{k}.
#' @return An object of class "cc.test", which is a list with the following components:
#' \item{lag}{The lags used.}
#' \item{cc}{The sample cross-correlations.}
#' \item{scb}{The lower and upper limit of the confidence bands based on the standard test statistics.}
#' \item{rcb}{The lower and upper limit of the confidence bands based on the robust test statistics.}
#' \item{t}{The \eqn{t}{t} test statistics.}
#' \item{pvt}{The p-values for the \eqn{t}{t} test statistics.}
#' \item{ttilde}{The \eqn{\widetilde{t}}{t-tilde} test statistics.}
#' \item{pvtttilde}{The p-values for the \eqn{\widetilde{t}}{t-tilde} test statistics.}
#' \item{hb}{The \eqn{HB} test statistics.}
#' \item{pvhb}{The p-values for the \eqn{HB} test statistics.}
#' \item{qtilde}{The \eqn{\widetilde{Q}}{Q-tilde} test statistics.}
#' \item{pvqtilde}{The p-values for the \eqn{\widetilde{Q}}{Q-tilde} test statistics.}
#' @note
#' Missing values are not allowed.
#' @author
#' Violetta Dalla, Liudas Giraitis and Peter C. B. Phillips
#' @references
#' Dalla, V., Giraitis, L. and Phillips, P. C. B. (2019). "Robust Tests for White Noise and Cross-Correlation". Cowles Foundation, Discussion Paper No. 2194, \url{https://cowles.yale.edu/sites/default/files/files/pub/d21/d2194.pdf}.
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' cc.test(x, y, max.lag = 10)
#' @importFrom stats acf ccf pchisq pnorm qchisq qnorm
#' @importFrom assertthat is.string
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom knitr kable
#' @importFrom methods show
#' @export

cc.test <- function(x, y, max.lag, alpha = 0.05, lambda = 2.576, plot = TRUE, table = TRUE, var.names = NULL, scale.font = 1) {

  if (is.null(x)) stop()
  if (is.null(y)) stop()
  if (is.ts(x) | is.numeric(x) | is.data.frame(x)) {}
  else{
    stop('argument "x" must be numeric time series object or numeric vector or data frame')
  }
  if (is.ts(y) | is.numeric(y) | is.data.frame(y)) {}
  else{
    stop('argument "y" must be numeric time series object or numeric vector or data frame')
  }
  if (!is.numeric(max.lag)) stop('argument "max.lag" must be numeric')
  if (!is.numeric(alpha)) stop('argument "alpha" must be numeric')
  if (!is.numeric(lambda)) stop('argument "lambda" must be numeric')
  if (!is.logical(plot)) stop('argument "plot" must be logical')
  if (!is.logical(table)) stop('argument "table" must be logical')
  if (!is.null(var.names) & NROW(var.names) == 2 & (!is.string(var.names[1]) | !is.string(var.names[2])))  stop('argument "var.names" must be NULL or string')
  if (!is.numeric(scale.font)) stop('argument "scale.font" must be numeric')
  if (NCOL(x) != 1) stop('argument "x" must be univariate')
  if (NCOL(y) != 1) stop('argument "y" must be univariate')
  if (any(is.na(x))) stop('argument "x" must not contain missing values')
  if (any(is.na(y))) stop('argument "y" must not contain missing values')
  if (NCOL(x) != NCOL(y)) stop('arguments "x" and "y" must have the same length')
  if (max.lag != ceiling(max.lag) | max.lag < 1 | max.lag > (NROW(x) - 1)) stop('argument "max.lag" must be an integer value greater or equal to one and less than the sample size')
  if (alpha < 0 | alpha > 1) stop('argument "alpha" must be between 0 and 1')
  if (lambda < 0) stop('argument "lambda" must be non-negative')
  if (is.null(var.names) & !is.null(colnames(x))) {
    if (NROW(colnames(x)) != NCOL(x)) stop('argument "x" must have two names')
  }
  if (!is.null(var.names)) {
    if (NROW(var.names) != 2) stop('argument "var.names" must contain two names')
  }
  if (scale.font <= 0) stop('argument "scale.font" must be positive')

  if (is.null(var.names)) {
    if (!is.null(colnames(x))) {
      my.names <- colnames(x)
    }
    if (is.null(colnames(x))) {
      my.names <- c("x","y")
    }
  }
  if (!is.null(var.names)) {
    my.names <- var.names
  }

  x <- as.matrix(x)
  x <- as.vector(x)
  y <- as.matrix(y)
  y <- as.vector(y)

  n <- NROW(x)

  results.tq.xyk <- cc.test.t.tmink(x, y, max.lag, alpha, lambda)
  results.tq.yxk <- cc.test.t.tmink(y, x, max.lag, alpha, lambda)

  results.tq <- rbind(results.tq.yxk[(max.lag + 1):2, ], results.tq.xyk)

  seqmaxlagmaxlag <- seq(-max.lag, max.lag, 1)

  cc <- results.tq[, 1]
  t.tilde <- results.tq[, 2]
  q.tilde <- results.tq[, 3]
  t <- results.tq[, 4]
  hb <- results.tq[, 5]

  z.cv <- qnorm(1 - alpha / 2)
  vec1 <- matrix(1, nrow = 2 * max.lag + 1, ncol = 1)
  s.cb <- z.cv / sqrt(n) * vec1
  r.cb <- z.cv * cc / t.tilde

  lu.s.cb <- cbind(-s.cb, s.cb)
  lu.r.cb <- cbind(-r.cb, r.cb)
  rownames(lu.r.cb) <- rownames(lu.s.cb)
  colnames(lu.r.cb) <- colnames(lu.s.cb)

  if (plot == TRUE) {
    chisq.cv <- qchisq(alpha, abs(seqmaxlagmaxlag) + 1, lower.tail = FALSE)

    plotcorr(max.lag, seqmaxlagmaxlag, cc, s.cb, r.cb, alpha, n, 2, my.names, scale.font)
    plotstat(max.lag, seqmaxlagmaxlag, hb, q.tilde, alpha, chisq.cv, n, 2, my.names, "HB", expression(tilde(Q)), 1, scale.font)
  }

  if (table == TRUE) {
    options(max.print = 1000000)

    results.cb <- cbind(cc, s.cb, r.cb)
    results.cb <- format(round(results.cb, 3), nsmall = 3)
    results.t <- cbind(t, 2 * (1 - pnorm(abs(t))), t.tilde, 2 * (1 - pnorm(abs(t.tilde))))
    results.t <- format(round(results.t, 3), nsmall = 3)
    results.q <- cbind(hb, 1 - pchisq(hb, abs(seqmaxlagmaxlag) + 1), q.tilde, 1 - pchisq(q.tilde, abs(seqmaxlagmaxlag) + 1))
    results.q <- format(round(results.q, 3), nsmall = 3)

    results.tq <- cbind(seqmaxlagmaxlag, results.cb, seqmaxlagmaxlag, results.t, seqmaxlagmaxlag, results.q)
    results.tq <- data.frame(results.tq, row.names = NULL)
    names(results.tq) <- c("Lag", "CC", paste("Stand. CB(", 100 * (1 - alpha), "%)", sep = ""), paste("Robust CB(", 100 * (1 - alpha), "%)", sep = ""), "Lag", "t", "p-value", "t-tilde", "p-value", "Lag", "HB", "p-value", "Q-tilde", "p-value")

    results.cb.lu <- cbind(s.cb, -s.cb, r.cb, -r.cb)
    results.cb.lu <- format(round(results.cb.lu, 3), nsmall = 3)
    results.tq[, 3] <- paste(results.cb.lu[, 2], paste(",", results.cb.lu[, 1], sep = ""), sep = "")
    results.tq[, 4] <- paste(results.cb.lu[, 4], paste(",", results.cb.lu[, 3], sep = ""), sep = "")
    results.tq[, 3] <- paste("(", paste(results.tq[, 3], ")", sep = ""), sep = "")
    results.tq[, 4] <- paste("(", paste(results.tq[, 4], ")", sep = ""), sep = "")

    table.tq <- kable(results.tq, align = "r")
    cat("\n")
    cat("Tests for zero cross-correlation of ", my.names[1], " and ", my.names[2], sep = "")
    print(table.tq)
    cat("\n")
  }

  colnames(seqmaxlagmaxlag) <- NULL
  rownames(seqmaxlagmaxlag) <- NULL
  colnames(cc) <- NULL
  rownames(cc) <- NULL
  colnames(lu.s.cb) <- NULL
  rownames(lu.s.cb) <- NULL
  colnames(lu.r.cb) <- NULL
  rownames(lu.r.cb) <- NULL
  colnames(t) <- NULL
  rownames(t) <- NULL
  colnames(t.tilde) <- NULL
  rownames(t.tilde) <- NULL
  colnames(hb) <- NULL
  rownames(hb) <- NULL
  colnames(q.tilde) <- NULL
  rownames(q.tilde) <- NULL

  invisible(structure(list(lag = as.numeric(seqmaxlagmaxlag), cc = as.numeric(cc), scb = as.matrix(lu.s.cb), rcb = as.matrix(lu.r.cb), t = as.numeric(t), pvt = as.numeric(2 * (1 - pnorm(abs(t)))), ttilde = as.numeric(t.tilde), pvttilde = as.numeric(2 * (1 - pnorm(abs(t.tilde)))), hb = as.numeric(hb), pvhb = as.numeric(1 - pchisq(hb, abs(seqmaxlagmaxlag) + 1)), qtilde = as.numeric(q.tilde), pvqtilde = as.numeric(1 - pchisq(q.tilde, abs(seqmaxlagmaxlag) + 1))), class = "cc.test"))

}

cc.test.t.tmink <- function(x, y, max.lag, alpha, lambda) {

  n <- NROW(x)

  x.tilde <- x - mean(x)
  y.tilde <- y - mean(y)

  matrix.e <- x.tilde * y.tilde
  for (kk in 1:max.lag) {
    matrix.e <- cbind(matrix.e, c(matrix(0, kk, 1), x.tilde[(kk + 1):n] * y.tilde[1:(n - kk)]))
  }

  matrix.esq <- matrix.e ^ 2

  t.tilde <- colSums(matrix.e) / sqrt(colSums(matrix.esq))

  big.rstar <- diag(max.lag + 1) / 2
  for (jj in 1:max.lag) {
    for (kk in (jj + 1):(max.lag + 1)) {
      max.jk <- max(jj, kk)
      num.r <- t(matrix.e[, jj]) %*% matrix.e[, kk]
      rjkstar.ind <- num.r / sqrt(t(matrix.esq[max.jk:n, jj]) %*% matrix.esq[max.jk:n, kk])
      if (abs(rjkstar.ind) > lambda) {
        big.rstar[jj, kk] <- num.r / sqrt(sum(matrix.esq[max.jk:n, jj]) * sum(matrix.esq[max.jk:n, kk]))
      }
    }
  }

  big.rstar <- big.rstar + t(big.rstar)

  q.tilde <- matrix(nrow = max.lag + 1, ncol = 1)
  for (mm in 1:(max.lag + 1)) {
    q.tilde[mm] <- t(t.tilde[1:mm]) %*% solve(big.rstar[1:mm, 1:mm]) %*% t.tilde[1:mm]
  }

  is.na(q.tilde) <- q.tilde <= 0

  cc <- ccf(x, y, max.lag, plot = FALSE)$acf[(max.lag + 1):(2 * max.lag + 1)]

  t <- sqrt(n) * cc

  seq0maxlag <- seq(0, max.lag, 1)
  sq.cc <- cc ^ 2
  seq.overnmink <- (n - seq0maxlag) ^ (-1)
  hb <- (n ^ 2) * cumsum(sq.cc * seq.overnmink)

  cbind(cc, t.tilde, q.tilde, t, hb)

}

########################################################################################################################
### BIVARIATE STATISTICS: special case k=0
########################################################################################################################

#' Testing zero Pearson correlation
#'
#' The function rcorr.test computes the test statistics for examining the null hypothesis
#' of zero Pearson correlation for multivariate series in Dalla, Giraitis and Phillips (2019).
#'
#' @usage rcorr.test(x, plot = TRUE, table = TRUE, var.names = NULL,
#'            scale.font = 1)
#' @param x A numeric matrix or a multivariate numeric time series object or a data frame.
#' @param plot Logical. If TRUE the sample Pearson correlations and the p-values for significance are plotted. Default is TRUE.
#' @param table Logical. If TRUE the sample Pearson correlations and the p-values for significance are printed out. Default is TRUE.
#' @param var.names NULL or a character string specifying the variable names. If NULL and x has names, the names of x are used. If NULL and x has no names, the string c("x[1]","x[2]",...) is used. Default is NULL.
#' @param scale.font A positive number indicating the scaling of the font size in the plots. Default is 1.
#' @details
#' The p-value of the robust \eqn{\widetilde{t}}{t-tilde} statistic is for testing the null hypothesis \eqn{H_0:\rho_{i,j}=0}{H[0]:\rho[i,j]=0},
#' where \eqn{\rho_{i,j}}{\rho[i,j]} denotes the correlation of \eqn{x_{i}}{x[i]} and \eqn{x_{j}}{x[j]}.
#' @return An object of class "rcorr.test", which is a list with the following components:
#' \item{pc}{The sample Pearson correlations.}
#' \item{pv}{The p-values for the \eqn{\widetilde{t}}{t-tilde} test statistics.}
#' @note
#' Missing values are not allowed.
#' @author
#' Violetta Dalla, Liudas Giraitis and Peter C. B. Phillips
#' @references
#' Dalla, V., Giraitis, L. and Phillips, P. C. B. (2019). "Robust Tests for White Noise and Cross-Correlation". Cowles Foundation, Discussion Paper No. 2194, \url{https://cowles.yale.edu/sites/default/files/files/pub/d21/d2194.pdf}.
#' @examples
#' x <- matrix(rnorm(400),100)
#' rcorr.test(x)
#' @importFrom stats cor pnorm
#' @importFrom assertthat is.string
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom forcats fct_rev
#' @importFrom methods show
#' @export

rcorr.test <- function(x, plot = TRUE, table = TRUE, var.names = NULL, scale.font = 1) {

  if (is.null(x)) stop()
  if (is.ts(x) | is.numeric(x) | is.data.frame(x)) {}
  else{
    stop('argument "x" must be numeric time series object or numeric matrix or data frame')
  }
  if (!is.logical(plot)) stop('argument "plot" must be logical')
  if (!is.logical(table)) stop('argument "table" must be logical')
  if (!is.null(var.names) & NROW(var.names) > 1 & !is.string(var.names)) stop('argument "var.names" must be NULL or string')
  if (!is.numeric(scale.font)) stop('argument "scale.font" must be numeric')
  if (NCOL(x) == 1) stop('argument "x" must be multivariate')
  if (any(is.na(x))) stop('argument "x" must not contain missing values')
  if (is.null(var.names) & !is.null(colnames(x))) {
    if (NROW(colnames(x)) != NCOL(x)) stop(paste('argument "x" must have', NCOL(x), 'names', sep = " "))
  }
  if (!is.null(var.names)) {
    if (NROW(var.names) != NCOL(x)) stop(paste('argument "var.names" must contain', NCOL(x), 'names', sep = " "))
  }
  if (scale.font <= 0) stop('argument "scale.font" must be positive')

  nv <- NCOL(x)

  if (is.null(var.names)) {
    if (!is.null(colnames(x))) {
      my.names <- colnames(x)
      label.xi <- FALSE
    }
    if (is.null(colnames(x))) {
      my.names <- paste0("x", 1:nv)
      label.xi <- TRUE
    }
  }
  if (!is.null(var.names)) {
    my.names <- var.names
    label.xi <- FALSE
  }

  x <- as.matrix(x)

  n <- NROW(x)

  x.tilde <- x - t(colMeans(x) * matrix(1, nv, n))

  pv <- matrix(0, nv, nv)
  for (ii in 1:(nv - 1)) {
    for (jj in (ii + 1):nv) {
      t.tilde <- (t(x.tilde[, ii]) %*% x.tilde[, jj]) / sqrt(t(x.tilde[, ii] ^ 2) %*% (x.tilde[, jj] ^ 2))
      pv[ii, jj] <- 2 * (1 - pnorm(abs(t.tilde)))
    }
  }

  pv <- pv + t(pv)

  diag(pv) <- NA

  pc <- cor(x)

  if (plot == TRUE) {
    plotcorrmat(pc, pv, my.names, scale.font, label.xi)
  }

  if (table == TRUE) {
    options(max.print = 1000000)

  df.pc <- as.data.frame(cbind(my.names, format(round(pc, 3), nsmall = 3)), stringsAsFactors = FALSE)

  for (ii in 1:nv) {
    df.pc[ii, (ii + 1)] <- 1
  }

  colnames(df.pc) <- c(" ", my.names)

  df.pv <- as.data.frame(cbind(my.names, format(round(pv, 3), nsmall = 3)), stringsAsFactors = FALSE)

  for (ii in 1:nv) {
    df.pv[ii, ii + 1] <- " "
  }

  colnames(df.pv) <- c(" ", my.names)

  cat("\n")
  cat("Matrix of Pearson correlations\n")
  cat("\n")
  print(df.pc, row.names = FALSE)
  cat("\n")
  cat("Matrix of p-values\n")
  cat("\n")
  print(df.pv, row.names = FALSE)
  cat("\n")
  }

  colnames(pc) <- NULL
  rownames(pc) <- NULL
  colnames(pv) <- NULL
  rownames(pv) <- NULL

  invisible(structure(list(pc = as.matrix(pc), pv = as.matrix(pv)), class = "rcorr.test"))

}

########################################################################################################################
### PLOT STATISTICS
########################################################################################################################

plotcorr <- function(max.lag, seq.max.lag, ac.cc, s.cb, r.cb, alpha, n, uni.biv, my.names, scale.font) {

  options(scipen = 999)

  results <- cbind(seq.max.lag, ac.cc, s.cb, -s.cb, r.cb, -r.cb)

  if (uni.biv == 1) {
    colnames(results) <- c("lag", "AC", "scb", "scbl", "rcb", "rcbl")
  }
  if (uni.biv == 2) {
    colnames(results) <- c("lag", "CC", "scb", "scbl", "rcb", "rcbl")
  }
  df.results <- data.frame(results)

  my.ticks.x <- seq(seq.max.lag[1], max.lag, ceiling(0.1 * max.lag))

  max.y.tick <- max(pretty_breaks()(c(0, results[, 2:6])))
  min.y.tick <- min(pretty_breaks()(c(0, results[, 2:6])))
  my.ticks.y <- pretty_breaks()(c(0, results[, 2:6]))

  label.scb <- paste("Standard CB(", 100 * (1 - alpha), "%)", sep = "")
  label.rcb <- paste("Robust CB(", 100 * (1 - alpha), "%)", sep = "")

  if (scale.font != 1) {
    scale.fontaxis <- 0.9 * scale.font
  } else {
    scale.fontaxis <- scale.font
  }

  scale.key <- 15 + 5 * scale.font

  g.corr <- ggplot(df.results, aes_(x = ~lag)) +
    theme_classic() +
    theme(plot.title = element_text(size = 13 * scale.font), legend.text = element_text(size = 12 * scale.font), axis.text = element_text(size = 10 * scale.fontaxis), axis.title = element_text(size = 10 * scale.font)) +
    theme(legend.position = "top", legend.key.size = unit(scale.key, "pt"), legend.text = element_text(margin = margin(r = 10, unit = "pt"))) +
    theme(axis.title.x = element_text(margin = margin(t = 10))) +
    theme(axis.text.x = element_text(margin = margin(t = 4))) +
    guides(fill = guide_legend(order = 1)) +
    theme(axis.title.y = element_blank()) +
    labs(x = "Lag") +
    scale_x_continuous(breaks = my.ticks.x) +
    scale_y_continuous(breaks = my.ticks.y, limits = c(min.y.tick, max.y.tick), expand = c(0, 0)) +
    geom_blank(aes(y = max.y.tick)) +
    geom_blank(aes(y = min.y.tick)) +
    geom_hline(yintercept = 0)

  if (uni.biv == 1) {
    g.corr <- g.corr + geom_col(aes_(y = ~AC, fill = "AC"), width = 0.2) + scale_fill_manual("", breaks = "AC", values = "#4F91CD")
  }
  if (uni.biv == 2) {
    g.corr <- g.corr + geom_col(aes_(y = ~CC, fill = "CC"), width = 0.2) + scale_fill_manual("", breaks = "CC", values = "#4F91CD")
  }

  g.corr <- g.corr + geom_line(aes_(y = ~scb, colour = "scb"), linetype = "dashed", size = 1) +
    geom_line(aes_(y = ~scbl), colour = "gray50", linetype = "dashed", size = 1) +
    geom_line(aes_(y = ~rcb, colour = "rcb"), linetype = "dashed", size = 1) +
    geom_line(aes_(y = ~rcbl), colour = "red", linetype = "dashed", size = 1) +
    scale_colour_manual("", breaks = c("scb", "rcb"), values = c("red", "gray50"), labels = c(label.scb, label.rcb))

  if (uni.biv == 1) {
    g.corr <- g.corr + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(bquote(Autocorrelation~of~.(my.names)[t]))
  }

  if (uni.biv == 2) {
    g.corr <- g.corr + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(bquote(Cross-correlation~of~.(my.names[1])[t]~and~.(my.names[2])[t-k]))
  }

  show(g.corr)

}

plotstat <- function(max.lag, seq.max.lag, stat1.val, stat2.val, alpha, stat.cv, n, uni.biv, my.names, stat1.name, stat2.name, stat.id, scale.font) {

  options(scipen = 999)

  results <- cbind(seq.max.lag, stat1.val, stat2.val, stat.cv)
  colnames(results) <- c("lag", "stat1", "stat2", "cv")
  df.results <- data.frame(results)

  my.ticks <- seq(seq.max.lag[1], max.lag, ceiling(0.1 * max.lag))

  max.y.tick <- max(pretty_breaks()(c(0, results[, 2:4])))
  min.y.tick <- min(pretty_breaks()(c(0, results[, 2:4])))
  my.ticks.y <- pretty_breaks()(c(0, results[, 2:4]))

  label.stat1 <- stat1.name
  label.stat2 <- stat2.name
  label.cv <- paste("cv(", 100 * alpha, "%)", sep = "")

  if (scale.font != 1) {
    scale.fontaxis <- 0.9 * scale.font
  } else {
    scale.fontaxis <- scale.font
  }

  scale.key <- 15 + 5 * scale.font

  g.stat <- ggplot(df.results, aes_(x = ~lag)) +
    theme_classic() +
    theme(plot.title = element_text(size = 13 * scale.font), legend.text = element_text(size = 12 * scale.font), axis.text = element_text(size = 10 * scale.fontaxis), axis.title = element_text(size = 10 * scale.font)) +
    theme(legend.position = "top", legend.key.size = unit(scale.key, "pt"), legend.text = element_text(margin = margin(r = 10, unit = "pt"))) +
    theme(axis.title.x = element_text(margin = margin(t = 10))) +
    theme(axis.text.x = element_text(margin = margin(t = 4))) +
    guides(fill = guide_legend(order = 1)) +
    theme(axis.title.y = element_blank()) +
    labs(x = "Lag") +
    scale_x_continuous(breaks = my.ticks) +
    scale_y_continuous(breaks = my.ticks.y, limits = c(min.y.tick, max.y.tick), expand = c(0, 0)) +
    geom_blank(aes(y = max.y.tick)) +
    geom_blank(aes(y = min.y.tick)) +
    geom_line(aes_(y = ~stat1, colour = "stat1", linetype = "stat1"), size = 1, na.rm = TRUE) +
    geom_line(aes_(y = ~stat2, colour = "stat2", linetype = "stat2"), size = 1, na.rm = TRUE) +
    geom_line(aes_(y = ~cv, colour = "cv", linetype = "cv"), size = 1, na.rm = TRUE)

  if (stat.id == 1) {
    g.stat <- g.stat + scale_colour_manual("", breaks = c("stat1", "stat2", "cv"), values = c("stat1" = "gray50", "stat2" = "red", "cv" = "black"), labels = c(label.stat1, label.stat2, label.cv)) +
    scale_linetype_manual("", breaks = c("stat1", "stat2", "cv"), values = c("stat1" = "solid", "stat2" = "solid", "cv" = "dashed"), labels = c(label.stat1, label.stat2, label.cv))
  }

  if (stat.id == 2 | stat.id == 3) {
    g.stat <- g.stat + scale_colour_manual("", breaks = c("stat1", "stat2", "cv"), values = c("stat1" = "red", "stat2" = "gray50", "cv" = "black"), labels = c(label.stat1, label.stat2, label.cv)) +
    scale_linetype_manual("", breaks = c("stat1", "stat2", "cv"), values = c("stat1" = "solid", "stat2" = "solid", "cv" = "dashed"), labels = c(label.stat1, label.stat2, label.cv))
  }

  if (uni.biv == 1 && stat.id == 1) {
    g.stat <- g.stat + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(bquote(Cumulative~tests~"for"~zero~autocorrelation~of~.(my.names)[t]))
  }

  if (uni.biv == 1 && stat.id == 2) {
    g.stat <- g.stat + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(bquote(Tests~"for"~iid~property~of~.(my.names)[t]))
  }

  if (uni.biv == 1 && stat.id == 3) {
    g.stat <- g.stat + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(bquote(Cumulative~tests~"for"~iid~property~of~.(my.names)[t]))
  }

  if (uni.biv == 2 && stat.id == 1) {
    g.stat <- g.stat + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(bquote(Cumulative~tests~"for"~zero~cross-correlation~of~.(my.names[1])[t]~and~.(my.names[2])[t-k]))
  }

  show(g.stat)

}

plotcorrmat <- function(cc, pv, my.names, scale.font, label.xi) {

  nv <- NCOL(cc)

  colnames(pv) <- my.names
  rownames(pv) <- my.names
  melted_pv <- melt(pv)
  colnames(melted_pv) <- c("Var1", "Var2", "pv")

  colnames(cc) <- my.names
  rownames(cc) <- my.names
  melted_cc <- melt(cc)
  colnames(melted_cc) <- c("Var1", "Var2", "cc")

  ccpv <- cbind(melted_cc, melted_pv[, 3])
  colnames(ccpv) <- c("Var1", "Var2", "cc", "pv")

  labelcc <- format(round(ccpv[, 3], 3), nsmall = 3)
  labelcc[is.na(ccpv[, 4])] <- 1
  labelpv <- paste("(", format(round(ccpv[, 4], 3), nsmall = 3), ")", sep = "")
  labelpv[is.na(ccpv[, 4])] <- ""
  mylabel <- paste(labelcc, "\n", labelpv, sep = "")

  if (label.xi == TRUE) {
    labelticks <- matrix(NA, nv, 1)
      for (ii in 1:nv) {
        labelticks[ii] <- as.expression(bquote(x[.(ii)]))
      }
  } else {
    labelticks <- my.names
  }

  ccpv[is.na(ccpv)] <- 1.3

  mycolors <- c("#FF0000", "#FF3F3F", "#FF7F7F", "#FFBFBF", "#FFFFFF", "#BEBEBE")

  colcode <- cut(ccpv[, 4], breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), include.lowest = TRUE)
  brk <- levels(colcode)

  colcode2 <- cut(ccpv[, 4], breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1, 1.5), include.lowest = TRUE)

  mylabelspv <- c("[0,0.001]", "(0.001,0.01]", "(0.01,0.05]", "(0.05,0.1]", "(0.1,1]")

  scale.key <- 5 + 10 * scale.font

  g <- ggplot(ccpv, aes_(~Var1, ~rev(Var2))) +
    geom_tile(aes(fill = colcode2), colour = "black") +
    scale_fill_manual(name = "p-value", breaks = brk, values = mycolors, labels = mylabelspv, drop = FALSE) +
    scale_x_discrete(position = "top", labels = labelticks) +
    scale_y_discrete(labels = rev(labelticks)) +
    theme(plot.title = element_text(size = 13 * scale.font), legend.title = element_text(size = 12 * scale.font), legend.text = element_text(size = 11 * scale.font), axis.text = element_text(size = 11 * scale.font)) +
    theme(legend.key.size = unit(scale.key, "pt")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(panel.background = element_blank(), axis.ticks = element_blank()) +
    theme(axis.ticks.length = unit(0, "pt")) +
    theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) +
    theme(plot.title = element_text(margin = margin(b = 20))) +
    geom_text(aes_(~Var1, ~rev(Var2), label = mylabel), color = "black", size = 3.5 * scale.font) +
    ggtitle("Pearson correlations and p-values for significance")

  show(g)

}
