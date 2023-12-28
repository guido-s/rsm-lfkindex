#point_lines <- function(x, test, xvals, nsim, pch, col) {
#  j <- 0
#  for (i in c(10, 20, 50, 100)) {
#    n.i <- length(xvals) / 4
#    dat.i <- subset(x, x$N_studies == i)
#    xvals.i <- xvals[j * n.i + seq_len(n.i)]
#    #
#    points(xvals.i, by(x[[test]], x$ia, sum)[xvals.i] / nsim, pch = pch,
#           type = "b", col = col)
#    #
#    j <- j + 1
#  }
#  #
#  invisible(NULL)
#}

point_lines <- function(x, test, xvals, nsim, pch, col) {
  k <- 0
  for (i in c(10, 20, 50, 100)) {
    for (j in c(0, 0.25, 0.5)) {
      n.ij <- length(xvals) / 12
      dat.ij <- subset(x, x$N_studies == i & x$tau == j)
      xvals.ij <- xvals[k * n.ij + seq_len(n.ij)]
      #
      points(xvals.ij, by(x[[test]], x$ia, sum)[xvals.ij] / nsim, pch = pch,
           type = "b", col = col, lwd = 1.5)
      #
      k <- k + 1
    }
  }
  #
  invisible(NULL)
}

nestedloop <- function(x,
                       ylim = c(0, 1),
                       createpdf = FALSE,
                       xlab = "Scenarios", ylab = "Proportion",
                       legend = "topleft",
                       main = NULL, sub = FALSE,
                       egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
                       type = "unadjusted",
                       width = 10, height = 8,
                       supplement = FALSE,
                       filename = NULL) {
  
  n.levs <- length(levels(x$ia))
  nsim <- attributes(x)$nsim
  #
  lev1 <- unique(x$N_studies)
  lev2 <- unique(x$tau)
  lev3 <- unique(x$theta_treat)
  lev4 <- unique(x$sd_treat)
  #
  n1 <- length(lev1)
  n2 <- length(lev2)
  n3 <- length(lev3)
  n4 <- length(lev4)
  #
  xvals <- seq_along(levels(x$ia))
  #
  rname <- deparse(substitute(x))
  #
  if (grepl("Schwarzer2002", rname))
    scenario <- "Schwarzer2002"
  else if (grepl("unif", rname))
    scenario <- "unif"
  else
    stop("Scenario unclear (name of argument 'x' must contain either ",
         "'Schwarzer2002' or 'unif'.")
  #
  if (grepl("type1", rname))
    absrho <- 0
  else if (grepl("rho0.5", rname))
    absrho <- 0.5
  else if (grepl("rho0.9", rname))
    absrho <- 0.9
  else
    stop("Unclear whether type-1 error or power provided.")
  #
  if (type == "inflated")
    x$lfk.sign.infl <- x$lfk.sign
  #
  if (absrho != 0 & type == "inflated") {
    text.type <- "_infl"
    ext <- ".infl"
  }
  else if (absrho != 0 & type == "adjusted") {
    text.type <- "_adj"
    ext <- ".adj"
  }
  else {
    text.type <- NULL
    ext <- NULL
  }
  #
  if (sub) {
    sub <- if (scenario == "unif")
      paste("Group sample sizes between 51 and 500",
            "drawn from uniform distribution")
    else
      "Group sample sizes according to Schwarzer et al. (2002)"
  }
  else
    sub <- NULL
  #
  dir <- if (supplement) "supplement/" else "graphics/"
  if (is.null(filename))
    filename <- paste0(dir,
                       if (absrho == 0) "type1_" else "power_",
                       scenario,
                       if (absrho != 0) paste0("_", absrho),
                       text.type,
                       ".pdf")
  else
    filename <- paste0(dir, filename)
  #
  if (createpdf)
    pdf(filename, width = width, height = height)
  #
  plot(xvals, by(x$lfk.sign, x$ia, sum) / nsim,
       axes = FALSE, ylim = ylim,
       xlab = "Scenarios", ylab = "Proportion",
       type = "n")
  #
  if (absrho == 0)
    abline(h = seq(0, 30, by = 5) / 100, col = "gray", lwd = 1)
  else
    abline(h = seq(0, 100, by = 10) / 100, col = "gray", lwd = 1)
  #
  if (absrho == 0)
    abline(h = 0.05, col = "black", lwd = 3, lty = 2)
  #
  abline(v = c(0, n1 * n2 * n3 * n4) + 0.5, col = "gray", lwd = 3)
  abline(v = seq_len(n1 - 1) * n.levs / n1 + 0.5, col = "gray", lwd = 3)
  abline(v = seq_len(n1 * n2 - 1) * n.levs / (n1 * n2) + 0.5,
         col = "gray", lwd = 1.5)
  abline(v = seq_len(n1 * n2 * n3 - 1) * n.levs / (n1 * n2 * n3) + 0.5,
         col = "gray", lwd = 1, lty = 2)
  #
  if (lfk)
    point_lines(x, paste0("lfk.sign", ext), xvals, nsim, 16, "black")
  if (egger)
    point_lines(x, paste0("egger.sign", ext), xvals, nsim, 17, "blue")
  if (rank)
    point_lines(x, paste0("rank.sign", ext), xvals, nsim, 4, "red")
  if (thompson)
    point_lines(x, paste0("thompson.sign", ext), xvals, nsim, 18, "green")
  #
  title(main = main, sub = sub)
  #
  axis(3, tick = FALSE, line = -1,
       at = 0:(n1 - 1) * (n1 * n2 * n3 * n4) / n1 + (n2 * n3 * n4) / 2,
       label = paste0("k=", lev1))
  lab2 <- rep(lev2, n1)
  lab2[1] <- paste0("tau=", lab2[1])
  axis(1, tick = FALSE, line = -1,
       at = 0:(n1 * n2 - 1) * (n1 * n2 * n3 * n4) / (n1 * n2) + (n3 * n4) / 2,
       label = lab2)
  #
  axis(2)
  #
  legend(legend, col = c("black", "blue", "red", "green"),
         bg = "lightgray",
         pch = c(16, 17, 4, 18), lty = rep(1, 4),
         lwd = 1.5,
         legend = c("LFK index test",
                    "Egger test",
                    "Rank test",
                    "Thompson-Sharp test"))
  #
  box()
  #
  if (createpdf)
    invisible(dev.off())
  #
  invisible(NULL)
}


ecdfs <- function(x, test, text = "", sub = "",
                  xlim = c(-2.5, 2.5),
                  step = 1) {
  var <- deparse(substitute(test))
  #
  x10 <- sort(x[[var]][x[["N_studies"]] == 10])
  x20 <- sort(x[[var]][x[["N_studies"]] == 20])
  x50 <- sort(x[[var]][x[["N_studies"]] == 50])
  x100 <- sort(x[[var]][x[["N_studies"]] == 100])
  #
  x10 <- x10[seq.int(1, length(x10), step)]
  x20 <- x20[seq.int(1, length(x20), step)]
  x50 <- x50[seq.int(1, length(x50), step)]
  x100 <- x100[seq.int(1, length(x100), step)]
  #
  plot(ecdf(x10),
       do.points = FALSE,
       main = text, xlab = sub, xlim = xlim)
  #
  plot(ecdf(x20),
       do.points = FALSE,
       add = TRUE, col = "blue")
  #
  plot(ecdf(x50),
       do.points = FALSE,
       add = TRUE, col = "red")
  #
  plot(ecdf(x100),
       do.points = FALSE,
       add = TRUE, col = "green")
  #
  abline(h = 0.5, v = 0, col = "gray")
  #
  legend("topleft", col = c("black", "blue", "red", "green"),
         lty = rep(1, 4),
         legend = c("k = 10",
                    "k = 20",
                    "k = 50",
                    "k = 100"),
         bg = "lightgray")
  #
  invisible(NULL)
}


ecdfs_tau <- function(x, test, text = "", sub = "",
                      xlim = c(-2.5, 2.5),
                      step = 1) {
  var <- deparse(substitute(test))
  #
  tau0 <- sort(x[[var]][x[["tau"]] == 0])
  tau1 <- sort(x[[var]][x[["tau"]] == 0.25])
  tau2 <- sort(x[[var]][x[["tau"]] == 0.5])
  #
  tau0 <- tau0[seq.int(1, length(tau0), step)]
  tau1 <- tau1[seq.int(1, length(tau1), step)]
  tau2 <- tau2[seq.int(1, length(tau2), step)]
  #
  plot(ecdf(tau0),
       do.points = FALSE,
       main = text, xlab = sub, xlim = xlim)
  #
  plot(ecdf(tau1),
       do.points = FALSE,
       add = TRUE, col = "blue")
  #
  plot(ecdf(tau2),
       do.points = FALSE,
       add = TRUE, col = "red")
  #
  abline(h = 0.5, v = 0, col = "gray")
  #
  legend("topleft", col = c("black", "blue", "red"),
         lty = rep(1, 3),
         legend = c("tau = 0",
                    "tau = 0.25",
                    "tau = 0.5"),
         bg = "lightgray")
  #
  invisible(NULL)
}


ecdfs_distr <- function(x, test, text = "", sub = "",
                      xlim = c(-2.5, 2.5),
                      step = 1) {
  var <- deparse(substitute(test))
  #
  black <- sort(x[[var]][x[["distr"]] == "Schwarzer2002"])
  unif <- sort(x[[var]][x[["distr"]] == "unif"])
  #
  black <- black[seq.int(1, length(black), step)]
  unif <- unif[seq.int(1, length(unif), step)]
  #
  plot(ecdf(black),
       do.points = FALSE,
       main = text, xlab = sub, xlim = xlim, col = "blue")
  #
  plot(ecdf(unif),
       do.points = FALSE,
       add = TRUE, col = "red")
  #
  abline(h = 0.5, v = 0, col = "gray")
  #
  legend("topleft", col = c("blue", "red"),
         lty = rep(1, 2),
         legend = c("small", "large"),
         bg = "lightgray")
  #
  invisible(NULL)
}


nestedloop_I2 <- function(x,
                          ylim = c(0, 1),
                          createpdf = FALSE,
                          xlab = "Scenarios", ylab = "Proportion",
                          main = NULL, sub = FALSE,
                          width = 10, height = 8,
                          supplement = FALSE,
                          filename = NULL) {
  
  n.levs <- length(levels(x$ia))
  nsim <- attributes(x)$nsim
  #
  lev1 <- unique(x$N_studies)
  lev2 <- unique(x$tau)
  lev3 <- unique(x$theta_treat)
  lev4 <- unique(x$sd_treat)
  #
  n1 <- length(lev1)
  n2 <- length(lev2)
  n3 <- length(lev3)
  n4 <- length(lev4)
  #
  xvals <- seq_along(levels(x$ia))
  #
  rname <- deparse(substitute(x))
  #
  if (grepl("Schwarzer2002", rname))
    scenario <- "Schwarzer2002"
  else if (grepl("unif", rname))
    scenario <- "unif"
  else
    stop("Scenario unclear (name of argument 'x' must contain either ",
         "'Schwarzer2002' or 'unif'.")
  #
  if (grepl("type1", rname))
    absrho <- "0.0"
  else if (grepl("rho0.5", rname))
    absrho <- "0.5"
  else if (grepl("rho0.9", rname))
    absrho <- "0.9"
  else
    stop("Unclear whether type-1 error or power provided.")
  #
  if (sub) {
    sub <- if (scenario == "unif")
      paste("Group sample sizes between 51 and 500",
            "drawn from uniform distribution")
    else
      "Group sample sizes according to Schwarzer et al. (2002)"
  }
  else
    sub <- NULL
  #
  dir <- if (supplement) "supplement/" else "graphics/"
  if (is.null(filename))
    filename <- paste0(dir,
                       "I2_",
                       scenario,
                       paste0("_rho", absrho),
                       ".pdf")
  else
    filename <- paste0(dir, filename)
  #
  if (createpdf)
    pdf(filename, width = width, height = height)
  #
  plot(xvals, by(x$I2, x$ia, median),
       axes = FALSE, ylim = ylim,
       xlab = "Scenarios", ylab = "Proportion",
       type = "n")
  #
  #
  abline(v = c(0, n1 * n2 * n3 * n4) + 0.5, col = "gray", lwd = 3)
  abline(v = seq_len(n1 - 1) * n.levs / n1 + 0.5, col = "gray", lwd = 3)
  abline(v = seq_len(n1 * n2 - 1) * n.levs / (n1 * n2) + 0.5,
         col = "gray", lwd = 1.5)
  abline(v = seq_len(n1 * n2 * n3 - 1) * n.levs / (n1 * n2 * n3) + 0.5,
         col = "gray", lwd = 1, lty = 2)
  #
  #points(xvals, by(x$I2, x$ia, median), pch = 16, type = "b", col = "black")
  k <- 0
  for (i in c(10, 20, 50, 100)) {
    for (j in c(0, 0.25, 0.5)) {
      n.ij <- length(xvals) / 12
      dat.ij <- subset(x, x$N_studies == i & x$tau == j)
      xvals.ij <- xvals[k * n.ij + seq_len(n.ij)]
      #
      points(xvals.ij, by(x$I2, x$ia, median)[xvals.ij], pch = 16,
             type = "b", col = "black", lwd = 1.5)
      #
      k <- k + 1
    }
  }
  #
  title(main = main, sub = sub)
  #
  axis(3, tick = FALSE, line = -1,
       at = 0:(n1 - 1) * (n1 * n2 * n3 * n4) / n1 + (n2 * n3 * n4) / 2,
       label = paste0("k=", lev1))
  lab2 <- rep(lev2, n1)
  lab2[1] <- paste0("tau=", lab2[1])
  axis(1, tick = FALSE, line = -1,
       at = 0:(n1 * n2 - 1) * (n1 * n2 * n3 * n4) / (n1 * n2) + (n3 * n4) / 2,
       label = lab2)
  #
  axis(2)
  #
  box()
  #
  if (createpdf)
    invisible(dev.off())
  #
  invisible(NULL)
}
