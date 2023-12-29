#
#
# (1) Calculations for woody plants example
#
#

library("metasens")

data(woodyplants)
m.woody <- metacont(n.elev, mean.elev, sd.elev, n.amb, mean.amb, sd.amb,
  data = woodyplants, sm = "ROM")

# LFK index test
#
lfkindex(m.woody)

# Two-sided tests for funnel plot asymmetry
#
metabias(m.woody)
metabias(m.woody, method = "Thompson")
metabias(m.woody, method = "Begg")

# One-sided tests
#
pt(metabias(m.woody)$statistic, df = m.woody$k - 2, lower = FALSE)
pt(metabias(m.woody, method = "Thompson")$statistic,
  df = m.woody$k - 2, lower = FALSE)
pnorm(metabias(m.woody, method = "Begg")$statistic, lower = FALSE)


#
#
# (2) Use of REML or DerSimonian-Laird estimator
#
#

#
# (2a) Simulations based on sample size distribution by Schwarzer et al. (2002)
#

load("type1_Schwarzer2002.rda")
load("power_Schwarzer2002_rho0.5.rda")
load("power_Schwarzer2002_rho0.9.rda")

#
# Use of REML or DerSimonian-Laird estimator in meta-analysis
#

with(type1_Schwarzer2002, 100 * table(method.tau) / length(method.tau))
with(power_Schwarzer2002_rho0.5, 100 * table(method.tau) / length(method.tau))
with(power_Schwarzer2002_rho0.9, 100 * table(method.tau) / length(method.tau))

#
# Use of REML or DerSimonian-Laird estimator in Thompson-Sharp test
#

with(type1_Schwarzer2002,
     100 * table(method.tau.thompson) / length(method.tau))
with(power_Schwarzer2002_rho0.5,
     100 * table(method.tau.thompson) / length(method.tau))
with(power_Schwarzer2002_rho0.9,
     100 * table(method.tau.thompson) / length(method.tau))

#
# Cross-tabulation (REML and DerSimonian-Laird estimator)
#

with(type1_Schwarzer2002, table(method.tau, method.tau.thompson))
with(power_Schwarzer2002_rho0.5, table(method.tau, method.tau.thompson))
with(power_Schwarzer2002_rho0.9, table(method.tau, method.tau.thompson))

#
# Summary of false positive rates (type I error)
#

nsim <- attr(type1_Schwarzer2002, "nsim")
#
round(100 * with(type1_Schwarzer2002,
                 summary(by(lfk.sign, ia, sum))) / nsim, 1)

#
# Ten smallest numbers of significant results for LFK index
#

with(type1_Schwarzer2002, table(by(lfk.sign, ia, sum)))[1:10]

#
# (2b) Simulations based on uniform distribution
#

load("type1_unif.rda")
load("power_unif_rho0.5.rda")
load("power_unif_rho0.9.rda")

#
# Use of REML or DerSimonian-Laird estimator in meta-analysis
#

with(type1_unif, 100 * table(method.tau) / length(method.tau))
with(power_unif_rho0.5, 100 * table(method.tau) / length(method.tau))
with(power_unif_rho0.9, 100 * table(method.tau) / length(method.tau))

#
# Use of REML or DerSimonian-Laird estimator in Thompson-Sharp test
#

with(type1_unif,
     100 * table(method.tau.thompson) / length(method.tau))
with(power_unif_rho0.5,
     100 * table(method.tau.thompson) / length(method.tau))
with(power_unif_rho0.9,
     100 * table(method.tau.thompson) / length(method.tau))

#
# Cross-tabulation (REML and DerSimonian-Laird estimator)
#

with(type1_unif, table(method.tau, method.tau.thompson))
with(power_unif_rho0.5, table(method.tau, method.tau.thompson))
with(power_unif_rho0.9, table(method.tau, method.tau.thompson))

#
# Summary of false positive rates (type I error)
#

nsim <- attr(type1_unif, "nsim")
#
round(100 * with(type1_unif, summary(by(lfk.sign, ia, sum))) / nsim, 1)

#
# Ten smallest numbers of significant results for LFK index
#

with(type1_unif, table(by(lfk.sign, ia, sum)))[1:10]


#
#
# (3) I2 statistic (median of medians)
#
#

#
# (3a) Simulations based on sample size distribution by Schwarzer et al. (2002)
#

meds <- vector("numeric", 0)
#
for (i in c(0, 0.25, 0.5))
  meds <-
    c(meds,
      median(with(subset(type1_Schwarzer2002, tau == i),
                  by(I2, list(N_studies, sd_treat, theta_treat), median))),
      median(with(subset(power_Schwarzer2002_rho0.5, tau == i),
                  by(I2, list(N_studies, sd_treat, theta_treat), median))),
      median(with(subset(power_Schwarzer2002_rho0.9, tau == i),
                  by(I2, list(N_studies, sd_treat, theta_treat), median))))
#
data.frame(tau = rep(c(0, 0.25, 0.5), rep(3, 3)),
           rho = rep(c(0, -0.5, -0.9), 3),
           median.I2 = paste0(round(100 * meds, 1), "%"))

#
# (3b) Simulations based on uniform distribution
#

meds <- vector("numeric", 0)
#
for (i in c(0, 0.25, 0.5))
  meds <-
  c(meds,
    median(with(subset(type1_unif, tau == i),
                by(I2, list(N_studies, sd_treat, theta_treat), median))),
    median(with(subset(power_unif_rho0.5, tau == i),
                by(I2, list(N_studies, sd_treat, theta_treat), median))),
    median(with(subset(power_unif_rho0.9, tau == i),
                by(I2, list(N_studies, sd_treat, theta_treat), median))))
#
data.frame(tau = rep(c(0, 0.25, 0.5), rep(3, 3)),
           rho = rep(c(0, -0.5, -0.9), 3),
           median.I2 = paste0(round(100 * meds, 1), "%"))


#
#
# (4) 
#
#

#
# (4a) Simulations based on sample size distribution by Schwarzer et al. (2002)
#

with(subset(type1_Schwarzer2002, tau == 0),
     by(lfkindex, N_studies, quantile, probs = c(0.05, 0.95)))
#
with(subset(type1_Schwarzer2002, tau == 0.25),
     by(lfkindex, N_studies, quantile, probs = c(0.05, 0.95)))
#
with(subset(type1_Schwarzer2002, tau == 0.5),
     by(lfkindex, N_studies, quantile, probs = c(0.05, 0.95)))

#
# (4b) Simulations based on uniform distribution
#

with(subset(type1_unif, tau == 0),
     by(lfkindex, N_studies, quantile, probs = c(0.05, 0.95)))
#
with(subset(type1_unif, tau == 0.25),
     by(lfkindex, N_studies, quantile, probs = c(0.05, 0.95)))
#
with(subset(type1_unif, tau == 0.5),
     by(lfkindex, N_studies, quantile, probs = c(0.05, 0.95)))
