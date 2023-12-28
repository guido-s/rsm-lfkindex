#
# Figure 1
#

library("metasens")
data(woodyplants)
m.woody <- metacont(n.elev, mean.elev, sd.elev, n.amb, mean.amb, sd.amb,
  data = woodyplants, sm = "ROM")
#
pdf("graphics/Figure1.pdf", width = 10, height = 6)
#
par(mfrow = c(1, 2))
funnel(m.woody, xlim = c(0.5, 4), axes = FALSE)
axis(1, at = c(0.5, 1, 2, 4))
axis(2)
box()
#
doiplot(m.woody, lfkindex = FALSE, axes = FALSE,
  xlab = "Ratio of Means", xlim = c(log(0.5), log(4)))
axis(1, at = log(c(0.5, 1, 2, 4)), labels = c(0.5, 1, 2, 4))
axis(2)
box()
#
dev.off()


#
#
# Simulations based on sample size distribution by Schwarzer et al. (2002)
#
#

source("funcs_plot.R")

load("type1_power_Schwarzer2002.rda")

#
# Figures 2, 4 and 5
#

nestedloop(type1_Schwarzer2002, createpdf = TRUE,
           ylim = c(0, 0.3),
           legend = "topright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           filename = "Figure2.pdf")

nestedloop(power_Schwarzer2002_rho0.5, createpdf = TRUE,
           legend = "bottomright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           filename = "Figure4.pdf")

nestedloop(power_Schwarzer2002_rho0.5, createpdf = TRUE,
           legend = "bottomright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           type = "inflated",
           filename = "Figure5.pdf")

#
# Supplement
#

nestedloop(power_Schwarzer2002_rho0.9, createpdf = TRUE,
           legend = "bottomright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           supplement = TRUE)

nestedloop(power_Schwarzer2002_rho0.9, createpdf = TRUE,
           legend = "bottomright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           type = "inflated",
           supplement = TRUE)

nestedloop(power_Schwarzer2002_rho0.5, createpdf = TRUE,
           legend = "bottomright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           type = "adjusted",
           supplement = TRUE)

nestedloop(power_Schwarzer2002_rho0.9, createpdf = TRUE,
           legend = "bottomright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           type = "adjusted",
           supplement = TRUE)

d0_Schwarzer2002 <- subset(type1_Schwarzer2002, tau == 0)
d1_Schwarzer2002 <- subset(type1_Schwarzer2002, tau == 0.25)
d2_Schwarzer2002 <- subset(type1_Schwarzer2002, tau == 0.5)

pdf("supplement/ecdf_Schwarzer2002.pdf", width = 10, height = 8)
#
par(mfrow = c(2, 2))
#
ecdfs(d0_Schwarzer2002, lfkindex, text = "LFK index", step = 50,
      xlim = c(-1.0, 1.0))
#
ecdfs(d0_Schwarzer2002, stat.egger, text = "Egger test", step = 50,
      xlim = c(qt(0.05, df = 8), qt(0.95, df = 8)))
#
ecdfs(d0_Schwarzer2002, stat.rank, text = "Rank test", step = 50,
      xlim = c(qnorm(0.05), qnorm(0.95)))
#
ecdfs(d0_Schwarzer2002, stat.thompson, text = "Thompson-Sharp test", step = 50,
      xlim = c(qt(0.05, df = 8), qt(0.95, df = 8)))
#
dev.off()

nestedloop_I2(type1_Schwarzer2002, createpdf = TRUE, supplement = TRUE)
nestedloop_I2(power_Schwarzer2002_rho0.5, createpdf = TRUE, supplement = TRUE)
nestedloop_I2(power_Schwarzer2002_rho0.9, createpdf = TRUE, supplement = TRUE)


#
#
# Simulations based on uniform distribution
#
#

load("type1_power_unif.rda")

#
# Figure 3
#

nestedloop(type1_unif, createpdf = TRUE,
           ylim = c(0, 0.3),
           legend = "topright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           filename = "Figure3.pdf")

#
# Figure 6
#

d0_unif <- subset(type1_unif, tau == 0)
d1_unif <- subset(type1_unif, tau == 0.25)
d2_unif <- subset(type1_unif, tau == 0.5)

pdf("graphics/Figure6.pdf", width = 10, height = 8)
#
par(mfrow = c(2, 2))
#
ecdfs(d0_unif, lfkindex, text = "LFK index", step = 50,
      xlim = c(-1.0, 1.0))
#
ecdfs(d0_unif, stat.egger, text = "Egger test", step = 50,
      xlim = c(qt(0.05, df = 8), qt(0.95, df = 8)))
#
ecdfs(d0_unif, stat.rank, text = "Rank test", step = 50,
      xlim = c(qnorm(0.05), qnorm(0.95)))
#
ecdfs(d0_unif, stat.thompson, text = "Thompson-Sharp test", step = 50,
      xlim = c(qt(0.05, df = 8), qt(0.95, df = 8)))
#
dev.off()

#
# Supplement
#

nestedloop(power_unif_rho0.5, createpdf = TRUE,
           legend = "topleft",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           supplement = TRUE)

nestedloop(power_unif_rho0.9, createpdf = TRUE,
           legend = "topleft",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           supplement = TRUE)

nestedloop(power_unif_rho0.5, createpdf = TRUE,
           legend = "topleft",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           type = "inflated",
           supplement = TRUE)

nestedloop(power_unif_rho0.9, createpdf = TRUE,
           legend = "topleft",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           type = "inflated",
           supplement = TRUE)

nestedloop(power_unif_rho0.5, createpdf = TRUE,
           legend = "bottomright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           type = "adjusted",
           supplement = TRUE)

nestedloop(power_unif_rho0.9, createpdf = TRUE,
           legend = "bottomright",
           egger = TRUE, rank = TRUE, thompson = TRUE, lfk = TRUE,
           type = "adjusted",
           supplement = TRUE)

nestedloop_I2(type1_unif, createpdf = TRUE, supplement = TRUE)
nestedloop_I2(power_unif_rho0.5, createpdf = TRUE, supplement = TRUE)
nestedloop_I2(power_unif_rho0.9, createpdf = TRUE, supplement = TRUE)
