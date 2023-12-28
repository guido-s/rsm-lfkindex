#
#
# Calculate power for sample sizes from Schwarzer et al. (2002)
#
#

library("magrittr")
library("tidyverse")

load("simu_LFK_unif.rda")
#
dat <- do.call("rbind", sim.res)
nsim <- attributes(sim.res[[1]])$nsim
#
rm(param.grid, param.names, params, sim.res)


#
#
# (1) Type I error rate (rho = 0)
#
#

d0 <- subset(dat, rho == 0)
#
d0$ia <-
  interaction(paste0("theta=", d0$theta_treat),
              paste0("sd=", d0$sd_treat),
              paste0("tau=", d0$tau),
              paste0("k=", d0$N_studies),
              sep = ",\n")
# Otherwise k=100 between k=10 and k=20 ...
d0 <- d0[order(d0$N_studies, d0$tau, d0$sd_treat, d0$theta_treat), ]
d0$ia <- factor(d0$ia, level = unique(as.character(d0$ia)))
#
d0$lfk.sign <- 1L * (d0$lfkindex < -1)
#
d0$egger.sign <- 1L * (d0$stat.egger < qt(0.05, d0$df.egger))
#
d0$rank.sign <- 1L * (d0$stat.rank < qnorm(0.05))
#
d0$thompson.sign <- 1L * (d0$stat.thompson < qt(0.05, d0$df.thompson))
#
d0 %<>% select(-(TE.common:pval.random)) %>% select(-(Q:pval.Q)) %>%
  select(-pval.egger, -pval.rank, -pval.thompson)
d0 %<>% relocate(ia, .after = last_col())
#
attr(d0, "nsim") <- nsim


#
#
# (2) Unadjusted Power (rho = -0.5)
#
#

d1 <- subset(dat, rho == -0.5)
#
d1$ia <-
  interaction(paste0("theta=", d1$theta_treat),
              paste0("sd=", d1$sd_treat),
              paste0("tau=", d1$tau),
              paste0("k=", d1$N_studies),
              sep = ",\n")
# Otherwise k=100 between k=10 and k=20 ...
d1 <- d1[order(d1$N_studies, d1$tau,
               d1$sd_treat, d1$theta_treat), ]
d1$ia <- factor(d1$ia, level = unique(as.character(d1$ia)))
#
d1$lfk.sign <- 1L * (d1$lfkindex < -1)
#
d1$egger.sign <- 1L * (d1$stat.egger < qt(0.05, d1$df.egger))
#
d1$rank.sign <- 1L * (d1$stat.rank < qnorm(0.05))
#
d1$thompson.sign <- 1L * (d1$stat.thompson < qt(0.05, d1$df.thompson))


#
#
# (3) Inflated power (rho = -0.5)
#
#

#
# (3a) Type-1 error rate for LFK index test
#

lfkindex_type1 <- tapply(d0$lfkindex < -1, d0$ia, sum) / nsim
#
lfkindex_type1 <- data.frame(ia = names(lfkindex_type1),
                             lfkindex_type1 = as.vector(lfkindex_type1))

#
# (3b) Inflated power
#

d1 <- merge(d1, lfkindex_type1, by = "ia")
#
d1$lfkindex_crit_z <- qnorm(d1$lfkindex_type1)
d1$lfkindex_crit_t <- qt(d1$lfkindex_type1, d1$N_studies - 1)
#
# Critical values based on type I error of LFK index test for k = 10
#
sel <- d1$N_studies <= 10
#
d1$egger.sign.infl <- 1L * (d1$stat.egger < qt(0.05, d1$df.egger))
d1$egger.sign.infl[sel] <- 1L * (d1$stat.egger < d1$lfkindex_crit_t)[sel]
#
d1$rank.sign.infl <- 1L * (d1$stat.rank < qnorm(0.05))
d1$rank.sign.infl[sel] <- 1L * (d1$stat.rank < d1$lfkindex_crit_z)[sel]
#
d1$thompson.sign.infl <- 1L * (d1$stat.thompson < qt(0.05, d1$df.thompson))
d1$thompson.sign.infl[sel] <- 1L * (d1$stat.thompson < d1$lfkindex_crit_t)[sel]


#
#
# (4) Adjusted power (rho = -0.5)
#
#

#
# (4a) Empirical 5% quantiles under null hypothesis
#

lfkindex_q0.05 <- by(d0$lfkindex, d0$ia, quantile, 0.05)
egger_q0.05 <- by(d0$stat.egger, d0$ia, quantile, 0.05)
rank_q0.05 <- by(d0$stat.rank, d0$ia, quantile, 0.05)
thompson_q0.05 <- by(d0$stat.thompson, d0$ia, quantile, 0.05)
#
lfkindex_q0.05 <- data.frame(ia = names(lfkindex_q0.05),
                             lfkindex_q0.05 = as.vector(lfkindex_q0.05))
#
egger_q0.05 <- data.frame(ia = names(egger_q0.05),
                          egger_q0.05 = as.vector(egger_q0.05))
#
rank_q0.05 <- data.frame(ia = names(rank_q0.05),
                         rank_q0.05 = as.vector(rank_q0.05))
#
thompson_q0.05 <- data.frame(ia = names(thompson_q0.05),
                             thompson_q0.05 = as.vector(thompson_q0.05))
#
d1 <- merge(d1, lfkindex_q0.05, by = "ia")
d1 <- merge(d1, egger_q0.05, by = "ia")
d1 <- merge(d1, rank_q0.05, by = "ia")
d1 <- merge(d1, thompson_q0.05, by = "ia")

#
# (4b) Adjusted power
#

d1$lfk.sign.adj <- 1L * (d1$lfkindex < d1$lfkindex_q0.05)
#
d1$egger.sign.adj <- 1L * (d1$stat.egger < d1$egger_q0.05)
#
d1$rank.sign.adj <- 1L * (d1$stat.rank < d1$rank_q0.05)
#
d1$thompson.sign.adj <- 1L * (d1$stat.thompson < d1$thompson_q0.05)
#
d1 %<>% select(-(TE.common:pval.random)) %>% select(-(Q:pval.Q)) %>%
  select(-pval.egger, -pval.rank, -pval.thompson)
d1 %<>% relocate(ia, .after = last_col())
#
attr(d1, "nsim") <- nsim


#
#
# (5) Unadjusted Power (rho = -0.9)
#
#

d2 <- subset(dat, rho == -0.9)
#
d2$ia <-
  interaction(paste0("theta=", d2$theta_treat),
              paste0("sd=", d2$sd_treat),
              paste0("tau=", d2$tau),
              paste0("k=", d2$N_studies),
              sep = ",\n")
# Otherwise k=100 between k=10 and k=20 ...
d2 <- d2[order(d2$N_studies, d2$tau,
               d2$sd_treat, d2$theta_treat), ]
d2$ia <- factor(d2$ia, level = unique(as.character(d2$ia)))
#
d2$lfk.sign <- 1L * (d2$lfkindex < -1)
#
d2$egger.sign <- 1L * (d2$stat.egger < qt(0.05, d2$df.egger))
#
d2$rank.sign <- 1L * (d2$stat.rank < qnorm(0.05))
#
d2$thompson.sign <- 1L * (d2$stat.thompson < qt(0.05, d2$df.thompson))


#
#
# (6) Inflated power (rho = -0.9)
#
#

d2 <- merge(d2, lfkindex_type1, by = "ia")
#
d2$lfkindex_crit_z <- qnorm(d2$lfkindex_type1)
d2$lfkindex_crit_t <- qt(d2$lfkindex_type1, d2$N_studies - 1)
#
# Critical values based on type I error of LFK index test for k = 10
#
sel <- d2$N_studies <= 10
#
d2$egger.sign.infl <- 1L * (d2$stat.egger < qt(0.05, d2$df.egger))
d2$egger.sign.infl[sel] <- 1L * (d2$stat.egger < d2$lfkindex_crit_t)[sel]
#
d2$rank.sign.infl <- 1L * (d2$stat.rank < qnorm(0.05))
d2$rank.sign.infl[sel] <- 1L * (d2$stat.rank < d2$lfkindex_crit_z)[sel]
#
d2$thompson.sign.infl <- 1L * (d2$stat.thompson < qt(0.05, d2$df.thompson))
d2$thompson.sign.infl[sel] <- 1L * (d2$stat.thompson < d2$lfkindex_crit_t)[sel]


#
#
# (7) Adjusted power (rho = -0.5)
#
#

#
# (7a) Empirical 5% quantiles under null hypothesis
#

d2 <- merge(d2, lfkindex_q0.05, by = "ia")
d2 <- merge(d2, egger_q0.05, by = "ia")
d2 <- merge(d2, rank_q0.05, by = "ia")
d2 <- merge(d2, thompson_q0.05, by = "ia")

#
# (7b) Adjusted power
#

d2$lfk.sign.adj <- 1L * (d2$lfkindex < d2$lfkindex_q0.05)
#
d2$egger.sign.adj <- 1L * (d2$stat.egger < d2$egger_q0.05)
#
d2$rank.sign.adj <- 1L * (d2$stat.rank < d2$rank_q0.05)
#
d2$thompson.sign.adj <- 1L * (d2$stat.thompson < d2$thompson_q0.05)
#
d2 %<>% select(-(TE.common:pval.random)) %>% select(-(Q:pval.Q)) %>%
  select(-pval.egger, -pval.rank, -pval.thompson)
d2 %<>% relocate(ia, .after = last_col())
#
attr(d2, "nsim") <- nsim


#
#
# (8) Save results
#
#

type1_unif <- d0
power_unif_rho0.5 <- d1
power_unif_rho0.9 <- d2
#
save(type1_unif,
     power_unif_rho0.5,
     power_unif_rho0.9,
     file = "results/type1_power_unif.rda")
