# Analyzing the spread of SARS-CoV-2 variants Alpha in Switzerland and Beta in South Africa
# Christian L. Althaus, 7 June 2021

# Load libraries
library(lubridate)
library(plotrix)
library(bbmle)
library(mvtnorm)
library(RColorBrewer)

# Set seed for random number generator
set.seed(6687695)

# Define colors
cols <- brewer.pal(4, "Set1")
t.cols <- cols
for(i in 1:length(cols)) {
  x <- col2rgb(cols[i])
  t.cols[i] <- rgb(x[1, ], x[2, ], x[3, ], alpha = 125, maxColorValue = 255)
}

# Sample size for parameter sets and bootstrapping
n_sim <- 1e4

# Generation time
gen_mean <- 5.2
gen_sd <- (6.78 - 3.78)/2/qnorm(0.975)
generation <- rnorm(n_sim, gen_mean, gen_sd)
generation <- ifelse(generation > 0, generation, rnorm(1, gen_mean, gen_sd))
par(mfrow = c(1, 1))
hist(generation)

# Reproduction number
#R_ge <- read.csv("https://raw.githubusercontent.com/covid-19-Re/dailyRe-Data/master/CHE-estimates.csv")
R_ge <- read.csv("../data/CHE-estimates.csv")
R_ge$date <- ymd(R_ge$date)
R_ge <- subset(R_ge,
               region == "GE"
               & data_type == "Confirmed cases"
               & estimate_type == "Cori_slidingWindow"
               & date >= ymd(20201101) # November 2020 to January 2021
               & date <= ymd(20210131))
R_ge <- R_ge$median_R_mean
R_ge <- sample(R_ge, n_sim, replace = TRUE)
par(mfrow = c(1, 2))
hist(R_ge)

#R_sa <- read.csv("https://raw.githubusercontent.com/covid-19-Re/dailyRe-Data/master/ZAF-estimates.csv")
R_sa <- read.csv("../data/ZAF-estimates.csv")
R_sa$date <- ymd(R_sa$date)
R_sa <- subset(R_sa,
               region == "ZAF"
               & data_type == "Confirmed cases"
               & estimate_type == "Cori_slidingWindow"
               & date >= ymd(20200901) # September and October 2020
               & date <= ymd(20201031))
R_sa <- R_sa$median_R_mean
R_sa <- sample(R_sa, n_sim, replace = TRUE)
hist(R_sa)

# Seroprevalence
# Geneva: Stringhini et al. (2021)
sero_ge_mean <- 0.211
sero_ge_sd <- (0.231 - 0.192)/2/qnorm(0.975)
sero_ge <- rnorm(n_sim, sero_ge_mean, sero_ge_sd)
sero_ge <- ifelse(sero_ge > 0, sero_ge, rnorm(1, sero_ge_mean, sero_ge_sd))
hist(sero_ge)

# South Africa: https://www.medrxiv.org/content/10.1101/2021.02.25.21252477v1
sero_sa <- rbinom(n_sim, 4387, 1324/4387)/4387
hist(sero_sa)

# Plot parameter distributions
pdf("../figures/parameters.pdf", width = 9, height = 3)
par(mfrow = c(1, 3))

# Generation time
d <- density(generation)
plot(d,
     xlim = c(0, 10), col = cols[4],
     xlab = "Generation time (days)", ylab = "Density", main = NA, frame = FALSE)
polygon(d, border = NA, col = t.cols[4])
mtext("A", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

# Reproduction number
d <- density(R_ge, adjust = 1)
plot(d,
     xlim = c(0, 1.5), ylim = c(0, 8), col = cols[3],
     xlab = "Effective reproduction number", ylab = "Density", main = NA, frame = FALSE)
polygon(d, border = NA, col = t.cols[3])
d <- density(R_sa, adjust = 1)
lines(d, col = cols[2])
polygon(d, border = NA, col = t.cols[2])
#abline(v = 1, lty = 3)
legend("topleft", inset = 0, c("Geneva, Switzerland", "South Africa"), col = cols[3:2], lty = 1, bty = "n")
mtext("B", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

# Seroprevalence
d <- density(sero_ge, adjust = 1)
plot(d,
     xlim = c(0, 0.4), ylim = c(0, 60), col = cols[3],
     xlab = "Seroprevalence", ylab = "Density", main = NA, axes = FALSE, frame = FALSE)
polygon(d, border = NA, col = t.cols[3])
d <- density(sero_sa, adjust = 1)
lines(d, col = cols[2])
polygon(d, border = NA, col = t.cols[2])
axis(1, seq(0, 0.4, 0.1), paste0(seq(0, 0.4, 0.1)*1e2, "%"))
axis(2)
legend("topleft", inset = 0, c("Geneva, Switzerland", "South Africa"), col = cols[3:2], lty = 1, bty = "n")
mtext("C", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

dev.off()

# Estimate logistic growth rate of variants

# Logistic growth function
logistic <- function(y_0, rho, x) {
  y <- 1/((1 - y_0)/y_0*exp(- rho*x) + 1)
  return(y)
}

# Negative log-likelihood function
nll <- function(p_0, rho) {
  p_0 <- plogis(p_0)
  rho <- rho
  times <- as.numeric(variant_data$date - min(variant_data$date))
  p <- logistic(p_0, rho, times)
  ll <- sum(dbinom(variant_data$variant, variant_data$total, p, log = TRUE))
  return(-ll)
}

# Load data from Geneva, Switzerland
variant_data <- read.csv("../data/variants_GE.csv")

# Clean data
variant_data$date <- ymd(variant_data$date)
names(variant_data)[3] <- "variant"
variant_data$frequency <- variant_data$variant/variant_data$total

# Add binomial confidence intervals
lower <- numeric()
upper <- numeric()
for(i in 1:length(variant_data$variant)) {
  int <- binom.test(variant_data$variant[i], variant_data$total[i])
  lower[i] <- int$conf.int[1]
  upper[i] <- int$conf.int[2]
}
variant_data <- cbind(variant_data, lower = lower, upper = upper)

# Fit logistic growth curve (using parameter transformation)
free <- c(p_0 = qlogis(0.01), rho = 0.1)
fit <- mle2(nll, start = as.list(free), method = "Nelder-Mead")
fit.ci <- confint(fit)

# Define time window
time_window <- c(ymd(20201101), ymd(20210401))

# Bootstrap sampling to compute confidence interval for logistic growth
m <- coef(fit)
rho_ge <- m[2]
sigma <- vcov(fit)
sim_coef <- data.frame(rmvnorm(n_sim, mean = m, sigma = sigma))
rho_sample_ge <- sim_coef[, 2]
sim_coef$p_0 <- plogis(sim_coef$p_0)
times <- as.numeric(as_date(time_window[1]:time_window[2]) - min(variant_data$date))
proportion <- array(NA, c(n_sim, length(times)))

for(i in 1:n_sim) {
  p <- logistic(sim_coef$p_0[i], sim_coef$rho[i], times)
  proportion[i, ] <- p
}
interval_95 <- apply(proportion, MAR = 2, FUN = quantile, probs = c(0.025, 0.975))

pdf("../figures/spread.pdf", width = 12, height = 5)
par(mfrow = c(1, 2))
plot(variant_data$date, variant_data$frequency,
     type = "n", xlim = time_window, ylim = c(0, 1),
     xlab = NA, ylab = "Proportion Alpha (501Y)", main = "Geneva, Switzerland", axes = FALSE, frame = FALSE)
label <- time_window[1] + months(0:5)
axis(1, label, paste(day(label), month(label, label = TRUE)))
axis(2)
mtext("A", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)
abline(v = label, h = seq(0, 1, 0.2), col = "lightgray", lty = 3)
polygon(c(time_window[1]:time_window[2], rev(time_window[1]:time_window[2])),
        c(interval_95[1, ], rev(interval_95[2, ])),
        col = t.cols[1], border = NA)
times <- as.numeric(as_date(time_window[1]:time_window[2]) - min(variant_data$date))
p <- logistic(plogis(coef(fit)[1]), coef(fit)[2], times)
lines(time_window[1]:time_window[2], p, col = cols[1])
plotCI(variant_data$date, variant_data$frequency, ui = upper, li = lower,
       pch = 21, pt.bg = "white", cex = 0.75, sfrac = 0.0025, col = cols[2], add = TRUE)

# Load data from South Africa
variant_data <- read.csv("../data/variants_SA.csv")

# Clean data
variant_data$date <- ymd(variant_data$date)
names(variant_data)[3] <- "variant"
variant_data$frequency <- variant_data$variant/variant_data$total

# Add binomial confidence intervals
lower <- numeric()
upper <- numeric()
for(i in 1:length(variant_data$variant)) {
  int <- binom.test(variant_data$variant[i], variant_data$total[i])
  lower[i] <- int$conf.int[1]
  upper[i] <- int$conf.int[2]
}
variant_data <- cbind(variant_data, lower = lower, upper = upper)

# Fit logistic growth curve (using parameter transformation)
free <- c(p_0 = qlogis(0.01), rho = 0.1)
fit <- mle2(nll, start = as.list(free), method = "Nelder-Mead")
fit.ci <- confint(fit)

# Define time window
time_window <- c(ymd(20200901), ymd(20210201))

# Bootstrap sampling to compute confidence interval for logistic growth
m <- coef(fit)
rho_sa <- m[2]
sigma <- vcov(fit)
sim_coef <- data.frame(rmvnorm(n_sim, mean = m, sigma = sigma))
rho_sample_sa <- sim_coef[, 2]
sim_coef$p_0 <- plogis(sim_coef$p_0)
times <- as.numeric(as_date(time_window[1]:time_window[2]) - min(variant_data$date))
proportion <- array(NA, c(n_sim, length(times)))

for(i in 1:n_sim) {
  p <- logistic(sim_coef$p_0[i], sim_coef$rho[i], times)
  proportion[i, ] <- p
}
interval_95 <- apply(proportion, MAR = 2, FUN = quantile, probs = c(0.025, 0.975))

plot(variant_data$date, variant_data$frequency,
     type = "n", xlim = time_window, ylim = c(0, 1),
     xlab = NA, ylab = "Proportion Beta", main = "South Africa", axes = FALSE, frame = FALSE)
label <- time_window[1] + months(0:5)
axis(1, label, paste(day(label), month(label, label = TRUE)))
axis(2)
mtext("B", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)
abline(v = label, h = seq(0, 1, 0.2), col = "lightgray", lty = 3)
polygon(c(time_window[1]:time_window[2], rev(time_window[1]:time_window[2])),
        c(interval_95[1, ], rev(interval_95[2, ])),
        col = t.cols[1], border = NA)
times <- as.numeric(as_date(time_window[1]:time_window[2]) - min(variant_data$date))
p <- logistic(plogis(coef(fit)[1]), coef(fit)[2], times)
lines(time_window[1]:time_window[2], p, col = cols[1])
plotCI(variant_data$date, variant_data$frequency, ui = upper, li = lower,
       pch = 21, pt.bg = "white", cex = 0.75, sfrac = 0.0025, col = cols[2], add = TRUE)
dev.off()

# Compute parameters
output <- data.frame(matrix(NA, nrow = 8, ncol = 6))
names(output) <- c("parameter", "region", "median", "lower", "upper", "p")

# Rho (logistic growth rate)
output[1, ] <- c("rho", "GE", round(rho_ge, 3), round(quantile(rho_sample_ge, probs = c(0.025, 0.975)), 3), NA)
output[2, ] <- c("rho", "SA", round(rho_sa, 3), round(quantile(rho_sample_sa, probs = c(0.025, 0.975)), 3), NA)

# Tau (increased transmissibility)
output[3, ] <- c("tau", "GE", round(quantile(rho_sample_ge*generation/R_ge, probs = c(0.5, 0.025, 0.975)), 2), NA)
output[4, ] <- c("tau", "SA", round(quantile(rho_sample_sa*generation/R_sa, probs = c(0.5, 0.025, 0.975)), 2), NA)

# Kappa (longer generation time)
output[5, ] <- c("kappa", "GE", round(quantile(rho_sample_ge/(1/generation - rho_sample_ge), probs = c(0.5, 0.025, 0.975)), 2), NA)
output[6, ] <- c("kappa", "SA", round(quantile(rho_sample_sa/(1/generation - rho_sample_sa), probs = c(0.5, 0.025, 0.975)), 2), NA)

# Epsilon (immune evasion)
s <- rho_sample_ge*(1 - sero_ge)*generation/(sero_ge*R_ge)
s <- s[s <= 1]
s <- length(s)/n_sim
output[7, ] <- c("epsilon", "GE", round(quantile(rho_sample_ge*(1 - sero_ge)*generation/(sero_ge*R_ge), probs = c(0.5, 0.025, 0.975)), 2), round(s, 2))
s <- rho_sample_sa*(1 - sero_sa)*generation/(sero_sa*R_sa)
s <- s[s <= 1]
s <- length(s)/n_sim
output[8, ] <- c("epsilon", "SA", round(quantile(rho_sample_sa*(1 - sero_sa)*generation/(sero_sa*R_sa), probs = c(0.5, 0.025, 0.975)), 2), round(s, 2))

write.csv(output, "../out/results.csv", quote = FALSE, row.names = FALSE)

# Sensitivity analysis

# Define functions
tau_func <- function(rho, kappa, epsilon, omega, D, R_w) {
  beta <- R_w/((1 - omega)*D)
  S <- 1 - omega
  gamma <- 1/D
  return(((kappa + 1)*(rho - beta*omega*epsilon) - gamma*kappa)/(beta*(kappa + 1)*(S + omega*epsilon)))
}

kappa_func <- function(rho, tau, epsilon, omega, D, R_w) {
  beta <- R_w/((1 - omega)*D)
  S <- 1 - omega
  gamma <- 1/D
  return((rho - beta*(S*tau + (tau + 1)*omega*epsilon))/(gamma - rho + beta*(S*tau + tau*omega*epsilon + omega*epsilon)))
}

epsilon_func <- function(rho, tau, kappa, omega, D, R_w) {
  beta <- R_w/((1 - omega)*D)
  S <- 1 - omega
  gamma <- 1/D
  return(((kappa + 1)*(rho - beta*S*tau) - gamma*kappa)/(beta*(kappa + 1)*(tau + 1)*omega))
}

# Estimates for South Africa for epsilon = 0.5
quantile(tau_func(rho_sample_sa, 0, 0.5, sero_sa, generation, R_sa), probs = c(0.5, 0.025, 0.975))
quantile(kappa_func(rho_sample_sa, 0, 0.5, sero_sa, generation, R_sa), probs = c(0.5, 0.025, 0.975))

pdf("../figures/relation.pdf", width = 9, height = 6)
par(mfrow = c(2, 3))
# Tau vs. kappa
plot(NA,
     xlim = c(0, 2), ylim = c(0, 1),
     xlab = expression(paste("Increased infectious duration (", kappa, ")")), ylab = expression(paste("Increased transmissibility (", tau, ")")),
     axes = FALSE, frame = FALSE)
axis(1, seq(0, 2, 0.5), paste0(seq(0, 2, 0.5)*1e2, "%"))
axis(2, seq(0, 1, 0.25), paste0(seq(0, 1, 0.25)*1e2, "%"))
abline(v = seq(0, 2, 0.5), h = seq(0, 1, 0.25), col = "lightgray", lty = 3)
mtext("A", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)
legend("topleft", inset = 0.05, "Alpha in Geneva, Switzerland", col = cols[3], lty = 1, bty = "n")

x <- seq(0, 2, 0.01)
y <- matrix(NA, nrow = length(x), ncol = 3)
for(i in 1:length(x)) {
  y[i, 1:3] <- quantile(tau_func(rho_sample_ge, x[i], 0, sero_ge, generation, R_ge), probs = c(0.5, 0.025, 0.975))
}

y[, 2][y[, 2] > 1] <- 1
y[, 3][y[, 3] > 1] <- 1
y[, 2][y[, 2] < 0] <- 0
y[, 3][y[, 3] < 0] <- 0
polygon(c(x, rev(x)),
        c(y[, 2], rev(y[, 3])),
        col = t.cols[3], border = NA)
lines(x[y[, 1] > 0 & y[, 1] < 1], y[, 1][y[, 1] > 0 & y[, 1] < 1], col = cols[3])

# Epsilon vs. tau
# Tau vs. kappa
plot(NA,
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = expression(paste("Increased transmissibility (", tau, ")")), ylab = expression(paste("Immune evasion (", epsilon, ")")),
     axes = FALSE, frame = FALSE)
axis(1, seq(0, 1, 0.25), paste0(seq(0, 1, 0.25)*1e2, "%"))
axis(2, seq(0, 1, 0.25), paste0(seq(0, 1, 0.25)*1e2, "%"))
abline(v = seq(0, 1, 0.25), h = seq(0, 1, 0.25), col = "lightgray", lty = 3)
mtext("B", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

x <- seq(0, 1, 0.001)
y <- matrix(NA, nrow = length(x), ncol = 3)
for(i in 1:length(x)) {
  y[i, 1:3] <- quantile(epsilon_func(rho_sample_ge, x[i], 0, sero_ge, generation, R_ge), probs = c(0.5, 0.025, 0.975))
}

y[, 2][y[, 2] > 1] <- 1
y[, 3][y[, 3] > 1] <- 1
y[, 2][y[, 2] < 0] <- 0
y[, 3][y[, 3] < 0] <- 0
polygon(c(x, rev(x)),
        c(y[, 2], rev(y[, 3])),
        col = t.cols[3], border = NA)
lines(x[y[, 1] > 0 & y[, 1] < 1], y[, 1][y[, 1] > 0 & y[, 1] < 1], col = cols[3])

# Epsilon vs. kappa
plot(NA,
     xlim = c(0, 2), ylim = c(0, 1),
     xlab = expression(paste("Increased infectious duration (", kappa, ")")), ylab = expression(paste("Immune evasion (", epsilon, ")")),
     axes = FALSE, frame = FALSE)
axis(1, seq(0, 2, 0.5), paste0(seq(0, 2, 0.5)*1e2, "%"))
axis(2, seq(0, 1, 0.25), paste0(seq(0, 1, 0.25)*1e2, "%"))
abline(v = seq(0, 3, 0.5), h = seq(0, 1, 0.25), col = "lightgray", lty = 3)
mtext("C", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

x <- seq(0, 2, 0.01)
y <- matrix(NA, nrow = length(x), ncol = 3)
for(i in 1:length(x)) {
  y[i, 1:3] <- quantile(epsilon_func(rho_sample_ge, 0, x[i], sero_ge, generation, R_ge), probs = c(0.5, 0.025, 0.975))
}

y[, 2][y[, 2] > 1] <- 1
y[, 3][y[, 3] > 1] <- 1
y[, 2][y[, 2] < 0] <- 0
y[, 3][y[, 3] < 0] <- 0
polygon(c(x, rev(x)),
        c(y[, 2], rev(y[, 3])),
        col = t.cols[3], border = NA)
lines(x[y[, 1] > 0 & y[, 1] < 1], y[, 1][y[, 1] > 0 & y[, 1] < 1], col = cols[3])

# Tau vs. kappa
plot(NA,
     xlim = c(0, 2), ylim = c(0, 1),
     xlab = expression(paste("Increased infectious duration (", kappa, ")")), ylab = expression(paste("Increased transmissibility (", tau, ")")),
     axes = FALSE, frame = FALSE)
axis(1, seq(0, 2, 0.5), paste0(seq(0, 2, 0.5)*1e2, "%"))
axis(2, seq(0, 1, 0.25), paste0(seq(0, 1, 0.25)*1e2, "%"))
abline(v = seq(0, 2, 0.5), h = seq(0, 1, 0.25), col = "lightgray", lty = 3)
mtext("D", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)
legend("topleft", inset = 0.05, "Beta in South Africa", col = cols[2], lty = 1, bty = "n")

x <- seq(0, 2, 0.01)
y <- matrix(NA, nrow = length(x), ncol = 3)
for(i in 1:length(x)) {
  y[i, 1:3] <- quantile(tau_func(rho_sample_sa, x[i], 0, sero_sa, generation, R_sa), probs = c(0.5, 0.025, 0.975))
}

y[, 2][y[, 2] > 1] <- 1
y[, 3][y[, 3] > 1] <- 1
y[, 2][y[, 2] < 0] <- 0
y[, 3][y[, 3] < 0] <- 0
polygon(c(x, rev(x)),
        c(y[, 2], rev(y[, 3])),
        col = t.cols[2], border = NA)
lines(x[y[, 1] > 0 & y[, 1] < 1], y[, 1][y[, 1] > 0 & y[, 1] < 1], col = cols[2])

# Epsilon vs. tau
# Tau vs. kappa
plot(NA,
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = expression(paste("Increased transmissibility (", tau, ")")), ylab = expression(paste("Immune evasion (", epsilon, ")")),
     axes = FALSE, frame = FALSE)
axis(1, seq(0, 1, 0.25), paste0(seq(0, 1, 0.25)*1e2, "%"))
axis(2, seq(0, 1, 0.25), paste0(seq(0, 1, 0.25)*1e2, "%"))
abline(v = seq(0, 1, 0.25), h = seq(0, 1, 0.25), col = "lightgray", lty = 3)
mtext("E", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

x <- seq(0, 1, 0.001)
y <- matrix(NA, nrow = length(x), ncol = 3)
for(i in 1:length(x)) {
  y[i, 1:3] <- quantile(epsilon_func(rho_sample_sa, x[i], 0, sero_sa, generation, R_sa), probs = c(0.5, 0.025, 0.975))
}

y[, 2][y[, 2] > 1] <- 1
y[, 3][y[, 3] > 1] <- 1
y[, 2][y[, 2] < 0] <- 0
y[, 3][y[, 3] < 0] <- 0
polygon(c(x, rev(x)),
        c(y[, 2], rev(y[, 3])),
        col = t.cols[2], border = NA)
lines(x[y[, 1] > 0 & y[, 1] < 1], y[, 1][y[, 1] > 0 & y[, 1] < 1], col = cols[2])

# Epsilon vs. kappa
plot(NA,
     xlim = c(0, 2), ylim = c(0, 1),
     xlab = expression(paste("Increased infectious duration (", kappa, ")")), ylab = expression(paste("Immune evasion (", epsilon, ")")),
     axes = FALSE, frame = FALSE)
axis(1, seq(0, 2, 0.5), paste0(seq(0, 2, 0.5)*1e2, "%"))
axis(2, seq(0, 1, 0.25), paste0(seq(0, 1, 0.25)*1e2, "%"))
abline(v = seq(0, 3, 0.5), h = seq(0, 1, 0.25), col = "lightgray", lty = 3)
mtext("F", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

x <- seq(0, 2, 0.01)
y <- matrix(NA, nrow = length(x), ncol = 3)
for(i in 1:length(x)) {
  y[i, 1:3] <- quantile(epsilon_func(rho_sample_sa, 0, x[i], sero_sa, generation, R_sa), probs = c(0.5, 0.025, 0.975))
}

y[, 2][y[, 2] > 1] <- 1
y[, 3][y[, 3] > 1] <- 1
y[, 2][y[, 2] < 0] <- 0
y[, 3][y[, 3] < 0] <- 0
polygon(c(x, rev(x)),
        c(y[, 2], rev(y[, 3])),
        col = t.cols[2], border = NA)
lines(x[y[, 1] > 0 & y[, 1] < 1], y[, 1][y[, 1] > 0 & y[, 1] < 1], col = cols[2])

dev.off()

# Comparison between Alpha and Beta
rho_func <- function(tau, kappa, epsilon, omega, D, R_w) {
  beta <- R_w/((1 - omega)*D)
  S <- 1 - omega
  gamma <- 1/D
  rho <- (1 + tau)*beta*(S + omega*epsilon) - gamma/(1 + kappa) - beta*S + gamma
  return(quantile(rho, probs = c(0.5, 0.025, 0.975)))
}

# Plot various combinations of R_e and epsilon
pdf("../figures/comparison.pdf", width = 9, height = 3)
par(mfrow = c(1, 3))
lab <- c("A", "B", "C")
omega_range <- seq(0, 0.999, 0.001)
R_e_range <- 1 #c(0.5, 1, 1.5)
epsilon_range <- c(0.25, 0.5, 0.75)

for(j in 1:length(R_e_range)) {
  R_e <- R_e_range[j]
  growth_ge <- data.frame(matrix(NA, nrow = length(omega_range), ncol = 3))
  tau_sample <- tau_func(rho_sample_ge, 0, 0, sero_ge, generation, R_ge)
  for(i in 1:length(omega_range)) {
    growth_ge[i, 1:3] <- rho_func(tau_sample, 0, 0, omega_range[i], generation, R_e)
  }
  for(k in 1:length(epsilon_range)) {
    epsilon_sa <- epsilon_range[k]
    growth_sa <- data.frame(matrix(NA, nrow = length(omega_range), ncol = 3))
    tau_sample <- tau_func(rho_sample_sa, 0, epsilon_sa, sero_sa, generation, R_sa)
    for(i in 1:length(omega_range)) {
      growth_sa[i, 1:3] <- rho_func(tau_sample, 0, epsilon_sa, omega_range[i], generation, R_e)
    }
    plot(NA,
         xlim = c(0, 1), ylim = c(0, 0.4),
         xlab = "Seroprevalence ('wild-type')", ylab = "Growth advantage (per day)", frame = FALSE)
    polygon(c(omega_range, rev(omega_range)), c(growth_ge[, 2], rev(growth_ge[, 3])), col = t.cols[3], border = NA)
    polygon(c(omega_range, rev(omega_range)), c(growth_sa[, 2], rev(growth_sa[, 3])), col = t.cols[2], border = NA)
    lines(omega_range, growth_ge[, 1], col = cols[3])
    lines(omega_range, growth_sa[, 1], col = cols[2])
    if(k == 1) legend("topleft", inset = 0, c(expression(paste("Alpha: ", epsilon, " = 0")), expression(paste("Beta: ", epsilon, " = ", 0.25))), col = cols[3:2], lty = 1, bty = "n")
    if(k == 2) legend("topleft", inset = 0, c(expression(paste("Alpha: ", epsilon, " = 0")), expression(paste("Beta: ", epsilon, " = ", 0.5))), col = cols[3:2], lty = 1, bty = "n")
    if(k == 3) legend("topleft", inset = 0, c(expression(paste("Alpha: ", epsilon, " = 0")), expression(paste("Beta: ", epsilon, " = ", 0.75))), col = cols[3:2], lty = 1, bty = "n")
    mtext(lab[k], side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)
    #w <- which(round(growth_sa[, 1], 3) == round(growth_ge[1, 1], 3))
    #w <- w[round(length(w)/2)]
    #points(omega_range[w], growth_ge[1, 1], pch = 16)
    #lines(c(omega_range[w], omega_range[w]), c(- 1, growth_ge[1, 1]), lty = 2)
  }
}
dev.off()
