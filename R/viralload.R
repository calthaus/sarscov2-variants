# Analyze viral load data from Geneva
# Christian L. Althaus, 7 June 2021

library(readxl)
library(plotrix)
library(lubridate)
library(beeswarm)
library(RColorBrewer)

# Define colors
cols <- brewer.pal(4, "Set1")
t.cols <- cols
for(i in 1:length(cols)) {
  x <- col2rgb(cols[i])
  t.cols[i] <- rgb(x[1, ], x[2, ], x[3, ], alpha = 125, maxColorValue = 255)
}

# Load data
load <- read_excel("../data/variants_final_20210401.xlsx")
load <- as.data.frame(load)
#load$date_prel <- ymd_hms(load$date_prel)

nvoc <- load$vl[load$mutation_n501y == 0]
voc <- load$vl[load$mutation_n501y == 1]

mean(log10(nvoc))
mean(log10(voc))

pdf("../figures/viralload.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
boxplot(log10(vl) ~ mutation_n501y, data = load,
        xlab = NA, ylab = "SARS-CoV-2 viral load (RNA copies/ml)",
        names = c("non-VOC", "Alpha"), col = t.cols[2:1], notch = TRUE,
        axes = FALSE, frame = FALSE)
axis(1, 1:2, c("non-VOC", "Alpha"))
axis(2, 4:10, expression(10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10), las = 2)
beeswarm(log10(vl) ~ mutation_n501y, data = load, cex = 0.5, col = cols[2:1], pch = 16, add = TRUE)
mtext("A", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

nvoc <- subset(load, mutation_n501y == 0)
voc <- subset(load, mutation_n501y == 1)
nvoc_m <- numeric(3)
nvoc_u <- numeric(3)
nvoc_l <- numeric(3)
voc_m <- numeric(3)
voc_u <- numeric(3)
voc_l <- numeric(3)
for(i in 1:3) {
  x <- nvoc$vl[nvoc$onset3 == i - 1]
  nvoc_m[i] <- mean(log10(x))
  nvoc_u[i] <- mean(log10(x)) + 1.96*sd(log10(x))/sqrt(length(x))
  nvoc_l[i] <- mean(log10(x)) - 1.96*sd(log10(x))/sqrt(length(x))
  
  x <- voc$vl[voc$onset3 == i - 1]
  voc_m[i] <- mean(log10(x))
  voc_u[i] <- mean(log10(x)) + 1.96*sd(log10(x))/sqrt(length(x))
  voc_l[i] <- mean(log10(x)) - 1.96*sd(log10(x))/sqrt(length(x))
}

plotCI(1:3, nvoc_m, ui = nvoc_u, li = nvoc_l,
       xlim = c(0.5, 3.5), ylim = c(5, 8), col = cols[2], pch = 16,
       xlab = "Post onset of symptoms (days)", ylab = "SARS-CoV-2 viral load (RNA copies/ml)",
       axes = FALSE, frame = FALSE)
lines(1:3, nvoc_m, ty = "b", col = cols[2], pch = 20)
plotCI(1:3, voc_m, ui = voc_u, li = voc_l,
       col = cols[1], pch = 16, add = TRUE)
lines(1:3, voc_m, ty = "b", col = cols[1], pch = 20)
axis(1, 1:3, c("0-2", "3-5", "6-11"))
axis(2, 5:8, expression(10^5, 10^6, 10^7, 10^8), las = 2)
abline(h = 6, lty = 2)
legend("bottomleft", inset = 0.05, c("non-VOC", "Alpha"), col = cols[2:1], pch = 16, bty = "n")
mtext("B", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)
dev.off()
