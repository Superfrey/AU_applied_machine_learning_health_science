# main6a.R
# Principal component analysis - body density data
# PMR Dec 2022

# ---- Import libraries and functions ----
library("stats")
library("corrplot")
library("ggplot2")
library(here)

# ---- Import data ----
X = read.table(here('data-raw/bodyMeasurements.txt'),header = TRUE, sep = ",")

# Make each feature zero mean (center data)
X = scale(X, center = TRUE, scale = FALSE)

# ---- Plot correlation matrix ----
C = cor(X)
corrplot(C)

# ---- Plot standard deviation of individual features ----
df = data.frame(std = apply(X,2,sd))
df = cbind(feature = rownames(df),df)

(p = ggplot(data =df, aes(x  =feature,y = std)) +
    geom_bar(stat = "identity", width = .75) +
    theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust=0))
  )

# ---- Make standardised data, where every column has standard deviation 1 and mean 0.
X_standardized = scale(X, center = TRUE, scale = TRUE)

# ---- Compute PCA on a pair of features ----

featureIdx = c(6, 13) # select a pair of features to analyse, try to select other combinations

# Compute PCA on non-scaled and scaled data
fit1 = prcomp(X[,featureIdx])
fit2 = prcomp(X_standardized[,featureIdx])

# Compute variance explained
fit1$pve = fit1$sdev^2/sum(fit1$sdev^2)
fit2$pve = fit2$sdev^2/sum(fit2$sdev^2)

# Compute scores
fit1$scores = X[,featureIdx]%*%fit1$rotation
fit2$scores = X_standardized[,featureIdx]%*%fit2$rotation

# Plot data and PC projections
par(mfrow = c(2,2))
# Non-standardized data
plot(X[,featureIdx[1]], X[,featureIdx[2]],
     xlab = colnames(X)[featureIdx[1]],
     ylab = colnames(X)[featureIdx[2]],
     asp = 1)
xx = seq(min(X[,featureIdx[1]]), max(X[,featureIdx[1]]), length.out = 10)
yy1 = xx*fit1$rotation[2,1]/fit1$rotation[1,1]
yy2 = xx*fit1$rotation[2,2]/fit1$rotation[1,2]
lines(xx,yy1, col="red")
lines(xx,yy2, col="blue")
#
# Standardized data
plot(X_standardized[,featureIdx[1]], X_standardized[,featureIdx[2]],
     xlab = colnames(X_standardized)[featureIdx[1]],
     ylab = colnames(X_standardized)[featureIdx[2]],
     asp = 1)
xx = seq(min(X_standardized[,featureIdx[1]]), max(X_standardized[,featureIdx[1]]), length.out = 10)
yy1 = xx*fit2$rotation[2,1]/fit2$rotation[1,1]
yy2 = xx*fit2$rotation[2,2]/fit2$rotation[1,2]
lines(xx,yy1, col="red")
lines(xx,yy2, col="blue")
#
# Scores, non-standardized data
plot(fit1$scores[,1],fit1$scores[,2],
     xlab = colnames(fit1$scores)[1],
     ylab = colnames(fit1$scores)[2],
     asp = 1)
#
# Scores, standardized data
plot(fit2$scores[,1],fit2$scores[,2],
     xlab = colnames(fit2$scores)[1],
     ylab = colnames(fit2$scores)[2],
     asp = 1)


# ---- Compute PCA from all features  ----
# Compute PCA on non-scaled and scaled data
fit1 = prcomp(X)
fit2 = prcomp(X_standardized)

# Compute variance explained
fit1$pve = fit1$sdev^2/sum(fit1$sdev^2)
fit2$pve = fit2$sdev^2/sum(fit2$sdev^2)

# Compute scores
fit1$scores = X%*%fit1$rotation
fit2$scores = X_standardized%*%fit2$rotation

# Plot data and PC projections
par(mfrow = c(1,3))
#
# Scores, non-standardized data
plot(fit1$scores[,1],fit1$scores[,2],
     col = "red",
     xlab = colnames(fit1$scores)[1],
     ylab = colnames(fit1$scores)[2])
#
# Scores, standardized data
plot(fit2$scores[,1],fit2$scores[,2],
     col = "blue",
     xlab = colnames(fit2$scores)[1],
     ylab = colnames(fit2$scores)[2])
#
# Plot pve
plot(fit1$pve, col = "red", type = "b",
     xlab = "principal component",
     ylab = "pve")
lines(fit2$pve, col = "blue", type = "b")

