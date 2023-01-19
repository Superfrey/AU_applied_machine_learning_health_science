# main6a.R
# Hierarchical clustering - CSF biomarker data
# PMR Dec 2022

# ---- Import libraries and functions ----
library("stats")
library("ggplot2")
library("ggdendro")

# ---- Import data etc. ----
X = read.table(here('data-raw/csfBiomarkers.txt',header = TRUE, sep = ",")

# Discard group variable
X = X[,!colnames(X)=="group"]

# ---- Plot standard deviation of individual features ----
df = data.frame(std = apply(X,2,sd))
df = cbind(feature = rownames(df),df)

(p = ggplot(data =df, aes(x  =feature,y = std)) +
    geom_bar(stat = "identity", width = .75) +
    theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust=0))
)

# ---- Standardize ----
X = scale(X, center = TRUE, scale = TRUE)

# ---- Hierarchical clustering of features ----
D = dist(t(X))
fit = hclust(D, method = "average")

# Plot dendrogram
ggdendrogram(fit)

