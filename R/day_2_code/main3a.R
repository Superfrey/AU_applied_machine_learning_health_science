# main3a.m
# Regularized linear regression on body density data.
# Cross-validation
#
# PMR Dec 2022

# ---- Import libraries and functions ----
library("stats")
library("ggplot2")
library("caret")
library("glmnet")
library("reshape2")
library("here")

# ---- Import data ----
data = read.table(here('data-raw/bodyMeasurements.txt'),header = TRUE, sep = ",")

# Extract predictors and response variables from tables into X and y
I = colnames(data)=='Body_Density' # Identify column with response variable (body density)
y = subset(data, select = I)
X = subset(data, select = !I)

p = ncol(X)

# ---- Create random partition of data ----
set.seed(0) # For reproducibility
K = 10
c = createFolds(y = as.matrix(y), k = K, returnTrain = TRUE) # Returns training indices
idxAll = 1:nrow(y)

# ---- Define sequence for regularization strengths ----
lambdas = 2^seq(5, -15, by = -0.5)

# ---- Run analysis ----

# Initialize array to store coefficients across CV-iterations
B = array(NaN, dim = c(p, length(lambdas),K))
errTrain = errTest = matrix(NaN, K, length(lambdas))

# Iterate over partitions
for (idx1 in 1:K) {

  # Print status
  print(sprintf('cross-validation iteration %d of %d',idx1,K))

  # Get training- and test sets
  I_train = is.element(idxAll, c[[idx1]])
  I_test = !is.element(idxAll, c[[idx1]])
  #
  Xtrain = X[I_train,, drop = FALSE]
  ytrain = y[I_train,, drop = FALSE]
  Xtest = X[I_test,, drop = FALSE]
  ytest = y[I_test,, drop = FALSE]

  # Standardize data
  scl = apply(Xtrain,2,sd)
  mn = apply(Xtrain,2,mean)
  Xtrain = as.data.frame(scale(Xtrain, scale = scl, center = mn))
  Xtest = as.data.frame(scale(Xtest, scale = scl, center = mn))

  # Fit regularized linear regression model
  # (glmnet function fits the models for all regularization strengths in lambdas in a single step)
  fit =glmnet(as.matrix(Xtrain), as.matrix(ytrain), family = "gaussian", alpha = 0, lambda = lambdas, standardize = FALSE)

  #Keep coefficients for plot
  B[,,idx1]= as.matrix(fit$beta)

  # Predict
  yhatTrain = predict(fit, newx = as.matrix(Xtrain))
  yhatTest = predict(fit, newx = as.matrix(Xtest))


  # Iterate over regularization strengths to compute training- and test
  # errors for individual regularization strengths.
  for (idx2 in 1:length(lambdas)){
    errTrain[idx1,idx2] = mean(as.matrix((ytrain-yhatTrain[,idx2])^2))
    errTest[idx1,idx2] = mean(as.matrix((ytest-yhatTest[,idx2])^2))
  }
}

# ---- Plot error vs. regularization strength ----
error_regular_plot <- (ggplot() +
   labs(x = "log2(lambda)", y = "mse")+
   geom_point(aes(x = log2(lambdas), y = apply(errTrain,2,mean), colour = "blue")) +
   geom_line (aes(x = log2(lambdas), y = apply(errTrain,2,mean), colour = "blue")) +
   geom_point(aes(x = log2(lambdas), y = apply(errTest,2,mean), colour = "red")) +
   geom_line (aes(x = log2(lambdas), y = apply(errTest,2,mean), colour = "red")) +
   scale_y_continuous(trans='log2') +
   scale_color_discrete(labels = c("train","test"))
)

## ---- Plot coefficient traces ----
Bmedian = t(apply(B,c(1,2),median));
Bmean = t(apply(B,c(1,2),mean));

# Convert to data frames for easy plotting
colnames(Bmedian) = names(X)
colnames(Bmean) = names(X)
dfBmedian = melt(data.frame(Bmedian, lambda = lambdas), id.vars = "lambda", variable.name = "beta")
dfBmean = melt(data.frame(Bmean, lambda = lambdas), id.vars = "lambda", variable.name = "beta")

# Plot
coef_regular_plot <- (ggplot(dfBmedian, aes(x = log2(lambda), y = value, colour = beta)) +
    geom_point()+
    geom_line()
)

# ---- Explore best model ----
# Look at best model
idxLambda = which.min(apply(errTest,2,mean))
bestLambda = lambdas[idxLambda]
print(sprintf('best lambda: log2(lambda)=%.3f',log2(lambdas[idxLambda])))

# Look at error rates
cat(sprintf('training error: %.12f',mean(errTrain[,idxLambda])))
cat(sprintf('validation error: %.12f',mean(errTest[,idxLambda])))

# Look at coefficients for chosen model (mean and median across CV-iterations)
# Sort features according to absolute median weight value for the chosen model
betaMed = abs(Bmedian[idxLambda,])
betaMean = abs(Bmean[idxLambda,])
idx = order(betaMed,decreasing = TRUE)
cat(sprintf('\n\nrank\t%35s\tbeta_med\t\tbeta_mean\n','feature'))
for (idx1 in 1:length(idx)){
   cat(sprintf('%d\t%35s\t%.16f\t%.16f\n',idx1,names(betaMean)[idx[idx1]],betaMed[idx[idx1]],betaMean[idx[idx1]]))
}

