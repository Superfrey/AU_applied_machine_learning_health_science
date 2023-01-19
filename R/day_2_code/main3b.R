# main3b.R
# Regularized logistic regression on CSF biomarker data.
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

# ---- Import data etc. ----
data = read.table(here('data-raw/csfBiomarkers.txt'),header = TRUE, sep = ",")

# Extract predictors and response variables from tables into X and y
Iresponse = names(data)=='group' # Identify column with response variable (group)
y = data[, Iresponse, drop = FALSE]
y$group = as.factor(y$group); # Convert to factor variable
X = data[,!Iresponse, drop = FALSE]

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
errTrain = errTest = accTrain = accTest = matrix(NaN, K, length(lambdas))
cMatTrain = cMatTest = list()

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

  # Get group categorical info for later use
  catInfo = levels(ytrain$group)

  # Standardize data
  scl = apply(Xtrain,2,sd)
  mn = apply(Xtrain,2,mean)
  Xtrain = as.data.frame(scale(Xtrain, scale = scl, center = mn))
  Xtest = as.data.frame(scale(Xtest, scale = scl, center = mn))

  # Fit regularized linear regression model
  # (glmnet function fits the models for all regularization strengths in lambdas in a single step)
  fit =glmnet(as.matrix(Xtrain), as.matrix(ytrain), family = "binomial", alpha = 1, lambda = lambdas, standardize = FALSE)

  #Keep coefficients for plot
  B[,,idx1]= as.matrix(fit$beta)

  # Predict
  yhatTrainProb = predict(fit, newx = as.matrix(Xtrain), type = "response")
  yhatTestProb = predict(fit, newx = as.matrix(Xtest), type = "response")

  yhatTrain = round(yhatTrainProb)
  yhatTest = round(yhatTestProb)

  # Iterate over regularization strengths to compute training- and test
  # errors for individual regularization strengths.
  for (idx2 in 1:length(lambdas)){

    # Make predictions categorical again (instead of 0/1 coding)

    yhatTrainCat = factor(yhatTrain[,idx2], levels = c("0","1"),labels = catInfo)
    yhatTestCat = factor(yhatTest[,idx2], levels = c("0","1"),labels = catInfo)

    #Evaluate classifier performance
    # Accuracy
    accTrain[idx1,idx2] = sum(yhatTrainCat==ytrain$group)/length(ytrain$group);
    accTest[idx1,idx2] = sum(yhatTestCat==ytest$group)/length(ytest$group);

    # Error rate
    errTrain[idx1,idx2] = 1 - accTrain[idx1,idx2];
    errTest[idx1,idx2] = 1 - accTest[idx1,idx2];

    # Compute confusion matrices
    cmTrain = confusionMatrix(data = yhatTrainCat, reference = ytrain$group, positive = "Impaired")
    cmTest = confusionMatrix(data = yhatTestCat, reference = ytest$group, positive = "Impaired")
    cmTrain = as.data.frame(cmTrain$table)
    cmTest = as.data.frame(cmTest$table)

    if (idx1==1){
      cMatTrain[[idx2]] = cmTrain
      cMatTest[[idx2]] = cmTest
    } else {
      cMatTrain[[idx2]][,3] = cMatTrain[[idx2]][,3] + cmTrain[,3]
      cMatTest[[idx2]][,3] = cMatTest[[idx2]][,3] + cmTest[,3]
    }
  }
}

# ---- Plot error vs. regularization strength ----
(ggplot() +
   labs(x = "log2(lambda)", y = "error rate")+
   geom_point(aes(x = log2(lambdas), y = apply(errTrain,2,mean), colour = "blue")) +
   geom_line (aes(x = log2(lambdas), y = apply(errTrain,2,mean), colour = "blue")) +
   geom_point(aes(x = log2(lambdas), y = apply(errTest,2,mean), colour = "red")) +
   geom_line (aes(x = log2(lambdas), y = apply(errTest,2,mean), colour = "red")) +
   scale_color_discrete(labels = c("train","test"))
)

# ---- Plot coefficient traces ----
Bmedian = t(apply(B,c(1,2),median));
Bmean = t(apply(B,c(1,2),mean));

# Convert to data frames for easy plotting
colnames(Bmedian) = names(X)
colnames(Bmean) = names(X)
dfBmedian = melt(data.frame(Bmedian, lambda = lambdas), id.vars = "lambda", variable.name = "beta")
dfBmean = melt(data.frame(Bmean, lambda = lambdas), id.vars = "lambda", variable.name = "beta")

# Plot
(ggplot(dfBmedian, aes(x = log2(lambda), y = value, colour = beta)) +
    geom_point() +
    geom_line() +
    theme(legend.position="none")
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

# Plot confusion matrix
ggplot(cMatTrain[[idxLambda]], aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("training")

ggplot(cMatTest[[idxLambda]], aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("test")

