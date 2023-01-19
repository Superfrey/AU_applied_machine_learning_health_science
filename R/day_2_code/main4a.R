# main4a.R
# Soft margin SVM (support vector classifier in ISL book) on CSF biomarker data.
# Cross-validation
#
# PMR Dec 2022

# ---- Import libraries and functions ----
library("stats")
library("ggplot2")
library("caret")
library("e1071")
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
Cvalues = 2^seq(15, 0, by = -0.5)

# ---- Run analysis ----

# Initialize array to store coefficients across CV-iterations
B = array(NaN, dim = c(p, length(Cvalues),K))
errTrain = errTest = accTrain = accTest = nSV = matrix(NaN, K, length(Cvalues))
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

  # Iterate over regularization strengths fit the svm, and to compute
  # training- and test errors for individual regularization strengths.
  for (idx2 in 1:length(Cvalues)){
    # Fit soft-margin svm
    fit =svm(Xtrain, ytrain, cost = 1/Cvalues[idx2], kernel = "linear", type = "C-classification", scale = FALSE) # Note: 1/C since R and ISL have different definitions for C.

    # Store the number of support vectors
    nSV[idx1,idx2] = fit$tot.nSV

    # Keep coefficients for plot
    B[,idx2,idx1]= as.matrix(t(fit$SV) %*% fit$coefs)

    # Predict
    yhatTrainCat = predict(fit, Xtrain)
    yhatTestCat = predict(fit, Xtest)

    # Evaluate classifier performance
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
error_plot_svm <- (ggplot() +
   labs(x = "log2(C)", y = "error rate")+
   geom_point(aes(x = log2(Cvalues), y = apply(errTrain,2,mean), colour = "blue")) +
   geom_line (aes(x = log2(Cvalues), y = apply(errTrain,2,mean), colour = "blue")) +
   geom_point(aes(x = log2(Cvalues), y = apply(errTest,2,mean), colour = "red")) +
   geom_line (aes(x = log2(Cvalues), y = apply(errTest,2,mean), colour = "red")) +
   scale_color_discrete(labels = c("train","test"))
)
error_plot_svm

# ---- Plot coefficient traces ----
Bmedian = t(apply(B,c(1,2),median));
Bmean = t(apply(B,c(1,2),mean));

# Convert to data frames for easy plotting
colnames(Bmedian) = names(X)
colnames(Bmean) = names(X)
dfBmedian = melt(data.frame(Bmedian, C = Cvalues), id.vars = "C", variable.name = "beta")
dfBmean = melt(data.frame(Bmean, C = Cvalues), id.vars = "C", variable.name = "beta")

# Plot
(ggplot(dfBmedian, aes(x = log2(C), y = value, colour = beta)) +
    geom_point() +
    geom_line() +
    theme(legend.position="none")
)

# ---- Explore best model ----
# Look at best model
idxC = which.min(apply(errTest,2,mean))
bestC = Cvalues[idxC]
print(sprintf('best C: log2(C)=%.3f',log2(Cvalues[idxC])))

# Look at error rates
cat(sprintf('training error: %.12f',mean(errTrain[,idxC])))
cat(sprintf('validation error: %.12f',mean(errTest[,idxC])))

# Look at coefficients for chosen model (mean and median across CV-iterations)
# Sort features according to absolute median weight value for the chosen model
betaMed = abs(Bmedian[idxC,])
betaMean = abs(Bmean[idxC,])
idx = order(betaMed,decreasing = TRUE)
cat(sprintf('\n\nrank\t%35s\tbeta_med\t\tbeta_mean\n','feature'))
for (idx1 in 1:length(idx)){
  cat(sprintf('%d\t%35s\t%.16f\t%.16f\n',idx1,names(betaMean)[idx[idx1]],betaMed[idx[idx1]],betaMean[idx[idx1]]))
}

evaluation_matrix <- function(TP,TN,FP,FN){
true_pr    <- TP / (TP+FP)
false_pr    <- 1-true_pr

true_nr    <- TN / (TN+FN)
false_nr    <- 1-true_nr

sens <- TP/ (TP+FN)
spec <- TN/ (TN+FP)


est <- cbind(true_pr,false_pr,true_nr,false_nr,sens, spec)
return(est)
}

# Plot confusion matrix
ggplot(cMatTrain[[idxC]], aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("training")

m_rates_train <- matrix(c(52,645,200,3),nrow = 2, ncol = 2, dimnames = list(c("imparied", "control"), c("control","impaired")))

ggplot(cMatTest[[idxC]], aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("test")
m_rates_test <- matrix(c(9,69,19,3),nrow = 2, ncol = 2, dimnames = list(c("imparied", "control"), c("control","impaired")))


# ---- Plot number of support vectors vs. regularization strength----
plot_sv <- (ggplot() +
   labs(x = "log2(C)", y = "# support vectors")+
   geom_point(aes(x = log2(Cvalues), y = apply(nSV,2,mean), colour = "blue")) +
   geom_line (aes(x = log2(Cvalues), y = apply(nSV,2,mean), colour = "blue"))
)

