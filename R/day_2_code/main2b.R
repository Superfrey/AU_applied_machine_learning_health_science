# main2b.R
# Logistic regression - CSF biomarker, prediction with a single input
# feature
# Cross-validation
#
# PMR Dec 2022

# ---- Import libraries and functions ----
library("stats")
library("caret")
library("ggplot2")

# ---- Import data etc. ----
data = read.table(here('data-raw/csfBiomarkers.txt'),header = TRUE, sep = ",")

Iresponse = names(data)=='group' # Identify column with response variable (group)
y = data[, Iresponse, drop = FALSE]
y$group = as.factor(y$group); # Convert to factor variable

Ipredictor = names(data)=='tau' # Only use a single predictor variable (tau)
X = data[,Ipredictor, drop = FALSE]


# ---- Create random partition of data ----
set.seed(0) # For reproducibility
K = 10
c = createFolds(y = as.matrix(y), k = K, returnTrain = TRUE) # Returns training indices
idxAll = 1:nrow(y)

## ---- Initialize array  ----
accTrain = accTest = errTrain = errTest = matrix(NaN, K, 1)

# ---- Run analysis ----

# Iterate over partitions
for (idx1 in 1:K){
  I_train = is.element(idxAll, c[[idx1]])
  I_test = !is.element(idxAll, c[[idx1]])
  #
  Xtrain = X[I_train,, drop = FALSE]
  ytrain = y[I_train,, drop = FALSE]
  Xtest = X[I_test,, drop = FALSE]
  ytest = y[I_test,, drop = FALSE]

  # Get group categorical info for later use
  catInfo = levels(ytrain$group)

  # Fit model
  fit = glm(group ~ ., data = cbind(ytrain,Xtrain), family = 'binomial')

  # Predict
  yhatTrainProb = predict(fit, Xtrain, type = "response")
  yhatTestProb = predict(fit, Xtest, type = "response")

  yhatTrain = round(yhatTrainProb)
  yhatTest = round(yhatTestProb)

  ## ---- Make predictions categorical again (instead of 0/1 coding) ----
  yhatTrainCat = factor(yhatTrain, levels = c("0","1"),labels = catInfo)
  yhatTestCat = factor(yhatTest, levels = c("0","1"),labels = catInfo)

  #Evaluate classifier performance
  # Accuracy
  accTrain[idx1] = sum(yhatTrainCat==ytrain$group)/length(ytrain$group);
  accTest[idx1] = sum(yhatTestCat==ytest$group)/length(ytest$group);

  # Error rate
  errTrain[idx1] = 1 - accTrain[idx1];
  errTest[idx1] = 1 - accTest[idx1];

  # Compute confusion matrices
  cmTrain = confusionMatrix(data = yhatTrainCat, reference = ytrain$group, positive = "Impaired")
  cmTest = confusionMatrix(data = yhatTestCat, reference = ytest$group, positive = "Impaired")
  cmTrain = as.data.frame(cmTrain$table)
  cmTest = as.data.frame(cmTest$table)

  if (idx1==1){
    cMatTrain = cmTrain
    cMatTest = cmTest
  } else {
    cMatTrain[,3] = cMatTrain[,3] + cmTrain[,3]
    cMatTest[,3] = cMatTest[,3] + cmTest[,3]
  }
}

# ---- Plot confusion matrix ----
cm_train_plot <- ggplot(cMatTrain, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("training")

cm_test_plot <- ggplot(cMatTest, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("test")


evaluation_matrix <- function(TP,TN,FP,FN){
    true_pr    <- TP / (TP+FP)
    false_pr    <- 1-true_pr

    true_nr    <- TN / (TN+FN)
    false_nr    <- 1-true_nr

    sens <- TP/ (TP+FN)
    spec <- TN/ (TN+FP)


    est <- cbind(true_pr,false_pr,true_nr,false_nr,sens, spec)
    return(est)}


cmTrain
report_prec_train <- evaluation_matrix(cmTrain[4,3],cmTrain[1,3],cmTrain[2,3],cmTrain[3,3])
report_prec_test <- evaluation_matrix(cmTest[4,3],cmTest[1,3],cmTest[2,3],cmTest[3,3])

report_prec_test
report_prec_train
