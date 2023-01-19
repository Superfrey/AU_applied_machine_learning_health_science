# main2a.R
# Logistic regression - CSF biomarker, prediction with a single input
# feature.
# Training- and test set.
#
# PMR Dec 2022

# ---- Import libraries and functions ----
library("stats")
library("caret")
library("ggplot2")
library("here")

# ----- Import data ----
data = read.table(here('data-raw/csfBiomarkers.txt'),header = TRUE, sep = ",")

Iresponse = names(data)=='group' # Identify column with response variable (group)
y = data[, Iresponse, drop = FALSE]
y$group = as.factor(y$group); # Convert to factor variable

Ipredictor = names(data)=='tau' # Only use a single predictor variable (tau)
X = data[,Ipredictor, drop = FALSE]

# Make high resolution x axis for plotting
xHighRes = data.frame(seq(0.5*min(X), 1.5*max(X), length.out = 1000))
names(xHighRes) = names(X)

# ---- Divide into training and test sets ----
set.seed(0) #For reproducibility
c = createDataPartition(y = y$group, p = 0.9) ## Cross-validation in Logistic regression

idxAll = 1:nrow(y)

I_train = is.element(idxAll, c[[1]])
I_test = !is.element(idxAll, c[[1]])
#
Xtrain = X[I_train,, drop = FALSE]
ytrain = y[I_train,, drop = FALSE]
Xtest = X[I_test,, drop = FALSE]
ytest = y[I_test,, drop = FALSE]

# ---- Train model, predict, and plot model ----
# Get group categorical info for later use
catInfo = levels(ytrain$group)

# Fit model
fit = glm(group ~ ., data = cbind(ytrain,Xtrain), family = 'binomial')

# Predict
yhatTrainProb = predict(fit, Xtrain, type = "response")
yhatTestProb = predict(fit, Xtest, type = "response")

yhatTrain = round(yhatTrainProb)
yhatTest = round(yhatTestProb)

# Plot model output
plot_response <- (ggplot() +
   labs(x = colnames(X), y = paste(c("P(Y=Impaired|", colnames(X), ")"), collapse = ""))+
   geom_line(aes(x = xHighRes[,1], y = predict(fit, newdata = data.frame(xHighRes), type = "response")))
)

## ---- Make predictions categorical again (instead of 0/1 coding) ----
yhatTrainCat = factor(yhatTrain, levels = c("0","1"),labels = catInfo)
yhatTestCat = factor(yhatTest, levels = c("0","1"),labels = catInfo)

# Plot confusion matrix
cmTrain = confusionMatrix(data = yhatTrainCat, reference = ytrain$group, positive = "Impaired")
cmTest = confusionMatrix(data = yhatTestCat, reference = ytest$group, positive = "Impaired")
cmTrain = as.data.frame(cmTrain$table)
cmTest = as.data.frame(cmTest$table)

ggplot(cmTrain, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("training")

ggplot(cmTest, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("test")

## ---- Evaluate classifier performance ----
# Accuracy
accTrain = sum(yhatTrainCat==ytrain$group)/length(ytrain$group);
accTest = sum(yhatTestCat==ytest$group)/length(ytest$group);
accTrain
accTest
# Error rate
errTrain = 1 - accTrain;
errTest = 1 - accTest;
errTest
