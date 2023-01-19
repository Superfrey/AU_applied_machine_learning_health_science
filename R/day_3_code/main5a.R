# main5a.R
# Feed-forward neural network on CSF biomarker data.
# Cross-validation
#
# PMR Dec 2022

# ---- Import libraries and functions ----
library("stats")
library("ggplot2")
library("caret")
library("torch")
library("luz")
library("reshape2")
library("here")

# ---- Import data etc. ----
data = read.table(here('data-raw/csfBiomarkers.txt'),header = TRUE, sep = ",")

# Extract predictors and response variables from tables into X and y
Iresponse = names(data)=='group' # Identify column with response variable (group)
y = data[, Iresponse, drop = FALSE]
y$group = as.factor(y$group); # Convert to factor variable
X = data[,!Iresponse, drop = FALSE]

# Create numeric version of y
yCat = y
yLevels = levels(yCat$group)
y = matrix(NaN,nrow(y),1)
for (idx1 in 1:length(yLevels)){
  y[yCat$group == yLevels[idx1]] = idx1 - 1
}

p = ncol(X)

# ---- Create random partition of data ----
set.seed(0) # For reproducibility
K = 10
c = createFolds(y = as.matrix(y), k = K, list = TRUE, returnTrain = TRUE) # Returns training indices
idxAll = 1:nrow(y)

# ---- Configure neural network ----
nnModel <- nn_module(
  initialize = function(inputSize = p) {
    self$linear1 <- nn_linear(inputSize, 8) # define the first hidden layer
    self$linear2 <- nn_linear(8,8) # second hidden layer
    self$output <- nn_linear(8,1) # output layer
    self$activation <- nn_relu() # final prediction
  },
  forward = function(x) {
    x %>%
      self$linear1() %>%
      self$activation() %>%
      self$linear2() %>%
      self$activation() %>%
      self$output()
  }
)

# Bias node is added automatically
nnModel <- nnModel %>%
  setup(
    loss = nn_bce_with_logits_loss(),
    optimizer = optim_adam,  # op
    metrics = list(luz_metric_binary_accuracy_with_logits() #meassurer error
  )
  )

# NN options
nEpochs = 50 # Number of epochs in NN fitting
validFrac = 0.1 # Fraction of training data to use as validation during NN fitting
learningRate = 0.01 # Learning rate in NN fitting

# ---- Define sequence for regularization strengths ----
lambdas = 2^seq(-12, 5, by = 1)

# ---- Run analysis ----

# Initialize array to store coefficients across CV-iterations
B = array(NaN, dim = c(p, length(lambdas),K))
errTrain = errTest = accTrain = accTest = nSV = matrix(NaN, K, length(lambdas))
cMatTrain = cMatTest = list()

# Iterate over partitions
for (idx1 in 1:K) {

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

  # Iterate over regularization strengths fit the svm, and to compute
  # training- and test errors for individual regularization strengths.
  for (idx2 in 1:length(lambdas)){

    # Print status
    cat(sprintf('\n\ncross-validation iteration %d of %d, inner iteration %d of %d\n\n',idx1,K, idx2, length(lambdas)))
    Sys.sleep(0.5)

    # Assemble training data
    dataTrain = list(
      torch_tensor(as.matrix(Xtrain)),
      torch_tensor(as.matrix(ytrain))
    )

    # Update regularization strength
    nnModelConfigured = nnModel %>%
      set_opt_hparams(lr = learningRate, weight_decay = lambdas[idx2])

    # Fit NN
    fitted <- nnModelConfigured %>% fit(dataTrain, valid_data = validFrac, epochs = nEpochs)

    # Plot fitting errors vs. epoch
    plot(fitted)

    # Predict
    yhatTrainProb = as.matrix(torch_sigmoid(predict(fitted, torch_tensor(as.matrix(Xtrain), dtype = torch_float()))))
    yhatTestProb = as.matrix(torch_sigmoid(predict(fitted, torch_tensor(as.matrix(Xtest)))))

    yhatTrain = round(yhatTrainProb)
    yhatTest = round(yhatTestProb)

    ## ---- Make predictions categorical again (instead of 0/1 coding) ----
    ytrainCat = factor(ytrain, levels = c(0:(length(yLevels)-1)),labels = yLevels)
    ytestCat = factor(ytest, levels = c(0:(length(yLevels)-1)),labels = yLevels)
    #
    yhatTrainCat = factor(yhatTrain, levels = c(0:(length(yLevels)-1)),labels = yLevels)
    yhatTestCat = factor(yhatTest, levels = c(0:(length(yLevels)-1)),labels = yLevels)

    # Evaluate classifier performance
    # Accuracy
    accTrain[idx1,idx2] = sum(yhatTrainCat==ytrainCat)/length(ytrainCat);
    accTest[idx1,idx2] = sum(yhatTestCat==ytestCat)/length(ytestCat);

    # Error rate
    errTrain[idx1,idx2] = 1 - accTrain[idx1,idx2];
    errTest[idx1,idx2] = 1 - accTest[idx1,idx2];

    # Compute confusion matrices
    cmTrain = confusionMatrix(data = yhatTrainCat, reference = ytrainCat, positive = "Impaired")
    cmTest = confusionMatrix(data = yhatTestCat, reference = ytestCat, positive = "Impaired")
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
error_plot<- (ggplot() +
   labs(x = "log2(lambda)", y = "error rate")+
   geom_point(aes(x = log2(lambdas), y = apply(errTrain,2,mean), colour = "blue")) +
   geom_line (aes(x = log2(lambdas), y = apply(errTrain,2,mean), colour = "blue")) +
   geom_point(aes(x = log2(lambdas), y = apply(errTest,2,mean), colour = "red")) +
   geom_line (aes(x = log2(lambdas), y = apply(errTest,2,mean), colour = "red")) +
   scale_color_discrete(labels = c("train","test"))
)

# ---- Explore best model ----
# Look at best model
idxLambda = which.min(apply(errTest,2,mean))
bestLambda = lambdas[idxLambda]
print(sprintf('best lambda: log2(lambda)=%.3f',log2(lambdas[idxLambda])))

# Look at error rates
cat(sprintf('training error: %.12f',mean(errTrain[,idxLambda])))
cat(sprintf('validation error: %.12f',mean(errTest[,idxLambda])))


# Plot confusion matrix
ggplot(cMatTrain[[idxLambda]], aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("training")

ggplot(cMatTest[[idxLambda]], aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="darkolivegreen", high="burlywood2") +
  ggtitle("test")

m_rates_train <- matrix(c(12,644,240,4),nrow = 2, ncol = 2, dimnames = list(c("imparied", "control"), c("control","impaired")))
m_rates_test <- matrix(c(10,69,18,3),nrow = 2, ncol = 2, dimnames = list(c("imparied", "control"), c("control","impaired")))

