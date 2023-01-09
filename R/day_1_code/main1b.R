# main1b.R
# Linear regression (polynomial regression) on body density data.
# Cross-validation
#
# PMR Dec 2022

# ---- Import libraries and functions ----
library("here")
library("stats")
library("ggplot2")
library("caret")
source("polyExpand.R")
source(here("R/library_packages.R"))

# ---- Import data ----
data = read.table(here('data-raw/bodyMeasurementsSingleCV.txt'),header = TRUE, sep = ",")

# Extract predictors and response variables from tables into X and y
I = colnames(data)=='Body_Density' # Identify column with response variable (body density)
y = subset(data, select = I)
X = subset(data, select = !I)

# ---- Run analysis ----

# Create random partition of data
set.seed(0) # For reproducibility
K = 10
c = createFolds(y = as.matrix(y), k = K, returnTrain = TRUE) # Returns training indices
idxAll = 1:nrow(y)

# Define polynomial orders
M = 1:7

# Iterate over partitions
errTrain = matrix(NaN, K, length(M))
errTest = errTrain
for (idx1 in 1:K){

    # Get training- and test sets
    I_train = is.element(idxAll, c[[idx1]])
    I_test = !is.element(idxAll, c[[idx1]])
    #
    Xtrain = X[I_train,, drop = FALSE]
    ytrain = y[I_train,, drop = FALSE]
    Xtest = X[I_test,, drop = FALSE]
    ytest = y[I_test,, drop = FALSE]

    # Iterate over polynomial orders
    for (m in M){

        # Create polynomial expansion of predictor variables
        XtrainPol = polyExpand(Xtrain,m);
        XtestPol = polyExpand(Xtest,m);

        # Scale predictor variables
        scl = apply(XtrainPol,2,sd)
        XtrainPol = as.data.frame(scale(XtrainPol, scale = scl, center = FALSE))
        XtestPol = as.data.frame(scale(XtestPol, scale = scl, center = FALSE))

        # Train linear regression model
        fit = glm( Body_Density ~ ., data = cbind(ytrain, XtrainPol), family = gaussian())

        # Use model to predict
        yhatTrain = predict(fit, newdata = XtrainPol)
        yhatTest = predict(fit, newdata = XtestPol)

        # Compute training and test error
        errTrain[idx1,m] = mean(as.matrix((ytrain-yhatTrain)^2))
        errTest[idx1,m] = mean(as.matrix((ytest-yhatTest)^2))
    }
}

# ---- Plot the data, and the model predictions ----
(ggplot() +
   labs(x = "model order", y = "mse")+
   geom_point(aes(x = M, y = apply(errTrain,2,mean), colour = "blue")) +
   geom_line (aes(x = M, y = apply(errTrain,2,mean), colour = "blue")) +
   geom_point(aes(x = M, y = apply(errTest,2,mean), colour = "red")) +
   geom_line (aes(x = M, y = apply(errTest,2,mean), colour = "red")) +
   scale_y_continuous(trans='log10') +
   scale_color_discrete(labels = c("train","test"))
)
