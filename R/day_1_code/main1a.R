# main1a.R
# Linear regression (polynomial regression) on body density data.
# Training- and test set.
#
# PMR Dec 2022

# ---- Import libraries and functions ----
library("stats")
library("ggplot2")
library("here")
source(here("R/day_1_code/polyExpand.R"))
source(here("R/library_packages.R"))



# ---- Import data ----
dataTrain = read.table(here('data-raw/bodyMeasurementsSingleTrain.txt'),header = TRUE, sep = ",")
dataTest = read.table(here('data-raw/bodyMeasurementsSingleTest.txt'), header = TRUE, sep = ",")

# Extract predictors and response variables from tables into X and y
I = colnames(dataTrain)=='Body_Density' # Identify column with response variable (body density)
ytrain = dataTrain[,I, drop = FALSE]
Xtrain = dataTrain[,!I, drop = FALSE]

I = colnames(dataTest)=='Body_Density'; # Identify column with response variable (body density)
ytest = dataTest[,I, drop = FALSE]
Xtest = dataTest[,!I, drop = FALSE]

# ---- Run analysis ----

# Create an x axis with high resolution for plotting the trained model
XHighRes = data.frame(seq(min(rbind(Xtrain,Xtest)), max(rbind(Xtrain,Xtest)),length.out = 1000))
colnames(XHighRes) = colnames(Xtrain)

# Define polynomial order
m = 1

# Create polynomial expansion of predictor variables
XtrainPol = polyExpand(Xtrain,m)
XtestPol = polyExpand(Xtest,m)

# Also create polynomial expansion of data used for plotting the trained model
XhighResPol = polyExpand(XHighRes,m)

# Scale predictor variables
scl = apply(XtrainPol,2,sd)

XtrainPol = as.data.frame(scale(XtrainPol, scale = scl, center = FALSE))
XtestPol = as.data.frame(scale(XtestPol, scale = scl, center = FALSE))
XhighResPol = as.data.frame(scale(XhighResPol, scale = scl, center = FALSE))

# Train linear regression model
head(ytrain)
head(XtrainPol)
fit = glm( Body_Density ~ ., data = cbind(ytrain, XtrainPol), family = gaussian())
summary(fit)

# Use model to predict
yhatTrain = predict(fit, newdata = XtrainPol)
yhatTest = predict(fit, newdata = XtestPol)

# Compute training and test error
errTrain = mean(data.matrix((ytrain-yhatTrain)^2))
errTest = mean(data.matrix((ytest-yhatTest)^2))
errTrain
errTest

# ---- Plot the data, and the model predictions ----
(ggplot() +
  labs(x = colnames(dataTrain)[!I], y = colnames(dataTrain)[I])+
  geom_point(aes(x = XtrainPol[,1], y = ytrain[,1])) +
  geom_line(aes(x = XhighResPol[,1], y = predict(fit, newdata = data.frame(XhighResPol)))) +
  geom_point(aes(x = XtestPol[,1], y = ytest[,1]), colour = "red")
)
################################################# Loop function ###################

pol <- 1:7
results <- for (i in pol) {


    # Create an x axis with high resolution for plotting the trained model
    XHighRes = data.frame(seq(min(rbind(Xtrain,Xtest)), max(rbind(Xtrain,Xtest)),length.out = 1000))
colnames(XHighRes) = colnames(Xtrain)

# Define polynomial order
m = i

# Create polynomial expansion of predictor variables
XtrainPol = polyExpand(Xtrain,m)
XtestPol = polyExpand(Xtest,m)

# Also create polynomial expansion of data used for plotting the trained model
XhighResPol = polyExpand(XHighRes,m)

# Scale predictor variables
scl = apply(XtrainPol,2,sd)

XtrainPol = as.data.frame(scale(XtrainPol, scale = scl, center = FALSE))
XtestPol = as.data.frame(scale(XtestPol, scale = scl, center = FALSE))
XhighResPol = as.data.frame(scale(XhighResPol, scale = scl, center = FALSE))

# Train linear regression model
head(ytrain)
head(XtrainPol)
fit = glm( Body_Density ~ ., data = cbind(ytrain, XtrainPol), family = gaussian())
summary(fit)

# Use model to predict
yhatTrain = predict(fit, newdata = XtrainPol)
yhatTest = predict(fit, newdata = XtestPol)


# Compute training and test error
errTrain = mean(data.matrix((ytrain-yhatTrain)^2))
errTest = mean(data.matrix((ytest-yhatTest)^2))
MSE <- cbind(m, errTrain,errTest)
# ---- Plot the data, and the model predictions ----
plots_lm <- (ggplot() +
     labs(x = colnames(dataTrain)[!I], y = colnames(dataTrain)[I])+
     geom_point(aes(x = XtrainPol[,1], y = ytrain[,1])) +
     geom_line(aes(x = XhighResPol[,1], y = predict(fit, newdata = data.frame(XhighResPol)))) +
     geom_point(aes(x = XtestPol[,1], y = ytest[,1]), colour = "red")+
    labs(title = i)
)

plot(plots_lm)
print(MSE)
}


#include multiple plot function from loop

datalist = list()

results <- for (i in pol) {


    # Create an x axis with high resolution for plotting the trained model
    XHighRes = data.frame(seq(min(rbind(Xtrain,Xtest)), max(rbind(Xtrain,Xtest)),length.out = 1000))
    colnames(XHighRes) = colnames(Xtrain)

    # Define polynomial order
    m = i

    # Create polynomial expansion of predictor variables
    XtrainPol = polyExpand(Xtrain,m)
    XtestPol = polyExpand(Xtest,m)

    # Also create polynomial expansion of data used for plotting the trained model
    XhighResPol = polyExpand(XHighRes,m)

    # Scale predictor variables
    scl = apply(XtrainPol,2,sd)

    XtrainPol = as.data.frame(scale(XtrainPol, scale = scl, center = FALSE))
    XtestPol = as.data.frame(scale(XtestPol, scale = scl, center = FALSE))
    XhighResPol = as.data.frame(scale(XhighResPol, scale = scl, center = FALSE))

    # Train linear regression model
    head(ytrain)
    head(XtrainPol)
    fit = glm( Body_Density ~ ., data = cbind(ytrain, XtrainPol), family = gaussian())
    summary(fit)

    # Use model to predict
    yhatTrain = predict(fit, newdata = XtrainPol)
    yhatTest = predict(fit, newdata = XtestPol)


    # Compute training and test error
    errTrain = mean(data.matrix((ytrain-yhatTrain)^2))
    errTest = mean(data.matrix((ytest-yhatTest)^2))
    MSE <- cbind(m, errTrain,errTest)
    MSE <- data.frame(MSE)
    MSE$i <- i
    datalist[[i]] <- MSE
}

MSE_poly <- do.call(rbind, datalist)

MSE_poly <- MSE_poly %>%
    mutate(ratioError=  errTest / errTrain,
           logerrTest = log(errTest),
           logerrTrain = log(errTrain))


MSE_error_plot <- MSE_poly %>%
    ggplot(aes(m,ratioError), label=Name)+
    geom_point(colour = "green") +
    geom_point(aes(m,logerrTest), colour= "red")+
    geom_point(aes(m,logerrTrain), colour= "blue") +
    ylim(c(-1,0.5))



# Set up the plot
bv_to<- plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="Flexibility", ylab="Error")

# Generate fake data for the curves
x <- seq(0, 1, length.out=100)
bias <- x^2
variance <- (1 - x)^2
training_error <- bias + variance
test_error <- training_error + 0.1*x
bayes <- 0.25*x

# Add the curves to the plot
lines(x, bias, col="red")
lines(x, variance, col="blue")
lines(x, training_error, col="green")
lines(x, test_error, col="orange")
lines(x, bayes, col="purple")

# Add a legend
legend("topright", c("Bias", "Variance", "Training Error", "Test Error", "Bayes Error"), col=c("red", "blue", "green", "orange", "purple"), lty=1)
