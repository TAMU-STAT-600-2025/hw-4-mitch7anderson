# Load the riboflavin data

# Uncomment below to install hdi package if you don't have it already; 
install.packages("hdi") 
library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene expression
dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
?riboflavin # this gives you more information on the dataset

# This is to make sure riboflavin$x can be converted and treated as matrix for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]


# Get matrix X and response vector Y
X = as.matrix(riboflavin$x)
Y = riboflavin$y

# Source your lasso functions
source("LassoFunctions.R")

# [ToDo] Use your fitLASSO function on the riboflavin data with 60 tuning parameters
lasso_fit <- fitLASSO(X = X, Y = Y, n_lambda = 60)

# [ToDo] Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter
nonzero_counts <- colSums(lasso_fit$beta_mat != 0)

plot(lasso_fit$lambda_seq, nonzero_counts, type = "b", log = "x",
     xlab = expression(lambda),
     ylab = "Number of nonzero coefficients",
     main = "LASSO Path on Riboflavin Data")

# [ToDo] Use microbenchmark 10 times to check the timing of your fitLASSO function above with 60 tuning parameters
library(microbenchmark)
microbenchmark(fitLASSO(X = X, Y = Y, n_lambda = 60), times = 10, unit = "s")

# [ToDo] Report your median timing in the comments here: (~5.8 sec for Irina on her laptop)
# Median runtime: 2.25 seconds

# [ToDo] Use cvLASSO function on the riboflavin data with 30 tuning parameters (just 30 to make it faster)
cvlasso_fit <- cvLASSO(X = X, Y = Y, n_lambda = 30) 

# [ToDo] Based on the above output, plot the value of CV(lambda) versus tuning parameter. Note that this will change with each run since the folds are random, this is ok.
plot(cvlasso_fit$lambda_seq, cvlasso_fit$cvm, type = "b", log = "x",
     xlab = expression(lambda),
     ylab = "CV MSE",
     main = "Cross-validated LASSO on Riboflavin")
