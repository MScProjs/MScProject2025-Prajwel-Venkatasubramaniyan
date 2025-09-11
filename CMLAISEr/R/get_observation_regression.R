#' Run a Machine Learning algorithm on a dataset.
#'
#' Runs cross-validation multiple times to evaluate a machine learning model's performance on a dataset.
#'
#'  @section References:
#' - F. Campelo
#'    CAISEr: Comparing Algorithms with Iterative Sample-size Estimation in R
#'
#' @param data The data should be a single data frame where:
#'    All feature columns (X) come first.
#'    The last column is the target variable (y), named as target.
#'    See [calc_nreps()] for details.
#'
#'
#' @param model_name Character string specifying the predictive model to fit.
#'   Supported options are:
#'   \describe{
#'     \item{"linear"}{Linear regression model using \code{lm()}, for regression tasks only.}
#'     \item{"random_forest"}{Random Forest model using \code{randomForest()}, supports both regression and classification.}
#'     \item{"svm"}{Support Vector Machine model using \code{svm()} from \code{e1071}, supports both regression and classification.}
#'     \item{"decision_tree"}{Decision tree model using \code{rpart()}, supports both regression and classification.}
#'     \item{"boosting"}{Gradient Boosting model using \code{gbm()}, with distribution set to "gaussian" for regression and "bernoulli" for classification.}
#'   }
#' @details
#'   Additional model-specific parameters can be passed through the \code{model_params} list argument.
#'   For classification tasks, the \code{target} variable must be a factor; for regression, it must be numeric.
#'
#' @param k_foldcv Integer specifying the number of folds to use in k-fold cross-validation.
#'
#' @param n number of observations to generate.
#'
#' @param performance_metric A function to evaluate model predictions.
#'   It should take two arguments: actual values and predicted values, and return a numeric performance score.
#'   If \code{NULL}, defaults to mean squared error for regression and accuracy for classification.
#'
#' @return vector of observed performance values
#'
#' @seealso [calc_nreps_regression()]
#'
#' @author Prajwel Venkatasubramaniyan (\email{prajwel16@gmail.com})
#'
#' @export
#'
#' @examples
#'
#' Regression algorithm
#' # Loading dataset
#' library(MASS)
#' data("Boston")
#'
#' # Feature matrix and target
#' X <- Boston[, -14]
#' y <- Boston$medv
#' boston <- as.data.frame(X)
#' boston$target <- y # Rename the target column to 'target' as expected by the function
#'
#' mse_results <- get_observation_regression(boston, k_foldcv = 5, model_name = "linear", n = 3)
#' print(mse_results)
#'
#' # Using another performance metric
#'
#' r2_metric <- function(actual, predicted) {
#' 1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)}
#' r2_results <- get_observation_regression(boston, k_foldcv = 5, model_name = "random_forest", n = 3, performance_metric = r2_metric)
#' print(r2_results)
#'
#'
# TESTED:OK
get_observation_regression <- function(data, k_foldcv, model_name, n, performance_metric = mse) {
  print(model_name)
  if(n==0){return(c())}
  # make sure that data_target would be the last column and named as target

  results <- numeric(n)

  for (i in 1:n) {
    folds <- createFolds(data$target, k = k_foldcv, list = TRUE, returnTrain = FALSE)
    fold_performances <- numeric(length(folds))

    for (j in seq_along(folds)) {
      test_indices <- folds[[j]]
      train_data <- data[-test_indices, ]
      test_data <- data[test_indices, ]

      # Fit the model
      model <- switch(model_name,
                      "linear" = lm(target ~ ., data = train_data),
                      "random_forest" = randomForest(target ~ ., data = train_data),
                      "svm" = svm(target ~ ., data = train_data),
                      "decision_tree" = rpart(target ~ ., data = train_data),
                      "boosting" = gbm(target ~ ., data = train_data,
                                       distribution = "gaussian",
                                       n.trees = 100, interaction.depth = 3,
                                       shrinkage = 0.1, n.minobsinnode = 10, verbose = FALSE),

                      stop("Unsupported model_name"))

      # Predict on test data
      preds <- predict(model, newdata = test_data)

      # Calculate performance
      fold_performances[j] <- performance_metric(test_data$target, preds)
    }

    results[i] <- mean(fold_performances)
  }

  return(results)
}

library(caret)
library(e1071)
library(randomForest)
library(rpart)
library(gbm)
library(Metrics)
