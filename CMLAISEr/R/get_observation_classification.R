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
#' @seealso [calc_nreps_classification()]
#' @seealso [calc_nreps_regression()]
#'
#' @author Prajwel Venkatasubramaniyan (\email{prajwel16@gmail.com})
#'
#' @export
#'
#'
#'
#'
# TESTED:OK

get_observation_classification <- function(data, k_foldcv, model_name, n, performance_metric = f1) {
  print(model_name)
  if (n == 0) { return(c()) }

  # make sure target is a factor
  data$target <- as.factor(data$target)

  results <- numeric(n)

  for (i in 1:n) {
    folds <- createFolds(data$Info_group, k = k_foldcv, list = TRUE, returnTrain = FALSE)
    fold_performances <- numeric(length(folds))

    for (j in seq_along(folds)) {
      test_indices <- folds[[j]]
      train_data <- data[-test_indices, ]
      test_data <- data[test_indices, ]

      # Fit the model
      model <- switch(model_name,
                      "random_forest" = randomForest(target ~ ., data = train_data),
                      "svm" = svm(target ~ ., data = train_data, probability = TRUE),
                      "decision_tree" = rpart(target ~ ., data = train_data, method = "class"),

                      stop("Unsupported model_name"))

      # ---- Prediction step (clean version) ----
      # if (model_name %in% c("random_forest", "decision_tree", "svm")) {
      #   # these models directly return class labels
      #   preds <- predict(model, newdata = test_data)
      #
      # } else if (model_name == "logistic") {
      #   # multinom from nnet returns classes when type="class"
      #   preds <- predict(model, newdata = test_data, type = "class")
      #
      # } else if (model_name == "boosting") {
      #   # gbm returns probabilities by default
      #   prob_preds <- predict(model, newdata = test_data, n.trees = 100, type = "response")
      #   preds <- colnames(prob_preds)[apply(prob_preds, 1, which.max)]
      # }
      #
      # # Make sure preds is a factor with the same levels as target
      # preds <- factor(preds, levels = levels(test_data$target))


      # ---- Predict class labels ----
      if (model_name %in% c("random_forest", "decision_tree", "svm")) {
        # these already return class labels
        preds <- predict(model, newdata = test_data)

      }

      # ensure predictions are factors with correct levels
      preds <- factor(preds, levels = levels(test_data$target))


      # Calculate F1 score (from Metrics package)
      fold_performances[j] <- accuracy(as.numeric(test_data$target), as.numeric(preds))
    }

    results[i] <- mean(fold_performances)
  }

  return(results)
}

library(randomForest)
library(e1071)
library(rpart)

library(caret)
library(e1071)
library(randomForest)
library(Metrics)
