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
#' @param task Character specifying the type of learning task:
#' either \code{"regression"} for predicting continuous outcomes or
#' \code{"classification"} for predicting categorical classes.
#'
#' @param model_name Character string specifying the predictive model to fit.
#'   Supported options are:
#'   \describe{
#'     \item{"linear"}{Linear regression model using \code{lm()}, for regression tasks only.}
#'     \item{"logistic"}{Logistic regression model using \code{glm()} with binomial family, for classification tasks only.}
#'     \item{"random_forest"}{Random Forest model using \code{randomForest()}, supports both regression and classification.}
#'     \item{"svm"}{Support Vector Machine model using \code{svm()} from \code{e1071}, supports both regression and classification.}
#'     \item{"decision_tree"}{Decision tree model using \code{rpart()}, supports both regression and classification.}
#'     \item{"knn"}{k-Nearest Neighbors using \code{knn.reg()} (regression) or \code{knn()} (classification), supports both tasks.}
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
#' @seealso [calc_nreps()]
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
#' result <- get_observation(data = boston,
#'                           task = "regression",
#'                           k_foldcv = 5,
#'                           model_name = "random_forest",
#'                           n = 3,
#'                           model_params = list(ntree = 200))
#' print(result)
#'
#' Example for Classification Algorithm
#' data(iris)
#' iris_data <- iris
#' names(iris_data)[which(names(iris_data) == "Species")] <- "target"
#' iris_data$target <- as.factor(iris_data$target) # Ensure target is a factor for classification
#'
#' result <- get_observation(data = iris_data,
#'                           task = "classification",
#'                           k_foldcv = 5,
#'                           model_name = "svm",
#'                           n = 5,
#'                           model_params = list(kernel = "radial"))  # additional model parameters
#'
#' print(result)
#'
#' result <- get_observation(data = iris_data,
#'                           task = "classification",
#'                           k_foldcv = 5,
#'                           model_name = "svm",
#'                           n = 5,
#'                           model_params = list(kernel = "radial"),
#'                           performance_metric = f1)  # using different performance metric
#'
#' print(result)
#'
#'
# TESTED
get_performance <- function(data,
                            task = c("regression", "classification"),
                            model_name,
                            k_foldcv = 5,
                            n = 1,
                            model_params = list(),
                            performance_metric = NULL) {

  task <- match.arg(task)
  # print(paste("Model:", model_name, "| Task:", task))

  if (n == 0) return(c())

  # Default performance metric
  if (is.null(performance_metric)) {
    performance_metric <- if (task == "regression") {
      function(actual, predicted) mean((actual - predicted)^2)  # MSE
    } else {
      function(actual, predicted) mean(actual == predicted)     # Accuracy
    }
  }

  results <- numeric(n)

  for (i in 1:n) {
    # performing cross validation
    folds <- createFolds(data$target, k = k_foldcv, list = TRUE, returnTrain = FALSE)
    fold_performances <- numeric(length(folds))

    for (j in seq_along(folds)) {
      # splitting the dataset into train and test set
      test_indices <- folds[[j]]
      train_data <- data[-test_indices, ]
      test_data <- data[test_indices, ]

      # fitting the model and predicting on test set
      preds <- switch(model_name,

                      "linear" = {
                        if (task != "regression") stop("Linear model only supports regression.")
                        model <- do.call(lm, c(list(formula = target ~ ., data = train_data), model_params))
                        predict(model, newdata = test_data)
                      },

                      "logistic" = {
                        if (task != "classification") stop("Logistic model only supports classification.")
                        model <- do.call(glm, c(list(formula = target ~ ., data = train_data, family = binomial), model_params))
                        pred_probs <- predict(model, newdata = test_data, type = "response")
                        as.factor(ifelse(pred_probs > 0.5, levels(train_data$target)[2], levels(train_data$target)[1]))
                      },

                      "random_forest" = {
                        model <- do.call(randomForest, c(list(formula = target ~ ., data = train_data), model_params))
                        predict(model, newdata = test_data)
                      },

                      "svm" = {
                        model <- do.call(svm, c(list(formula = target ~ ., data = train_data), model_params))
                        predict(model, newdata = test_data)
                      },

                      "decision_tree" = {
                        model <- do.call(rpart, c(list(formula = target ~ ., data = train_data), model_params))
                        predict(model, newdata = test_data, type = ifelse(task == "classification", "class", "vector"))
                      },

                      "knn" = {
                        train_x <- train_data[, setdiff(names(train_data), "target")]
                        test_x <- test_data[, setdiff(names(test_data), "target")]
                        train_y <- train_data$target
                        if (task == "regression") {
                          preds_knn <- do.call(knn.reg, c(list(train = train_x, test = test_x, y = train_y), model_params))
                          preds_knn$pred
                        } else {
                          preds_knn <- do.call(knn, c(list(train = train_x, test = test_x, cl = train_y), model_params))
                          preds_knn
                        }
                      },

                      "boosting" = {
                        distribution_type <- ifelse(task == "regression", "gaussian", "bernoulli")
                        model <- do.call(gbm, c(list(formula = target ~ ., data = train_data,
                                                     distribution = distribution_type,
                                                     n.trees = 100, interaction.depth = 3,
                                                     shrinkage = 0.1, n.minobsinnode = 10, verbose = FALSE),
                                                model_params))
                        if (task == "regression") {
                          predict(model, newdata = test_data, n.trees = model$n.trees)
                        } else {
                          probs <- predict(model, newdata = test_data, n.trees = model$n.trees, type = "response")
                          as.factor(ifelse(probs > 0.5, levels(train_data$target)[2], levels(train_data$target)[1]))
                        }
                      },

                      stop("Unsupported model_name")
      )
      # calculating performance metric
      fold_performances[j] <- performance_metric(test_data$target, preds)
    }

    # considering the average of performances
    results[i] <- mean(fold_performances)
  }

  return(results)
}

library(caret)
library(e1071)
library(randomForest)
library(Metrics)
library(rpart)
library(FNN)
library(gbm)
