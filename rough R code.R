# Load necessary libraries
library(CAISEr)
library(tidyverse)
library(caret)
library(e1071)
library(randomForest)
library(glmnet)
library(doParallel)


cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


dataset_path <- "~/Downloads/datasets"
dataset_files <- list.files(dataset_path, pattern = "*.csv", full.names = TRUE)


preprocess_data <- function(data) {
  # Remove metadata
  data <- data %>% select(-Info_PepID, -Info_pos)
  
  # Ensure target and grouping are factors
  data$Class <- as.factor(data$Class)
  data$Info_group <- as.factor(data$Info_group)
  
  # Normalize and apply PCA
  predictors <- data %>% select(starts_with("feat"))
  pca_model <- prcomp(predictors, scale. = TRUE)
  
  # Keep enough PCs to explain 90% variance
  var_explained <- summary(pca_model)$importance[3,]
  n_components <- which(cumsum(var_explained) >= 0.9)[1]
  pca_data <- as.data.frame(pca_model$x[, 1:n_components])
  
  # Combine with group and target
  out_data <- bind_cols(pca_data,
                        Class = data$Class,
                        Info_group = data$Info_group)
  return(out_data)
}


get_observations <- function(data, classifier = "rf") {
  folds <- groupKFold(data$Info_group, k = 5)
  accs <- c()
  
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    train_data <- data[train_idx,]
    test_data <- data[-train_idx,]
    
    # Drop Info_group before training
    train_data <- train_data %>% select(-Info_group)
    test_data <- test_data %>% select(-Info_group)
    
    # Ensure Class has the same levels
    train_data$Class <- factor(train_data$Class)
    test_data$Class <- factor(test_data$Class, levels = levels(train_data$Class))
    
    if (classifier == "rf") {
      class_freq <- table(train_data$Class)
      
      if (length(class_freq) < 2) {
        message("Skipping fold with only one class.")
        next
      }
      
      classwt <- as.numeric(sum(class_freq) / (length(class_freq) * class_freq))
      names(classwt) <- names(class_freq)
      
      model <- randomForest(Class ~ ., data = train_data, classwt = classwt)
    } else if (classifier == "logreg") {
      model <- glm(Class ~ ., data = train_data, family = "binomial")
    }
    
    preds <- predict(model, test_data, type = "response")
    if (classifier == "logreg") {
      preds <- as.factor(ifelse(preds > 0.5, levels(train_data$Class)[2], levels(train_data$Class)[1]))
    }
    
    acc <- mean(preds == test_data$Class)
    accs <- c(accs, acc)
  }
  return(accs)
}





# Function to run one experiment on a dataset
run_experiment_on_dataset <- function(file) {
  data <- read.csv(file)
  processed_data <- preprocess_data(data)
  
  # Get observations for both classifiers
  obs_rf <- get_observations(processed_data, "rf")
  obs_lr <- get_observations(processed_data, "logreg")
  
  # Run CAISEr adaptive experiment
  results <- calc_nreps(
    x1 = obs_rf,
    x2 = obs_lr,
    se.max = 0.01,      # precision control
    power = 0.8,
    d = 1,              # minimal relevant effect size
    sig.level = 0.05,
    test.type = "t.test",
    alternative = "two.sided"
  )
  
  list(
    dataset = basename(file),
    n_required = results$n,
    mean_rf = mean(obs_rf),
    mean_lr = mean(obs_lr),
    se = results$se,
    effect_size = results$d,
    p_value = results$p.value
  )
}

library(CAISEr)

# Define function that returns the difference in accuracy
diff_fun <- function() {
  # Draw random sample / fold
  folds <- groupKFold(data$Info_group, k = 5)
  fold_idx <- sample(seq_along(folds), 1)
  
  train_idx <- folds[[fold_idx]]
  train_data <- data[train_idx,] %>% select(-Info_group)
  test_data <- data[-train_idx,] %>% select(-Info_group)
  
  # Make sure class levels are consistent
  train_data$Class <- factor(train_data$Class)
  test_data$Class <- factor(test_data$Class, levels = levels(train_data$Class))
  
  # Random Forest
  rf_class_freq <- table(train_data$Class)
  rf_classwt <- as.numeric(sum(rf_class_freq) / (length(rf_class_freq) * rf_class_freq))
  names(rf_classwt) <- names(rf_class_freq)
  rf_model <- randomForest(Class ~ ., data = train_data, classwt = rf_classwt)
  rf_preds <- predict(rf_model, test_data)
  acc_rf <- mean(rf_preds == test_data$Class)
  
  # Logistic Regression
  logreg_model <- glm(Class ~ ., data = train_data, family = "binomial")
  logreg_preds <- predict(logreg_model, test_data, type = "response")
  logreg_preds <- as.factor(ifelse(logreg_preds > 0.5, levels(train_data$Class)[2], levels(train_data$Class)[1]))
  acc_logreg <- mean(logreg_preds == test_data$Class)
  
  # Return difference in performance
  return(acc_rf - acc_logreg)
}

diff_fun <- function() {
  acc1 <- runif(1, 0.7, 0.9)  # simulate accuracy
  acc2 <- runif(1, 0.65, 0.85)
  return(acc1 - acc2)
}

result <- calc_nreps(
  FUN = diff_fun,
  se.max = 0.01,
  nstart = 5,
  nmax = 1000
)

print(result)

# Run adaptive sampling
result <- calc_nreps(
  FUN = diff_fun,
  se.max = 0.01,
  init.n = 5,
  max.n = 1000
)

print(result)

# Run experiments over all datasets
results_list <- lapply(dataset_files, run_experiment_on_dataset)
results_df <- bind_rows(results_list)

# Save results
write.csv(results_df, "caise_experiment_results.csv", row.names = FALSE)

# Stop parallel backend
stopCluster(cl)
