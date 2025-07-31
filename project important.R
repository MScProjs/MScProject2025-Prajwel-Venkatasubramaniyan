
boot_sdm <- function(x,             # vector of observations
                     boot.R = 999,  # number of bootstrap resamples
                     ncpus  = 1,    # number of cores to use
                     seed   = NULL) # PRNG seed
{
  
  # ========== Error catching ========== #
  assertthat::assert_that(
    is.numeric(x), length(x) > 1,
    assertthat::is.count(boot.R), boot.R > 1,
    assertthat::is.count(ncpus))
  # ==================================== #
  
  # set PRNG seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Perform bootstrap
  if(ncpus > 1){
    x.boot <- parallel::mclapply(1:boot.R,
                                 function(i){
                                   mean(sample(x,
                                               size = length(x),
                                               replace = TRUE))},
                                 mc.cores = ncpus)
  } else {
    x.boot <- lapply(1:boot.R,
                     function(i){
                       mean(sample(x,
                                   size = length(x),
                                   replace = TRUE))})
  }
  return(unlist(x.boot))
}

################################################################################

calc_instances <- function(ncomparisons,               # number of comparisons
                           d,                          # MRES
                           ninstances   = NULL,        # number of instances
                           power        = NULL,        # power
                           sig.level    = 0.05,        # significance level
                           alternative.side = "two.sided", # type of H1
                           test         = "t.test",    # type of test
                           power.target = "mean")      # target power design
{
  
  test    <- match.arg(tolower(test),
                       c("t.test", "wilcoxon", "binomial"))
  alternative.side  <- match.arg(tolower(alternative.side),
                                 c("one.sided", "two.sided"))
  power.target <- match.arg(tolower(power.target),
                            c("worst.case", "mean", "median"))
  
  
  # ========== Error catching ========== #
  assertthat::assert_that(
    assertthat::is.count(ncomparisons),
    is.null(ninstances) || (assertthat::is.count(ninstances) && ninstances > 1),
    is.null(power) || (is.numeric(power) && power > 0 && power < 1),
    is.null(d) || (is.numeric(d) && d > 0),
    sum(c(is.null(ninstances), is.null(power), is.null(d))) == 1,
    is.numeric(sig.level) && sig.level > 0 && sig.level < 1,
    alternative.side %in% c("one.sided", "two.sided"),
    test %in% c("t.test", "wilcoxon", "binomial"),
    power.target %in% c("worst.case", "mean", "median"))
  # ==================================== #
  
  # Calculate correction multiplier depending on test type
  # Based on the ARE of the tests (See Sheskin 1996)
  corr.factor <- switch(test,
                        t.test   = 1,
                        wilcoxon = 1 / 0.86,
                        binomial = 1 / 0.637)
  
  
  if (is.null(ninstances)){ # Estimate sample size
    if (power.target == "mean"){
      # Start by calculating N without any correction:
      N <- stats::power.t.test(delta = d, sd = 1,
                               sig.level = sig.level, power = power,
                               type = "paired",
                               alternative = alternative.side)$n - 1
      
      p.mean <- 0 # mean power
      while (p.mean < power){
        N <- ceiling(N) + 1
        p <- numeric(ncomparisons)
        a <- numeric(ncomparisons)
        for (i in seq_along(p)){
          ss <- stats::power.t.test(delta = d, sd = 1, n = N,
                                    sig.level = sig.level / (ncomparisons - i + 1),
                                    type = "paired",
                                    alternative = alternative.side)
          p[i] <- ss$power
          a[i] <- ss$sig.level
        }
        p.mean <- mean(p)
      }
      p.median <- stats::median(p)
      #
    } else {
      if (power.target == "worst.case") r <- 1 # Bonferroni correction
      if (power.target == "median")     r <- floor(ncomparisons / 2)
      
      N <- stats::power.t.test(delta = d, sd = 1,
                               sig.level = sig.level / (ncomparisons - r + 1),
                               power = power,
                               type = "paired", alternative = alternative.side)$n
      
      # Calculate individual power of each comparison, + mean and median power
      p <- numeric(ncomparisons)
      a <- numeric(ncomparisons)
      for (i in seq_along(p)){
        ss <- stats::power.t.test(delta = d, sd = 1, n = ceiling(N),
                                  sig.level = sig.level / (ncomparisons - i + 1),
                                  type = "paired",
                                  alternative = alternative.side)
        p[i] <- ss$power
        a[i] <- ss$sig.level
      }
      p.mean   <- mean(p)
      p.median <- stats::median(p)
    }
    
    # Adjust sample size depending on the type of test to be performed
    N <- ceiling(N * corr.factor)
    
  } else if (is.null(power)){ # calculate power profile of experiment
    N <- ninstances
    p <- numeric(ncomparisons)
    a <- numeric(ncomparisons)
    for (i in seq_along(p)){
      ss <- stats::power.t.test(delta = d, sd = 1, n = ceiling(N / corr.factor),
                                sig.level = sig.level / (ncomparisons - i + 1),
                                type = "paired",
                                alternative = alternative.side)
      p[i] <- ss$power
      a[i] <- ss$sig.level
    }
    p.mean <- mean(p)
    p.median <- stats::median(p)
  }
  
  
  output <- list(ninstances   = N,
                 power        = p,
                 mean.power   = p.mean,
                 median.power = p.median,
                 d            = d,
                 sig.level    = a,
                 alternative  = alternative.side,
                 test         = test,
                 power.target = power.target)
  
  return(output)
}

################################################################################

se_boot <- function(Xk,                  # vector of observations
                    dif = "simple",      # type of difference
                    comparisons = "all.vs.all", # standard errors to calculate
                    boot.R = 999, # number of bootstrap resamples
                    ...)
{
  
  # ========== Error catching ========== #
  assertthat::assert_that(
    is.list(Xk),
    all(sapply(Xk, is.numeric)),
    all(sapply(Xk, function(x){length(x) >= 2})),
    dif %in% c('simple', 'perc'),
    comparisons %in% c("all.vs.all", "all.vs.first"),
    assertthat::is.count(boot.R), boot.R > 1)
  # ==================================== #
  
  nalgs <- length(Xk)
  Nk    <- sapply(Xk, length)
  
  # Get pairs for calculation
  algo.pairs <- t(utils::combn(1:length(Xk), 2))
  if (comparisons == "all.vs.first") algo.pairs <- algo.pairs[1:(nalgs - 1), , drop = FALSE]
  
  # Calculate point estimates and standard errors for all required pairs
  Phik  <- numeric(nrow(algo.pairs))
  SEk   <- numeric(nrow(algo.pairs))
  Roptk <- numeric(nrow(algo.pairs))
  
  
  for (k in 1:nrow(algo.pairs)){
    ind      <- as.numeric(algo.pairs[k, ])
    phi.hat  <- numeric(boot.R)
    ropt.hat <- numeric(boot.R)
    
    for(i in 1:boot.R){
      # Resample everyone with replacement
      Xk.b <- mapply(FUN = sample,
                     Xk, lapply(Xk, length),
                     MoreArgs = list(replace = TRUE),
                     SIMPLIFY = FALSE)
      
      # Calculate relevant statistics for this bootstrap replicate
      Vark     <- sapply(Xk.b, stats::var)
      Xbark    <- sapply(Xk.b, mean)
      Xbar.all <- mean(Xbark)
      
      if (dif == "simple") {
        # mu1 - mu2
        phi.hat[i] <- Xbark[ind[1]] - Xbark[ind[2]]
        # s1 / s2
        ropt.hat[i] <- sqrt(Vark[ind[1]] / Vark[ind[2]])
        #
      } else if (dif == "perc"){
        if (comparisons == "all.vs.all"){
          # (mu1 - mu2) / mu
          phi.hat[i] <- (Xbark[ind[1]] - Xbark[ind[2]]) / Xbar.all
          # r = s1 / s2
          ropt.hat[i] <- sqrt(Vark[ind[1]] / Vark[ind[2]])
          #
        } else if (comparisons == "all.vs.first"){
          # (mu1 - mu2) / mu1
          phi.hat[i] <- 1 - Xbark[ind[2]] / Xbark[ind[1]]
          # r = (s1 / s2) * (mu2 / mu1)
          ropt.hat[i] <- sqrt(Vark[ind[1]] / Vark[ind[2]]) * (Xbark[ind[2]] / Xbark[ind[1]])
          #
        } else stop("comparisons option *", comparisons, "* not recognized.")
        #
      } else stop ("dif option *", dif, "* not recognized.")
    }
    # Estimate quantities of interest
    Phik[k]  <- mean(phi.hat)
    SEk[k]   <- stats::sd(phi.hat)
    Roptk[k] <- mean(ropt.hat)
  }
  
  # Assemble data frame with results
  output <- data.frame(Alg1 = algo.pairs[, 1],
                       Alg2 = algo.pairs[, 2],
                       N1   = Nk[algo.pairs[, 1]],
                       N2   = Nk[algo.pairs[, 2]],
                       Phi  = Phik,
                       SE   = SEk,
                       r    = Nk[algo.pairs[, 1]] / Nk[algo.pairs[, 2]],
                       ropt = Roptk)
  return(output)
  return(output)
}

################################################################################

se_param <- function(Xk,                  # vector of observations
                     dif = "simple",      # type of difference
                     comparisons = "all.vs.all", # standard errors to calculate
                     ...)
{
  
  # ========== Error catching ========== #
  assertthat::assert_that(
    is.list(Xk),
    all(sapply(Xk, is.numeric)),
    all(sapply(Xk, function(x){length(x) >= 2})),
    dif %in% c('simple', 'perc'),
    comparisons %in% c("all.vs.all", "all.vs.first"))
  # ==================================== #
  
  # Estimates
  nalgs    <- length(Xk)
  Vark     <- sapply(Xk, stats::var)
  Xbark    <- sapply(Xk, mean)
  Nk       <- sapply(Xk, length)
  Xbar.all <- mean(Xbark)
  
  # Get pairs for comparison
  algo.pairs <- t(utils::combn(1:length(Xk), 2))
  if (comparisons == "all.vs.first") algo.pairs <- algo.pairs[1:(nalgs - 1), , drop = FALSE]
  
  # Calculate point estimates and standard errors for all required pairs
  Phik  <- numeric(nrow(algo.pairs))
  SEk   <- numeric(nrow(algo.pairs))
  Roptk <- numeric(nrow(algo.pairs))
  for (i in 1:nrow(algo.pairs)){
    ind <- as.numeric(algo.pairs[i, ])
    if (dif == "simple") {
      # mu1 - mu2
      Phik[i]  <- Xbark[ind[1]] - Xbark[ind[2]]
      # se = s1/sqrt(n1) + s2/sqrt(n2)
      SEk[i]   <- sqrt(sum(Vark[ind] / Nk[ind]))
      # r = s1 / s2
      Roptk[i] <- sqrt(Vark[ind[1]] / Vark[ind[2]])
      #
    } else if (dif == "perc"){
      if (comparisons == "all.vs.all"){
        # (mu1 - mu2) / mu
        Phik[i] <- (Xbark[ind[1]] - Xbark[ind[2]]) / Xbar.all
        # c1 = 1/mu^2 + (mu1 - mu2)^2 / (A * mu^2)^2
        #    = (1 + phi^2 / A^2) / mu^2
        C1 <- (1 + Phik[i] ^ 2 / nalgs ^ 2) / Xbar.all ^ 2
        # c2 = sum_{k!=ind}(s_k^2/n_k^2) * phi^2 / (A^2 * mu^2)
        C2 <- sum(Vark[-ind] / Nk[-ind]) * Phik[i] ^ 2 / (nalgs ^ 2 * Xbar.all ^ 2)
        # se = sqrt(c1 (s1^2/n1 + s2^2/n2) + c2)
        SEk[i] <- sqrt(C1 * (sum(Vark[ind] / Nk[ind])) + C2)
        # r = s1 / s2
        Roptk[i] <- sqrt(Vark[ind[1]] / Vark[ind[2]])
        #
      } else if (comparisons == "all.vs.first"){
        # 1 - mu2/mu1
        Phik[i] <- 1 - Xbark[ind[2]] / Xbark[ind[1]]
        # c1 = s1^2 * (mu_2 / mu_1^2)^2
        C1 <- Vark[ind[1]] * (Xbark[ind[2]] / (Xbark[ind[1]] ^ 2)) ^2
        # c2 = s2^2 / mu_1^2
        C2 <- Vark[ind[2]] / (Xbark[ind[1]] ^ 2)
        # se = sqrt(c1 / n1 + c2 / n2)
        SEk[i] <- sqrt(C1 / Nk[ind[1]] + C2 / Nk[ind[2]])
        # r* = s1/s2 * mu2/mu1
        Roptk[i] <- sqrt(Vark[ind[1]] / Vark[ind[2]]) * (Xbark[ind[2]] / Xbark[ind[1]])
        #
      } else stop("comparisons option *", comparisons, "* not recognized.")
    } else stop ("dif option *", dif, "* not recognized.")
  }
  
  # Assemble data frame with results
  output <- data.frame(Alg1 = algo.pairs[, 1],
                       Alg2 = algo.pairs[, 2],
                       N1   = Nk[algo.pairs[, 1]],
                       N2   = Nk[algo.pairs[, 2]],
                       Phi  = Phik,
                       SE   = SEk,
                       r    = Nk[algo.pairs[, 1]] / Nk[algo.pairs[, 2]],
                       ropt = Roptk)
  return(output)
}

################################################################################

calc_se <- function(Xk,                  # vector of observations
                    dif = "simple",      # type of difference
                    comparisons = "all.vs.all", # standard errors to calculate
                    method = "param",    # method for calculating CI
                    boot.R = 999)        # number of bootstrap resamples
  
{
  
  # ========== Error catching ========== #
  assertthat::assert_that(
    is.list(Xk),
    all(sapply(Xk, is.numeric)),
    all(sapply(Xk, function(x){length(x) >= 2})),
    dif %in% c('simple', 'perc'),
    comparisons %in% c("all.vs.all", "all.vs.first"),
    method %in% c('param', 'boot'),
    assertthat::is.count(boot.R))
  # ==================================== #
  
  # Calculate point estimates and standard errors
  if (method == "param"){
    Diffk <- se_param(Xk = Xk, dif = dif, comparisons = comparisons)
  } else if (method == "boot"){
    Diffk <- se_boot(Xk = Xk, dif = dif, comparisons = comparisons,
                     boot.R = boot.R)
  }
  
  # Fix NaN problem that happens if some variance = 0
  Diffk$SE[is.nan(Diffk$SE)] <- 0
  
  return(Diffk)
}

################################################################################
library(caret)
library(e1071)
library(randomForest)
get_observation_regression <- function(data, k_foldcv, model_name, n, performance_metric = function(actual, predicted) mean((actual - predicted)^2)) {
  print(model_name)
  if(n==0){return(c())}
  # Combine X and y into one data frame
  # data <- as.data.frame(X)
  # data$target <- y
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


# Load data
library(MASS)
data("Boston")

# Feature matrix and target
X <- Boston[, -14]
y <- Boston$medv  # Median value of owner-occupied homes
boston <- as.data.frame(X)
boston$target <- y

install.packages("mlbench")

# Load dataset
library(mlbench)
data("BostonHousing")

# Feature matrix and target
X1 <- BostonHousing[, -14]
y1 <- BostonHousing$medv  # Median value of owner-occupied homes

# Combine into a single data frame
california <- as.data.frame(X1)
california$target <- y1

# Example 1: MSE with linear regression
mse_results <- get_observation_regression(boston, k_foldcv = 5, model_name = "linear", n = 3)
print(mse_results)

# Example 2: R-squared with random forest
r2_metric <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}
r2_results <- get_observation_regression(boston, k_foldcv = 5, model_name = "random_forest", n = 3, performance_metric = r2_metric)
print(r2_results)

##################################################################################

instance <- list(alias = "datainstance", Data = boston)


instances <- list(list(alias = "instance1", Data = boston), list(alias = "instance2", Data = california))


#algorithms <- list(svm="svm",random_forest="random_forest")
algorithms <- c("svm", "random_forest")

calc_nreps <- function(instance,            # instance parameters
                       algorithms,          # algorithm parameters
                       se.max = 0.1,              # desired (max) standard error
                       dif = "simple",      # type of difference
                       comparisons = "all.vs.all", # differences to consider
                       method = "param",    # method ("param", "boot")
                       nstart = 15,         # initial number of samples
                       nmax   = 500,       # maximum allowed sample size
                       seed   = NULL,       # seed for PRNG
                       boot.R = 499,        # number of bootstrap resamples
                       ncpus  = 1,          # number of cores to use
                       force.balanced = FALSE, # force balanced sampling?
                       load.folder = NA,    # folder to load results from
                       save.folder = NA)    # folder to save results to
{
  
  # ========== Error catching ========== #
  # assertthat::assert_that(
  #   is.list(instance),
  #   assertthat::has_name(instance, "alias"),
  #   is.list(algorithms),
  #   all(sapply(X = algorithms, FUN = is.list)),
  #   all(sapply(X = algorithms,
  #              FUN = function(x){assertthat::has_name(x, "FUN")})),
  #   is.numeric(se.max) && length(se.max) == 1,
  #   dif %in% c("simple", "perc"),
  #   comparisons %in% c("all.vs.all", "all.vs.first"),
  #   method %in% c("param", "boot"),
  #   assertthat::is.count(nstart),
  #   is.infinite(nmax) || assertthat::is.count(nmax),
  #   nmax >= length(algorithms) * nstart,
  #   is.null(seed) || seed == seed %/% 1,
  #   assertthat::is.count(boot.R), boot.R > 1,
  #   is.logical(force.balanced), length(force.balanced) == 1,
  #   is.na(save.folder) || (length(save.folder) == 1 && is.character(save.folder)),
  #   is.na(load.folder) || (length(load.folder) == 1 && is.character(load.folder)))
  # # ==================================== #
  
  # set PRNG seed
  if (is.null(seed)) seed <- as.numeric(Sys.time())
  set.seed(seed)
  
  # Set instance alias if needed
  # if (!("alias" %in% names(instance))) {
  #   instance$alias <- instance$FUN
  # }
  
  # Set algorithm aliases if needed
  # for (i in seq_along(algorithms)){
  #   if (!("alias" %in% names(algorithms[[i]]))) {
  #     algorithms[[i]]$alias <- algorithms[[i]]$FUN
  #   }
  # }
  # 
  # Initialize vectors
  Xk <- vector(mode = "list", length = length(algorithms))
  Nk <- numeric(length(algorithms))
  names(Xk) <- algorithms
  names(Nk) <- names(Xk)
  
  # # Load results (if required)
  # if (!is.na(load.folder)){
  #   if(load.folder == "") load.folder <- "./"
  #   # Check that folder exists
  #   if (dir.exists(load.folder)){
  #     filepath <- paste0(normalizePath(load.folder),
  #                        "/", instance$alias, ".rds")
  #     # Check that a file for this instance exists in the folder
  #     if (file.exists(filepath)){
  #       data.in.file  <- readRDS(filepath)
  #       algos.in.file <- names(data.in.file$Nk)
  #       cat("\nExisting data loaded for instance:", instance$alias)
  #       # Extract relevant observations from loaded results
  #       for (i in seq_along(algos.in.file)){
  #         if (algos.in.file[i] %in% names(Xk)){
  #           indx <- which(algos.in.file[i] == names(Xk))
  #           Xk[[indx]] <- data.in.file$Xk[[i]]
  #           Nk[[indx]] <- data.in.file$Nk[[i]]
  #           cat("\n", Nk[[indx]],
  #               "observations retrieved for algorithm:", algos.in.file[i])
  #         }
  #       }
  #     } else {
  #       cat("\nNOTE: Instance file '", filepath, "' not found.")
  #     }
  #   } else {
  #     cat("\nNOTE: folder '", normalizePath(load.folder), "' not found.")
  #   }
  # }
  n.loaded <- Nk
  
  cat("\nSampling algorithms on instance", instance$alias, ": ")
  
  # generate initial samples (if required)
  n0 <- ifelse(rep(force.balanced, length(Nk)),
               yes = max(c(Nk, nstart)) - Nk,
               no  = nstart - pmin(nstart, Nk))
  
  # from here 
  # newX <- parallel::mcmapply(FUN      = get_observation_regression,
  #                            model_name = algorithms,
  #                            data = instance$Data,
  #                            k_foldcv = 2,
  #                            n=n0)
                             
  # newX <- mcmapply(
  #   FUN        = get_observation_regression,
  #   model_name = algorithms,
  #   MoreArgs   = list(data = instance$Data, k_foldcv = 2, n = n0),
  #   mc.cores   = ncups  # or as many cores as needed
  # )
  
  newX <- list()
  k <- 1
  for (algo in algorithms) {
    newX[[algo]] <- get_observation_regression(
      data = instance$Data,
      k_foldcv = 2,
      model_name = algo,
      n = n0[k]
    )
    k <- k+1
  }
  
  #                            algo     = algorithms,
  #                            n        = n0,
  #                            MoreArgs = list(instance = instance),
  #                            mc.cores = ncpus,
  #                            SIMPLIFY = FALSE)
  # boston, k_foldcv = 5, model_name = "linear", n = 3
  
  # Append new observation to each algo list and update sample size counters
  Xk <- mapply(FUN = c, Xk, newX,
               SIMPLIFY = FALSE)
  Nk <- sapply(Xk, length)
  
  # Calculate point estimates, SEs, and sample size ratios (current x optimal)
  Diffk <- calc_se(Xk     = Xk,
                   dif    = dif,
                   comparisons = comparisons,
                   method = method,
                   boot.R = boot.R)
  
  while(any(Diffk$SE > se.max) & (sum(Nk) - sum(n.loaded) < nmax)){
    # Echo something for the user
    if (!(sum(Nk) %% nstart)) cat(".")
    
    # Determine which algorithm(s) should get new observation
    n <- numeric(length(algorithms))
    if(force.balanced){
      ind <- 1:length(algorithms)
    } else {
      # Get pair that has the worst SE
      worst.se <- Diffk[which.max(Diffk$SE), ]
      
      # Determine algorithm from worst.se that should receive a new observation
      if (worst.se$r <= worst.se$ropt){
        ind <- worst.se[1, 1]
      } else {
        ind <- worst.se[1, 2]
      }
    }
    n[ind] <- 1
    
    # Generate new observation(s)
    newX <- list()
    k <- 1
    for (algo in algorithms) {
      newX[[algo]] <- get_observation_regression(
        data = instance$Data,
        k_foldcv = 2,
        model_name = algo,
        n = n[k]
      )
      k <- k+1
    }
    
    # Append new observation(s) and update sample size counters
    # Xk <- mapply(FUN = c, Xk, newX,
    #              SIMPLIFY = FALSE)
    for (model in names(newX)) {
      Xk[[model]] <- c(Xk[[model]], newX[[model]])
    }
    Nk[ind] <- Nk[ind] + 1
    
    # Recalculate point estimates, SEs, and sample size ratios
    Diffk <- calc_se(Xk     = Xk,
                     dif    = dif,
                     comparisons   = comparisons,
                     method = method,
                     boot.R = boot.R)
  }
  
  # Assemble output list
  output    <- list(instance    = instance$alias,
                    Xk          = Xk,
                    Nk          = Nk,
                    n.loaded    = n.loaded,
                    Diffk       = Diffk,
                    dif         = dif,
                    method      = method,
                    comparisons = comparisons,
                    seed        = seed)
  class(output) <- c("nreps", "list")
  
  # Save to file if required
  # if (!is.na(save.folder)){
  #   # Check save folder
  #   if(save.folder == "") save.folder <- "./"
  #   save.folder <- normalizePath(save.folder)
  #   if(!dir.exists(save.folder)) dir.create(save.folder)
  #   
  #   # Prepare save filename
  #   save.file <- paste0(save.folder, "/", instance$alias, ".rds")
  #   
  #   # save output to file
  #   cat("\nWriting file", basename(save.file))
  #   saveRDS(output, file = save.file)
  # }
  
  # Return output
  return(output)
}

result_nreps <- calc_nreps(instance, algorithms)

summary.nreps <- function(object, ...)
{
  # Print summary
  cat("#====================================")
  cat("\nInstance:", object$instance)
  cat("\nNumber of algorithms:", length(object$Nk))
  for (i in seq_along(object$Nk)){
    cat(paste0("\n", names(object$Nk)[i], ": ", object$Nk[i], " runs"))
  }
  cat("\n --------------------")
  cat("\nTotal runs:", sum(object$Nk))
  cat("\nComparisons:", object$comparisons)
  cat("\n#====================================\n\n")
  print(signif(object$Diffk, 3))
  cat("\n#====================================")
  
}
summary.nreps(result_nreps)

plot.nreps <- function(x, y = NULL, ...,
                       instance.name = NULL,
                       latex = FALSE,
                       show.SE = TRUE,
                       show.CI = TRUE,
                       sig.level = 0.05,
                       show.text = TRUE)
{
  
  object <- x
  # Extract a single instance if plotting from CAISEr object
  if ("CAISEr" %in% class(object)){
    assertthat::assert_that(is.character(instance.name),
                            length(instance.name) == 1,
                            instance.name %in% unique(object$data.summary$Instance))
    obj <- list()
    obj$Diffk <- object$data.summary[object$data.summary$Instance == instance.name, ]
    nk <- table(object$data.raw$Algorithm[object$data.raw$Instance == instance.name])
    obj$Nk <- as.numeric(nk)
    names(obj$Nk) <- names(nk)
    obj$instance <- instance.name
    object <- obj
    class(object) <- "nreps"
  }
  
  assertthat::assert_that(all(c("Diffk", "Nk", "instance") %in% names(object)),
                          is.logical(latex), length(latex) == 1,
                          is.logical(show.SE), length(show.SE) == 1,
                          is.logical(show.CI), length(show.CI) == 1,
                          is.logical(show.text), length(show.text) == 1,
                          is.numeric(sig.level), length(sig.level) == 1,
                          sig.level > 0, sig.level < 1,
                          any(c("CAISEr", "nreps") %in% class(object)))
  
  
  df      <- object$Diffk
  algs    <- names(object$Nk)
  df$Alg1 <- algs[df$Alg1]
  df$Alg2 <- algs[df$Alg2]
  df$CIHW <- df$SE * stats::qt(p = 1 - sig.level / 2,
                               df = df$N1 + df$N2)
  
  if (latex){
    pairx <- " $\\times$ "
    ylabtxt <- "$\\phi_{ij}$"
    setxt <- paste0("$SE_{ij} = ", signif(df$SE, 2), "$")
    citxt <- paste0("$CI_{",
                    1 - sig.level,
                    "} = [", signif(df$Phi - df$CIHW, 2), ", ",
                    signif(df$Phi + df$CIHW, 3), "]$")
  } else {
    pairx <- " x "
    ylabtxt <- "diff"
    setxt <- paste0("SE = ", signif(df$SE, 2))
    citxt <- paste0("CI(",
                    1 - sig.level,
                    ") = [", signif(df$Phi - df$CIHW, 2), ", ",
                    signif(df$Phi + df$CIHW, 2), "]")
  }
  
  df$pair <- paste0(df$Alg1, pairx, df$Alg2)
  
  mp <- ggplot2::ggplot(df,
                        ggplot2::aes(x = "pair",
                                            y = "Phi",
                                            ymin = "Phi - SE",
                                            ymax = "Phi + SE")) +
    ggplot2::theme_minimal() +
    ggplot2::geom_abline(slope = 0, intercept = 0,
                         lty = 3, col = "red", lwd = 1.4,
                         alpha = .5)
  if (show.CI){
    mp <- mp +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = "Phi - CIHW",
                                                 ymax = "Phi + CIHW"),
                             alpha = .5, col = 2,
                             width = .12, size = 1.2)
  }
  
  if (show.SE){
    mp <- mp +
      ggplot2::geom_linerange(size = 1.8)
  }
  
  mp <- mp + ggplot2::geom_point(size = 2.5) +
    ggplot2::coord_flip() +
    ggplot2::xlab("Pair") +
    ggplot2::ylab(ylabtxt) +
    ggplot2::labs(caption = paste0("Instance: ", object$instance))
  
  if(show.text & show.SE){
    mp <- mp +
      ggplot2::geom_text(ggplot2::aes(label = setxt),
                         nudge_x = .2, size = 2.5)
  }
  
  if(show.text & show.CI){
    mp <- mp +
      ggplot2::geom_text(ggplot2::aes(label = citxt),
                         nudge_x = -.2, size = 2.5)
  }
  
  return(mp)
}
plot.nreps(result_nreps)

################################################################################

run_experiment <- function(instances, algorithms, d=0.5, se.max=0.1,
                           power = 0.8, sig.level = 0.1,
                           power.target = "mean",
                           dif = "simple", comparisons = "all.vs.all",
                           alternative = "two.sided", test = "t.test",
                           method = "param",
                           nstart = 20, nmax = 100 * length(algorithms),
                           force.balanced = FALSE,
                           ncpus = 2, boot.R = 499, seed = NULL,
                           save.partial.results = NA,
                           load.partial.results = NA,
                           save.final.result    = NA)
{
  
  # ================ Preliminary bureaucracies ================ #
  
  # one-sided tests only make sense for all-vs-one experiments
  if (alternative %in% c("less", "greater")){
    assertthat::assert_that(comparisons == "all.vs.first")
    alternative.side <- "one.sided"
  } else {
    alternative.side <- "two.sided"
  }
  
  # Fix a common mistake
  if (tolower(dif) == "percent") dif <- "perc"
  
  # set PRNG seed
  assertthat::assert_that(is.null(seed) || seed == seed %/% 1)
  if (is.null(seed)) seed <- as.numeric(Sys.time())
  set.seed(seed)
  
  # Capture input parameters to be returned later
  var.input.pars <- as.list(environment())
  
  # Set up parallel processing
  # assertthat::assert_that(assertthat::is.count(ncpus))
  # if ((.Platform$OS.type == "windows") & (ncpus > 1)){
  #   cat("\nAttention: multicore not currently available for Windows.\n
  #       Forcing ncpus = 1.")
  #   ncpus <- 1
  # } else {
  #   available.cores <- parallel::detectCores()
  #   if (ncpus >= available.cores){
  #     cat("\nAttention: ncpus too large, we only have ", available.cores,
  #         " cores.\nUsing ", available.cores - 1,
  #         " cores for run_experiment().")
  #     ncpus <- available.cores - 1
  #   }
  # }
  
  # Fill up instance aliases if needed
  # assertthat::assert_that(is.list(instances), length(instances) > 1)
  # for (i in 1:length(instances)){
  #   if (!("alias" %in% names(instances[[i]]))) {
  #     instances[[i]]$alias <- instances[[i]]$FUN
  #   }
  # }
  
  # Fill up algorithm aliases if needed
  # assertthat::assert_that(is.list(algorithms), length(algorithms) > 1)
  # for (i in 1:length(algorithms)){
  #   if (!("alias" %in% names(algorithms[[i]]))) {
  #     algorithms[[i]]$alias <- algorithms[[i]]$FUN
  #   }
  # }
  
  # ================ Start the actual method ================ #
  
  # Calculate N*
  n.available   <- length(instances)
  n.algs        <- length(algorithms)
  n.comparisons <- switch(comparisons,
                          all.vs.all   = n.algs * (n.algs - 1) / 2,
                          all.vs.first = n.algs - 1)
  
  if (power >= 1) {
    ss.calc <- calc_instances(ncomparisons = n.comparisons,
                              d            = d,
                              ninstances   = n.available,
                              sig.level    = sig.level,
                              alternative.side  = alternative.side,
                              test         = test,
                              power.target = power.target)
    N.star <- n.available
  } else {
    ss.calc <- calc_instances(ncomparisons = n.comparisons,
                              d            = d,
                              power        = power,
                              sig.level    = sig.level,
                              alternative.side  = alternative.side,
                              test         = test,
                              power.target = power.target)
    
    N.star <- ceiling(ss.calc$ninstances)
    if (N.star < n.available){
      # Randomize order of presentation for available instances
      instances <- instances[sample.int(n.available)]
    }
  }
  inst.to.use <- min(N.star, n.available)
  
  # Echo some information for the user
  cat("CAISEr running")
  cat("\n-----------------------------")
  cat("\nRequired number of instances:", N.star)
  cat("\nAvailable number of instances:", n.available)
  cat("\nUsing", ncpus, "cores.")
  cat("\n-----------------------------")
  
  # Sample algorithms on instances
  if(ncpus > 1){
    my.results <- pbmcapply::pbmclapply(X   = instances[1:inst.to.use],
                                        FUN = calc_nreps,
                                        # Arguments for calc_nreps:
                                        algorithms     = algorithms,
                                        se.max         = se.max,
                                        dif            = dif,
                                        comparisons    = comparisons,
                                        method         = method,
                                        nstart         = nstart,
                                        nmax           = nmax,
                                        boot.R         = boot.R,
                                        force.balanced = force.balanced,
                                        save.folder    = save.partial.results,
                                        load.folder    = load.partial.results,
                                        # other arguments for pbmclapply:
                                        mc.cores       = ncpus,
                                        mc.preschedule = FALSE)
  } else {
    my.results <- lapply(X   = instances[1:inst.to.use],
                         FUN = calc_nreps,
                         # Arguments for calc_nreps:
                         algorithms     = algorithms,
                         se.max         = se.max,
                         dif            = dif,
                         comparisons    = comparisons,
                         method         = method,
                         nstart         = nstart,
                         nmax           = nmax,
                         boot.R         = boot.R,
                         force.balanced = force.balanced,
                         save.folder    = save.partial.results,
                         load.folder    = load.partial.results)
  }
  
  # Consolidate raw data
  data.raw <- lapply(X   = my.results,
                     FUN = function(x){
                       inst  <- x$instance
                       nj    <- sum(x$Nk)
                       data.frame(Algorithm = do.call(what = c,
                                                      mapply(rep,
                                                             names(x$Nk),
                                                             x$Nk,
                                                             SIMPLIFY = FALSE)),
                                  Instance    = rep(inst, nj),
                                  Observation = do.call(c, x$Xk))})
  
  data.raw <- do.call(rbind, data.raw)
  rownames(data.raw) <- NULL
  
  # Consolidate summary data
  data.summary <- lapply(X   = my.results,
                         FUN = function(x){
                           cbind(Instance = rep(x$instance, nrow(x$Diffk)),
                                 x$Diffk)})
  
  data.summary <- do.call(rbind, data.summary)
  algonames <- algorithms
  rownames(data.summary) <- NULL
  data.summary$Alg1 <- as.factor(algonames[data.summary$Alg1])
  data.summary$Alg2 <- as.factor(algonames[data.summary$Alg2])
  
  
  # Assemble output
  output <- list(Configuration     = var.input.pars,
                 data.raw          = data.raw,
                 data.summary      = data.summary,
                 N                 = min(N.star, n.available),
                 N.star            = N.star,
                 total.runs        = nrow(data.raw),
                 instances.sampled = unique(data.raw$Instance),
                 Underpowered      = (N.star > n.available),
                 samplesize.calc   = ss.calc)
  
  class(output) <- c("CAISEr", "nreps", "list")
  
  # Save output (if required)
  # assertthat::assert_that(is.na(save.final.result) ||
  #                           (is.character(save.final.result) &
  #                              length(save.final.result) == 1))
  # if(!is.na(save.final.result)){
  #   # Check save folder
  #   if(save.final.result == "") save.final.result <- "./"
  #   save.folder <- normalizePath(save.final.result)
  #   if(!dir.exists(save.folder)) dir.create(save.folder)
  #   
  #   # Prepare save filename
  #   filepath <- paste0(save.folder, "/CAISEr_results_",
  #                      gsub("[- ::alpha::]", "", Sys.time()),
  #                      ".rds")
  #   
  #   # save output to file
  #   cat("\nWriting file", basename(filepath))
  #   saveRDS(output, file = filepath)
  # }
  
  return(output)
}

run_experiment(instances, algorithms)
summary.CAISEr <- function(object, test = NULL,
                           alternative = NULL,
                           sig.level = NULL,
                           ...)
{
  
  # Standard value assignment and error checking
  if (is.null(test)) test <- object$Configuration$test
  if (is.null(alternative)) alternative <- object$Configuration$alternative
  if (is.null(sig.level)) sig.level <- object$Configuration$sig.level
  
  assertthat::assert_that("CAISEr" %in% class(object),
                          is.character(test), length(test) == 1,
                          test %in% c("t.test", "wilcoxon", "binomial"),
                          is.character(alternative), length(alternative) == 1,
                          alternative %in% c("less", "greater", "two.sided"),
                          is.numeric(sig.level), length(sig.level) == 1,
                          sig.level > 0, sig.level < 1)
  
  
  # ===========================================================================
  algonames <- as.character(unique(object$data.raw$Algorithm))
  algoruns  <- as.numeric(table(object$data.raw$Algorithm))
  algopairs <- paste(object$data.summary$Alg1,
                     object$data.summary$Alg2,
                     sep = " x ")
  
  # perform initial tests just to calculate p-values
  # (ignoring significance correction)
  my.tests <- vector(mode = "list", length = length(unique(algopairs)))
  for (i in seq_along(unique(algopairs))){
    tmp <- object$data.summary[algopairs == unique(algopairs)[i], ]
    my.tests[[i]]$comparison <- unique(algopairs)[i]
    my.tests[[i]]$data <- tmp
    my.tests[[i]]$d <- mean(tmp$Phi) / stats::sd(tmp$Phi)
    
    if (test == "t.test"){
      my.tests[[i]]$test <- stats::t.test(tmp$Phi,
                                          conf.level = 1 - sig.level,
                                          alternative = alternative)
      
    } else if (test == "wilcoxon"){
      my.tests[[i]]$test <- stats::wilcox.test(tmp$Phi,
                                               conf.level = 1 - sig.level,
                                               alternative = alternative)
      
    } else if  (test == "binomial"){
      x <- tmp$Phi[tmp$Phi != 0]
      n <- length(x)
      x <- sum(x > 0)
      my.tests[[i]]$test <- stats::binom.test(x, n,
                                              conf.level = 1 - sig.level,
                                              alternative = alternative)
      
    } else stop("Test", test, "not recognised in function summary.CAISEr")
    
    my.tests[[i]]$pval <- my.tests[[i]]$test$p.value
  }
  
  # Reorder the tests in increasing order of p-values
  my.tests <- my.tests[order(sapply(my.tests, function(x) x$pval))]
  
  # Re-evaluate tests based on corrected significance values
  alpha <- object$samplesize.calc$sig.level
  for (i in seq_along(my.tests)){
    if (test == "t.test"){
      my.tests[[i]]$test <- stats::t.test(my.tests[[i]]$data$Phi,
                                          conf.level = 1 - alpha[i],
                                          alternative = alternative)
      
    } else if (test == "wilcoxon"){
      my.tests[[i]]$test <- stats::wilcox.test(my.tests[[i]]$data$Phi,
                                               conf.level = 1 - alpha[i],
                                               alternative = alternative,
                                               conf.int = TRUE)
      
    } else if  (test == "binomial"){
      x <- my.tests[[i]]$data$Phi[my.tests[[i]]$data$Phi != 0]
      n <- length(x)
      x <- sum(x > 0)
      my.tests[[i]]$test <- stats::binom.test(x, n,
                                              conf.level = 1 - alpha[i],
                                              alternative = alternative)
    }
    
    my.tests[[i]]$pval <- my.tests[[i]]$test$p.value
  }
  
  
  # Print summary
  cat("#====================================")
  cat("\n CAISEr object:")
  cat("\n Number of instances sampled:", object$N)
  cat("\n Number of instances required:", object$N.star)
  cat("\n Adequate power:", !object$Underpowered)
  for (i in seq_along(algonames)){
    cat("\n Total runs of", algonames[i], ":", algoruns[i])
  }
  cat("\n#====================================")
  cat("\n Pairwise comparisons of interest:")
  cat("\n Test:", test)
  cat("\n H1:", alternative)
  cat("\n Comparisons:", object$Configuration$comparisons)
  cat("\n Alpha (FWER):", sig.level)
  cat("\n Power target:", object$Configuration$power.target)
  cat("\n Desired power:", object$Configuration$power)
  cat("\n#====================================")
  cat("\nTests using Holm's step-down procedure:")
  stflag <- FALSE
  for (i in seq_along(my.tests)){
    if (!stflag && (my.tests[[i]]$pval > alpha[i])){
      cat("\n\n ----- Stop rejecting H0 at this point -----\n")
      stflag <- TRUE
    } else cat("\n")
    cat("\n Test", i, ":", my.tests[[i]]$comparison)
    cat("\n H0:", switch(test,
                         t.test = "mean = 0",
                         wilcoxon = "median = 0",
                         binomial = "prob = 0.5"))
    cat("\n alpha\t\t=", signif(alpha[i], 4),
        "\n p-value\t=", signif(my.tests[[i]]$pval, 4))
    cat(paste0("\n Est. ", switch(test,
                                  t.test = "mean",
                                  wilcoxon = "median",
                                  binomial = "prob"),
               "\t="), signif(my.tests[[i]]$test$estimate, 4))
    cat("\n CI{1-alpha}\t= [", signif(my.tests[[i]]$test$conf.int, 4), "]")
    cat("\n d\t\t=", my.tests[[i]]$d)
  }
  cat("\n#====================================")
  
  # Return invisibly
  invisible(list(test.info = my.tests,
                 algoruns  = algoruns,
                 algonames = algonames,
                 algopairs = algopairs))
}