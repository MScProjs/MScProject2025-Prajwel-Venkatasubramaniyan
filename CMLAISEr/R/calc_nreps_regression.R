#' Determine sample sizes for a set of ML algorithms on a single data instance
#'
#' Iteratively calculates the required sample sizes for K ML algorithms
#' on a given problem instance, so that the standard errors of the estimates of
#' the pairwise differences in performance is controlled at a predefined level.
#'
#' @section Instance:
#' @param instance A list containing information about the problem instance to evaluate.
#'   Must include:
#'   \describe{
#'     \item{\code{alias}}{A character string naming the instance (used for saving/loading results).}
#'     \item{\code{Data}}{A data frame with predictors and a target column named \code{target}.}
#'   }
#' @section Algorithms:
#' @param algorithms A character vector of model names to compare.
#' In what follows `algorithm[[k]]`
#' refers to any algorithm specified in the `algorithms` list.
#'
#' @section Initial Number of Observations:
#' In the **general case** the initial number of observations per algorithm
#' (`nstart`) should be relatively high. For the parametric case
#' we recommend between 10 and 20 if outliers are not expected, or between 30
#' and 50 if that assumption cannot be made. For the bootstrap approach we
#' recommend using at least 20. However, if some distributional assumptions can
#' be made - particularly low skewness of the population of algorithm results on
#' the test instances), then `nstart` can in principle be as small as 5 (if the
#' output of the algorithms were known to be normal, it could be 1).
#'
#' In general, higher sample sizes are the price to pay for abandoning
#' distributional assumptions. Use lower values of `nstart` with caution.
#'
#' @section Pairwise Differences:
#' Parameter `dif` informs the type of difference in performance to be used
#' for the estimation (\eqn{\mu_a} and \eqn{\mu_b} represent the mean
#' performance of any two algorithms on the test instance, and \eqn{mu}
#' represents the grand mean of all algorithms given in `algorithms`):
#'
#' - If `dif == "perc"` and `comparisons == "all.vs.first"`, the estimated quantity is
#'    \eqn{\phi_{1b} = (\mu_1 - \mu_b) / \mu_1 = 1 - (\mu_b / \mu_1)}.
#'
#' - If `dif == "perc"` and `comparisons == "all.vs.all"`, the estimated quantity is
#'    \eqn{\phi_{ab} = (\mu_a - \mu_b) / \mu}.
#'
#' - If `dif == "simple"` it estimates \eqn{\mu_a - \mu_b}.
#'
#' @param instance a list object containing the definitions of the problem
#'    instance.
#'    See Section `Instance` for details.
#' @param algorithms a list object containing the definitions of all algorithms.
#'    See Section `Algorithms` for details.
#' @param se.max desired upper limit for the standard error of the estimated
#'        difference between pairs of algorithms. See Section
#'        `Pairwise Differences` for details.
#' @param dif type of difference to be used. Accepts "perc" (for percent
#'          differences) or "simple" (for simple differences)
#' @param comparisons type of comparisons being performed. Accepts "all.vs.first"
#'          (in which cases the first object in `algorithms` is considered to be
#'          the reference algorithm) or "all.vs.all" (if there is no reference
#'          and all pairwise comparisons are desired).
#' @param method method to use for estimating the standard errors. Accepts
#'          "param" (for parametric) or "boot" (for bootstrap)
#' @param nstart initial number of algorithm runs for each algorithm.
#'      See Section `Initial Number of Observations` for details.
#' @param nmax maximum **total** allowed number of runs to execute. Loaded
#'             results (see `load.folder` below) do not count towards this
#'             total.
#' @param seed seed for the random number generator
#' @param boot.R number of bootstrap resamples to use (if `method == "boot"`)
#' @param ncpus number of cores to use
#' @param force.balanced logical flag to force the use of balanced sampling for
#'        the algorithms on each instance
#' @param load.folder name of folder to load results from. Use either "" or
#'        "./" for the current working directory. Accepts relative paths.
#'        Use `NA` for not saving. `calc_nreps()` will look for a .RDS file
#'        with the same name
#' @param save.folder name of folder to save the results. Use either "" or
#'        "./" for the current working directory. Accepts relative paths.
#'        Use `NA` for not saving.
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{instance} - alias for the problem instance considered
#'    \item \code{Xk} - list of observed performance values for all `algorithms`
#'    \item \code{Nk} - vector of sample sizes generated for each algorithm
#'    \item \code{Diffk} - data frame with point estimates, standard errors and
#'    other information for all algorithm pairs of interest
#'    \item \code{seed} - seed used for the PRNG
#'    \item \code{dif} - type of difference used
#'    \item \code{method} - method used ("param" / "boot")
#'    \item \code{comparisons} - type of pairings ("all.vs.all" / "all.vs.first")
#' }
#'
#' @author Prajwel Venkatasubramaniyan (\email{prajwel16@gmail.com})
#'
#' @export
#'
#' @references
#' - F. Campelo
#'    CAISEr: Comparing Algorithms with Iterative Sample-size Estimation in R
#' - F. Campelo, F. Takahashi:
#'    Sample size estimation for power and accuracy in the experimental
#'    comparison of algorithms. Journal of Heuristics 25(2):305-338, 2019.
#' - P. Mathews.
#'    Sample size calculations: Practical methods for engineers and scientists.
#'    Mathews Malnar and Bailey, 2010.
#' -  A.C. Davison, D.V. Hinkley:
#'    Bootstrap methods and their application. Cambridge University Press (1997)
#' -  E.C. Fieller:
#'     Some problems in interval estimation. Journal of the Royal Statistical
#'     Society. Series B (Methodological) 16(2), 175–185 (1954)
#' - V. Franz:
#'    Ratios: A short guide to confidence limits and proper use (2007).
#'    https://arxiv.org/pdf/0710.2024v1.pdf
#' - D.C. Montgomery, C.G. Runger:
#'    Applied Statistics and Probability for Engineers, 6th ed. Wiley (2013)
#'
#' @examples
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
#' instance <- list(alias = "boston_data", Data = boston)
#' algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
#' result_nreps <- calc_nreps_regression(instance, algorithms)
#' summary(result_nreps)
#' plot(result_nreps)

calc_nreps_regression <- function(instance,            # instance parameters
                       algorithms,          # algorithm parameters
                       se.max = 0.05,              # desired (max) standard error
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
      k_foldcv = 3,
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
        k_foldcv = 3,
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

library(caret)
library(e1071)
library(randomForest)
library(Metrics)
library(rpart)
library(gbm)
