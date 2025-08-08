instances <- list(list(alias = "instance1", Data = boston), list(alias = "instance2", Data = california))
algorithms <- c("svm", "random_forest")
result_nreps <- calc_nreps(instance, algorithms)

run_MLmodel <- function(instances, algorithms, d=0.5, se.max=0.1,
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
