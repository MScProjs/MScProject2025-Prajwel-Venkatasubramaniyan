library(CAISEr)
algorithms <- mapply(FUN = function(i, m, s){
                            list(FUN   = "dummyalgo",
                                 alias = paste0("algo", i),
                                 distribution.fun  = "rnorm",
                                 distribution.pars = list(mean = m, sd = s))},
                       i = c(alg1 = 1, alg2 = 2, alg3 = 3, alg4 = 4, alg5 = 5),
                       m = c(15, 10, 30, 15, 20),
                       s = c(2, 4, 6, 8, 10),
                       SIMPLIFY = FALSE)

# Make a dummy instance with a centered (zero-mean) exponential distribution:
instance = list(FUN = "dummyinstance", distr = "rexp", rate = 5, bias = -1/5)

# Explicitate all other parameters (just this one time:
# most have reasonable default values)
myreps <- calc_nreps(instance   = instance,
                        algorithms = algorithms,
                        se.max     = 0.05,          # desired (max) standard error
                        dif        = "perc",        # type of difference
                        comparisons = "all.vs.all", # differences to consider
                        method     = "param",       # method ("param", "boot")
                        nstart     = 15,            # initial number of samples
                        nmax       = 1000,          # maximum allowed sample size
                        seed       = 1234,          # seed for PRNG
                        boot.R     = 499,           # number of bootstrap resamples (unused)
                        ncpus      = 1,             # number of cores to use
                        force.balanced = FALSE,     # force balanced sampling?
                        load.folder   = NA,         # file to load results from
                        save.folder = NA)         # folder to save results
summary(myreps)
plot(myreps)
  