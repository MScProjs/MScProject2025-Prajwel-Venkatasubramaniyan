algorithms <- mapply(FUN = function(i, m, s){
    list(FUN   = "dummyalgo",
         alias = paste0("algo", i),
         distribution.fun  = "rnorm",
         distribution.pars = list(mean = m, sd = s))},
    i = c(alg1 = 1, alg2 = 2, alg3 = 3, alg4 = 4),
    m = c(15, 10, 30, 15),
    s = c(2, 4, 6, 8),
    SIMPLIFY = FALSE)

# plot(my.results)Generate 100 dummy instances with centered exponential distributions
instances <- lapply(1:100,
                      function(i) {rate <- runif(1, 1, 10)
                                   list(FUN   = "dummyinstance",
                                        alias = paste0("Inst.", i),
                                        distr = "rexp", rate = rate,
                                        bias  = -1 / rate)})

my.results <- run_experiment(instances, algorithms,
                               d = .5, se.max = .1,
                               power = .9, sig.level = .05,
                               power.target = "mean",
                               dif = "perc", comparisons = "all.vs.all",
                               ncpus = 1, seed = 1234)

# Take a look at the results
summary(my.results)
plot(my.results)
