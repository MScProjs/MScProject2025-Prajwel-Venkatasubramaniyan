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
                            c("worst.case", "mean", "median", "var", "skewness", "kurtosis"))
  
  
  # ========== Error catching ========== #
  assertthat::assert_that(
    assertthat::is.count(ncomparisons),
    is.null(ninstances) || (assertthat::is.count(ninstances) && ninstances > 1),
    is.null(power) || (is.numeric(power) && power > 0),
    is.null(d) || (is.numeric(d) && d > 0),
    sum(c(is.null(ninstances), is.null(power), is.null(d))) == 1,
    is.numeric(sig.level) && sig.level > 0 && sig.level < 1,
    alternative.side %in% c("one.sided", "two.sided"),
    test %in% c("t.test", "wilcoxon", "binomial"),
    power.target %in% c("worst.case", "mean", "median", "var", "skewness", "kurtosis"))
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
      p.median   <- stats::median(p)
      p.var      <- stats::var(p)
      p.skewness <- moments::skewness(p)
      p.kurtosis <- moments::kurtosis(p)
      #
    } else {
      if (power.target == "worst.case") {
        r <- 1
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
        p.mean     <- mean(p)
        p.median   <- stats::median(p)
        p.var      <- stats::var(p)
        p.skewness <- moments::skewness(p)
        p.kurtosis <- moments::kurtosis(p)
      } # Bonferroni correction
      
      if (power.target == "median")     {
        r <- floor(ncomparisons / 2)
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
        p.mean     <- mean(p)
        p.median   <- stats::median(p)
        p.var      <- stats::var(p)
        p.skewness <- moments::skewness(p)
        p.kurtosis <- moments::kurtosis(p) 
      }
      if (power.target == "var") 
      {
        N <- stats::power.t.test(delta = d, sd = 1,
                                 sig.level = sig.level, power = power,
                                 type = "paired",
                                 alternative = alternative.side)$n - 1
        
        p.var <- 0 # mean power
        while (p.var < power){
          
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
          p.var <- var(p)
        }
        p.mean <- mean(p)
        p.median      <- stats::median(p)
        p.skewness <- moments::skewness(p)
        p.kurtosis <- moments::kurtosis(p)
      }
      if (power.target == "skewness")   {
        N <- stats::power.t.test(delta = d, sd = 1,
                                 sig.level = sig.level, power = power,
                                 type = "paired",
                                 alternative = alternative.side)$n - 1
        
        p.skewness <- 0 # mean power
        while (p.skewness < power){
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
          p.skewness <- moments::skewness(p)
        }
        p.mean     <- mean(p)
        p.median   <- stats::median(p)
        p.var      <- stats::var(p)
        p.kurtosis <- moments::kurtosis(p)
      }
      if (power.target == "kurtosis")   {
        N <- stats::power.t.test(delta = d, sd = 1,
                                 sig.level = sig.level, power = power/10,
                                 type = "paired",
                                 alternative = alternative.side)$n - 1
        
        p.kurtosis <- 1 # mean power
        while (p.kurtosis/10 < power/10){
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
          p.kurtosis <- moments::kurtosis(p)
        }
        p.mean     <- mean(p)
        p.median   <- stats::median(p)
        p.var      <- stats::var(p)
        p.skewness <- moments::skewness(p)
      }
      
      
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
    p.mean     <- mean(p)
    p.median   <- stats::median(p)
    p.var      <- stats::var(p)
    p.skewness <- moments::skewness(p)
    p.kurtosis <- moments::kurtosis(p) 
  }
  
  
  output <- list(ninstances   = N,
                 power        = p,
                 mean.power   = p.mean,
                 median.power = p.median,
                 var.power    = p.var,
                 skew.power   = p.skewness,
                 kurt.power   = p.kurtosis,
                 d            = d,
                 sig.level    = a,
                 alternative  = alternative.side,
                 test         = test,
                 power.target = power.target)
  
  return(output)
}
