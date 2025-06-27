# mean 

K      <- 10   # number of comparisons
alpha  <- 0.05 # significance level
power  <- 0.5  # desired power
d      <- 0.5  # MRES

out <- calc_instances(K, d,
                      power     = power,
                      sig.level = alpha, power.target = "mean")

# Plot power of each comparison to detect differences of magnitude d
plot(1:K, out$power,
     type = "b", pch = 20, las = 1, ylim = c(0, 1), xlab = "comparison",
     ylab = "power", xaxs = "i", xlim = c(0, 11))
#########################################################################################

# variance
K      <- 10   # number of comparisons
alpha  <- 0.05 # significance level
power  <- 0.01  # desired power
d      <- 0.5  # MRES

out <- calc_instances(K, d,
                      power     = power,
                      sig.level = alpha, power.target = "var")
# Plot power of each comparison to detect differences of magnitude d
plot(1:K, out$power,
     type = "b", pch = 20, las = 1, ylim = c(0, 1), xlab = "comparison",
     ylab = "power", xaxs = "i", xlim = c(0, 11))
#############################################################################################

# skewness
K      <- 10   # number of comparisons
alpha  <- 0.05 # significance level
power  <- 0.5  # desired power
d      <- 0.5  # MRES

out <- calc_instances(K, d,
                      power     = power,
                      sig.level = alpha, power.target = "skewness")

# Plot power of each comparison to detect differences of magnitude d
plot(1:K, out$power,
     type = "b", pch = 20, las = 1, ylim = c(0, 1), xlab = "comparison",
     ylab = "power", xaxs = "i", xlim = c(0, 11))
#########################################################################################

# kurtosis

K      <- 10   # number of comparisons
alpha  <- 0.05 # significance level
power  <- 4  # desired power
d      <- 0.5  # MRES

out <- calc_instances(K, d,
                      power     = power,
                      sig.level = alpha, power.target = "kurtosis")

# Plot power of each comparison to detect differences of magnitude d
plot(1:K, out$power/10,
     type = "b", pch = 20, las = 1, xlab = "comparison",
     ylab = "power", xaxs = "i", xlim = c(0, 11))

##########################################################################################
