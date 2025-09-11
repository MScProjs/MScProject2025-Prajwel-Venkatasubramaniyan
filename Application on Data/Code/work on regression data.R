library(MASS)
data("Boston")

# Feature matrix and target
X <- Boston[, -14]
y <- Boston$medv  # Median value of owner-occupied homes
boston <- as.data.frame(X)
boston$target <- y

instance <- list(alias = "boston", Data = boston)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps1 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps1)
plot.nreps(result_nreps1)

####################################################################################

# 2 california_housing
# library(mlbench)
# data(BostonHousing)
# df <- as.data.frame(BostonHousing)
# names(df)[which(names(df) == "medv")] <- "target"
# instance <- list(alias = "california_housing", Data = df)
# algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
# result_nreps2 <- calc_nreps(instance, algorithms)
# summary.nreps(result_nreps2)
# plot.nreps(result_nreps2)

california <- read.csv("housing[1].csv")
california <- na.omit(california)
head(california)
X <- california[, -9]
X$ocean_proximity <- as.integer (as.factor(X$ocean_proximity))
y <- california$median_house_value  # Median value of owner-occupied homes
cal <- as.data.frame(X)
cal$target <- y
head(cal)

instance <- list(alias = "California Housing", Data = cal)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps2 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps2)

#####################################################################################

library(datasets)
head(airquality) #3
airquality <- na.omit(airquality)
X <- airquality[, -1]
y <- airquality$Ozone  # Median value of owner-occupied homes
air <- as.data.frame(X)
air$target <- y
head(air)

instance <- list(alias = "airquality", Data = air)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps3 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps3)
plot.nreps(result_nreps3)

######################################################################################

df <- read.csv("forestfires.csv")
tail(df)
names(df)[which(names(df) == "area")] <- "target"
df$month <- as.integer(df$month)
df <- df[,-4]
df <- df[,-3]
fire <- df
instance <- list(alias = "forestfires", Data = df) #4
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps4 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps4)
plot.nreps(result_nreps4)

######################################################################################

head(attenu) #5
attenu <- na.omit(attenu)
X <- attenu[, -5]
y <- attenu$accel  # Median value of owner-occupied homes
atten <- as.data.frame(X)
atten$target <- y
head(atten)

df <- as.data.frame(BostonHousing)
names(df)[which(names(df) == "medv")] <- "target"
atten <- df

instance <- list(alias = "attenu", Data = atten)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps5 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps5)
plot.nreps(result_nreps5)

####################################################################################

head(cars) # 6

head(ChickWeight) # 7

head(CO2) #8

head(DNase) # 9

library(ISLR2)

head(Wage)
X <- Wage[, -11]
y <- Wage$wage  # Median value of owner-occupied homes
w <- as.data.frame(X)

w$maritl <- as.integer(w$maritl)
w$race <- as.integer(w$race)
w$education <- as.integer(w$education)
w$region <- as.integer(w$region)
w$jobclass <- as.integer(w$jobclass)
w$health <- as.integer(w$health)
w$health_ins <- as.integer(w$health_ins)

w$target <- y
head(w)

instance <- list(alias = "Wage", Data = w)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps13 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps13)
plot.nreps(result_nreps13)

######################################################################################

head(Auto)
X <- Auto[, -c(1,9)]
y <- Auto$mpg  # Median value of owner-occupied homes
auto <- as.data.frame(X)
auto$target <- y
head(auto)

instance <- list(alias = "Auto", Data = auto)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps14 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps14)
plot.nreps(result_nreps14)
#########################################################################################

insurance <- read.csv("insurance[1].csv")
head(insurance)
X <- insurance[, -7]
y <- insurance$charges  # Median value of owner-occupied homes
ins <- as.data.frame(X)
ins$target <- y
ins$sex <- as.integer(as.factor(ins$sex))
ins$smoker <- as.integer(as.factor(ins$smoker))
ins$region <- as.integer(as.factor(ins$region))
head(ins)

instance <- list(alias = "Insurance", Data = ins)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps15 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps15)

#########################################################################################
head(Bikeshare)

X <- Bikeshare[,-15]
y <- Bikeshare[,15]
X$mnth <- as.integer(X$mnth)
X$hr <- as.integer(X$hr)
X$weathersit <- as.integer(X$weathersit)


bike <- as.data.frame(X)
bike$target <- y
head(bike)

instance <- list(alias = "Bikeshares", Data = bike)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps12 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps12)
plot.nreps(result_nreps12)



###############################################################################
head(Carseats)
X <- Carseats[,-1]
y <- Carseats[,1]
X$ShelveLoc <- as.integer(X$ShelveLoc)
X$Urban <- as.integer(X$Urban)
X$US <- as.integer(X$US)

car <- as.data.frame(X)
car$target <- y
head(car)

instance <- list(alias = "Carseats", Data = car)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps11 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps11)
plot.nreps(result_nreps11)

###################################################################################

head(College)
X <- College[,-18]
y <- College[,18]
X$Private <- as.integer(X$Private)

college <- as.data.frame(X)
college$target <- y
head(college)

instance <- list(alias = "College", Data = college)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps10 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps10)
plot.nreps(result_nreps10)

####################################################################################

head(Credit)
dim(Credit)
X <- Credit[,-11]
y <- Credit[,11]
X$Own <- as.integer(X$Own)
X$Student <- as.integer(X$Student)
X$Married <- as.integer(X$Married)
X$Region <- as.integer(X$Region)

credit <- as.data.frame(X)
credit$target <- y
head(credit)

instance <- list(alias = "Credit", Data = credit)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps9 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps9)
plot.nreps(result_nreps8)

#################################################################################

head(Default)
X <- Default[,-4]
y <- Default[,4]
X$default <- as.integer(X$default)
X$student <- as.integer(X$student)

default <- as.data.frame(X)
default$target <- y
head(default)

instance <- list(alias = "Default", Data = default)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps8 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps8)
plot.nreps(result_nreps8)

##################################################################################
head(ISLR2::Hitters)
Hitters <- na.omit(ISLR2::Hitters)
X <- Hitters[,-19]
y <- Hitters$Salary
X$League <- as.integer(X$League)
X$Division <- as.integer(X$Division)
X$NewLeague <- as.integer(X$NewLeague)

Hitters <- as.data.frame(X)
Hitters$target <- y
head(Hitters)

instance <- list(alias = "Hitters", Data = Hitters)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps7 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps7)
plot.nreps(result_nreps7)

#################################################################################
head(OJ)
dim(OJ)


####################################################################################
head(USArrests)
X <- USArrests[,-1]
y <- USArrests[,1]
Usa <- as.data.frame(X)
Usa$target <- y
head(Usa)

set.seed(123)

Usa100 <- rbind(Usa, Usa)
Usa100 <- Usa100[sample(1:100), ]   # shuffle rows

head(Usa100)
Usa100 <- df

instance <- list(alias = "USAarrests", Data = Usa100)
algorithms <- c("svm", "random_forest", "linear", "decision_tree", "boosting")
result_nreps6 <- calc_nreps(instance, algorithms)
summary.nreps(result_nreps6)
plot.nreps(result_nreps6)

##################################################################################





algorithms <- c("svm", "random_forest")


library(MASS)
data("Aids2")
head(Aids2)

head(Animals)

head(anorexia) # 16
X <- anorexia[,-3]
y <- anorexia[,3]
ano <- as.data.frame(X)
ano$Treat <- as.integer(as.factor(ano$Treat))
ano$target <- y
head(ano)
list(alias = "anorexia", Data = ano)

head(beav1) # 17
X <- beav1[,-3]
y <- beav1[,3]
bea1 <- as.data.frame(X)
bea1$target <- y
head(bea1)
list(alias = "beav1", Data = bea1)


head(beav2) # 18
X <- beav2[,-3]
y <- beav2[,3]
bea2 <- as.data.frame(X)
bea2$target <- y
head(bea2)
list(alias = "beav2", Data = bea2)


head(birthwt) # 19
X <- birthwt[,-10]
y <- birthwt[,10]
bw <- as.data.frame(X)
bw$target <- y
head(bw)
list(alias = "birthwt", Data = bw)


head(cabbages) # 20
X <- cabbages[,-4]
y <- cabbages[,4]
X$Cult <- as.integer(X$Cult)
X$Date <- as.integer(X$Date)
cb <- as.data.frame(X)
cb$target <- y
head(cb)
set.seed(123)

cb100 <- rbind(cb, cb)
cb100 <- cb100[sample(1:120), ]   # shuffle rows

head(cb100)
list(alias = "cabbages", Data = cb100)


head(Cars93) # 21
Cars93 <- na.omit(Cars93)
X <- Cars93[,-c(1,2,3,5,9,10,16,26,27)]
y <- Cars93[,5]
car93 <- as.data.frame(X)
car93$target <- y
head(car93)
list(alias = "Cars93", Data = car93)


head(cats) # 22
cats <- na.omit(cats)
X <- cats[,-3]
y <- cats[,3]
X$Sex <- as.integer(X$Sex)
cat <- as.data.frame(X)
cat$target <- y
head(cat)
list(alias = "Cats", Data = cat)


head(cement) # 23
cement <- na.omit(cement)
X <- cement[,-5]
y <- cement[,5]
cem <- as.data.frame(X)
cem$target <- y
head(cem)
list(alias = "cement", Data = cem)


head(coop) # 24
coop <- na.omit(coop)
X <- coop[,-4]
y <- coop[,4]
co <- as.data.frame(X)
co$target <- y
co$Lab <- as.integer(co$Lab)
co$Spc <- as.integer(co$Spc)
co$Bat <- as.integer(co$Bat)
head(co)
list(alias = "coop", Data = co)

head(cpus) # 25
cpus <- na.omit(MASS::cpus)
cpus <- cpus[,-1]
X <- cpus[,-8]
y <- cpus[,8]
cp <- as.data.frame(X)
cp$target <- y

head(cp)
list(alias = "cpus", Data = cp)

head(crabs) # 26
crabs <- na.omit(MASS::crabs)
crabs <- crabs[,-c(1,2,3)]
X <- crabs[,-5]
y <- crabs[,5]
crab <- as.data.frame(X)
crab$target <- y
head(crab)
list(alias = "crab", Data = crab)

head(eagles) # 27
eagles <- na.omit(eagles)
X <- eagles[,-1]
y <- eagles[,1]
eagle <- as.data.frame(X)
eagle$target <- y
eagle$P <- as.integer(eagle$P)
eagle$A <- as.integer(eagle$A)
eagle$V <- as.integer(eagle$V)
eagle100 <- rbind(eagle, eagle, eagle,eagle,eagle,eagle,eagle,eagle,eagle,eagle,eagle,eagle,eagle)
cb100 <- cb100[sample(1:120), ]   # shuffle rows

head(cb100)
head(eagle)
list(alias = "eagles", Data = eagle100)

head(epil) # 28
epil <- na.omit(epil)
X <- epil[,-9]
y <- epil[,9]
ep <- as.data.frame(X)
ep$target <- y
ep$trt <- as.integer(ep$trt)
head(ep)
list(alias = "epil", Data = ep)

head(gilgais) # 29
gilgais <- na.omit(gilgais)
X <- gilgais[,-9]
y <- gilgais[,9]
gil <- as.data.frame(X)
gil$target <- y
head(gil)
list(alias = "gilgais", Data = gil)


head(housing) # 30
housing <- na.omit(housing)
X <- housing[,-5]
y <- housing[,5]
X$Sat <- as.integer(X$Sat)
X$Infl <- as.integer(X$Infl)
X$Type <- as.integer(X$Type)
X$Cont <- as.integer(X$Cont)
hou <- as.data.frame(X)
hou$target <- y
head(hou)
list(alias = "Housing", Data = hou)


n.algs <- 8
n.comparisons <- switch(comparisons,
                        all.vs.all   = n.algs * (n.algs - 1) / 2,
                        all.vs.first = n.algs - 1)
n.comparisons

ss.calc <- calc_instances(ncomparisons = n.comparisons,
                          d            = d,
                          power        = power,
                          sig.level    = sig.level,
                          alternative.side  = alternative.side,
                          test         = test,
                          power.target = power.target)

N.star <- ceiling(ss.calc$ninstances)
N.star

# Load library
library(ggplot2)

# Create the data
data <- data.frame(
  Algorithms = c(2, 3, 4, 5, 6, 7, 8),
  Datasets = c(27, 33, 39, 43, 47, 50, 53)
)

# Plot
ggplot(data, aes(x = Algorithms, y = Datasets)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 3, color = "red") +
  labs(
    x = "Number of Algorithms",
    y = "Number of Datasets Required"
  ) +
  theme_minimal()


########################################################################################

instances <- list(list(alias = "attenu", Data = atten),
                  list(alias = "boston", Data = boston),
                  list(alias = "california_housing", Data = cal),
                  list(alias = "airquality", Data = air),
                  list(alias = "forestfires", Data = fire),
                  list(alias = "Wage", Data = w),
                  list(alias = "Auto", Data = auto),
                  list(alias = "Insurance", Data = ins),
                  list(alias = "Bikeshares", Data = bike),
                  list(alias = "Carseats", Data = car),
                  list(alias = "College", Data = college),
                  list(alias = "Credit", Data = credit),
                  list(alias = "Default", Data = default),
                  list(alias = "Hitters", Data = Hitters),
                  list(alias = "USAarrests", Data = Usa100),
                  list(alias = "anorexia", Data = ano),
                  list(alias = "beav1", Data = bea1),
                  list(alias = "beav2", Data = bea2),
                  list(alias = "birthwt", Data = bw),
                  list(alias = "cabbages", Data = cb100),
                  list(alias = "Cars93", Data = car93),
                  list(alias = "Cats", Data = cat),
                  list(alias = "cement", Data = cem),
                  list(alias = "coop", Data = co),
                  list(alias = "cpus", Data = cp),
                  list(alias = "crab", Data = crab),
                  list(alias = "eagles", Data = eagle100),
                  list(alias = "epil", Data = ep),
                  list(alias = "gilgais", Data = gil),
                  list(alias = "Housing", Data = hou)

)

algorithms <- c("svm", "random_forest")
final_result <- run_experiment_regression(instances, algorithms)
summary(final_result)
plot(final_result)
