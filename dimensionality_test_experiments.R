# First we generate the process with arbitrary number of memory variables 
require(MASS)
require(Rcpp)
# require(inline)
require(futile.logger) # setting up logging
require(parallel)
sourceCpp("example.cpp")

rate_function <- function(x, const=1){
  # here for simplicity we use the same family of functions
  if (x<log(20)){
    value <- const*exp(x)
  } else {
    value <- 40*const/(1+400*exp(-2*x))
  }
  return(value)
}

hawkes_approximation <- function(N, delta, n_pop, n_neur, eta, nu, c_rate, K){
  n_total = n_pop+sum(eta)
  Z = matrix(nrow = n_total, ncol = N) # slot for the process
  ind_rough = numeric(n_pop)
  # now we are checking the indexes of the rough variables
  for (i in 1:n_pop){ind_rough[i] = sum(eta[1:i])+i}
  Z[,1] = rep(0,n_total)
  for (n in 1:(N-1)){
    for (k in 1:n_total){
      i = min(which(ind_rough >= k))
      if (k %in% ind_rough) {
        if (k+1 <= n_total) {j = k + 1} else {j = 1}
        if (i+1 <= n_pop) {i_ = i + 1} else {i_ = 1}
        p = n_neur[i_] #/sum(n_neur)
        Z[k, n + 1] = Z[k,n] + delta*(-nu[i]*Z[k, n] + c_rate[i]*rate_function(x = Z[j,n], const = K[i_])) + c_rate[i]*sqrt(delta)*rnorm(1, mean = 0, sd = 1)*sqrt(rate_function(x = Z[j,n], const = K[i_])/p)
      } else { 
        Z[k, n + 1] = Z[k,n] + delta*(-nu[i]*Z[k, n] + Z[k + 1,n])
      }
    }
  }
  return(Z)
}

FHN.simulate <- function(parameters, N, delta, class = 1){
  # Explanations for classes:
  # 1 = hypoelliptic
  # 2 = elliptic
  # 3 = deterministic
  
  gamma <- parameters[1]
  beta <- parameters[2]
  eps <- parameters[3]
  s <- parameters[4]
  sigma <- parameters[5]
  
  X <- numeric()
  Y <- numeric()
  
  Bx <- rnorm(N, mean = 0, sd = 1) # Generate increments of Brownian motion
  By <- rnorm(N, mean = 0, sd = 1) 
  X[1] <- 0  # Starting points
  Y[1] <- 0
  if (class == 1){
    for (i in 1:(N-1)){
      X[i+1] <- X[i] + delta/eps*(X[i]-X[i]^3-Y[i] + s) + delta^2/(2*eps)*((1-3*X[i]^2)/eps*(X[i]-X[i]^3-Y[i] + s) - (gamma*X[i]-Y[i]+beta))  + (delta^(3/2)*By[i] + delta^(3/2)*Bx[i]/sqrt(3))/eps*sigma/2
      Y[i+1] <- Y[i] + delta*(gamma*X[i]-Y[i]+beta) + delta^2/2*(gamma/eps*(X[i]-X[i]^3-Y[i] + s) - (gamma*X[i]-Y[i]+beta)) + sqrt(delta)*By[i]*sigma
    }} else if (class == 2) {
      for (i in 1:(N-1)){
        X[i+1] <- X[i] + delta/eps*(X[i]-X[i]^3-Y[i] + s) + delta^2/(2*eps)*((1-3*X[i]^2)/eps*(X[i]-X[i]^3-Y[i] + s) - (gamma*X[i]-Y[i]+beta))  + sqrt(delta)*Bx[i]*sigma
        Y[i+1] <- Y[i] + delta*(gamma*X[i]-Y[i]+beta) + delta^2/2*(gamma/eps*(X[i]-X[i]^3-Y[i] + s) - (gamma*X[i]-Y[i]+beta)) + sqrt(delta)*By[i]*sigma
      }} else {
        for (i in 1:(N-1)){
          X[i+1] <- X[i] + delta/eps*(X[i]-X[i]^3-Y[i] + s) + delta^2/(2*eps)*((1-3*X[i]^2)/eps*(X[i]-X[i]^3-Y[i] + s) - (gamma*X[i]-Y[i]+beta))
          Y[i+1] <- Y[i] + delta*(gamma*X[i]-Y[i]+beta) + delta^2/2*(gamma/eps*(X[i]-X[i]^3-Y[i] + s) - (gamma*X[i]-Y[i]+beta))
        }
      }
  return(rbind(X, Y))
}

BM.simulate <- function(N, delta, class = 1){
  X <- numeric(N)
  Y <- numeric(N)
  Bx <- rnorm(N, mean = 0, sd = 1) # Generate increments of Brownian motion
  By <- rnorm(N, mean = 0, sd = 1) 
  X[1] <- 0  # Starting points
  Y[1] <- 0
  if (class == 1){
    for (i in 1:(N-1)){
      X[i+1] <- X[i] # + (delta^(3/2)*By[i] + delta^(3/2)*Bx[i]/sqrt(3))*sigma/2
      Y[i+1] <- Y[i] + sqrt(delta)*By[i]
    }} else if (class == 2) {
      for (i in 1:(N-1)){
        X[i+1] <- X[i] + sqrt(delta)*Bx[i]
        Y[i+1] <- Y[i] +  sqrt(delta)*By[i]
      }} 
  return(rbind(X, Y))
}

ML.simulate <- function(parameters, N, delta, class = 1){
  X <- numeric()
  Y <- numeric()
  I <- parameters[1]   # applied current
  C <- parameters[2]   # membrane capacitance
  g <- parameters[3:5] # conductances
  V <- parameters[6:12] # equilibrium potentials + tuning parameters
  phi <- parameters[13] # reference frequency
  for (i in 1:N){
    Mss = 0.5*(1+tanh((Y - V[4])/V[5]))
    alpha = 0.5
  }
}

build_plot <- function(Z, ind_rough){
  dimZ <- dim(Z)
  ylim_max <- max(Z)
  ylim_min <- min(Z)
  plot(Z[1,], type = "l", xlab = "", ylab = "", ylim = c(ylim_min, ylim_max), col = "grey")
  for (i in 2:dimZ[1]){
    if (i %in% ind_rough) {col_i = "black"} else {col_i = "grey"}
    lines(Z[i,], col = col_i)
  }
}

construct_test <- function(Z, delta, j){
  # As an indput we take a matrix dxN, where d is the dimensionality of the process, N --- available observations
  dimZ <- dim(Z)
  j = 1
  sigma = 1
  
  if (is.null(dimZ)) {
    d = 1
    N = length(Z)
    Z = t(as.matrix(Z))
  } else {
    d = dimZ[1]
    N = dimZ[2]
  }
  
  # modify the process, generating another Brownian motion
  if (d==1) {
    W_d <- t(as.matrix(cumsum(rnorm(n = N, mean = 0, sd = sigma*sqrt(delta)))))
    } else {
      W <- t(mvrnorm(n = N, mu = rep(0, times = d), Sigma = sigma*sqrt(delta)*diag(d)))
      W_d <- matrix(unlist(lapply(c(1:d), function(i) cumsum(W[i,]))), nrow = d, byrow = TRUE)
      }
  
  Z1 = add_rcpp(Z,(delta)^(1/2)*W_d)
  Z2 = add_rcpp(Z,(2*delta)^(1/2)*W_d)
  
  # Now the calculations:
  # Here S1 and S2 correspond to (2.13) on page 2396 from paper, 
  # V - to (3.17) on page 2401 

  i_lim = floor(N/(2*d)) - 1
  
  Z_mat1 <- lapply(c(1:i_lim), function(i) (Z1[,(2*d*i+1):(2*d*i+d)]-Z1[,(2*d*i):(2*d*i+d-1)])/sqrt(delta))
  Z_mat2 <- lapply(c(1:i_lim), function(i) (Z2[,seq(from = (2*d*i+2), to = (2*d*i+2*d), by = 2)]-Z2[,seq(from = (2*d*i), to = (2*d*i+2*d-2), by = 2)])/sqrt(2*delta))
  if (d==1){
    S1_vec <- unlist(lapply(Z_mat1, function(x) x^2))
    S2_vec <- unlist(lapply(Z_mat2, function(x) x^2))
  } else {
    S1_vec <- unlist(lapply(Z_mat1, edet))
    S2_vec <- unlist(lapply(Z_mat2, edet))
  }

  S1 <- 2*d*delta*sum_rcpp(S1_vec)
  S2 <- 2*d*delta*sum_rcpp(S2_vec) 
  
  R = d - log(abs(S2/S1))/log(2)  # Value of the "estimator" R hat, formula (3.10)
  V <- sum_rcpp((add_rcpp_vec(S1_vec, - S2_vec*2^(R-d)))^2)/(delta*(sum_rcpp(S1_vec)*log(2))^2) # the same value, as below, but computed faster
  # V <- (V11 + 2^(2*(R-d))*V22 - 2^(1+R-d)*V12)/(log(2)*S1)^2
  return(c(R, V))
}

########################################################
#############  Executable part    ######################
########################################################
setwd("University/neuroscience/Dimensionality_tests") # comment/uncomment if necessary
date <- Sys.Date()          # Create a subdirectory with the current date
wd <- getwd()               # Save the name of the working directory 
path_to_logs <- file.path(wd,date)
dir.create(path_to_logs)
file <- paste(path_to_logs, "logfile", sep = "/")
# flog.appender(appender.file(file))
flog.appender(appender.tee(file)) # write both to console and to file 
flog.threshold(DEBUG)        # By default set threshold to INFO (because I can)
flog.debug("Debugging is on!") 
my_colors <- c(adjustcolor("skyblue", alpha.f = 0.5),  adjustcolor("chartreuse", alpha.f = 0.5), adjustcolor("coral", alpha.f = 0.5), adjustcolor("darkgoldenrod1", alpha.f = 0.5))

# Symmetric quantile values 
q_95 <- 1.959964
q_99 <- 2.575829
q_999 <- 3.290527

##### Experimental part: Hawkes process  ######

n_pop = 2 # number of populations
n_neur = c(20, 20) # number of neurons in population, vector of integers of length N
eta = c(3,3) # number of memory variables, vector of integers of length N
nu = c(1,1) # auxilliary constants
c_rate = c(-1,1) # rates of population
K = c(1, 10) # constants for the rate functions

# N_set = c(500, 5000, 50000, 500000)     # number of observations
delta_set = c(0.1, 0.01, 0.001, 0.0001, 0.00001)    # discretization step 
N_trials <- 100
delta_gen = 0.00001
N_gen = 5000000
true_decisions <- numeric()
true_rejection <- numeric()

flog.debug("We are working with Hawkes model, number of population is %s, number of neurons is %s, eta = %s, nu = %s, c_rate = %s", n_pop, toString(n_neur), toString(eta), toString(nu), toString(c_rate))
R_noint <- matrix(nrow = length(delta_set), ncol = N_trials)
test_normalized <- matrix(nrow = length(delta_set), ncol = N_trials)

# Setting up cluster
# cl <- makeCluster(3)

for (k in 1:length(delta_set)){
  Z = hawkes_approximation(N = N_gen, delta = delta_gen, n_pop = n_pop, n_neur = n_neur, eta = eta, nu = nu, c_rate = c_rate, K = K)
  Z_sub = Z[,seq(1,N_gen,as.integer(delta_set[k]/delta_gen))]
  TEST <- lapply(c(1:N_trials), construct_test, Z = Z_sub, delta = delta_set[k])
  cl <- makeCluster(4)
  clusterExport(cl = cl, varlist = c("TEST", "delta_set", "n_pop", "N_trials", "k"))
  R_noint[k,] <- unlist(parLapply(cl = cl, c(1:N_trials), function(i) TEST[[i]][1]))
  test_normalized[k,] <- unlist(parLapply(cl = cl, c(1:N_trials), function(i) abs(TEST[[i]][1] - n_pop)/sqrt(delta_set[k]*TEST[[i]][2])))
  # pdf(paste(path_to_logs, paste("R_noint_Hawkes_", format(Sys.time(), format = "%H:%M:%S"), ".pdf", sep=""), sep = "/"))
  par(mfrow = c(2,1))
  plot(density(R_noint[k,]), main = "Density of R_hat")
  polygon(density(R_noint[k,]), col = my_colors[3])
  plot(density(test_normalized[k,]), col = my_colors[1], main = "Density of normalized test statistics", ylim = c(0,0.8))
  polygon(density(test_normalized[k,]), col = my_colors[1])
  polygon(density(abs(rnorm(N_trials))), col = my_colors[2])
  # dev.off()
  true_decisions[k] = length(which(round(R_noint[k,]) == n_pop))/N_trials
  #true_rejection[k] = length(which(test_h0[k,] == TRUE))/N_trials
  true_rejection[k] = length(which(test_normalized[k,]>q_95))/N_trials
  flog.debug("Percent of true decisions for Delta = %s, N = %s, for %s of trials is: %s, null hypothesis is rejected in %s cases", delta_set[k], length(Z_sub[1,]), N_trials, true_decisions[k], true_rejection[k])
}
# stopCluster(cl)

pdf(paste(path_to_logs, paste("R_density_Hawkes_", format(Sys.time(), format = "%H:%M:%S"), ".pdf", sep=""), sep = "/"))
plot(density(R_noint[1,]), main = "", xlab = "", ylab = "", xlim = c(min(R_noint),max(R_noint)), ylim = c(0, mean(R_noint[2,])))
polygon(density(R_noint[1,]), col = my_colors[1])
polygon(density(R_noint[2,]),  col = my_colors[2])
polygon(density(R_noint[3,]), col = my_colors[3])
polygon(density(R_noint[4,]),  col = my_colors[4])
legend("topleft", inset=.02, title="Size of delta",
       c("0.1","0.01","0.001", "0.0001"), fill=my_colors, horiz=TRUE, cex=0.8)
dev.off()

write.csv(R_noint, file = paste(path_to_logs, paste("R_noint_Hawkes_", format(Sys.time(), format = "%H:%M")), sep = "/"), row.names = FALSE)

ind_rough = numeric(n_pop)
for (i in 1:n_pop){ind_rough[i] = sum(eta[1:i])+i}
build_plot(Z, ind_rough)

####### Experimental part: FHN model #######

real_parameters <- c(1.5, 0.3, 0.1, 0.01, 0.6) # first set
# real_parameters <- c(1.2, 1.3, 0.1, 0.01, 0.4) # second set

delta_set = c(0.1, 0.01, 0.001, 0.0001, 0.00001)    # discretization step 
N_trials <- 1000
true_decisions <- numeric()
true_rejection <- numeric()
delta_gen = 0.00001
N_gen = 1000000
dim_true <- 1

flog.debug("FitzHugh-Nagumo model, set: %s", toString(real_parameters))
R_noint <- matrix(nrow = length(delta_set), ncol = N_trials)
test_normalized <- matrix(nrow = length(delta_set), ncol = N_trials)

for (k in 1:length(delta_set)){
  Z = FHN.simulate(parameters = real_parameters, N = N_gen, delta = delta_gen, class = dim_true)
  Z_h = Z[,seq(1,N_gen,as.integer(delta_set[k]/delta_gen))]
  TEST <- lapply(c(1:N_trials), construct_test, Z = Z_h, delta = delta_set[k])
  cl <- makeCluster(4)
  clusterExport(cl = cl, varlist = c("TEST", "delta_set", "dim_true", "N_trials", "k"))
  R_noint[k,] <- unlist(parLapply(cl = cl, c(1:N_trials), function(i) TEST[[i]][1]))
  test_normalized[k,] <- unlist(parLapply(cl = cl, c(1:N_trials), function(i) abs(TEST[[i]][1] - dim_true)/sqrt(delta_set[k]*TEST[[i]][2])))
  stopCluster(cl)

  # pdf(paste(path_to_logs, paste("R_noint_FHN_", format(Sys.time(), format = "%H:%M:%S"),".pdf",sep=""), sep = "/"))
  par(mfrow = c(2,1))
  plot(density(R_noint[k,]), main = "Density of R_hat")
  polygon(density(R_noint[k,]), col = my_colors[3])
  plot(density(test_normalized[k,]), main = "Density of normalized test statistics")
  polygon(density(test_normalized[k,]), col = my_colors[1])
  polygon(density(abs(rnorm(N_trials))), col = my_colors[2])
  # dev.off()
  
  true_decisions[k] = length(which(round(R_noint[k,]) == dim_true))/N_trials
  true_rejection[k] = length(which(test_normalized[k,]>q_95))/N_trials
  
  flog.debug("Percent of true decisions for Delta = %s, N = %s, for %s of trials is: %s, null hypothesis is rejected in %s cases", delta_set[k], length(Z_h[1,]), N_trials, true_decisions[k], true_rejection[k])
}

write.csv(R_noint, file = paste(path_to_logs, paste("R_noint_FHN_", format(Sys.time(), format = "%H:%M")), sep = "/"), row.names = FALSE)

pdf(paste(path_to_logs, paste("R_density_FHN_", format(Sys.time(), format = "%H:%M:%S"), ".pdf", sep=""), sep = "/"))
my_colors <- c(adjustcolor("skyblue", alpha.f = 0.5),  adjustcolor("chartreuse", alpha.f = 0.5), adjustcolor("coral", alpha.f = 0.5), adjustcolor("darkgoldenrod1", alpha.f = 0.5))
plot(density(R_noint[3,]), main = "", xlab = "", ylab = "", xlim = c(min(R_noint),max(R_noint)))
polygon(density(R_noint[1,]), col = my_colors[1])
polygon(density(R_noint[2,]),  col = my_colors[2])
polygon(density(R_noint[3,]), col = my_colors[3])
legend("topleft", inset=.02, title="Size of delta",
       c("0.1","0.01","0.001"), fill=my_colors[1:3], horiz=TRUE, cex=0.8)
dev.off()

##### Experiments: "pure Brownian motion"  ######

# N_set = c(500, 5000, 50000, 500000)     # number of observations
delta_set = c(0.1, 0.01, 0.001, 0.0001, 0.00001)    # discretization step 
N_trials <- 1000
true_decisions <- numeric()
true_rejection <- numeric()
delta_gen = 0.00001
N_gen = 1000000
dim_true <- 1

flog.info("Brownian motion, dimension = %s", dim_true)
R_noint <- matrix(nrow = length(delta_set), ncol = N_trials)
test_normalized <- matrix(nrow = length(delta_set), ncol = N_trials)

Z = BM.simulate(N = N_gen, delta = delta_gen, class = 2)
for (k in 1:length(delta_set)){
  Z_h = Z[,seq(1,N_gen,as.integer(delta_set[k]/delta_gen))]
  TEST <- lapply(c(1:N_trials), construct_test, Z = Z_h[1,], delta = delta_set[k])
  cl <- makeCluster(4)
  clusterExport(cl = cl, varlist = c("TEST", "delta_set", "dim_true", "N_trials", "k"))
  R_noint[k,] <- unlist(parLapply(cl = cl, c(1:N_trials), function(i) TEST[[i]][1]))
  test_normalized[k,] <- unlist(parLapply(cl = cl, c(1:N_trials), function(i) abs(TEST[[i]][1] - dim_true)/sqrt(delta_set[k]*TEST[[i]][2])))
  stopCluster(cl)
  
  # pdf(paste(path_to_logs, paste("R_noint_FHN_", format(Sys.time(), format = "%H:%M:%S"),".pdf",sep=""), sep = "/"))
  par(mfrow = c(2,1))
  plot(density(R_noint[k,]), main = "Density of R_hat")
  polygon(density(R_noint[k,]), col = my_colors[3])
  plot(density(test_normalized[k,]), main = "Density of normalized test statistics")
  polygon(density(test_normalized[k,]), col = my_colors[1])
  polygon(density(abs(rnorm(N_trials))), col = my_colors[2])
  # dev.off()
  true_decisions[k] = length(which(round(R_noint[k,]) == dim_true))/N_trials
  true_rejection[k] = length(which(test_normalized[k,]>q_95))/N_trials
  flog.debug("Percent of true decisions for Delta = %s, N = %s, for %s of trials is: %s, null hypothesis is rejected in %s cases", delta_set[k], length(Z_h[1,]), N_trials, true_decisions[k], true_rejection[k])
  flog.debug("Standard deviation of R_hat is %s", (R_noint[k,]-dim_true)^2)
}

##### Another toy example #####

f_t <- function(t){
  return(1 + (2*t - 1)^2)
}
delta_gen = 0.0001
N_gen = as.integer(1/delta_gen)
S <- 0.5
dim_true <- 2-ifelse((S==0), 1, 0)
N_trials = 1000

R_noint <- numeric()
X <- numeric()
Y <- numeric()
X[1] <- 0
Y[1] <- 0

for (i in 1:(N_gen-1)){
  X[i+1] <- X[i] + 2*delta_gen + f_t(i*delta_gen)*rnorm(1, sd = sqrt(delta_gen))*ifelse((i*delta_gen < S), 1, 0)
  Y[i+1] <- Y[i] + 2*delta_gen + f_t(i*delta_gen)*rnorm(1, sd = sqrt(delta_gen))
}

Z_xy <- rbind(X, Y)

TEST <- lapply(c(1:N_trials), construct_test, Z = Z_xy, delta = delta_gen)
cl <- makeCluster(4)
clusterExport(cl = cl, varlist = c("TEST", "delta_gen", "dim_true", "N_trials"))
R_noint <- unlist(parLapply(cl = cl, c(1:N_trials), function(i) TEST[[i]][1]))
test_normalized <- unlist(parLapply(cl = cl, c(1:N_trials), function(i) abs(TEST[[i]][1] - dim_true)/sqrt(delta_gen*TEST[[i]][2])))
stopCluster(cl)

# pdf(paste(path_to_logs, paste("R_noint_FHN_", format(Sys.time(), format = "%H:%M:%S"),".pdf",sep=""), sep = "/"))
par(mfrow = c(2,1))
plot(density(R_noint), main = "Density of R_hat")
polygon(density(R_noint), col = my_colors[3])
plot(density(test_normalized), main = "Density of normalized test statistics")
polygon(density(test_normalized), col = my_colors[1])
polygon(density(abs(rnorm(N_trials))), col = my_colors[2])
# dev.off()

true_decisions = length(which(round(R_noint) == dim_true))/N_trials
true_rejection = length(which(test_normalized>q_95))/N_trials

flog.debug("Percent of true decisions for Delta = %s, for %s of trials is: %s, null hypothesis is rejected in %s cases", delta_gen, N_trials, true_decisions, true_rejection)
