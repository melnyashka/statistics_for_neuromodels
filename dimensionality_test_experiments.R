# First we generate the process with arbitrary number of memory variables 
require(MASS)

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

gather_statistics <- function(Z, delta){
  # As an input we take a matrix dxN, where d is the dimensionality of the process, N --- available observations
  dimZ <- dim(Z)
  d = dimZ[1]
  S1_vec <- numeric()
  S2_vec <- numeric()
  # modify the process, generating another Brownian motion
  W_d = mvrnorm(n = dimZ[2], mu = rep(0, times = d), Sigma = diag(d))
  Z1 = Z + delta*t(W_d)
  Z2 = Z + 2*delta*t(W_d)
  
  # Now the calculations:
  i_lim = floor(N/(2*d)) - 1
  for (i in 1:i_lim){
    Z_mat1 <- (Z1[,(2*d*i+1):(2*d*i+d)]-Z1[,(2*d*i):(2*d*i+d-1)])/sqrt(delta)
    Z_mat2 <- (Z2[,seq(from = (2*d*i+2), to = (2*d*i+2*d), by = 2)]-Z2[,seq(from = (2*d*i), to = (2*d*i+2*d-2), by = 2)])/sqrt(2*delta)
    S1_vec[i] <- det(Z_mat1%*%Z_mat1)
    S2_vec[i] <- det(Z_mat2%*%Z_mat2)
    # S1_vec[i] <- det(Z_mat1)^2
    # S2_vec[i] <- det(Z_mat2)^2
  }
  S1 <- 2*d*delta*sum(S1_vec)
  S2 <- 2*d*delta*sum(S2_vec)
  R = d - log(abs(S2/S1))/log(2) 
  return(R)
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

########################################################
#############  Executable part    ######################
########################################################

n_pop = 2 # number of populations
n_neur = c(20, 20) # number of neurons in population, vector of integers of length N
eta = c(3,2) # number of memory variables, vector of integers of length N
nu = c(0.8,0.8) # auxilliary constants
c_rate = c(-1, 1) # rates of population
K = c(1, 10) # constants for the rate functions

N = 10000       # number of observations
delta = 0.01     # discretization step 

Z = hawkes_approximation(N = N, delta = delta, n_pop = n_pop, n_neur = n_neur, eta = eta, nu = nu, c_rate = c_rate, K = K)
gather_statistics(Z, delta)

ind_rough = numeric(n_pop)
for (i in 1:n_pop){ind_rough[i] = sum(eta[1:i])+i}
build_plot(Z, ind_rough)

# now some experiments with FHN


real_parameters <- c(1.5, 0.3, 0.1, 0.01, 0.6) # first set

Z_h = FHN.simulate(parameters = real_parameters, N = N, delta = delta, class = 1)
gather_statistics(Z_h, delta)