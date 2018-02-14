#####################################
######     Auxilliary part     ######
#####################################

FHN.simulate <- function(parameters, N, delta, start = c(0,0)){
  gamma <- parameters[1]
  beta <- parameters[2]
  eps <- parameters[3]
  s <- parameters[4]
  sigma <- parameters[5]
  X <- numeric()
  Y <- numeric()
  Bx <- rnorm(N, mean = 0, sd = 1) # Generate increments of Brownian motion
  By <- rnorm(N, mean = 0, sd = 1) 
  X[1] <- start[1]  # Starting points
  Y[1] <- start[2] 
  for (i in 1:(N-1)){
    X[i+1] <- X[i] + delta/eps*(X[i]-X[i]^3-Y[i] + s) + delta^2/(2*eps)*((1-3*X[i]^2)/eps*(X[i]-X[i]^3-Y[i] + s) - (gamma*X[i]-Y[i]+beta))  + (delta^(3/2)*By[i] + delta^(3/2)*Bx[i]/sqrt(3))/eps*sigma/2
    Y[i+1] <- Y[i] + delta*(gamma*X[i]-Y[i]+beta) + delta^2/2*(gamma/eps*(X[i]-X[i]^3-Y[i] + s) - (gamma*X[i]-Y[i]+beta)) + sqrt(delta)*By[i]*sigma
    }
  return(rbind(X, Y))
}

FHN.sample_contrast <- function(theta, x, y, delta) {
  n = length(x)
  gamma = theta[1]
  beta = theta[2]
  eps = theta[3]
  # eps = 0.1
  # s = theta[4]
  s = 0.01
  x_diff <- x[2:n] - x[1:(n-1)]
  y_diff <- y[2:n] - y[1:(n-1)]
  var_y = mean(y_diff^2)
  var_x = mean(x_diff^2)
  cov_xy = mean(x_diff*y_diff)
  sigma2 = var_y/delta
  # eps = (delta*var_y)/(2*cov_xy)
  crit_x = x[2:n]- x[1:(n-1)] - delta/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - delta^2/(2*eps)*((1-3*x[1:(n-1)]^2)/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  crit_y = y[2:n] - y[1:(n-1)]- delta*(gamma*x[1:(n-1)]-y[1:(n-1)]+beta) - delta^2/2*(gamma/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  # crit = sum(crit_y^2/var_y - 3*crit_y*crit_x/cov_xy + 3*crit_x^2/var_x) # + (n-2)*log(sigma2/eps)
  crit = sum(crit_y^2/var_y - 2*crit_y*crit_x/cov_xy + crit_x^2/var_x) 
  return(crit)
}

FHN.classic_contrast <- function(theta, x, y, delta) {
  n = length(x)
  gamma = theta[1]
  beta = theta[2]
  eps = theta[3]
  sigma2 = theta[4]
  s = 0.01
  crit_x = x[2:n]- x[1:(n-1)] - delta/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - delta^2/(2*eps)*((1-3*x[1:(n-1)]^2)/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  crit_y = y[2:n] - y[1:(n-1)]- delta*(gamma*x[1:(n-1)]-y[1:(n-1)]+beta) - delta^2/2*(gamma/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  varx_c <- delta^3/(3*eps^2)
  covxy_c <- delta^2/(3*eps) 
  vary_c <- delta 
    
  crit = 0.5*sum(crit_y^2/vary_c - crit_y*crit_x/covxy_c + crit_x^2/varx_c)/sigma2  - (n-2)*log(eps) + (n-2)*log(sigma2)
  return(crit)
}

FHN.vf_contrast <- function(theta, x, y, delta) {
  n = length(x)
  gamma = theta[1]
  beta = theta[2]
  eps = theta[3]
  s = 0.01
  crit_x = x[2:n]- x[1:(n-1)] - delta/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - delta^2/(2*eps)*((1-3*x[1:(n-1)]^2)/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  crit_y = y[2:n] - y[1:(n-1)]- delta*(gamma*x[1:(n-1)]-y[1:(n-1)]+beta) - delta^2/2*(gamma/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  return(sum(crit_x^2 + crit_y^2))
}

FHN.compute_eps <- function(theta0, x, y, delta){
  n = length(x)
  Vm = x[1:(n-1)]; Vp = x[2:n]; 
  Hm = y[1:(n-1)]; Hp = y[2:n]; 
  gamma = theta0[1]
  beta = theta0[2]
  eps_true = theta0[3]
  s = theta0[4]
  sigma2 = theta0[5]^2
  Zm1 = Vm-Vm^3-Hm*(1-delta/2)+s-beta*delta/2-gamma*delta/2*Vm
  Zm2 = (1-3*Vm^2)*(Vm-Vm^3-Hm+s)
  S <- numeric()
  S[1] =  -2*sum((Vp-Vm)^2) #for 1/eps
  S[2] = delta*3*sum(Zm1*(Vp-Vm))
  S[3] = delta^2*sum(-Zm1^2+Zm2*(Vp-Vm)+2/3*delta*sigma2^2)
  S[4] = -delta^3/2*sum(Zm1*Zm2)
  root<-polyroot(S[1:4])
  whichReal <- (abs(Im(root))< 0.0001)
  root1 <- Re(root)
  root2 <- root1[whichReal]
  rootp<-root2[root2>0]
  eps = 1/min(rootp)
  # flog.info("True eps is %s, computed eps is %s", eps_true, eps)
}

FHN.compute_all <- function(theta, eps, x, y, delta){
  n = length(x)
  s = 0.01
  gamma = theta[1]
  beta = theta[2]
  sigma2 = theta[3]
  s = 0.01
  crit_y = y[2:n] - y[1:(n-1)]- delta*(gamma*x[1:(n-1)]-y[1:(n-1)]+beta) - delta^2/2*(gamma/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  vary_c <- delta 
  crit = sum(crit_y^2/(sigma2*delta)) + (n-2)*log(sigma2)
  return(crit)
}

FHN.sigma_explicit <- function(theta, x, y, delta){
  n = length(x)
  gamma = theta[1]
  beta = theta[2]
  eps = theta[3]
  s = 0.01
  crit_x = x[2:n]- x[1:(n-1)] - delta/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - delta^2/(2*eps)*((1-3*x[1:(n-1)]^2)/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  crit_y = y[2:n] - y[1:(n-1)]- delta*(gamma*x[1:(n-1)]-y[1:(n-1)]+beta) - delta^2/2*(gamma/eps*(x[1:(n-1)]-x[1:(n-1)]^3-y[1:(n-1)] + s) - (gamma*x[1:(n-1)]-y[1:(n-1)]+beta))
  varx_c <- delta^3/(3*eps^2)
  covxy_c <- delta^2/(3*eps) 
  vary_c <- delta 
  crit = 2*sum(crit_y^2/vary_c - crit_y*crit_x/covxy_c + crit_x^2/varx_c)/(n-2)
  return(sqrt(crit))
}

write_dataframe <- function(table, real_parameters, name){   # Make a dataframe with results given vector of estimated parameters for every trial, filter out unplausible values
  require(xtable)
  table_filter <- subset(table, table$V1 > 0 & table$V1 < 10 & table$V2 > 0 & table$V2 < 5& table$V3 > 0 & table$V3 < 1)
 # table_filter <- table
  real_parameters <- c(real_parameters[1:3], real_parameters[5])
  
  Values_mean <- round(apply(table_filter, 2, mean), digits = 3)
  Values_sd <- round(apply((table_filter - real_parameters)^2, 2, mean), digits = 3)
  Values <- rbind(Values_mean, Values_sd)
  
  print(xtable(Values, digits = 3), type.placement = "h!", file = paste(name,".tex"))
  
  mycol <- rgb(0, 0, 255, max = 255, alpha = 75)
  # mycol <- rgb(255, 0, 0, max = 255, alpha = 75)
  d1 <- density(table_filter[,1])
  d2 <- density(table_filter[,2])
  d3 <- density(table_filter[,3])
  d4 <- density(table_filter[,4])
  
  pdf(paste(name,".pdf"))
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  #hist(table_filter[,1], xlab = "", ylab = "", main = "Gamma")
  plot(d1, type = "l", main = "Gamma",  xlab = "", ylab = "")
  polygon(d1, col = mycol)
  abline(v = real_parameters[1], col = "red", lwd = 2)
  #hist(table_filter[,2], xlab = "", ylab = "", main = "Beta")
  plot(d2, type = "l", main = "Beta",  xlab = "", ylab = "")
  polygon(d2, col = mycol)
  abline(v = real_parameters[2], col = "red", lwd = 2)
  #hist(table_filter[,3],  xlab = "", ylab = "",  main = "Epsilon")
  plot(d3, type = "l", main = "Epsilon",  xlab = "", ylab = "")
  polygon(d3, col = mycol)
  abline(v = real_parameters[3], col = "red", lwd = 2)
  #hist(table_filter[,4],  xlab = "", ylab = "", main = "Sigma")
  plot(d4, type = "l", main = "Sigma",  xlab = "", ylab = "")
  polygon(d4, col = mycol)
  abline(v = real_parameters[4], col = "red", lwd = 2)
  dev.off()
}

###################################
#####    Executable part ##########
###################################

###  Setting logger information

require("futile.logger")
# require("ggplot2")
# require("Matrix")
date <- Sys.Date()          # Create a subdirectory with the current date
wd <- getwd()               # Save the name of the working directory 
path_to_logs <- file.path(wd,date)
dir.create(path_to_logs)
file <- paste(path_to_logs, "logfile", sep = "/")
#flog.appender(appender.file(file))
flog.appender(appender.tee(file))
flog.threshold(DEBUG)        # By default set threshold to INFO (because I can)
flog.debug("Debugging is on!") # Further instruction to use: use debugging ONLY for the long list of values (numerous iterations etc), for "formal" messages use 

###  Generate data to work with

# parameters are to be given as a vector in a form (gamma, beta, eps, s, sigma), then:
# FHN.simulate(parameters, N, delta) (with NO defaults!)
# We will have two "good" sets: c(1.5, 0.3, 0.1, 0.01, 0.6)
# and ... 

K = 100
real_parameters <- c(1.5, 0.3, 0.1, 0.01, 0.6) # first set
# real_parameters <- c(1.2, 1.3, 0.1, 0.01, 0.4) # second set (choose wisely!)

########
# Now we are going to look at the variance, because we can
flog.info("All the same contrasts [11.2017]")
N = 500000
delta = 0.001
X = matrix(nrow = K, ncol = N)
Y = matrix(nrow = K, ncol = N)
  for (k in 1:K){
    Z = FHN.simulate(parameters = real_parameters, N = N, delta = delta)
    X[k,] <- Z[1,]
    Y[k,] <- Z[2,]
  }
  
flog.info("Data is generated with N = %s, delta = %s", N, delta)

#### Let's have some plots ####
pdf("FHN1.pdf")
z = FHN.simulate(parameters = c(1.5, 0.3, 0.1, 0.01, 0.6), N = 5000, delta = 0.01)
x <- z[1,]
y <- z[2,]
layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE))
plot(x[100:length(x)], y[100:length(x)], xlab = "", ylab = "", main = "Phase portrait", type = "l")
plot(x, ylab = "", xlab = "", xaxt = "n", type = "l", main = "Membrane potential")
plot(y, ylab = "", xlab = "", xaxt = "n", type = "l", main = "Channel kinetic")
dev.off()
pdf("FHN2.pdf")
layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE))
z = FHN.simulate(parameters = c(1.2, 1.3, 0.1, 0.01, 0.4), N = 5000, delta = 0.01)
x <- z[1,]
y <- z[2,]
plot(x[100:length(x)], y[100:length(x)], xlab = "", ylab = "", main = "Phase portrait", type = "l")
plot(x, ylab = "", xlab = "", xaxt = "n", type = "l", main = "Membrane potential")
plot(y, ylab = "", xlab = "", xaxt = "n", type = "l", main = "Channel kinetic")
dev.off()

#### Let's compute some differences ######

DELTAs <- c(0.01, 0.01,  0.001, 0.001, 0.0001)
Ns <- c(500000, 5000000, 500000, 5000000, 5000000)

# We just want to illustrate how the variance is computed 
# from the observations of the first and the second coordinate: 

eps = 0.1 # real parameter, we must know that 
s = 0.01
sigma_bar <- matrix(ncol = 2, nrow = length(Ns))
for (n in c(1:length(Ns))){
  delta <- DELTAs[n]
  N <- Ns[n]
  Z = FHN.simulate(parameters = real_parameters, N = N, delta = delta)
  x <- Z[1,]
  y <- Z[2,]
  diff_x2 <- mean((x[2:N] - x[1:(N-1)])^2)
  diff_y2 <- mean((y[2:N] - y[1:(N-1)])^2)
  # Statistics from the second coordinate is obvious:
  sigma_bar[n, 1] = sqrt(diff_y2/delta)
  sigma_bar[n, 2] = sqrt(eps^2*3*diff_x2/(delta^3))
  
  drift <- delta/eps*(x-x^3-y+s)
  flog.debug("Sigma1 = %s, Sigma2 = %s, mean drift = %s", sqrt(diff_y2/delta), sqrt(eps^2*3*diff_x2/(delta^3)), mean(drift))
}

##### And here is where "normal" contrast starts ######

K = 100
N = 500000
delta = 0.001
X = matrix(nrow = K, ncol = N)
Y = matrix(nrow = K, ncol = N)
real_parameters <- c(1.5, 0.3, 0.1, 0.01, 0.6) # first set

for (k in 1:K){
  Z = FHN.simulate(parameters = real_parameters, N = N, delta = delta)
  X[k,] <- Z[1,]
  Y[k,] <- Z[2,]
  sigma_exp <- FHN.sigma_explicit(theta = real_parameters, x = X[k,], y = Y[k,], delta = delta)
  flog.debug("Here sigma goes %s", sigma_exp)
}

flog.info("Computing classical contrast for Set 1")
for (i in c(10)){
flog.info("For subsampling we take each %s-th point in the trajectory", i)
  in_values <- numeric()
  VarY <- numeric(K)
  optimals <- matrix(nrow = K, ncol = 4)
  
  X_sub = X[,seq(1,N,i)]
  Y_sub = Y[,seq(1,N,i)]
  
  N_sub = length(X_sub[1,])
  delta_sub = delta*i
  # for (k in 1:K){
  #   VarY[k] <- mean((Y_sub[k, 1:(N_sub-1)] - Y_sub[k, 2:N_sub])^2)/delta_sub
  # }
  # 
  # sigma2 = mean(VarY)
  flog.info("Trial: i = %s, Delta = %s, N = %s", i, delta_sub, N_sub)
  # flog.info("Sigma is equal to %s, eps is equal to %s", sqrt(sigma2), eps)
  flog.info("And now optimization on the same trial!")
  for (k in 1:K){
    in_values[1] <- real_parameters[1] + sample(c(-0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5),1)
    in_values[2] <- real_parameters[2] + sample(c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),1)
    in_values[3] <- real_parameters[3] + sample(c(-0.05, -0.01, 0.01, 0.05), 1)
    # in_values[4] <- real_parameters[5]^2 + sample(c(-0.3, -0.2, -0.1, 0, 0.1),1)
    in_values[4] <- real_parameters[5]^2 + sample(c(-0.1, 0, 0.1),1)
    opt <- optim(par = in_values, x = X_sub[k,], y = Y_sub[k,], delta = delta_sub, fn = FHN.classic_contrast, method = "CG")
    # flog.debug(paste(toString(opt$par), toString(sqrt(VarY[k]))))
    # optimals[k,] <- c(opt$par, sqrt(VarY[k]))
    optimals[k,] <- c(opt$par[1:3], sqrt(opt$par[4]))
    flog.debug(paste(toString(optimals[k,])))
  }
  write.csv(optimals, 
              file = paste(path_to_logs, format(Sys.time(), format="%H:%M:%S"), sep = "/"), 
              row.names = FALSE) # write csv file with ALL the parameters, estimated on each run
}
CC1 <- optimals

flog.info("Computing VF contrast for Set 1")
for (i in c(10)){
  flog.info("For subsampling we take each %s-th point in the trajectory", i)
  VarY <- numeric(K)
  in_values <- numeric()
  optimals <- matrix(nrow = K, ncol = 4)
  
  X_sub = X[,seq(1,N,i)]
  Y_sub = Y[,seq(1,N,i)]
  
  N_sub = length(X_sub[1,])
  delta_sub = delta*i
  for (k in 1:K){
     VarY[k] <- mean((Y_sub[k, 1:(N_sub-1)] - Y_sub[k, 2:N_sub])^2)/delta_sub
  }
  flog.info("Trial: i = %s, Delta = %s, N = %s", i, delta_sub, N_sub)
  flog.info("And now optimization on the same trial!")
  for (k in 1:K){
    in_values[1] <- real_parameters[1] + sample(c(-0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5),1)
    in_values[2] <- real_parameters[2] + sample(c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),1)
    in_values[3] <- real_parameters[3] + sample(c(-0.05, -0.01, 0.01, 0.05), 1)
    opt <- optim(par = in_values, x = X_sub[k,], y = Y_sub[k,], delta = delta_sub, fn = FHN.vf_contrast, method = "CG")
    flog.debug(paste(toString(opt$par), toString(sqrt(VarY[k]))))
    optimals[k,] <- c(opt$par, sqrt(VarY[k]))
  }
  write.csv(optimals, 
            file = paste(path_to_logs, format(Sys.time(), format="%H:%M:%S"), sep = "/"), 
            row.names = FALSE) # write csv file with ALL the parameters, estimated on each run
}
VF1 <- optimals

####### Now we will be trying to build density for both estimators: original and the "variance-free" #######

mycol1 <- rgb(0, 0, 255, max = 255, alpha = 75)
mycol2 <- rgb(255, 0, 0, max = 255, alpha = 75)
d1 <- density(CC1[,1])
d2 <- density(CC1[,2])
d3 <- density(CC1[,3])
d4 <- density(CC1[,4])

dv1 <- density(VF1[,1])
dv2 <- density(VF1[,2])
dv3 <- density(VF1[,3])
dv4 <- density(VF1[,4])

pdf(paste(path_to_logs, "densities1.pdf", sep = "/"))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(d1, type = "l", main = "Gamma",  xlab = "", ylab = "", ylim = c(0,1.5))
polygon(d1, col = mycol1)
polygon(dv1, col = mycol2)
abline(v = real_parameters[1], col = "red", lwd = 1)
plot(d2, type = "l", main = "Beta",  xlab = "", ylab = "", ylim = c(0, 2.4))
polygon(d2, col = mycol1)
polygon(dv2, col = mycol2)
abline(v = real_parameters[2], col = "red", lwd = 1)
plot(d3, type = "l", main = "Epsilon",  xlab = "", ylab = "", xlim = c(0.099, 0.10105))
polygon(d3, col = mycol1)
polygon(dv3, col = mycol2)
abline(v = real_parameters[3], col = "red", lwd = 1)
plot(d4, type = "l", main = "Sigma",  xlab = "", ylab = "", xlim = c(0.58, 0.75), ylim = c(0, 50))
polygon(d4, col = mycol1)
polygon(dv4, col = mycol2)
abline(v = real_parameters[5], col = "red", lwd = 1)
dev.off()
########################################################################################################

# Second set of the parameters

real_parameters <- c(1.2, 1.3, 0.1, 0.01, 0.4) # second set (choose wisely!)

for (k in 1:K){
  Z = FHN.simulate(parameters = real_parameters, N = N, delta = delta)
  X[k,] <- Z[1,]
  Y[k,] <- Z[2,]
}

flog.info("Computing classical contrast for Set 2")
for (i in c(10)){
  flog.info("For subsampling we take each %s-th point in the trajectory", i)
  in_values <- numeric()
  VarY <- numeric(K)
  optimals <- matrix(nrow = K, ncol = 4)
  
  X_sub = X[,seq(1,N,i)]
  Y_sub = Y[,seq(1,N,i)]
  
  N_sub = length(X_sub[1,])
  delta_sub = delta*i
  flog.info("Trial: i = %s, Delta = %s, N = %s", i, delta_sub, N_sub)
  # flog.info("Sigma is equal to %s, eps is equal to %s", sqrt(sigma2), eps)
  flog.info("And now optimization on the same trial!")
  for (k in 1:K){
    in_values[1] <- real_parameters[1] + sample(c(-0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5),1)
    in_values[2] <- real_parameters[2] + sample(c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),1)
    in_values[3] <- real_parameters[3] + sample(c(-0.05, -0.01, 0.01, 0.05), 1)
    # in_values[4] <- real_parameters[5]^2 + sample(c(-0.3, -0.2, -0.1, 0, 0.1),1)
    in_values[4] <- real_parameters[5]^2 + sample(c(-0.1, 0, 0.1),1)
    # opt <- optim(par = in_values, x = X_sub[k,], y = Y_sub[k,], delta = delta_sub, fn = FHN.classic_contrast, method = "CG")
    opt <- optim(par = c(in_values[1:2], in_values[4]), x = X_sub[k,], y = Y_sub[k,], eps = real_parameters[3], delta = delta_sub, fn = FHN.compute_all, method = "CG")
    # optimals[k,] <- c(opt$par[1:3], sqrt(opt$par[4]))
    optimals[k, ] <- c(opt$par[1:2], FHN.compute_eps(real_parameters, x = X_sub[k,], y = Y_sub[k,], delta = delta_sub), sqrt(opt$par[3]))
    flog.debug(paste(toString(optimals[k,])))
  }
  write.csv(optimals, 
            file = paste(path_to_logs, format(Sys.time(), format="%H:%M:%S"), sep = "/"), 
            row.names = FALSE) # write csv file with ALL the parameters, estimated on each run
}
# CC2 <- optimals
ad_su1 <- optimals

flog.info("Computing VF contrast for Set 2")
for (i in c(10)){
  flog.info("For subsampling we take each %s-th point in the trajectory", i)
  VarY <- numeric(K)
  in_values <- numeric()
  optimals <- matrix(nrow = K, ncol = 4)
  
  X_sub = X[,seq(1,N,i)]
  Y_sub = Y[,seq(1,N,i)]
  
  N_sub = length(X_sub[1,])
  delta_sub = delta*i
  for (k in 1:K){
    VarY[k] <- mean((Y_sub[k, 1:(N_sub-1)] - Y_sub[k, 2:N_sub])^2)/delta_sub
  }
  flog.info("Trial: i = %s, Delta = %s, N = %s", i, delta_sub, N_sub)
  flog.info("And now optimization on the same trial!")
  for (k in 1:K){
    in_values[1] <- real_parameters[1] + sample(c(-0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5),1)
    in_values[2] <- real_parameters[2] + sample(c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),1)
    in_values[3] <- real_parameters[3] + sample(c(-0.05, -0.01, 0.01, 0.05), 1)
    opt <- optim(par = in_values, x = X_sub[k,], y = Y_sub[k,], delta = delta_sub, fn = FHN.vf_contrast, method = "CG")
    flog.debug(paste(toString(opt$par), toString(sqrt(VarY[k]))))
    optimals[k,] <- c(opt$par, sqrt(VarY[k]))
  }
  write.csv(optimals, 
            file = paste(path_to_logs, format(Sys.time(), format="%H:%M:%S"), sep = "/"), 
            row.names = FALSE) # write csv file with ALL the parameters, estimated on each run
}
VF2 <- optimals

####### Write the results now ########

name <- paste(path_to_logs, "anyname", sep = "/") # insert any name of the saved table here
table <- read.csv(file = name)
write_dataframe(table = table, real_parameters = real_parameters, name = name) # write filtered dataframe
