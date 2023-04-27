setwd('~')

source("SupportingFunctions_ChenDunson.R")

# simulate the data scenario in Chen and Dunson's paper 

n_i <- 8
G <- 200;
ids <- rep(1:G, each = n_i)
N <- length(ids)
#X <- cbind(rep(1, N), rnorm(N), rnorm(N), rnorm(N))
X <- cbind(rep(1, N), runif(N, -2, 2), runif(N, -2, 2), runif(N, -2, 2))
q1 <- ncol(X)

# random intercept and slope 

Z <- X[,c(1,2,3,4)]
q2 <- ncol(Z)

lambda_true <- c(3, 1.2, .8,0); 
gamma_true <- c(1.33, .25, .71, rep(0,3)); trueGam <- diag(rep(1,q2), ncol = q2); 
trueGam[lower.tri(trueGam)] <- gamma_true
trueGam[4,1] <- 0; trueGam[3,2] <- .71 # lower tri goes by column not row 
cov_rand <- diag(lambda_true) %*% trueGam %*% t(trueGam) %*% diag(lambda_true)
true_d <- diag(cov_rand)

#cov2cor(cov_rand)

rand_effs <- mvrnorm(G, rep(0,q2),Sigma = cov_rand)

randbetas <- rand_effs[rep(1:nrow(rand_effs), each = n_i),]

y <- X%*%c(1, 1, 1, 1) + rowSums(Z*randbetas) + rnorm(N, mean = 0, sd = 1)

#####################
######## Analysis ###
#####################

trialme <- Bayesian_LMM_Selection(y=y, X=X, Z=Z, G=G, mprior = rep(0,4), sprior = rep(30,4),
                                  p_prior = rep(.2,4), ids = ids, NBurn = 5000, NSamps = 30000,
                                  thin = 6)


# looking over summary statistics and convergence. 

par(mfrow = c(2,4))

true_vals_f <- rep(1, 4)
true_vals_r <- c(3, 1.2, .8,0)

sapply(1:4, function(x){
  titlesd <- paste0('Fixed Effect ',x,'\nTrue Value = ',true_vals_f[x])
  xlablab <- paste0("Fixed Effect ", x)
  hist(trialme$alpha_samps[,x], main = titlesd, las = 1, breaks = 100, xlab = xlablab)})
sapply(1:4, function(x){
  titlesd <- paste0('Random Effect Cholesky Diag ',x,'\nTrue Value = ',lambda_true[x])
  xlablab <- paste0("Cholesky Diagonal ", x)
  hist(trialme$lambda_samps[,x], main = titlesd, las = 1, breaks = 100, xlab = xlablab)})

# summary of probability(lambda = 0)
summary(trialme$p_samps)  

# posterior probability of having all random effects be nonzero
mean(apply(trialme$lambda_samps, 1, function(x) sum(x!=0)==4))
mean(apply(trialme$p_samps, 1, function(x) prod(1-x)))

# posterior probability of having first 3 out of 4 random effects be nonzero
mean(apply(trialme$lambda_samps[,1:3], 1, function(x) sum(x!=0)==3)*trialme$lambda_samps[,4]==0)
mean(apply(trialme$p_samps[,1:3], 1, function(x) prod(1-x))*(trialme$p_samps[,4]))

# posterior probability of having first 2 be nonzero and last two be zero
mean(apply(trialme$lambda_samps[,1:2], 1, function(x) sum(x!=0)==2)*
       apply(trialme$lambda_samps[,3:4], 1, function(x) sum(x==0)==2))

# posterior probability of having first 1 be nonzero and last three be zero
mean((trialme$lambda_samps[,1]!=0)*
       apply(trialme$lambda_samps[,2:4], 1, function(x) sum(x==0)==3))

sapply(1:4, function(x){
  titlesd <- paste0('Random Effect Cholesky Diag ',x,'\nTrue Value = ',round(true_d[x],3))
  xlablab <- paste0("Cholesky Diagonal ", x)
  plot(trialme$DiagonalD[,x], main = titlesd, las = 1, xlab = xlablab)})
sapply(1:4, function(x){
  titlesd <- paste0('Fixed Effect ',x,'\nTrue Value = ',true_vals_f[x])
  xlablab <- paste0("Fixed Effect ", x)
  plot(trialme$alpha_samps[,x], main = titlesd, las = 1, xlab = xlablab)})

sapply(1:4, function(x){
  titlesd <- paste0('Random Effect Diagonal\n',x,' True Value = ',lambda_true[x])
  xlablab <- paste0("Diagonal ", x)
  plot(trialme$lambda_samps[,x], main = titlesd, las = 1, xlab = xlablab)})
sapply(1:4, function(x){
  titlesd <- paste0('Random Effect Off-Diagonal\n',x,' True Value = ',gamma_true[x])
  xlablab <- paste0("Off-Diagonal ", x)
  plot(trialme$gamma_samps[,x], main = titlesd, las = 1, xlab = xlablab)})
# gammas are swapped. may need to check the order 

cor_trues <- c(.8, .2, .5,0,0,0)
cor_labs <- c('[2,1]', '[3,1]', '[3,2]', '[4,1]', '[4,2]', '[4,3]')
sapply(1:6, function(x){
  titlesd <- paste0('Correlation Entry ',cor_labs[x],'\nTrue Value = ',cor_trues[x])
  xlablab <- paste0("Correlation ", cor_labs[x])
  plot(trialme$Correlations[,x], main = titlesd, las = 1, xlab = xlablab)})

summary(trialme$Sigma)

apply(trialme$Correlations, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
apply(trialme$DiagonalD, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
