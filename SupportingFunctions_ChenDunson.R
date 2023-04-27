# Wrapper functions for various substeps of Chen and Dunson 2003
# random effects selection model 
# finalized on April 27 2023, by Adam Kaplan 

# sub-functions for posterior sampler below: 

# sampling for fixed effects

#install.packages('truncnorm') <- for sampling from truncated normal
#install.packages('MASS') <- for sampling from multivariate normal
require(truncnorm)
require(MASS)

Sample_alpha <- function(BigGamma, lambda, b, y, XtX, Z, inv.priorA, sigma){
  A <- solve(sigma^(-1)*XtX + inv.priorA); # this part will be different with polya gamma 
  alpha <- A%*%(sigma^(-1)*t(X)%*%(y - rowSums(Z%*%diag(lambda)%*%BigGamma*b))) # no prior alpha
  alpha.samp <- mvrnorm(1, alpha, A)
  return(alpha.samp)
}

# sampling for error variance 
Sample_Sigma <- function(c_0, d_0, N, y, X, alpha, Z, BigGamma, lambda, b){
  a_hat <- c_0 + N/2
  b_hat <- d_0 + (sum((y - X%*%(alpha) - rowSums(Z%*%diag(lambda)%*%BigGamma*b))^2)/2)
  sigma <- 1/rgamma(1, shape = a_hat, rate = b_hat)
  return(sigma)
}

# sample vector b for each group, but 'long_b_i' creates
# an nrow(Z) matrix of bs (duplicating each b_i for the 
# observations that belong to that group id)
Sample_b <- function(Z, y,X,alpha,lambda, BigGamma, ids, sigma){
  id_lengths <- table(ids)
  v <-(Z%*%diag(lambda)%*%BigGamma)
  v_i <- lapply(unique(ids), function(x) v[ids == x,])
  H_i <- lapply(unique(ids), function(x) solve(sigma^(-1)*crossprod(v_i[[x]],v_i[[x]]) + diag(1, nrow = ncol(v_i[[x]]))))
  h_i <- lapply(unique(ids), function(x) sigma^(-1)*H_i[[x]]%*%(t(v_i[[x]]) %*% (y[ids == x] - c(X[ids == x,]%*%(alpha)))))
  b_i <- lapply(unique(ids), function(x) mvrnorm(1, h_i[[x]], H_i[[x]]))
  long_b_i <- do.call(rbind, b_i)
  long_b_i <- long_b_i[rep(1:nrow(long_b_i), times = id_lengths),]
  return(long_b_i)
}

# create the Us and Ts that relate to the effect of 
# gammas and lambdas, respectively. 
create_Us_Ts <- function(b, lambda, Z, BigGamma){
  u_ijs <- lapply(1:(q2-1), function(x) 
  {out1 <- lapply((x+1):q2, 
                  function(y){
                    outnum <- b[,x] * lambda[y] * Z[,y];         
                    labord <- c(y,x); return(list(outnum, labord));
                  }); 
  outord <- t(sapply(out1, function(x) x[[2]]))
  outnumb <- sapply(out1, function(x) x[[1]])
  return(list(outord, outnumb))
  })
  ordneed <- do.call(rbind,(sapply(u_ijs, function(x) x[[1]])))
  ordneed2 <- order(ordneed[,1], ordneed[,2]) # this gives the proper order for the 
  # u columns that will match the gammas' order
  
  u_cols <- do.call(cbind, lapply(u_ijs, function(x) x[[2]]))[,ordneed2]
  
  
  t_ijs <- lapply(1:q2, function(x){
    return(
      Z[,x] * (b[,x] + as.numeric(x > 1)*rowSums(t(as.matrix(t(b[,1:(x-1)])*BigGamma[x,1:(x-1)]))))
    )
  })
  return(list('us' = u_cols, 'ts' = do.call(cbind,t_ijs)))
}

Sample_Gammas <- function(us, inv.priorR, y, X,Z, alpha, b,gamma_0, lambda,q2, sigma){
  R <- solve(sigma^(-1) * crossprod(us,us) + inv.priorR)
  gamma_hat <- R%*%(sigma^(-1)*(t(us)%*%(y - X%*%(alpha) - rowSums(t(t(Z)*lambda)*b))) + inv.priorR%*%gamma_0)
  R_lambda_set <- which(lambda == 0)
  newGamma.b40 <- mvrnorm(1, gamma_hat, R); gamma_st <- newGamma.b40
  BigGamma <- diag(1, nrow = q2)
  # the lower.diag function does not match how Chen and Dunson aligned their 
  # gammas in the Gamma matrix. So had to do a for-loop
  if(q2 > 1){
    for(k in 2:q2){
      BigGamma[k,1:(k-1)] <- gamma_st[1:(k-1)]
      gamma_st <- gamma_st[-c(1:(k-1))]
    }
  }
  BigGamma[,R_lambda_set] <- BigGamma[R_lambda_set,] <- 0
  diag(BigGamma) <- 1

  # lower.tri matrix to match Dunson/Chen paper notation
  if(q2 > 1){
    gamma_back <- c()
    for(k in 2:q2){
      gamma_back <- c(gamma_back,BigGamma[k,1:(k-1)])
    }
  }else{gamma_back <- 1; BigGamma <- 1}
  return(list('BigGamma' = BigGamma, 'gammas' = gamma_back))
}

# collecting the errors in regressing on fixed effects
# and all but one lambda's effect
create_etas <- function(y, X, alpha, ts, lambda, q2){
  etas <- lapply(1:q2, function(x){
    errors_l <- y - X%*%(alpha) - ((ts[,-x])%*%lambda[-x])
    return(errors_l)
  })
  return(do.call(cbind,etas)) 
}

# computing probability of random effect being zero
# and hyperparameters for lambda's sampling step
prob_samps <- function(ts, etas, sprior, mprior, p_prior, q2, sigma){
  
  ps <- lapply(1:q2, function(l){
    etas_l <- etas[,l]
    sprior_l <- sprior[l]
    s_prior_sd <- sprior_l^(1/2); # standard deviation in prior variance 
    m_prior_l <- mprior[l]; # prior guess at eigenvalue-like mean
    p_prior_l <- p_prior[l]; # prior probability of random effect being excluded 
    
    # what the paper says to do
    omega_ls <- ts[,l]%*%ts[,l]/sigma;
    lambda_tilde <- (ts[,l] %*%etas[,l])/(ts[,l]%*%ts[,l]);
    sigma_hat_l <- (omega_ls + sprior_l^(-1))^(-1) # maximum likelihood estimate for lambda_l_hat
    lambda_hat <- sigma_hat_l*(omega_ls*lambda_tilde + sprior_l^(-1) * m_prior_l ) # questionable code because the paper gave 
    
    # either way you get same exact results
    # between the way above and below
    
    # double checked for posterior derivation below, results match
    # what was supplied in paper. 
    #post_var_lambdal <- (sprior_l^(-1) + ts[,l]%*%ts[,l]/sigma)^(-1)
    #regress_lambdal <- (sprior_l^(-1)*m_prior_l + sigma^(-1)*ts[,l]%*%etas[,l])
    #post_mean_lambdal <- post_var_lambdal * regress_lambdal
    
    # 3 extra steps to calculate the posterior value of lambda_l hat, powers and labels in paper are wrong. Fantastic. 
    a = (- (etas[,l] %*% etas[,l])/(2*sigma))
    
    b1 = log(sigma_hat_l^(1/2)) + pnorm(m_prior_l/s_prior_sd, log.p = T) - (log(s_prior_sd) - pnorm(lambda_hat/sigma_hat_l, log.p = T)) 
    #if(denom_prob > .99999){denom_prob <- .999}
    #b1 = log((sigma_hat_l^(1/2)/s_prior_sd) * ((1-pnorm(-m_prior_l /s_prior_sd))/(1-denom_prob)))
    b2 = -sum((etas_l - c(lambda_tilde) * ts[,l])^2)/(2*sigma)
    #if(b2 < 1e-7){b2 <- 0.00001}
    b3 = -((lambda_tilde^2 * omega_ls/2) + (m_prior_l^2/(2*sprior_l)) - (lambda_hat^2/(2*sigma_hat_l)))
    
    b = (b1+b2+b3)
    
    # some log and exponential tricks to avoid small digits rounding floating error stuff
    log_p_l <-  - log(1 + exp(-log(p_prior_l)- a +log(1-p_prior_l) + b))
    
    p_l = exp(log_p_l);
    if(p_l < 0.0001){p_l = 0}
    if(p_l > 1){p_l = 1}
    return(list('lambda_hats' = lambda_hat, 'ps' = p_l, 'sigma_hats' = sigma_hat_l))
  })
  
  p_1s <- sapply(ps, function(x) x$ps)
  lambda_hats <- sapply(ps, function(x) x$lambda_hats)
  sigma_hats <- sapply(ps, function(x) x$sigma_hats)
  return(list('ps' = p_1s, 'lambda_hats' = lambda_hats, 'Sigma_hats_lambda' = sigma_hats))
}

# sample for lambda from zero-inflated truncated normal density
Lambda_Samp <- function(ps, lambda_hats, sigma_hats,q2){
  Upper <- rep(Inf, q2)
  Lower <- rep(0, q2)
  lam_samp <- rtruncnorm(1, Lower, Upper, mean = lambda_hats, sd = sqrt(sigma_hats))
  indic_zero <- as.logical(rbinom(n = length(ps), size = 1, prob = ps))
  lam_samp[indic_zero] <- 0
  return(lam_samp)
}

########################################################################################################
# sampler below: 

Bayesian_LMM_Selection <- function(y, X, Z,G, mprior, sprior, p_prior, ids, NBurn, NSamps,c_0 = .05, d_0 = .05, thin){
  
  q2 <- dim(Z)[2]; q1 <- dim(X)[2]
  #store XtX before model for loop 
  XtX <- crossprod(X,X)
  N <- nrow(X)
  
  # initialize lambda, BigGamma, and b and other parameters
  # using the prior precision for R supplied in paper
  lambda <- rtruncnorm(q2, a = 0, b = Inf, mean = 0, sd = 1)
  n_lowertriag <- q2*(q2-1)/2
  gammas <- rnorm(n_lowertriag, 0, 1)
  BigGamma <- diag(1, nrow = q2)
  BigGamma[lower.tri(BigGamma)] <- gammas
  b <- sapply(1:q2, function(x) rep(rnorm(G), each = n_i))
  priorA <- diag(1000, nrow = q1)
  inv.priorA <- solve(priorA)
  priorR <- diag(.5, nrow = n_lowertriag)
  inv.priorR <- solve(priorR)
  gamma_0 <- rep(0, n_lowertriag)
  rand_mat_id_index <- unlist(sapply(1:q2, function(x) rep(x,x-1)))
  alpha <- rnorm(q1)
  sigma<-1
  
  # Storage 
  alpha_col <- matrix(NA, ncol = q1, nrow = NSamps)
  lambda_col <- matrix(NA, ncol = q2, nrow = NSamps)
  gamma_col <- matrix(NA, ncol = n_lowertriag, nrow = NSamps)
  sigma_col <- matrix(NA, ncol = 1, nrow = NSamps)
  rand_vars <- matrix(NA, ncol = q2, nrow = NSamps)
  rand_offdiags <- matrix(NA, ncol = n_lowertriag, nrow = NSamps)
  ps_col <- matrix(NA, ncol = q2, nrow = NSamps)
  corrs_stor <- matrix(NA, ncol = n_lowertriag, nrow = NSamps)
  
  
  # Begin Sampler 
  start_collect <- 0
  for(k in 1:(NSamps + NBurn)){
    
    alpha <-Sample_alpha(BigGamma, lambda, b, y, XtX, Z, inv.priorA, sigma)
    
    sigma <- Sample_Sigma(c_0, d_0,N, y, X, alpha, Z, BigGamma, lambda, b)
    
    b <- Sample_b(Z, y,X,alpha,lambda, BigGamma, ids, sigma)
    
    Ts_Us <- create_Us_Ts(b, lambda, Z, BigGamma)
    
    us <- Ts_Us$us; ts <- Ts_Us$ts
    
    gammas <- Sample_Gammas(us, inv.priorR, y, X, Z, alpha,b, gamma_0, lambda,q2, sigma)
    
    BigGamma <- gammas$BigGamma; gamma <- gammas$gammas
    
    etas <- create_etas(y, X, alpha, ts, lambda, q2)
    
    ps <- prob_samps(ts, etas, sprior, mprior, p_prior, q2, sigma)
    
    if(sum(is.na(ps$ps))>0){stop('stopped for missing p')}
    
    lambda_hats <- ps$lambda_hats; sigma_hats <- c(ps$Sigma_hats_lambda); ps <- ps$ps
    
    lambda <- Lambda_Samp(ps, lambda_hats, sigma_hats,q2)
    
    if(k > NBurn){
      start_collect <- start_collect+1
      alpha_col[start_collect,] <- alpha
      lambda_col[start_collect,] <- lambda
      gamma_col[start_collect,] <- gamma
      ps_col[start_collect,] <- ps
      sigma_col[start_collect] <- sigma
      
      # the lower.diag function does not match how Chen and Dunson aligned their 
      # gammas in the Gamma matrix. So had to do a for-loop
      # prepping D matrix 
      lambda_mat <- diag(lambda); 
      gammas_mat <- diag(1,nrow = q2)
      gammast <- gamma
      for(k in 2:q2){
        gammas_mat[k,1:(k-1)] <- gammast[1:(k-1)]
        gammast <- gammast[-c(1:(k-1))]
      }
      
      outD <- lambda_mat %*% gammas_mat %*% t(gammas_mat) %*% lambda_mat
      
      offdiags <- c()
      for(k in 2:q2){
        offdiags <- c(offdiags, outD[k,1:(k-1)])
      }
 
      variances_rand <- diag(outD)
      corr_mat <- suppressWarnings(cov2cor(outD))
      corr_mat[is.na(corr_mat)] <- 0
      corrs_vec <- c()
      for(k in 2:q2){
        corrs_vec <- c(corrs_vec, corr_mat[k,1:(k-1)])
      }
      
      rand_vars[start_collect,] <- variances_rand
      rand_offdiags[start_collect,] <- offdiags
      corrs_stor[start_collect,] <- corrs_vec
    }
  }
  
  for_out <- list('alpha_samps' = alpha_col, 'lambda_samps' = lambda_col, 'gamma_samps' = gamma_col, 'p_samps' = ps_col,
                  "DiagonalD" = rand_vars, "OffDiagonalD" = rand_offdiags, 'Sigma' = sigma_col,
                  "Correlations" = corrs_stor)
  if(thin > 1){for_out <- lapply(for_out, function(x) x[c(rep(FALSE, thin-1), TRUE),])}else{
    for_out <- for_out}
  
  return(for_out)
}